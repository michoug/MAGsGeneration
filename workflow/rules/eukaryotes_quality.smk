rule eukaryotes_quality:
    input:
        os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotes_filter_mags"),
        # os.path.join(
        #     RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_deepmicroclass.txt"
        # ),
        # os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_busco_done.txt"),
        # os.path.join(RESULTS_DIR, "Bins/Eukaryotic_mags_dereplicated/done.txt"),
    output:
        touch("status/eukaryotes_quality.done"),


rule get_all_mags:
    input:
        contigs=rules.filter_euk.output,
        metabat=rules.euk_metabat2.output,
        rosella=rules.euk_rosella.output,
        concoct=rules.euk_concoct.output,
        semibin=rules.euk_semibin.output,
        vamb=rules.euk_taxvamb.output,
        remag=rules.remag.output.bins,
    output:
        dir=directory(
            temp(
                os.path.join(
                    RESULTS_DIR, "Bins/intermediate/euk_potential_mags/{sample}"
                )
            )
        ),
        done=os.path.join(
            RESULTS_DIR, "Bins/intermediate/euk_potential_mags/{sample}/done.txt"
        ),
        concoct=temp(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/euk_potential_mags/{sample}/concoct_list.txt",
            )
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/euk_potential/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/euk_potential/{sample}.err.log"),
    params:
        value="{sample}",
    message:
        "Collecting all the MAGs from different binning tools"
    shell:
        """
        set -euo pipefail
        (
            date
            mkdir -p {output.dir}
            for file in {input.metabat}/*.fa; do
                [ -e "$file" ] || continue
                cp "$file" {output.dir}
            done
            # Skip Rosella's leftover contigs, they are not a bin
            for file in "$(dirname {input.rosella})"/*.fna; do
                [ -e "$file" ] || continue
                [[ "$(basename "$file")" == *_unbinned.fna ]] && continue
                cp "$file" {output.dir}
            done
            for file in "$(dirname {input.vamb})"/bins/*.fna; do
                [ -e "$file" ] || continue
                base=$(basename "$file" .fna)
                cp "$file" "{output.dir}/vamb_${{base}}.fa"
            done
            perl -pe 's/,(\\d+)/,concoct_$1/g' {input.concoct} >{output.concoct}
            pixi r -e concoct extract_fasta_bins.py {input.contigs} {output.concoct} --output_path {output.dir}
            for file in "$(dirname {input.semibin})"/output_bins/*.fa.gz; do
                [ -e "$file" ] || continue
                cp "$file" {output.dir}
            done
            for file in {output.dir}/*.gz; do
                [ -e "$file" ] || continue
                gunzip "$file"
            done
            for file in {output.dir}/*.fna; do
                [ -e "$file" ] || continue
                mv "$file" "${{file%.fna}}.fa"
            done
            for file in {input.remag}/*.fa; do
                [ -e "$file" ] || continue
                cp "$file" "{output.dir}/remag_$(basename "$file")"
            done
            touch {output.done}
            date
        ) >{log.out} 2>{log.err}
        """


rule rename_mags:
    input:
        rules.get_all_mags.output.dir,
    output:
        temp(
            directory(
                os.path.join(
                    RESULTS_DIR,
                    "Bins/intermediate/euk_potential_mags_renamed/{sample}_renamed",
                )
            )
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/rename_mags/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/rename_mags/{sample}.err.log"),
    localrule: True
    params:
        value="{sample}",
    message:
        "Renaming the MAGs"
    shell:
        """
        set -euo pipefail
        date
        mkdir -p {output}
        for file in {input}/*fa; do
            base=$(basename "$file")
            mv "$file" "{output}/{params}_$base"
        done
        date
        """


rule get_length:
    input:
        rules.rename_mags.output,
    output:
        txt=temp(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/euk_potential_mags_length/{sample}_stats_fasta_all.txt",
            )
        ),
        final=os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/euk_potential_mags_length/{sample}_length_filtered_list.txt",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/get_length/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/get_length/{sample}.err.log"),
    threads: 4
    message:
        "Get the length of all the MAGs"
    shell:
        """
        date
        pixi r -e seqkit seqkit stats -j {threads} {input}/*fa -T >{output.txt}
        awk '$5>=2500000 && $7 >= 2500' {output.txt} | cut -f1 | tail -n +2 >{output.final}
        date
        """


checkpoint get_filtered_mags:
    input:
        list=expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/euk_potential_mags_length/{sample}_length_filtered_list.txt",
            ),
            sample=ASSEMBLY,
        ),
        dirs=expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/euk_potential_mags_renamed/{sample}_renamed",
            ),
            sample=ASSEMBLY,
        ),
    output:
        files=directory(
            os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotes_filter_mags")
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/get_filter_euk.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/get_filter_euk.err.log"),
    shell:
        """
        set -euo pipefail
        (
            date
            mkdir -p {output.files}
            cat {input.list} | while read -r line; do
                base=$(basename "$line")
                cp $line {output.files}/$base
            done
            date
        ) 2>{log.err} >{log.out}
        """


# checkpoint get_eukaryotes_mags:
#     input:
#         folder=rules.get_filtered_mags.output.files,
#         deep=os.path.join(
#             RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_deepmicroclass.txt"
#         ),
#         remag=rules.get_length_for_remag.output.final,
#     output:
#         folder=directory(os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotic_mags")),
#         done=os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotic_mags/done.txt"),
#     log:
#         out=os.path.join(RESULTS_DIR, "logs/get_eukaryotic_mags.out.log"),
#         err=os.path.join(RESULTS_DIR, "logs/get_eukaryotic_mags.err.log"),
#     localrule: True
#     message:
#         "Getting the eukaryotic MAGs"
#     shell:
#         """
#         set -euo pipefail
#         (
#             date
#             rm -rf {output.folder}
#             mkdir -p {output.folder}
#             cat {input.deep} | while read -r line; do
#                 base=$(basename "$line")
#                 cp {input.folder}/$base {output.folder}/$base
#             done
#             while read -r line; do
#                 base=$(basename "$line")
#                 cp $line {output.folder}/$base
#             done <{input.remag}
#             touch {output.done}
#             date
#         ) 2>{log.err} >{log.out}
#         """


# def get_file_eukaryotes(wildcards):
#     ck_output = checkpoints.get_eukaryotes_mags.get(**wildcards).output[0]
#     return expand(
#         os.path.join(
#             RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_busco/{mags}/{mags}.done"
#         ),
#         mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
#     )


# rule busco:
#     input:
#         folder=os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotic_mags/{mags}.fa"),
#         db=rules.busco_db_euk.output.dir,
#     output:
#         folder=directory(
#             os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_busco/{mags}")
#         ),
#         done=os.path.join(
#             RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_busco/{mags}/{mags}.done"
#         ),
#     log:
#         out=os.path.join(RESULTS_DIR, "logs/busco/busco_{mags}.out.log"),
#         err=os.path.join(RESULTS_DIR, "logs/busco/busco_{mags}.err.log"),
#     threads: 10
#     resources:
#         runtime=2 * 60,
#         mem_mb=get_mem_mb,
#     message:
#         "Running BUSCO on the eukaryotic MAGs"
#     shell:
#         """
#         (
#             date
#             rm -rf {output.folder}
#             pixi r -e busco busco -i {input.folder} --out_path {output.folder} -m genome --auto-lineage-euk -c 10 -f --metaeuk --download_path {input.db} --offline || echo "BUSCO failed (possibly no genes found) — continuing anyway"
#             touch {output.done}
#             date
#         ) 2>{log.err} >{log.out}
#         """


# rule busco_done:
#     input:
#         get_file_eukaryotes,
#     output:
#         os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_busco_done.txt"),
#     localrule: True
#     message:
#         "All BUSCO runs on eukaryotic MAGs are done"
#     shell:
#         """
#         touch {output}
#         """


# rule busco_multiqc:
#     input:
#         rules.busco_done.output,
#     output:
#         html=os.path.join(RESULTS_DIR, "qc/busco/multiqc_report.html"),
#         stats=os.path.join(
#             RESULTS_DIR, "qc/busco/multiqc_data/multiqc_general_stats.txt"
#         ),
#         source=os.path.join(RESULTS_DIR, "qc/busco/multiqc_data/multiqc_sources.txt"),
#         folder=directory(os.path.join(RESULTS_DIR, "qc/busco/")),
#     log:
#         err=os.path.join(RESULTS_DIR, "logs/multiqc_busco.err.log"),
#         out=os.path.join(RESULTS_DIR, "logs/multiqc_busco.out.log"),
#     resources:
#         runtime=60,
#         mem_mb=get_mem_mb,
#     message:
#         "Running MultiQC on BUSCO results"
#     shell:
#         """
#         (
#             date
#             rm -rf {output.folder}
#             mkdir -p {output.folder}/multiqc_data

#             pixi r -e multiqc multiqc $(dirname {input[0]})/eukaryotic_mags_busco \
#                 --verbose --ignore auto_lineage/ -o {output.folder} --force || true

#             if [ ! -s {output.html} ]; then
#                 printf '%s\\n' \
#                     '<!doctype html>' \
#                     '<html lang="en">' \
#                     '<head><meta charset="utf-8"><title>BUSCO MultiQC</title></head>' \
#                     '<body><p>No BUSCO results were found for eukaryotic MAGs.</p></body>' \
#                     '</html>' >{output.html}
#             fi
#             if [ ! -s {output.stats} ]; then
#                 printf "Sample\\tbusco-complete\\tbusco-complete_single_copy\\tbusco-fragmented\\tbusco-missing\\tbusco-complete_duplicated\\n" >{output.stats}
#             fi
#             if [ ! -s {output.source} ]; then
#                 printf "Sample Name\\tSource\\n" >{output.source}
#             fi
#             date
#         ) 2>{log.err} >{log.out}
#         """


# rule quality_eukaryotes:
#     input:
#         stats=rules.busco_multiqc.output.stats,
#         source=rules.busco_multiqc.output.source,
#         mags_done=rules.get_eukaryotes_mags.output.done,
#     output:
#         os.path.join(RESULTS_DIR, "qc/busco_euk.csv"),
#     log:
#         out=os.path.join(RESULTS_DIR, "logs/quality_eukaryotes.out.log"),
#         err=os.path.join(RESULTS_DIR, "logs/quality_eukaryotes.err.log"),
#     localrule: True
#     message:
#         "Extracting quality metrics for eukaryotic MAGs from MultiQC report"
#     run:
#         from pathlib import Path
#         import pandas as pd

#         df = pd.read_csv(input.stats, sep="\t")
#         src = pd.read_csv(input.source, sep="\t")
#         required_cols = {
#             "Sample",
#             "busco-complete",
#             "busco-complete_single_copy",
#             "busco-fragmented",
#             "busco-missing",
#         }
#         missing_cols = required_cols - set(df.columns)
#         if missing_cols:
#             raise ValueError(f"Missing required columns: {sorted(missing_cols)}")
#         required_src_cols = {"Sample Name", "Source"}
#         missing_src_cols = required_src_cols - set(src.columns)
#         if missing_src_cols:
#             raise ValueError(
#                 f"Missing required columns in sources file: {sorted(missing_src_cols)}"
#             )
#         busco_total = (
#             df["busco-complete_single_copy"]
#             + df["busco-fragmented"]
#             + df["busco-missing"]
#         )
#         if "busco-complete_duplicated" in df.columns:
#             duplicated = df["busco-complete_duplicated"]
#         else:
#             duplicated = (
#                 df["busco-complete"] - df["busco-complete_single_copy"]
#             ).clip(lower=0)
#         src["sample_key"] = src["Sample Name"].str.strip()
#         src["genome"] = (
#             src["Source"]
#             .str.extract(r"/eukaryotic_mags_busco/([^/]+)/", expand=False)
#             .fillna(
#                 src["Source"].str.extract(
#                     r"/BUSCO_([^/]+)\.fa(?:/|$)", expand=False
#                 )
#             )
#         )
#         sample_to_genome = (
#             src.dropna(subset=["genome"])
#             .drop_duplicates(subset=["sample_key"])
#             .set_index("sample_key")["genome"]
#         )
#         df["sample_key"] = df["Sample"].str.strip()
#         df["genome"] = df["sample_key"].map(sample_to_genome)
#         fallback = (
#             df["sample_key"]
#             .str.replace(
#                 r"^short_summary\.specific\.eukaryota_odb12\.", "", regex=True
#             )
#             .str.replace(r"^BUSCO_", "", regex=True)
#             .str.replace(r"\.fa(?:\.txt)?$", "", regex=True)
#         )
#         df["genome"] = df["genome"].fillna(fallback)
#         df["completeness"] = (
#             (df["busco-complete"] + df["busco-fragmented"]) * 100 / busco_total
#         )
#         df["contamination"] = duplicated * 100 / busco_total
#         result = (
#             df[["genome", "completeness", "contamination"]]
#             .sort_values("completeness", ascending=False)
#             .drop_duplicates(subset="genome", keep="first")
#         )
#         mags_folder = Path(input.mags_done).parent
#         all_genomes = {p.stem for p in mags_folder.glob("*.fa")}
#         present_genomes = set(result["genome"].astype(str))
#         missing_genomes = sorted(all_genomes - present_genomes)
#         if missing_genomes:
#             missing_df = pd.DataFrame(
#                 {
#                     "genome": missing_genomes,
#                     "completeness": 0.0,
#                     "contamination": 0.0,
#                 }
#             )
#             result = pd.concat([result, missing_df], ignore_index=True)
#         result = result.sort_values("genome")
#         result.to_csv(output[0], index=False)

# rule galah_euk:
#     input:
#         bins=rules.get_eukaryotes_mags.output.folder,
#         busco=rules.quality_eukaryotes.output,
#     output:
#         final=directory(os.path.join(RESULTS_DIR, "Bins/Eukaryotic_mags_dereplicated")),
#         done=os.path.join(RESULTS_DIR, "Bins/Eukaryotic_mags_dereplicated/done.txt"),
#         cluster=os.path.join(RESULTS_DIR, "Bins/intermediate/cluster_galah_euk.txt"),
#     log:
#         out=os.path.join(RESULTS_DIR, "logs/galah/galah_euk.out.log"),
#         err=os.path.join(RESULTS_DIR, "logs/galah/galah_euk.err.log"),
#     threads: 10
#     resources:
#         runtime=60 * 2,
#         mem_mb=get_mem_mb,
#     params:
#         comp=70,
#         cont=10,
#     message:
#         "Running Galah on all Eukaryotic MAGs (skip clustering if BUSCO table is empty)"
#     shell:
#         """
#         (
#             date
#             if [[ $(wc -l <"{input.busco}") -le 1 ]]; then
#                 echo "No BUSCO quality rows available; skipping galah clustering."
#                 mkdir -p {output.final}
#                 : >{output.cluster}
#             else
#                 pixi r -e galah galah cluster --genome-fasta-directory {input.bins} --genome-info {input.busco} -x fa --output-representative-fasta-directory-copy {output.final} -t {threads} --output-cluster-definition {output.cluster}
#             fi
#             touch {output.done}
#             date
#         ) 2>{log.err} >{log.out}
#         """
