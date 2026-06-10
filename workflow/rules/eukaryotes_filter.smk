rule eukaryotes_filter:
    input:
        os.path.join(RESULTS_DIR, "Bins/intermediate/remag_all_bins/done"),
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/euk_potential_mags_length/{sample}_length_filtered_list.txt",
            ),
            sample=ASSEMBLY,
        ),
        os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotes_filter_mags"),
        os.path.join(
            RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_deepmicroclass.txt"
        ),
        os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_busco/done.txt"),
        os.path.join(RESULTS_DIR, "Bins/Eukaryotic_mags_dereplicated/done.txt"),
    output:
        touch("status/eukaryotes_filter.done"),


rule get_all_mags:
    input:
        contigs=rules.filter_length.output.out,
        metabat=rules.metabat2.output,
        rosella=rules.rosella.output,
        concoct=rules.concoct.output,
        semibin=rules.semibin.output,
        vamb=rules.taxvamb.output,
    output:
        dir=temp(
            directory(
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
        date
        mkdir -p {output.dir}
        cp {input.metabat}/*.fa {output.dir}
        cp $(dirname {input.rosella})/*.fna {output.dir}
        for file in "$(dirname {input.vamb})"/bins/*.fna; do
            [ -e "$file" ] || continue
            base=$(basename "$file" .fna)
            cp "$file" "{output.dir}/vamb_${{base}}.fa"
        done
        perl -pe 's/,(\\d+)/,concoct_$1/g' {input.concoct} > {output.concoct}
        pixi r -e concoct extract_fasta_bins.py {input.contigs} {output.concoct} --output_path {output.dir}
        for file in "$(dirname {input.semibin})"/output_bins/*.fa.gz; do
            [ -e "$file" ] || continue
            cp "$file" {output.dir}
        done
        for file in {output.dir}/*gz; do [[ -f "$file" && "$file" == *.gz ]] && gunzip "$file"; done
        for file in {output.dir}/*.fna; do [[ -f "$file" ]] && mv "$file" "${{file%.fna}}.fa"; done
        touch {output.done}
        date
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
        pixi r -e seqkit seqkit stats -j {threads} {input}/*fa -T > {output.txt}
        awk '$5>=2500000 && $7 >= 2500' {output.txt} | cut -f1 | tail -n +2 > {output.final}
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
        (date
        mkdir -p {output.files}
        cat {input.list} | while read -r line; do
            base=$(basename "$line")
            cp $line {output.files}/$base
        done
        date) 2> {log.err} > {log.out}
        """


def get_file_names_mm_deepmicroclass(wildcards):
    ck_output = checkpoints.get_filtered_mags.get(**wildcards).output[0]
    return expand(
        os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/deepmicroclass/collated_output/{mags}_percentage.txt",
        ),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )


rule deepmicroclass:
    input:
        os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotes_filter_mags/{mags}.fa"),
    output:
        os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/deepmicroclass/{mags}.fa_pred_one-hot_hybrid.tsv",
        ),
    log:
        out=os.path.join(
            RESULTS_DIR, "logs/deepmicroclass/deepmicroclass_{mags}.out.log"
        ),
        err=os.path.join(
            RESULTS_DIR, "logs/deepmicroclass/deepmicroclass_{mags}.err.log"
        ),
    group:
        "deepmicroclass"
    threads: 3
    resources:
        runtime=30,
        mem_mb=get_mem_mb,
    message:
        "Running DeepMicroClass on {wildcards.mags}"
    shell:
        "(date && pixi r -e deepmicroclass DeepMicroClass predict --input {input} --output_dir $(dirname {output}) --device cpu --cpu_thread {threads} && date) 2> {log.err} > {log.out}"


rule deepmicroclass_percentage:
    input:
        pred=os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/deepmicroclass/{mags}.fa_pred_one-hot_hybrid.tsv",
        ),
    output:
        percent=os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/deepmicroclass/collated_output/{mags}_percentage.txt",
        ),
    log:
        os.path.join(RESULTS_DIR, "logs/collate/deepmicroclass_{mags}.log"),
    localrule: True
    message:
        "Getting the percentage for contigs identified as 'eukarya' in {wildcards.mags}"
    run:
        import pandas as pd

        df = pd.read_csv(input.pred, sep="\t", header=0)

        choice_labels = {
            1: "Eukaryote",
            2: "Eukaryote",
            3: "Prokaryote",
            4: "Prokaryote",
            5: "Prokaryote",
        }

        df["pred"] = df["best_choice"].map(choice_labels)

        count = pd.DataFrame(df.groupby("pred")["pred"].count())

        count["percent"] = (count["pred"] / count["pred"].sum()) * 100

        count_trans = count.transpose()
        count_trans["bin"] = os.path.basename(input.pred)
        count_trans["bin"] = count_trans["bin"].str.replace(
            "_pred_one-hot_hybrid.tsv", ""
        )
        count_trans_edited = count_trans.drop(labels="pred")
        count_trans_edited.to_csv(output.percent, sep="\t", header=True)


rule concatenate_deepmicroclass:
    input:
        all_files=get_file_names_mm_deepmicroclass,
    output:
        collated=os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/deepmicroclass/collated_output/collated_percentages.txt",
        ),
        final_list=os.path.join(
            RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_deepmicroclass.txt"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/all_files_collation_deepmicroclass.log"),
    message:
        "Collecting all the percentages into one file"
    run:
        import pandas as pd
        import contextlib

        with open(log.out, "w") as f_log:
            with (
                contextlib.redirect_stdout(f_log),
                contextlib.redirect_stderr(f_log),
            ):
                dfs = []

                for file in input.all_files:
                    df = pd.read_csv(file, sep="\t", header=0)
                    print(f"Loaded file: {file}, shape: {df.shape}")
                    dfs.append(df)

                concatenated = pd.concat(dfs, ignore_index=True)
                print(f"Concatenated DataFrame shape: {concatenated.shape}")

                # Save the collated output
                concatenated.to_csv(output.collated, sep="\t", index=False)
                print(f"Collated output written to: {output.collated}")

                # Filter and save eukaryotic bins
                eukaryotic_bins = concatenated[concatenated["Eukaryote"] >= 80]
                eukaryotic_bins[["bin"]].to_csv(
                    output.final_list, sep="\t", index=False, header=False
                )
                print(f"Eukaryotic bins written to: {output.final_list}")


rule remag:
    input:
        contig=rules.filter_length.output.out,
        tsv=rules.metabat2_mapping.output,
    output:
        folder=directory(
            os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_remag/")
        ),
        file=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_remag/bins.csv"
        ),
        bins=directory(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_remag/bins/"
            )
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/remag/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/remag/{sample}.err.log"),
    threads: 10
    resources:
        runtime=60 * 10,
        mem_mb=get_mem_mb,
    message:
        "Running REMAG on {wildcards.sample}"
    shell:
        """
        (date
        pixi r -e remag remag --fasta {input.contig} --coverage {input.tsv} -o {output.folder} -t {threads}
        date) 2> {log.err} > {log.out}
        """


rule concat_remag:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_remag/bins.csv"
            ),
            sample=SAMPLES,
        ),
    output:
        folder=directory(os.path.join(RESULTS_DIR, "Bins/intermediate/remag_all_bins")),
        file=os.path.join(RESULTS_DIR, "Bins/intermediate/remag_all_bins/done"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/remag_concat.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/remag_concat.err.log"),
    localrule: True
    message:
        "Concatenating all REMAG bins into one folder"
    shell:
        r"""
        set -euo pipefail
        (date
        rm -rf {output.folder}
        mkdir -p {output.folder}
        for file in {input}; do
            sample=$(basename $(dirname $file))
            for bins in "$(dirname $file)/bins/"*fa; do
                [ -e "$bins" ] || continue
                cp "$bins" {output.folder}/$sample\_$(basename "$bins")
            done
        done
        touch {output.file}
        date) 2> {log.err} > {log.out}
        """


rule get_length_for_remag:
    input:
        rules.concat_remag.output.folder,
    output:
        txt=temp(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/remag_all_bins/remag_stats_fasta_all.txt",
            )
        ),
        final=os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/remag_all_bins/remag_length_filtered_list.txt",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/get_length_remag.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/get_length_remag.err.log"),
    threads: 4
    shell:
        """
        (date
        pixi r -e seqkit seqkit stats -j {threads} {input}/*fa -T > {output.txt}
        awk '$5>=2500000' {output.txt} | cut -f1 | tail -n +2 > {output.final}
        date) 2> {log.err} > {log.out}
        """


checkpoint get_eukaryotes_mags:
    input:
        folder=rules.get_filtered_mags.output.files,
        deep=os.path.join(
            RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_deepmicroclass.txt"
        ),
        remag=rules.get_length_for_remag.output.final,
    output:
        folder=directory(os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotic_mags")),
        done=os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotic_mags/done.txt"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/get_eukaryotic_mags.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/get_eukaryotic_mags.err.log"),
    localrule: True
    message:
        "Getting the eukaryotic MAGs"
    shell:
        """
        set -euo pipefail
        (date
        rm -rf {output.folder}
        mkdir -p {output.folder}
        cat {input.deep} | while read -r line; do
            base=$(basename "$line")
            cp {input.folder}/$base {output.folder}/$base
        done
        while read -r line; do
            base=$(basename "$line")
            cp $line {output.folder}/$base
        done < {input.remag}
        touch {output.done}
        date) 2> {log.err} > {log.out}
        """


def get_file_eukaryotes(wildcards):
    ck_output = checkpoints.get_eukaryotes_mags.get(**wildcards).output[0]
    return expand(
        os.path.join(
            RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_busco/{mags}/{mags}.done"
        ),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )


rule busco:
    input:
        folder=os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotic_mags/{mags}.fa"),
        db=rules.busco_db_euk.output.dir,
    output:
        folder=directory(
            os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_busco/{mags}")
        ),
        done=os.path.join(
            RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_busco/{mags}/{mags}.done"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/busco/busco_{mags}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/busco/busco_{mags}.err.log"),
    threads: 10
    resources:
        runtime=2 * 60,
        mem_mb=get_mem_mb,
    message:
        "Running BUSCO on the eukaryotic MAGs"
    shell:
        """
        (date
        rm -rf {output.folder}
        pixi r -e busco busco -i {input.folder} --out_path {output.folder} -m genome --auto-lineage-euk -c 10 -f --metaeuk --download_path {input.db} --offline || echo "BUSCO failed (possibly no genes found) — continuing anyway"
        touch {output.done}
        date) 2> {log.err} > {log.out}
        """


rule busco_done:
    input:
        get_file_eukaryotes,
    output:
        os.path.join(RESULTS_DIR, "Bins/intermediate/eukaryotic_mags_busco/done.txt"),
    localrule: True
    message:
        "All BUSCO runs on eukaryotic MAGs are done"
    shell:
        """
        touch {output}
        """


rule busco_multiqc:
    input:
        rules.busco_done.output,
    output:
        html=os.path.join(RESULTS_DIR, "qc/busco/multiqc_report.html"),
        stats=os.path.join(
            RESULTS_DIR, "qc/busco/multiqc_data/multiqc_general_stats.txt"
        ),
        source=os.path.join(RESULTS_DIR, "qc/busco/multiqc_data/multiqc_sources.txt"),
        folder=directory(os.path.join(RESULTS_DIR, "qc/busco/")),
    log:
        err=os.path.join(RESULTS_DIR, "logs/multiqc_busco.err.log"),
        out=os.path.join(RESULTS_DIR, "logs/multiqc_busco.out.log"),
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    message:
        "Running MultiQC on BUSCO results"
    shell:
        """
        (date
        rm -rf {output.folder}
        pixi r -e multiqc multiqc $(dirname {input}) --verbose --ignore auto_lineage/ -o {output.folder} || echo "Multiqc failed (possibly no eukaryotic bins found) — continuing anyway"
        date)   2> {log.err} > {log.out}
        """


rule quality_eukaryotes:
    input:
        stats=rules.busco_multiqc.output.stats,
        source=rules.busco_multiqc.output.source,
    output:
        os.path.join(RESULTS_DIR, "qc/busco_euk.csv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/quality_eukaryotes.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/quality_eukaryotes.err.log"),
    localrule: True
    message:
        "Extracting quality metrics for eukaryotic MAGs from MultiQC report"
    run:
        import pandas as pd

        df = pd.read_csv(input.stats, sep="\t")
        src = pd.read_csv(input.source, sep="\t")

        required_cols = {
            "Sample",
            "busco-complete",
            "busco-complete_single_copy",
            "busco-fragmented",
            "busco-missing",
        }
        missing_cols = required_cols - set(df.columns)
        if missing_cols:
            raise ValueError(f"Missing required columns: {sorted(missing_cols)}")

        required_src_cols = {"Sample Name", "Source"}
        missing_src_cols = required_src_cols - set(src.columns)
        if missing_src_cols:
            raise ValueError(
                f"Missing required columns in sources file: {sorted(missing_src_cols)}"
            )

        busco_total = (
            df["busco-complete_single_copy"]
            + df["busco-fragmented"]
            + df["busco-missing"]
        )

        if "busco-complete_duplicated" in df.columns:
            duplicated = df["busco-complete_duplicated"]
        else:
            duplicated = (
                df["busco-complete"] - df["busco-complete_single_copy"]
            ).clip(lower=0)

        src["sample_key"] = src["Sample Name"].str.strip()
        src["genome"] = (
            src["Source"]
            .str.extract(r"/eukaryotic_mags_busco/([^/]+)/", expand=False)
            .fillna(
                src["Source"].str.extract(
                    r"/BUSCO_([^/]+)\.fa(?:/|$)", expand=False
                )
            )
        )
        sample_to_genome = (
            src.dropna(subset=["genome"])
            .drop_duplicates(subset=["sample_key"])
            .set_index("sample_key")["genome"]
        )

        df["sample_key"] = df["Sample"].str.strip()
        df["genome"] = df["sample_key"].map(sample_to_genome)
        fallback = (
            df["sample_key"]
            .str.replace(
                r"^short_summary\.specific\.eukaryota_odb12\.", "", regex=True
            )
            .str.replace(r"^BUSCO_", "", regex=True)
            .str.replace(r"\.fa(?:\.txt)?$", "", regex=True)
        )
        df["genome"] = df["genome"].fillna(fallback)

        df["completeness"] = (
            (df["busco-complete"] + df["busco-fragmented"]) * 100 / busco_total
        )
        df["contamination"] = duplicated * 100 / busco_total

        df[["genome", "completeness", "contamination"]].drop_duplicates().to_csv(
            output[0], index=False
        )


rule galah_euk:
    input:
        bins=rules.get_eukaryotes_mags.output.folder,
        busco=rules.quality_eukaryotes.output,
    output:
        final=directory(os.path.join(RESULTS_DIR, "Bins/Eukaryotic_mags_dereplicated")),
        done=os.path.join(RESULTS_DIR, "Bins/Eukaryotic_mags_dereplicated/done.txt"),
        cluster=os.path.join(RESULTS_DIR, "Bins/intermediate/cluster_galah_euk.txt"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/galah/galah_euk.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/galah/galah_euk.err.log"),
    threads: 10
    resources:
        runtime=60 * 2,
        mem_mb=get_mem_mb,
    params:
        comp=70,
        cont=10,
    message:
        "Running Galah on all Eukaryotic MAGs"
    shell:
        """
        (date
        pixi r -e galah galah cluster --genome-fasta-directory {input.bins} --genome-info {input.busco} -x fa --output-representative-fasta-directory-copy {output.final} -t {threads}  --output-cluster-definition {output.cluster}
        touch {output.done}
        date) 2> {log.err} > {log.out}
        """
