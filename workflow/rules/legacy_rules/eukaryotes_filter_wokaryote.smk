rule eukaryotes_filter:
    input:
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/euk_potential_mags_length/{sample}_length_filtered_list.txt",
            ),
            sample=SAMPLES,
        ),
        os.path.join(RESULTS_DIR, "Bins/eukaryotes_filter_mags"),
        os.path.join(RESULTS_DIR, "Bins/eukaryotic_mags_whokaryote.txt"),
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
            directory(os.path.join(RESULTS_DIR, "Bins/euk_potential_mags/{sample}"))
        ),
        done=os.path.join(RESULTS_DIR, "Bins/euk_potential_mags/{sample}/done.txt"),
        concoct=temp(
            os.path.join(
                RESULTS_DIR, "Bins/euk_potential_mags/{sample}/concoct_list.txt"
            )
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/euk_potential/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/euk_potential/{sample}.err.log"),
    message:
        "Collecting all the MAGs from different binning tools"
    conda:
        os.path.join(WORK_DIR, "envs/concoct.yaml")
    params:
        value="{sample}",
    shell:
        """
        set -euo pipefail
        date
        mkdir -p {output.dir}
        cp {input.metabat}/*.fa {output.dir}
        cp $(dirname {input.rosella})/*.fna {output.dir}
        for file in "$(dirname {input.vamb})"/bins/*.fna; do
            base=$(basename "$file" .fna)
            cp "$file" "{output.dir}/vamb_${{base}}.fa"
        done
        perl -pe 's/,(\\d+)/,concoct_$1/g' {input.concoct} > {output.concoct}
        extract_fasta_bins.py {input.contigs} {output.concoct} --output_path {output.dir}
        cp $(dirname {input.semibin})/output_bins/*.fa.gz {output.dir}
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
                    RESULTS_DIR, "Bins/euk_potential_mags_renamed/{sample}_renamed"
                )
            )
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/rename_mags/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/rename_mags/{sample}.err.log"),
    message:
        "Renaming the MAGs"
    params:
        value="{sample}",
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
                "Bins/euk_potential_mags_length/{sample}_stats_fasta_all.txt",
            )
        ),
        final=os.path.join(
            RESULTS_DIR,
            "Bins/euk_potential_mags_length/{sample}_length_filtered_list.txt",
        ),
    conda:
        os.path.join(WORK_DIR, "envs/seqkit.yaml")
    threads: 4
    message:
        "Get the length of all the MAGs"
    log:
        out=os.path.join(RESULTS_DIR, "logs/get_length/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/get_length/{sample}.err.log"),
    shell:
        """
        date
        seqkit stats -j {threads} {input}/*fa -T > {output.txt}
        awk '$5>=2500000 && $7 >= 2500' {output.txt} | cut -f1 | tail -n +2 > {output.final}
        date
        """


checkpoint get_filtered_mags:
    input:
        list=expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/euk_potential_mags_length/{sample}_length_filtered_list.txt",
            ),
            sample=SAMPLES,
        ),
        dirs=expand(
            os.path.join(
                RESULTS_DIR, "Bins/euk_potential_mags_renamed/{sample}_renamed"
            ),
            sample=SAMPLES,
        ),
    output:
        files=directory(os.path.join(RESULTS_DIR, "Bins/eukaryotes_filter_mags")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/get_filter_euk.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/get_filter_euk.err.log"),
    shell:
        """
        set -euo pipefail
        date
        mkdir -p {output.files}
        cat {input.list} | while read -r line; do
            base=$(basename "$line")
            cp $line {output.files}/$base
        done
        date
        """


def get_file_names_mm(wildcards):
    ck_output = checkpoints.get_filtered_mags.get(**wildcards).output[0]
    return expand(
        os.path.join(
            RESULTS_DIR, "Bins/whokaryote/collated_output/{mags}_percentage.txt"
        ),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )


rule whokaryote:
    input:
        os.path.join(RESULTS_DIR, "Bins/eukaryotes_filter_mags/{mags}.fa"),
    output:
        os.path.join(RESULTS_DIR, "Bins/whokaryote/{mags}/whokaryote_predictions_T.tsv"),
    conda:
        os.path.join(WORK_DIR, "envs/whokaryote.yaml")
    threads: 3
    resources:
        runtime=30,
        mem_mb=get_mem_mb,
    log:
        out=os.path.join(RESULTS_DIR, "logs/whokaryote/whokaryote_{mags}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/whokaryote/whokaryote_{mags}.err.log"),
    message:
        "Running WHOKARYOTE on {wildcards.mags}"
    shell:
        """
        date
        whokaryote.py --contigs {input} --outdir $(dirname {output}) --minsize 2500 --threads {threads}
        date
        """


rule percentage_whokaryote:
    input:
        pred=os.path.join(
            RESULTS_DIR, "Bins/whokaryote/{mags}/whokaryote_predictions_T.tsv"
        ),
    output:
        percent=os.path.join(
            RESULTS_DIR, "Bins/whokaryote/collated_output/{mags}_percentage.txt"
        ),
    log:
        os.path.join(RESULTS_DIR, "logs/collate/{mags}.log"),
    message:
        "Getting the percentage for contigs identified as 'eukaryote' in {wildcards.mags}"
    localrule: True
    run:
        import pandas as pd

        df = pd.read_csv(input.pred, sep="\t", header=0)

        count = pd.DataFrame(df.groupby("predicted")["predicted"].count())

        count["percent"] = (count["predicted"] / count["predicted"].sum()) * 100

        count_trans = count.transpose()
        count_trans["bin"] = os.path.basename(os.path.dirname(input.pred))
        count_trans_edited = count_trans.drop(labels="predicted")
        count_trans_edited.to_csv(output.percent, sep="\t", header=True)


rule concatenate_whokaryote:
    input:
        all_files=get_file_names_mm,
    output:
        collated=os.path.join(
            RESULTS_DIR, "Bins/whokaryote/collated_output/collated_percentages.txt"
        ),
        final_list=os.path.join(RESULTS_DIR, "Bins/eukaryotic_mags_whokaryote.txt"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/all_files_collation_whokaryote.log"),
    message:
        "Collecting all the percentages into one file"
    run:
        import pandas as pd
        import contextlib

        with open(log.out, "w") as f_log:
            with contextlib.redirect_stdout(f_log), contextlib.redirect_stderr(f_log):
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
                eukaryotic_bins = concatenated[concatenated["eukaryote"] >= 80]
                eukaryotic_bins[["bin"]].to_csv(
                    output.final_list, sep="\t", index=False
                )
                print(f"Eukaryotic bins written to: {output.final_list}")
