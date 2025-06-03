rule eukaryotes_filter:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/euk_potential_mags_length/{sample}_length_filtered_list.txt"), sample=ASSEMBLY),
        os.path.join(RESULTS_DIR, "Bins/eukaryotes_filter_mags"),
        os.path.join(RESULTS_DIR, "Bins/eukaryotic_mags_whokaryote.txt"),
        os.path.join(RESULTS_DIR, "Bins/eukaryotic_mags_deepmicroclass.txt"),
    output:
        touch("status/eukaryotes_filter.done")

rule get_all_mags:
    input:
        contigs=rules.filter_length.output.out,
        metabat=rules.metabat2.output,
        rosella=rules.rosella.output,
        # taxvamb=rules.taxvamb_euk.output,
        concoct=rules.concoct.output,
        semibin=rules.semibin.output,
    output:
        dir=temp(directory(os.path.join(RESULTS_DIR, "Bins/euk_potential_mags/{sample}"))),
        done=os.path.join(RESULTS_DIR, "Bins/euk_potential_mags/{sample}/done.txt"),
        concoct=temp(os.path.join(RESULTS_DIR, "Bins/euk_potential_mags/{sample}/concoct_list.txt")),
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
        (
        set -euo pipefail
        date
        mkdir -p {output.dir}
        cp {input.metabat}/*.fa {output.dir}
        cp $(dirname {input.rosella})/*.fna {output.dir}
        perl -pe 's/,(\\d+)/,concoct_$1/g' {input.concoct} > {output.concoct}
        extract_fasta_bins.py {input.contigs} {output.concoct} --output_path {output.dir}
        cp $(dirname {input.semibin})/output_bins/*.fa.gz {output.dir}
        for file in {output.dir}/*gz; do
            if [[ -f "$file" && "$file" == *.gz ]]; then
                gunzip "$file"
            fi
        done
        for file in {output.dir}/*.fna; do
            if [[ -f "$file" ]]; then
                mv "$file" "${{file%.fna}}.fa"
            fi
        done
        touch {output.done}
        date) 2> {log.err} > {log.out}
        """

rule rename_mags:
    input:
        rules.get_all_mags.output.dir,
    output:
        temp(directory(os.path.join(RESULTS_DIR, "Bins/euk_potential_mags_renamed/{sample}_renamed"))),
    message:
        "Renaming the MAGs"
    log:
        out=os.path.join(RESULTS_DIR, "logs/rename_mags/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/rename_mags/{sample}.err.log"),
    params:
        value="{sample}",
    shell:
        """
        (
        set -euo pipefail
        date
        mkdir -p {output}
        for file in {input}/*fa; do
            base=$(basename "$file")
            mv "$file" "{output}/{params}_$base"
        done
        date) 2> {log.err} > {log.out}
        """

rule get_length:
    input:
        rules.rename_mags.output,
    output:
        txt=temp(os.path.join(RESULTS_DIR, "Bins/euk_potential_mags_length/{sample}_stats_fasta_all.txt")),
        final=os.path.join(RESULTS_DIR, "Bins/euk_potential_mags_length/{sample}_length_filtered_list.txt"),
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
        (date && seqkit stats -j {threads} {input}/*fa -T > {output.txt} &&
        awk '$5>=2500000 && $7 >= 2500' {output.txt} | cut -f1 | tail -n +2 > {output.final} && 
        date) 2> {log.err} > {log.out}
        """


checkpoint get_filtered_mags:
    input:
        list=expand(
            os.path.join(RESULTS_DIR, "Bins/euk_potential_mags_length/{sample}_length_filtered_list.txt"),
            sample=ASSEMBLY,
        ),
        dirs = expand(
            os.path.join(RESULTS_DIR, "Bins/euk_potential_mags_renamed/{sample}_renamed"),
            sample=ASSEMBLY,
        ),
    output:
        files=directory(os.path.join(RESULTS_DIR, "Bins/eukaryotes_filter_mags")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/get_filter_euk.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/get_filter_euk.err.log"),
    shell:
        """
        (
        set -euo pipefail
        date
        mkdir -p {output.files}
        cat {input.list} | while read -r line; do
            base=$(basename "$line")
            cp $line {output.files}/$base
        done
        date) 2> {log.err} > {log.out}
        """


def get_file_names_mm(wildcards):
    ck_output = checkpoints.get_filtered_mags.get(**wildcards).output[0]
    return expand(
        os.path.join(RESULTS_DIR, "Bins/whokaryote/collated_output/{mags}_percentage.txt"),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )

def get_file_names_mm_deepmicroclass(wildcards):
    ck_output = checkpoints.get_filtered_mags.get(**wildcards).output[0]
    return expand(
        os.path.join(RESULTS_DIR, "Bins/deepmicroclass/collated_output/{mags}_percentage.txt"),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )

# Running WHOKARYOTE
rule whokaryote:
    input:
        os.path.join(RESULTS_DIR, "Bins/eukaryotes_filter_mags/{mags}.fa"),
    output:
        os.path.join(RESULTS_DIR, "Bins/whokaryote/{mags}/tiara_pred.txt"),
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
        "(date && whokaryote.py --contigs {input} --outdir $(dirname {output}) --minsize 2500 --threads {threads} && date) 2> {log.err} > {log.out}"


rule percentage_whokaryote:
    input:
        pred=os.path.join(RESULTS_DIR, "Bins/whokaryote/{mags}/tiara_pred.txt"),
    output:
        percent=os.path.join(
            RESULTS_DIR, "Bins/whokaryote/collated_output/{mags}_percentage.txt"
        ),
    log:
        os.path.join(RESULTS_DIR, "logs/collate/{mags}.log"),
    message:
        "Getting the percentage for contigs identified as 'eukarya' in {wildcards.mags}"
    localrule: True
    run:
        import pandas as pd
        df = pd.read_csv(input.pred, sep="\t", header=0)

        count = pd.DataFrame(df.groupby("class_fst_stage")["class_fst_stage"].count())

        count["percent"] = (
            count["class_fst_stage"] / count["class_fst_stage"].sum()
        ) * 100

        count_trans = count.transpose()
        count_trans["bin"] = os.path.basename(os.path.dirname(input.pred))
        count_trans_edited = count_trans.drop(labels="class_fst_stage")
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
                eukaryotic_bins = concatenated[concatenated["eukarya"] >= 80]
                eukaryotic_bins[["bin"]].to_csv(output.final_list, sep="\t", index=False)
                print(f"Eukaryotic bins written to: {output.final_list}")



rule deepmicroclass:
    input:
        os.path.join(RESULTS_DIR, "Bins/eukaryotes_filter_mags/{mags}.fa"),
    output:
        os.path.join(RESULTS_DIR, "Bins/deepmicroclass/{mags}.fa_pred_one-hot_hybrid.tsv"),
    conda:
        os.path.join(WORK_DIR, "envs/deepmicroclass.yaml")
    threads: 3
    resources:
        runtime=30,
        mem_mb=get_mem_mb,
    log:
        out=os.path.join(RESULTS_DIR, "logs/deepmicroclass/deepmicroclass_{mags}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/deepmicroclass/deepmicroclass_{mags}.err.log"),
    message:
        "Running DeepMicroClass on {wildcards.mags}"
    shell:
        "(date && DeepMicroClass predict --input {input} --output_dir $(dirname {output}) --device cpu --cpu_thread {threads} && date) 2> {log.err} > {log.out}"


rule deepmicroclass_percentage:
    input:
        pred=os.path.join(RESULTS_DIR, "Bins/deepmicroclass/{mags}.fa_pred_one-hot_hybrid.tsv"),
    output:
        percent=os.path.join(
            RESULTS_DIR, "Bins/deepmicroclass/collated_output/{mags}_percentage.txt"
        ),
    log:
        os.path.join(RESULTS_DIR, "logs/collate/deepmicroclass_{mags}.log"),
    message:
        "Getting the percentage for contigs identified as 'eukarya' in {wildcards.mags}"
    localrule: True
    run:
        import pandas as pd

        df = pd.read_csv(input.pred, sep="\t", header=0)

        choice_labels = {
            1: "Eukaryote",
            2: "Eukaryote",
            3: "Prokaryote",
            4: "Prokaryote",
            5: "Prokaryote"
        }

        df["pred"] = df["best_choice"].map(choice_labels)

        count = pd.DataFrame(df.groupby("pred")["pred"].count())

        count["percent"] = (
                    count["pred"] / count["pred"].sum()
                ) * 100

        count_trans = count.transpose()
        count_trans["bin"] = os.path.basename(input.pred)
        count_trans_edited = count_trans.drop(labels="pred")
        count_trans_edited.to_csv(output.percent, sep="\t", header=True)

rule concatenate_deepmicroclass:
    input:
        all_files=get_file_names_mm_deepmicroclass,
    output:
        collated=os.path.join(
            RESULTS_DIR, "Bins/deepmicroclass/collated_output/collated_percentages.txt"
        ),
        final_list=os.path.join(RESULTS_DIR, "Bins/eukaryotic_mags_deepmicroclass.txt"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/all_files_collation_deepmicroclass.log"),
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
                eukaryotic_bins = concatenated[concatenated["Eukaryote"] >= 80]
                eukaryotic_bins[["bin"]].to_csv(output.final_list, sep="\t", index=False)
                print(f"Eukaryotic bins written to: {output.final_list}")