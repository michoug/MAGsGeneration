def euk_assemblies(wildcards):
    """Assemblies with at least one eukaryotic contig according to DeepMicroClass.

    Assemblies without eukaryotic contigs are skipped by the eukaryotic binning,
    as the mapping and binning tools fail on an empty fasta.
    """
    samples = []
    for sample in ASSEMBLY:
        lst = checkpoints.deepmicroclass_list.get(sample=sample).output.lst
        with open(lst) as f:
            if any(line.strip() for line in f):
                samples.append(sample)
    return samples


def euk_expand(path, **kwargs):
    """expand() over the assemblies with eukaryotic contigs, as an input function."""
    return lambda wildcards: expand(path, sample=euk_assemblies(wildcards), **kwargs)


rule euk_mags_generation:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR,
                "Assembly/intermediate/{sample}_filter.fasta_pred_one-hot_hybrid.tsv",
            ),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR,
                "Assembly/intermediate/{sample}_deepmicroclass_list_euk.txt",
            ),
            sample=ASSEMBLY,
        ),
        euk_expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter_euk.fasta"),
        ),
        euk_expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles_euk/{sample}/{sample}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        euk_expand(
            os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/Metabat_euk/"),
        ),
        euk_expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Rosella_euk/UMAP_projection_with_clusters_final.png",
            ),
        ),
        euk_expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Vamb_euk/Vamb_Result/vaevae_clusters_unsplit.tsv",
            ),
        ),
        euk_expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Concoct_euk/concoct_clustering_merged.csv",
            ),
        ),
        euk_expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/Semibin_euk/contig_bins.tsv"
            ),
        ),
        euk_expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_remag/bins.csv"
            ),
        ),
    output:
        touch("status/euk_mags_generation.done"),


rule deepmicroclass:
    input:
        rules.filter_length.output.out,
    output:
        os.path.join(
            RESULTS_DIR,
            "Assembly/intermediate/{sample}_filter.fasta_pred_one-hot_hybrid.tsv",
        ),
    log:
        out=os.path.join(
            RESULTS_DIR, "logs/deepmicroclass/deepmicroclass_{sample}.out.log"
        ),
        err=os.path.join(
            RESULTS_DIR, "logs/deepmicroclass/deepmicroclass_{sample}.err.log"
        ),
    threads: 10
    resources:
        runtime=60 * 4,
        mem_mb=get_mem_mb,
    message:
        "Running DeepMicroClass on {wildcards.sample} to predict eukaryotic contigs"
    shell:
        """
        (
            date
            pixi r -e deepmicroclass DeepMicroClass predict --input {input} --output_dir $(dirname {output}) --device cpu --cpu_thread {threads}
            date
        ) 2>{log.err} >{log.out}
        """


checkpoint deepmicroclass_list:
    input:
        deep=os.path.join(
            RESULTS_DIR,
            "Assembly/intermediate/{sample}_filter.fasta_pred_one-hot_hybrid.tsv",
        ),
    output:
        lst=os.path.join(
            RESULTS_DIR,
            "Assembly/intermediate/{sample}_deepmicroclass_list_euk.txt",
        ),
    log:
        os.path.join(
            RESULTS_DIR, "logs/deepmicroclass/deepmicroclass_{sample}_list.log"
        ),
    localrule: True
    message:
        "Getting the contigs identified as 'eukarya' in {wildcards.sample}"
    run:
        import pandas as pd

        df = pd.read_csv(input.deep, sep="\t", header=0)
        # best_choice is 1-based: 1 Eukaryote, 2 EukaryoteVirus, 3 Plasmid,
        # 4 Prokaryote, 5 ProkaryoteVirus. Eukaryotic viruses are kept with the
        # eukaryotes on purpose.
        choice_labels = {
            1: "Eukaryote",
            2: "Eukaryote",
            3: "Prokaryote",
            4: "Prokaryote",
            5: "Prokaryote",
        }
        df["pred"] = df["best_choice"].map(choice_labels)
        filtered_df = df[df["pred"] == "Eukaryote"]
        filtered_df_list = filtered_df["Sequence Name"]
        filtered_df_list.to_csv(output.lst, sep="\t", header=False, index=False)


rule filter_euk:
    input:
        assembly=rules.filter_length.output.out,
        lst=rules.deepmicroclass_list.output.lst,
    output:
        os.path.join(RESULTS_DIR, "Assembly/{sample}_filter_euk.fasta"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/filter_length_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/filter_length_euk/{sample}.err.log"),
    threads: 4
    resources:
        mem_mb=get_mem_mb,
    message:
        "Get Eukaryotic contigs from {wildcards.sample} based on DeepMicroClass results"
    shell:
        """
        (
            date
            pixi r -e seqkit seqkit faidx -U -l {input.lst} -j {threads} -o {output} {input.assembly}
            date
        ) 2>{log.err} >{log.out}
        """


use rule prok_indexAssembly as euk_indexAssembly with:
    input:
        fa=rules.filter_euk.output,
    output:
        temp(
            os.path.join(
                RESULTS_DIR, "Assembly/{sample}_filter_euk.fasta.r{number}.sti"
            )
        ),
    log:
        out=os.path.join(
            RESULTS_DIR,
            "logs/indexAssembly_euk/indexAssembly_{sample}_{number}.out.log",
        ),
        err=os.path.join(
            RESULTS_DIR,
            "logs/indexAssembly_euk/indexAssembly_{sample}_{number}.err.log",
        ),


rule euk_check_reads_assembly_for_mapping:
    input:
        expand(
            os.path.join(RESULTS_DIR, "TrimFinal/{sample}_1.fq.gz"),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(RESULTS_DIR, "TrimFinal/{sample}_2.fq.gz"),
            sample=SAMPLES,
        ),
        euk_expand(
            os.path.join(
                RESULTS_DIR, "Assembly/{sample}_filter_euk.fasta.r{number}.sti"
            ),
            number=[100, 125, 150, 250, 400],
        ),
    output:
        os.path.join(RESULTS_DIR, "Assembly/euk_check_reads_for_mapping.done"),
    localrule: True
    shell:
        "touch {output}"


use rule prok_mapping as euk_mapping with:
    input:
        check_reads=rules.euk_check_reads_assembly_for_mapping.output,
        read1=lambda wildcards: os.path.join(
            RESULTS_DIR, f"TrimFinal/{wildcards.sample_2}_1.fq.gz"
        ),
        read2=lambda wildcards: os.path.join(
            RESULTS_DIR, f"TrimFinal/{wildcards.sample_2}_2.fq.gz"
        ),
        fasta=rules.filter_euk.output,
        index=euk_expand(
            os.path.join(
                RESULTS_DIR, "Assembly/{sample}_filter_euk.fasta.r{number}.sti"
            ),
            number=[100, 125, 150, 250, 400],
        ),
    output:
        bam=os.path.join(
            RESULTS_DIR, "intermediate/BamFiles_euk/{sample}/{sample}_{sample_2}.bam"
        ),
        bai=os.path.join(
            RESULTS_DIR,
            "intermediate/BamFiles_euk/{sample}/{sample}_{sample_2}.bam.bai",
        ),
    # group:
    #     "mapping_group"
    log:
        out=os.path.join(RESULTS_DIR, "logs/mapping_euk/{sample}_{sample_2}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/mapping_euk/{sample}_{sample_2}.err.log"),


rule euk_check_bam_files_for_binning:
    input:
        euk_expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles_euk/{sample}/{sample}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
    output:
        os.path.join(
            RESULTS_DIR, "intermediate/BamFiles_euk/check_bam_files_for_binning.done"
        ),
    localrule: True
    shell:
        "touch {output}"


use rule prok_metabat2_mapping as euk_metabat2_mapping with:
    input:
        bam=expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles_euk/{{sample}}/{{sample}}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        check=rules.euk_check_bam_files_for_binning.output,
    output:
        os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_metabat_cov_euk.txt"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabat_mapping_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabat_mapping_euk/{sample}.err.log"),


use rule prok_metabat2 as euk_metabat2 with:
    input:
        contig=rules.filter_euk.output,
        cov=rules.euk_metabat2_mapping.output,
    output:
        directory(os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/Metabat_euk/")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabat2_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabat2_euk/{sample}.err.log"),


use rule prok_clean_coverage as euk_clean_coverage with:
    input:
        bam=expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles_euk/{{sample}}/{{sample}}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        check=rules.euk_check_bam_files_for_binning.output,
    output:
        temp=temp(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/{sample}_clean_coverage_temp_euk.txt",
            )
        ),
        final=os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/{sample}/{sample}_clean_coverage_cov_euk.txt",
        ),
    log:
        out=os.path.join(
            RESULTS_DIR, "logs/clean_coverage_mapping_euk/{sample}.out.log"
        ),
        err=os.path.join(
            RESULTS_DIR, "logs/clean_coverage_mapping_euk/{sample}.err.log"
        ),


use rule prok_rosella as euk_rosella with:
    input:
        cont=rules.filter_euk.output,
        cov=expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles_euk/{{sample}}/{{sample}}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        check=rules.euk_check_bam_files_for_binning.output,
    output:
        os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/{sample}/Rosella_euk/UMAP_projection_with_clusters_final.png",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/rosella_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/rosella_euk/{sample}.err.log"),


use rule prok_vamb_coverage as euk_vamb_coverage with:
    input:
        rules.euk_clean_coverage.output.final,
        rules.euk_check_bam_files_for_binning.output,
    output:
        os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Vamb_euk/{sample}_aemb_cov.txt"
        ),


use rule prok_vamb_length as euk_vamb_length with:
    input:
        rules.filter_euk.output,
    output:
        temp(
            os.path.join(
                RESULTS_DIR, "Assembly/intermediate/{sample}_filter_vamb_euk.fasta"
            )
        ),


use rule prok_vamb_filter as euk_vamb_filter with:
    input:
        rules.tax_add_missing.output,
        rules.euk_vamb_length.output,
    output:
        temp(os.path.join(RESULTS_DIR, "Assembly/intermediate/{sample}_tax_euk.tsv")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/vamb_filter_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/vamb_filter_euk/{sample}.err.log"),


use rule prok_taxvamb as euk_taxvamb with:
    input:
        contig=rules.euk_vamb_length.output,
        tax=rules.euk_vamb_filter.output,
        cov=rules.euk_vamb_coverage.output,
    output:
        os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/{sample}/Vamb_euk/Vamb_Result/vaevae_clusters_unsplit.tsv",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/vamb_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/vamb_euk/{sample}.err.log"),


use rule prok_concoct_prepare as euk_concoct_prepare with:
    input:
        contig=rules.filter_euk.output,
        bam=expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles_euk/{{sample}}/{{sample}}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        check=rules.euk_check_bam_files_for_binning.output,
    output:
        bed=os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_10k_euk.bed"),
        cont=os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_10k_euk.fa"),
        cov=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_concoct_cov_euk.txt"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/concoct_prepare_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/concoct_prepare_euk/{sample}.err.log"),


use rule prok_concoct as euk_concoct with:
    input:
        contig=rules.euk_concoct_prepare.output.cont,
        cov=rules.euk_concoct_prepare.output.cov,
    output:
        os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/{sample}/Concoct_euk/concoct_clustering_merged.csv",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/concoct_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/concoct_euk/{sample}.err.log"),


use rule prok_semibin as euk_semibin with:
    input:
        cont=rules.filter_euk.output,
        bam=expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles_euk/{{sample}}/{{sample}}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        check=rules.euk_check_bam_files_for_binning.output,
    output:
        os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Semibin_euk/contig_bins.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/semibin_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/semibin_euk/{sample}.err.log"),


rule remag:
    input:
        contig=rules.filter_euk.output,
        tsv=rules.euk_metabat2_mapping.output,
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
        (
            date
            pixi r -e remag remag --fasta {input.contig} --coverage {input.tsv} -o {output.folder} -t {threads}
            date
        ) 2>{log.err} >{log.out}
        """
