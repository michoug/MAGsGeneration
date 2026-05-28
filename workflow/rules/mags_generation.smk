
rule mags_generation:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR, "intermediate/BamFiles/{sample}/{sample}_{sample_2}.bam"
            ),
            sample=ASSEMBLY,
            sample_2=SAMPLES,
        ),
        expand(
            os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/Metabat"),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Rosella/UMAP_projection_with_clusters_final.png",
            ),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Vamb/Vamb_Results/vaevae_clusters_unsplit.tsv",
            ),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Concoct/concoct_clustering_merged.csv",
            ),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Comebin/comebin_res/comebin_res.tsv",
            ),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/Semibin/contig_bins.tsv"
            ),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/MetaCAT/{sample}_metacat.done",
            ),
            sample=ASSEMBLY,
        ),
    output:
        touch("status/mags_generation.done"),


rule indexAssembly:
    input:
        fa=rules.filter_length.output.out,
    output:
        temp(os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta.r{number}.sti")),
    log:
        out=os.path.join(
            RESULTS_DIR, "logs/indexAssembly/indexAssembly_{sample}_{number}.out.log"
        ),
        err=os.path.join(
            RESULTS_DIR, "logs/indexAssembly/indexAssembly_{sample}_{number}.err.log"
        ),
    threads: 10
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    params:
        numb="{number}",
    shell:
        """
        (date 
        pixi r -e coverm strobealign --create-index -t {threads} {input.fa} -r {params.numb}
        date ) 2> {log.err} > {log.out}
        """


rule check_reads_assembly_for_mapping:
    input:
        expand(
            os.path.join(RESULTS_DIR, "TrimFinal/{sample}_1.fq.gz"),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(RESULTS_DIR, "TrimFinal/{sample}_2.fq.gz"),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta.r{number}.sti"),
            sample=ASSEMBLY,
            number=[100, 125, 150, 250, 400],
        ),
    output:
        os.path.join(RESULTS_DIR, "Assembly/check_reads_for_mapping.done"),
    localrule: True
    shell:
        "touch {output}"


rule mapping:
    input:
        check_reads=rules.check_reads_assembly_for_mapping.output,
        read1=lambda wildcards: os.path.join(
            RESULTS_DIR, f"TrimFinal/{wildcards.sample_2}_1.fq.gz"
        ),
        read2=lambda wildcards: os.path.join(
            RESULTS_DIR, f"TrimFinal/{wildcards.sample_2}_2.fq.gz"
        ),
        fasta=rules.filter_length.output.out,
        index=expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta.r{number}.sti"),
            sample=ASSEMBLY,
            number=[100, 125, 150, 250, 400],
        ),
    output:
        bam=os.path.join(
            RESULTS_DIR, "intermediate/BamFiles/{sample}/{sample}_{sample_2}.bam"
        ),
        bai=os.path.join(
            RESULTS_DIR, "intermediate/BamFiles/{sample}/{sample}_{sample_2}.bam.bai"
        ),
    # group:
    #     "mapping_group"
    log:
        out=os.path.join(RESULTS_DIR, "logs/mapping/{sample}_{sample_2}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/mapping/{sample}_{sample_2}.err.log"),
    threads: 8
    resources:
        runtime=60 * 1,
        mem_mb=get_mem_mb,
    shell:
        """
        (date
        rm -f {output.bam}.tmp*
        pixi r -e coverm strobealign -t {threads} --use-index {input.fasta} {input.read1} {input.read2} | 
        pixi r -e coverm samtools sort -@ {threads} -o {output.bam}
        pixi r -e coverm samtools index {output.bam}
        touch -r {output.bam} -d '+1 second' {output.bai}
        date ) 2> {log.err} > {log.out}
        """


rule check_bam_files_for_binning:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "intermediate/BamFiles/{sample}/{sample}_{sample_2}.bam"
            ),
            sample=ASSEMBLY,
            sample_2=SAMPLES,
        ),
    output:
        os.path.join(
            RESULTS_DIR, "intermediate/BamFiles/check_bam_files_for_binning.done"
        ),
    localrule: True
    shell:
        "touch {output}"


rule metabat2_mapping:
    input:
        bam=expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles/{{sample}}/{{sample}}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        check=rules.check_bam_files_for_binning.output,
    output:
        os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_metabat_cov.txt"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabat_mapping/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabat_mapping/{sample}.err.log"),
    threads: 10
    resources:
        runtime=60 * 2,
        mem_mb=get_mem_mb,
    message:
        "Running the first part of metabat2 to estimate the coverage"
    shell:
        "(date && pixi r -e metabat jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam} && "
        "date) 2> {log.err} > {log.out}"


rule metabat2:
    input:
        contig=rules.filter_length.output.out,
        cov=rules.metabat2_mapping.output,
    output:
        directory(os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/Metabat/")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabat2/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabat2/{sample}.err.log"),
    priority: 20
    threads: 10
    resources:
        runtime=60 * 2,
        mem_mb=get_mem_mb,
    message:
        "Running metabat2 to obtain bins"
    shell:
        "(date && pixi r -e metabat metabat2 -i {input.contig} -a {input.cov} -o {output}/metabat -t {threads} && "
        "date) 2> {log.err} > {log.out}"


rule metabinner_coverage:
    input:
        bam=expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles/{{sample}}/{{sample}}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        check=rules.check_bam_files_for_binning.output,
    output:
        temp=temp(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_metabinner_temp.txt"
            )
        ),
        final=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_metabinner_cov.txt"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabinner_mapping/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabinner_mapping/{sample}.err.log"),
    threads: 14
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    params:
        length=1500,
    message:
        "Running the first part of metabinner to estimate the coverage"
    shell:
        """
        (date && pixi r -e metabat jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth {output.temp} {input.bam} && 
        cat {output.temp} | awk '{{if ($2>{params.length}) print $0 }}' | cut -f -1,4- > {output.final} && 
        date) 2> {log.err} > {log.out}
        """


rule rosella:
    input:
        cont=rules.filter_length.output.out,
        cov=expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles/{{sample}}/{{sample}}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        check=rules.check_bam_files_for_binning.output,
    output:
        os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/{sample}/Rosella/UMAP_projection_with_clusters_final.png",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/rosella/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/rosella/{sample}.err.log"),
    threads: 28
    resources:
        runtime=60 * 48,
        mem_mb=get_mem_mb,
    message:
        "Running Rosella"
    shell:
        "(date && rm -r $(dirname {output}) && "
        "pixi r -e rosella rosella recover -b {input.cov} -r {input.cont} --threads {threads} -o $(dirname {output}) "
        "&& date) 2> {log.err} > {log.out}"


rule vamb_coverage:
    input:
        rules.metabinner_coverage.output.final,
        rules.check_bam_files_for_binning.output,
    output:
        os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Vamb/{sample}_aemb_cov.txt"
        ),
    localrule: True
    params:
        name="{sample}",
    run:
        import pandas as pd

        df = pd.read_csv(input[0], sep="\t")
        df.rename(columns={"contigName": "contigname"}, inplace=True)
        df.to_csv(output[0], sep="\t", index=False, header=True)


rule vamb_length:
    input:
        rules.filter_length.output.out,
    output:
        temp(
            os.path.join(
                RESULTS_DIR, "Assembly/intermediate/{sample}_filter_vamb.fasta"
            )
        ),
    localrule: True
    threads: 4
    message:
        "Remove contigs less than 1501 bp for Vamb"
    shell:
        "pixi r -e seqkit seqkit seq -j {threads} --remove-gaps -o {output} -m 1501 {input}"


rule vamb_filter:
    input:
        rules.tax_add_missing.output,
        rules.vamb_length.output,
    output:
        temp(os.path.join(RESULTS_DIR, "Assembly/intermediate/{sample}_tax.tsv")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/vamb_filter/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/vamb_filter/{sample}.err.log"),
    localrule: True
    run:
        import pandas as pd

        ids = []
        with open(input[1], "r") as file:
            for line in file:
                if line.startswith(">"):
                    ids.append(line[1:].strip())
        df = pd.read_csv(input[0], sep="\t")
        filtered_df = df[df.iloc[:, 0].isin(ids)]
        filtered_df.to_csv(output[0], sep="\t", index=False, header=True)


rule taxvamb:
    input:
        contig=rules.vamb_length.output,
        tax=rules.vamb_filter.output,
        cov=rules.vamb_coverage.output,
    output:
        os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/{sample}/Vamb/Vamb_Results/vaevae_clusters_unsplit.tsv",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/vamb/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/vamb/{sample}.err.log"),
    threads: 14
    resources:
        runtime=60 * 24,
        mem_mb=get_mem_mb,
    shell:
        "(date && "
        "rm -r $(dirname {output}) && "
        "pixi r -e vamb vamb bin taxvamb --outdir $(dirname {output}) -m 1501 --fasta {input.contig} --taxonomy {input.tax} --minfasta 500000 --abundance_tsv {input.cov} -p {threads} -o && "
        "date) 2> {log.err} > {log.out}"


rule concoct_prepare:
    input:
        contig=rules.filter_length.output.out,
        bam=expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles/{{sample}}/{{sample}}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        check=rules.check_bam_files_for_binning.output,
    output:
        bed=os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_10k.bed"),
        cont=os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_10k.fa"),
        cov=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/{sample}_concoct_cov.txt"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/concoct_prepare/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/concoct_prepare/{sample}.err.log"),
    threads: 10
    resources:
        runtime=60 * 10,
        mem_mb=get_mem_mb,
    params:
        script=os.path.join(WORK_DIR, "scripts/concoct_coverage_table_threads.py"),
    message:
        "Preparing files for Concoct"
    shell:
        "(date && rm -f {output.bed} {output.cont} {output.cov} && pixi r -e concoct cut_up_fasta.py {input.contig} -c 10000 -o 0 --merge_last -b {output.bed} > {output.cont} && "
        "pixi r -e concoct {params.script} {output.bed} {input.bam} --threads {threads} > {output.cov} && [[ -s {output.cov} ]] && "
        "date) 2> {log.err} > {log.out}"


rule concoct:
    input:
        contig=rules.concoct_prepare.output.cont,
        cov=rules.concoct_prepare.output.cov,
    output:
        os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/{sample}/Concoct/concoct_clustering_merged.csv",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/concoct/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/concoct/{sample}.err.log"),
    threads: 20
    resources:
        runtime=60 * 8,
        mem_mb=get_mem_mb,
    message:
        "Running concoct to obtain bins"
    shell:
        "(date && pixi r -e concoct concoct --composition_file {input.contig} --coverage_file {input.cov} -b $(dirname {output})/concoct -t {threads} && "
        "pixi r -e concoct merge_cutup_clustering.py $(dirname {output})/concoct_clustering_gt1000.csv > {output} && "
        "date) 2> {log.err} > {log.out}"


rule comebin:
    input:
        cont=rules.filter_length.output.out,
        check=rules.check_bam_files_for_binning.output,
    output:
        os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/{sample}/Comebin/comebin_res/comebin_res.tsv",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/comebin/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/comebin/{sample}.err.log"),
    threads: 20
    resources:
        runtime=60 * 24,
        mem_mb=get_mem_mb,
    params:
        bam=os.path.join(RESULTS_DIR, "intermediate/BamFiles/{sample}"),
        folder=os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/Comebin/"),
    message:
        "Running Comebin"
    shell:
        "(date && "
        "rm -r {params.folder} && "
        "pixi r -e comebin run_comebin.sh -a {input.cont} -o {params.folder} -t {threads} -p {params.bam}"
        "&& date) 2> {log.err} > {log.out}"


rule semibin:
    input:
        cont=rules.filter_length.output.out,
        bam=expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles/{{sample}}/{{sample}}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        check=rules.check_bam_files_for_binning.output,
    output:
        os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/Semibin/contig_bins.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/semibin/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/semibin/{sample}.err.log"),
    threads: 20
    resources:
        runtime=60 * 48,
        mem_mb=get_mem_mb,
    message:
        "Running Semibin"
    shell:
        """
        (date
        rm -rf $(dirname {output})
        pixi r -e semibin SemiBin2 single_easy_bin -i {input.cont} -b {input.bam} -o $(dirname {output}) -t {threads} --no-recluster --minfasta 200000
        touch {output}
        date) 2> {log.err} > {log.out}
        """


rule metacat_cov:
    input:
        db=rules.metacat_install.output.file,
        bam=expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles/{{sample}}/{{sample}}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        check=rules.check_bam_files_for_binning.output,
    output:
        os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/MetaCAT/{sample}_metacat.cov"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/metacat_cov/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metacat_cov/{sample}.err.log"),
    threads: 20
    resources:
        runtime=60 * 4,
        mem_mb=get_mem_mb,
    message:
        "Running MetaCAT coverage"
    shell:
        """
        (date
        rm -rf {output}
        pixi r -e metacat MetaCAT coverage -b {input.bam} -ti 20 -tc 20 -o {output}
        date) 2> {log.err} > {log.out}
        """


rule metacat_seed:
    input:
        db=rules.metacat_install.output.file,
        cont=rules.filter_length.output.out,
    output:
        os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/MetaCAT/{sample}_metacat.seed"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/metacat_seed/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metacat_seed/{sample}.err.log"),
    threads: 20
    resources:
        runtime=60 * 2,
        mem_mb=get_mem_mb,
    message:
        "Running MetaCAT seed"
    shell:
        """
        (date
        rm -rf {output}
        pixi r -e metacat MetaCAT seed -f {input.cont} -o {output} -t 20
        date) 2> {log.err} > {log.out}
        """


rule metacat_cluster:
    input:
        db=rules.metacat_install.output.file,
        cont=rules.filter_length.output.out,
        seed=rules.metacat_seed.output,
        cov=rules.metacat_cov.output,
    output:
        folder=directory(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/MetaCAT/{sample}_metacat/Bins/",
            )
        ),
        done=os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/{sample}/MetaCAT/{sample}_metacat.done",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/metacat_cluster/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metacat_cluster/{sample}.err.log"),
    threads: 20
    resources:
        runtime=60 * 4,
        mem_mb=get_mem_mb,
    params:
        name="{sample}",
    message:
        "Running MetaCAT"
    shell:
        """
        (date
        rm -rf {output.folder}
        mkdir -p {output.folder}
        pixi r -e metacat MetaCAT cluster -t {threads} -f {input.cont} -s {input.seed} -c {input.cov} -o {output.folder}/{params.name}_metacat
        touch {output.done}
        date) 2> {log.err} > {log.out}
        """
