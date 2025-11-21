
rule mags_generation:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta.r{number}.sti"),
            sample=SAMPLES,
            number=[100, 125, 150, 250, 400],
        ),
        expand(
            os.path.join(RESULTS_DIR, "Assembly/index_{sample}.done"),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(RESULTS_DIR, "BamFiles/{sample}/{sample}_{sample_2}.bam"),
            sample=SAMPLES,
            sample_2=SAMPLES,
        ),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Metabat"), sample=SAMPLES),
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/{sample}/Rosella/UMAP_projection_with_clusters_final.png",
            ),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/{sample}/Vamb/Vamb_Results/vaevae_clusters_unsplit.tsv",
            ),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/{sample}/Concoct/concoct_clustering_merged.csv"
            ),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/{sample}/Comebin/comebin_res/comebin_res.tsv"
            ),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(RESULTS_DIR, "Bins/{sample}/Semibin/contig_bins.tsv"),
            sample=SAMPLES,
        ),
    output:
        touch("status/mags_generation.done"),


rule indexAssembly:
    input:
        fa=rules.filter_length.output.out,
    output:
        os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta.r{number}.sti"),
    threads: 10
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    params:
        numb="{number}",
    log:
        out=os.path.join(
            RESULTS_DIR, "logs/indexAssembly/indexAssembly_{sample}_{number}.out.log"
        ),
        err=os.path.join(
            RESULTS_DIR, "logs/indexAssembly/indexAssembly_{sample}_{number}.err.log"
        ),
    shell:
        """
        (date 
        pixi r -e coverm strobealign --create-index -t {threads} {input.fa} -r {params.numb}
        date ) 2> {log.err} > {log.out}
        """


rule checkIndexAssembyDone:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "Assembly/{{sample}}_filter.fasta.r{number}.sti"
            ),
            number=[100, 125, 150, 250, 400],
        ),
    output:
        os.path.join(RESULTS_DIR, "Assembly/index_{sample}.done"),
    localrule: True
    shell:
        "touch {output}"


rule mapping:
    input:
        read1=rules.getdedupreads.output.read1,
        read2=rules.getdedupreads.output.read2,
        fasta=rules.filter_length.output.out,
        index=rules.checkIndexAssembyDone.output,
    output:
        os.path.join(RESULTS_DIR, "BamFiles/{sample}/{sample}_{sample_2}.bam"),
    threads: 8
    resources:
        runtime=60 * 2,
        mem_mb=get_mem_mb,
    log:
        out=os.path.join(RESULTS_DIR, "logs/mapping/{sample}_{sample_2}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/mapping/{sample}_{sample_2}.err.log"),
    shell:
        """
        (date
        rm -f {output}.tmp*
        pixi r -e coverm strobealign -t {threads} --use-index {input.fasta} {input.read1} {input.read2} | 
        pixi r -e coverm samtools sort -@ {threads} -o {output}
        sleep 30
        pixi r -e coverm samtools index {output}
        date ) 2> {log.err} > {log.out}
        """


rule metabat2_mapping:
    input:
        expand(
            os.path.join(RESULTS_DIR, "BamFiles/{{sample}}/{{sample}}_{sample_2}.bam"),
            sample_2=SAMPLES,
        ),
    output:
        os.path.join(RESULTS_DIR, "Bins/{sample}/{sample}_metabat_cov.txt"),
    threads: 10
    resources:
        runtime=60 * 2,
        mem_mb=get_mem_mb,
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabat_mapping/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabat_mapping/{sample}.err.log"),
    message:
        "Running the first part of metabat2 to estimate the coverage"
    shell:
        "(date && pixi r -e metabat jgi_summarize_bam_contig_depths --outputDepth {output} {input} && "
        "date) 2> {log.err} > {log.out}"


rule metabat2:
    priority: 20
    input:
        contig=rules.filter_length.output.out,
        cov=rules.metabat2_mapping.output,
    output:
        directory(os.path.join(RESULTS_DIR, "Bins/{sample}/Metabat/")),
    threads: 10
    resources:
        runtime=60 * 2,
        mem_mb=get_mem_mb,
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabat2/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabat2/{sample}.err.log"),
    message:
        "Running metabat2 to obtain bins"
    shell:
        "(date && pixi r -e metabat metabat2 -i {input.contig} -a {input.cov} -o {output}/metabat -t {threads} && "
        "date) 2> {log.err} > {log.out}"


rule metabinner_coverage:
    input:
        expand(
            os.path.join(RESULTS_DIR, "BamFiles/{{sample}}/{{sample}}_{sample_2}.bam"),
            sample_2=SAMPLES,
        ),
    output:
        temp=temp(
            os.path.join(RESULTS_DIR, "Bins/{sample}/{sample}_metabinner_temp.txt")
        ),
        final=os.path.join(RESULTS_DIR, "Bins/{sample}/{sample}_metabinner_cov.txt"),
    threads: 14
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    params:
        length=1500,
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabinner_mapping/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabinner_mapping/{sample}.err.log"),
    message:
        "Running the first part of metabinner to estimate the coverage"
    shell:
        """
        (date && pixi r -e metabat jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth {output.temp} {input} && 
        cat {output.temp} | awk '{{if ($2>{params.length}) print $0 }}' | cut -f -1,4- > {output.final} && 
        date) 2> {log.err} > {log.out}
        """


rule rosella:
    input:
        cont=rules.filter_length.output.out,
        cov=expand(
            os.path.join(RESULTS_DIR, "BamFiles/{{sample}}/{{sample}}_{sample_2}.bam"),
            sample_2=SAMPLES,
        ),
    output:
        os.path.join(
            RESULTS_DIR,
            "Bins/{sample}/Rosella/UMAP_projection_with_clusters_final.png",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/rosella/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/rosella/{sample}.err.log"),
    threads: 14
    resources:
        runtime=60 * 16,
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
    output:
        os.path.join(RESULTS_DIR, "Bins/{sample}/Vamb/{sample}_aemb_cov.txt"),
    params:
        name="{sample}",
    localrule: True
    run:
        import pandas as pd

        df = pd.read_csv(input[0], sep="\t")
        df.rename(columns={"contigName": "contigname"}, inplace=True)
        df.to_csv(output[0], sep="\t", index=False, header=True)


rule vamb_length:
    input:
        rules.filter_length.output.out,
    output:
        temp(os.path.join(RESULTS_DIR, "Assembly/{sample}_filter_vamb.fasta")),
    threads: 4
    localrule: True
    message:
        "Remove contigs less than 1501 bp for Vamb"
    shell:
        "pixi r -e seqkit seqkit seq -j {threads} --remove-gaps -o {output} -m 1501 {input}"


rule vamb_filter:
    input:
        rules.tax_add_missing.output,
        rules.vamb_length.output,
    output:
        temp(os.path.join(RESULTS_DIR, "Assembly/{sample}_tax.tsv")),
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
            RESULTS_DIR, "Bins/{sample}/Vamb/Vamb_Results/vaevae_clusters_unsplit.tsv"
        ),
    threads: 14
    resources:
        runtime=60 * 24,
        mem_mb=get_mem_mb,
    log:
        out=os.path.join(RESULTS_DIR, "logs/vamb/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/vamb/{sample}.err.log"),
    shell:
        "(date && "
        "rm -r $(dirname {output}) && "
        "pixi r -e vamb vamb bin taxvamb --outdir $(dirname {output}) -m 1501 --fasta {input.contig} --taxonomy {input.tax} --minfasta 500000 --abundance_tsv {input.cov} -p {threads} -o && "
        "date) 2> {log.err} > {log.out}"


rule concoct_prepare:
    input:
        contig=rules.filter_length.output.out,
        bam=expand(
            os.path.join(RESULTS_DIR, "BamFiles/{{sample}}/{{sample}}_{sample_2}.bam"),
            sample_2=SAMPLES,
        ),
    output:
        bed=os.path.join(RESULTS_DIR, "Bins/{sample}/{sample}_10k.bed"),
        cont=os.path.join(RESULTS_DIR, "Bins/{sample}/{sample}_10k.fa"),
        cov=os.path.join(RESULTS_DIR, "Bins/{sample}/{sample}_concoct_cov.txt"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/concoct_prepare/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/concoct_prepare/{sample}.err.log"),
    message:
        "Preparing files for Concoct"
    resources:
        runtime=60 * 5,
    shell:
        "(date && rm -f {output.bed} {output.cont} {output.cov} && pixi r -e concoct cut_up_fasta.py {input.contig} -c 10000 -o 0 --merge_last -b {output.bed} > {output.cont} && "
        "pixi r -e concoct concoct_coverage_table.py {output.bed} {input.bam} > {output.cov} && [[ -s {output.cov} ]] && "
        "date) 2> {log.err} > {log.out}"


rule concoct:
    input:
        contig=rules.concoct_prepare.output.cont,
        cov=rules.concoct_prepare.output.cov,
    output:
        os.path.join(RESULTS_DIR, "Bins/{sample}/Concoct/concoct_clustering_merged.csv"),
    threads: 20
    resources:
        runtime=60 * 8,
        mem_mb=get_mem_mb,
    log:
        out=os.path.join(RESULTS_DIR, "logs/concoct/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/concoct/{sample}.err.log"),
    message:
        "Running concoct to obtain bins"
    shell:
        "(date && pixi r -e concoct concoct --composition_file {input.contig} --coverage_file {input.cov} -b $(dirname {output})/concoct -t {threads} && "
        "pixi r -e concoct merge_cutup_clustering.py $(dirname {output})/concoct_clustering_gt1000.csv > {output} && "
        "date) 2> {log.err} > {log.out}"


rule comebin:
    input:
        cont=rules.filter_length.output.out,
        tmp=os.path.join(RESULTS_DIR, "Bins/{sample}/{sample}_metabat_cov.txt"),  # To make sure that the bam files are there
    output:
        os.path.join(RESULTS_DIR, "Bins/{sample}/Comebin/comebin_res/comebin_res.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/comebin/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/comebin/{sample}.err.log"),
    params:
        bam=os.path.join(RESULTS_DIR, "BamFiles/{sample}"),
        folder=os.path.join(RESULTS_DIR, "Bins/{sample}/Comebin/"),
    threads: 20
    resources:
        runtime=60 * 24,
        mem_mb=get_mem_mb,
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
            os.path.join(RESULTS_DIR, "BamFiles/{{sample}}/{{sample}}_{sample_2}.bam"),
            sample_2=SAMPLES,
        ),
    output:
        os.path.join(RESULTS_DIR, "Bins/{sample}/Semibin/contig_bins.tsv"),
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
        pixi r -e semibin SemiBin2 single_easy_bin -i {input.cont} -b {input.bam} -o $(dirname {output}) -t {threads} --no-recluster
        touch {output}
        date) 2> {log.err} > {log.out}
        """
