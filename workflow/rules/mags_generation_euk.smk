
rule mags_generation_euk:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta.r{number}.sti"),
            sample=ASSEMBLY,
            number=[100, 125, 150, 250, 400],
        ),
        expand(
            os.path.join(RESULTS_DIR, "BamFiles/{sample}/{sample}_{sample_2}.bam"),
            sample=ASSEMBLY,
            sample_2=SAMPLES,
        ),
        expand(
            os.path.join(RESULTS_DIR, "Bins/{sample}/Metabat/metabat.done"),
            sample=ASSEMBLY,
        ),
        # expand(
        #     os.path.join(
        #         RESULTS_DIR,
        #         "Bins/{sample}/Vamb/Vamb_Results/vaevae_clusters_unsplit.tsv",
        #     ),
        #     sample=ASSEMBLY,
        # ),
    output:
        touch("status/mags_generation.done"),


rule indexAssembly:
    input:
        fa=rules.filter_length.output.out,
    output:
        os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta.r{number}.sti"),
    conda:
        os.path.join(WORK_DIR, "envs/coverm.yaml")
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
        strobealign --create-index -t {threads} {input.fa} -r {params.numb}
        date ) 2> {log.err} > {log.out}
        """


rule checkIndexAssembyDone:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta.r{number}.sti"),
            sample=ASSEMBLY,
            number=[100, 125, 150, 250, 400],
        ),
    output:
        os.path.join(RESULTS_DIR, "Assembly/index.done"),
    localrule: True
    shell:
        "touch {output}"


rule mapping:
    input:
        read1=os.path.join(RESULTS_DIR, "TrimFinal/{sample_2}_1.fq.gz"),
        read2=os.path.join(RESULTS_DIR, "TrimFinal/{sample_2}_2.fq.gz"),
        fasta=rules.filter_length.output.out,
        index=rules.checkIndexAssembyDone.output,
    output:
        os.path.join(RESULTS_DIR, "BamFiles/{sample}/{sample}_{sample_2}.bam"),
    conda:
        os.path.join(WORK_DIR, "envs/coverm.yaml")
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
        strobealign -t {threads} --use-index {input.fasta} {input.read1} {input.read2} | 
        samtools sort -@ {threads} -o {output}
        samtools index {output}
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
    conda:
        os.path.join(WORK_DIR, "envs/metabat.yaml")
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
        "(date && jgi_summarize_bam_contig_depths --outputDepth {output} {input} && "
        "date) 2> {log.err} > {log.out}"


rule metabat2:
    priority: 20
    input:
        contig=rules.filter_length.output.out,
        cov=rules.metabat2_mapping.output,
    output:
        folder=directory(os.path.join(RESULTS_DIR, "Bins/{sample}/Metabat/")),
        file=os.path.join(RESULTS_DIR, "Bins/{sample}/Metabat/metabat.done"),
    conda:
        os.path.join(WORK_DIR, "envs/metabat.yaml")
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
        """
        (date
        metabat2 -i {input.contig} -a {input.cov} -o {output.folder}/metabat -t {threads}
        touch {output.file}
        date) 2> {log.err} > {log.out}
        """


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
    conda:
        os.path.join(WORK_DIR, "envs/metabat.yaml")
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
        (date && jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth {output.temp} {input} && 
        cat {output.temp} | awk '{{if ($2>{params.length}) print $0 }}' | cut -f -1,4- > {output.final} && 
        date) 2> {log.err} > {log.out}
        """


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
    conda:
        os.path.join(WORK_DIR, "envs/seqkit.yaml")
    threads: 4
    localrule: True
    message:
        "Remove contigs less than 1501 bp for Vamb"
    shell:
        "seqkit seq -j {threads} --remove-gaps -o {output} -m 1501 {input}"


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
        file=rules.vamb_install.output,
        contig=rules.vamb_length.output,
        tax=rules.vamb_filter.output,
        cov=rules.vamb_coverage.output,
    output:
        os.path.join(
            RESULTS_DIR, "Bins/{sample}/Vamb/Vamb_Results/vaevae_clusters_unsplit.tsv"
        ),
    conda:
        os.path.join(WORK_DIR, "envs/vamb.yaml")
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
        "vamb bin taxvamb --outdir $(dirname {output}) -m 1501 --fasta {input.contig} --taxonomy {input.tax} --minfasta 500000 --abundance_tsv {input.cov} -p {threads} -o && "
        "date) 2> {log.err} > {log.out}"
