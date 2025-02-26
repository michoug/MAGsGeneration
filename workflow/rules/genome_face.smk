rule genome_face:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/GenomeFace/concat_assembly_cov.txt")),
    output:
        touch("status/genome_face.done"),


rule renameAssembly:
    input:
        fa=rules.filter_length.output.out,
    output:
        temp(os.path.join(RESULTS_DIR, "Assembly/{sample}_rename.fa")),
    localrule: True
    log:
        out=os.path.join(RESULTS_DIR, "logs/renameAssembly_{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/renameAssembly_{sample}.err.log"),
    shell:
        """
        (date 
        perl -pe "s/_k141/_C_k141/g" {input.fa} > {output}
        date ) 2> {log.err} > {log.out}
        """


rule concat_assembly:
    input:
        expand(os.path.join(RESULTS_DIR, "Assembly/{sample}_rename.fa"), sample=SAMPLES),
    output:
        os.path.join(RESULTS_DIR, "Bins/GenomeFace/concat_assembly.fa"),
    message:
        "Concatenating all renamed assemblies"
    localrule: True
    shell:
        "cat {input} > {output}"


rule indexAssemblyGenomeFace:
    input:
        fa=rules.concat_assembly.output,
    output:
        os.path.join(RESULTS_DIR, "Bins/GenomeFace/concat_assembly.fa.r{number}.sti"),
    conda:
        os.path.join(WORK_DIR, "envs/coverm.yaml")
    threads: 10
    resources:
        runtime=60,
        mem_mb=lambda wildcards, threads: threads * config["memory_per_thread"],
    params:
        numb="{number}",
    log:
        out=os.path.join(RESULTS_DIR, "logs/indexAssemblyGenomeFace_{number}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/indexAssemblyGenomeFace_{number}.err.log"),
    shell:
        """
        (date 
        strobealign --create-index -t {threads} {input.fa} -r {params.numb}
        date ) 2> {log.err} > {log.out}
        """


rule checkindexAssemblyGenomeFaceDone:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/GenomeFace/concat_assembly.fa.r{number}.sti"
            ),
            sample=ASSEMBLY,
            number=[125, 150],
        ),
    output:
        os.path.join(RESULTS_DIR, "Bins/GenomeFace/index.done"),
    localrule: True
    shell:
        "touch {output}"


rule coverage_genome_face_bam:
    input:
        read1=rules.getfinalreads.output.read1,
        read2=rules.getfinalreads.output.read2,
        cont=os.path.join(RESULTS_DIR, "Bins/GenomeFace/concat_assembly.fa"),
        index=os.path.join(RESULTS_DIR, "Bins/GenomeFace/index.done"),
    output:
        temp(
            os.path.join(
                RESULTS_DIR, "Bins/GenomeFace/Bam/concat_assembly.fa.{sample}.bam"
            )
        ),
    conda:
        os.path.join(WORK_DIR, "envs/coverm.yaml")
    threads: 10
    resources:
        runtime=60 * 2,
        mem_mb=lambda wildcards, threads: threads * config["memory_per_thread"],
    log:
        out=os.path.join(RESULTS_DIR, "logs/coverage_genome_face_{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/coverage_genome_face_{sample}.err.log"),
    message:
        "Running coverm to produce sorted bams for GenomeFace"
    shell:
        """
        (date
        rm -f {output}.tmp*
        strobealign -t {threads} --use-index {input.cont} {input.read1} {input.read2} | 
        samtools sort -@ {threads} -o {output}
        date) 2> {log.err} > {log.out}
        """


rule coverage_genome_face:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/GenomeFace/Bam/concat_assembly.fa.{sample}.bam"
            ),
            sample=SAMPLES,
        ),
    output:
        os.path.join(RESULTS_DIR, "Bins/GenomeFace/concat_assembly_cov.txt"),
    conda:
        os.path.join(WORK_DIR, "envs/coverm.yaml")
    threads: 10
    resources:
        runtime=60 * 2,
        mem_mb=lambda wildcards, threads: threads * config["memory_per_thread"],
    log:
        out=os.path.join(RESULTS_DIR, "logs/coverage_genome_face/coverage.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/coverage_genome_face/coverage.err.log"),
    message:
        "Calculating coverage for GenomeFace"
    shell:
        "(date && "
        "coverm contig -b {input} -o {output} -t {threads} -m metabat && "
        "date) 2> {log.err} > {log.out}"
