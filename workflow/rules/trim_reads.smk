rule trim_reads:
    input:
        # expand(
        #     os.path.join(RESULTS_DIR, "DB/hg38.fa.gz.r{number}.sti"),
        #     number=[100, 125, 150, 200, 250, 300],
        # ),
        # expand(
        #     os.path.join(RESULTS_DIR, "Trim/{sample}_R1_001_val_1.fq.gz"),
        #     sample=SAMPLES,
        # ),
        expand(os.path.join(RESULTS_DIR, "TrimDedup/{sample}_1.fq.gz"), sample=SAMPLES),
    output:
        touch("status/trimming.done"),


rule trimming:
    input:
        read1=os.path.join(READS_DIR, "{sample}_R1_001.fastq.gz"),
        read2=os.path.join(READS_DIR, "{sample}_R2_001.fastq.gz"),
    output:
        trim1=temp(os.path.join(RESULTS_DIR, "Trim/{sample}_R1_001_val_1.fq.gz")),
        trim2=temp(os.path.join(RESULTS_DIR, "Trim/{sample}_R2_001_val_2.fq.gz")),
    threads: 4
    resources:
        runtime=60 * 3,
        mem_mb=get_mem_mb,
    priority: 20
    log:
        out=os.path.join(RESULTS_DIR, "logs/trimming/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/trimming/{sample}.err.log"),
    shell:
        "(date && "
        "pixi r -e trimgalore trim_galore --paired -j {threads} --output_dir $(dirname {output.trim1}) {input.read1} {input.read2} && "
        "date) 2> {log.err} > {log.out}"


rule indexHuman:
    input:
        fa=rules.downloadHuman.output,
    output:
        os.path.join(RESULTS_DIR, "DB/hg38.fa.gz.r{number}.sti"),
    threads: 10
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    params:
        # numb = lambda wildcards, output: output[0][-7:-4]
        numb="{number}",
    log:
        out=os.path.join(RESULTS_DIR, "logs/indexHuman_{number}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/indexHuman_{number}.err.log"),
    shell:
        """
        (date
        pixi r -e coverm strobealign --create-index -t {threads} {input.fa} -r {params.numb}
        date ) 2> {log.err} > {log.out}
        """


rule checkIndexDone:
    input:
        expand(
            os.path.join(RESULTS_DIR, "DB/hg38.fa.gz.r{number}.sti"),
            number=[100, 125, 150, 250, 400],
        ),
    output:
        touch("status/index.done"),
    localrule: True


rule mappingHuman:
    input:
        read1=rules.trimming.output.trim1,
        read2=rules.trimming.output.trim2,
        fasta=rules.downloadHuman.output,
        index=rules.checkIndexDone.output,
    output:
        temp(os.path.join(RESULTS_DIR, "TrimHuman/{sample}.bam")),
    threads: 10
    priority: 20
    resources:
        runtime=60 * 5,
        mem_mb=get_mem_mb,
    log:
        out=os.path.join(RESULTS_DIR, "logs/mappingHuman/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/mappingHuman/{sample}.err.log"),
    shell:
        """
        (date
        pixi r -e coverm strobealign -t {threads} --use-index {input.fasta} {input.read1} {input.read2} | 
        pixi r -e coverm samtools view -@ {threads} -bS - > {output} &&
        date ) 2> {log.err} > {log.out}
        """


rule removeHuman:
    input:
        rules.mappingHuman.output,
    output:
        temp(os.path.join(RESULTS_DIR, "TrimHuman/{sample}.noHuman.bam")),
    threads: 10
    priority: 20
    resources:
        runtime=60 * 3,
        mem_mb=get_mem_mb,
    log:
        out=os.path.join(RESULTS_DIR, "logs/removeHuman/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/removeHuman/{sample}.err.log"),
    shell:
        """
        (date
        rm -f {output}.tmp*
        pixi r -e coverm samtools view -@ {threads} -b -f 12 -F 256 {input} | pixi r -e coverm samtools sort -@ {threads} - -o {output}
        date ) 2> {log.err} > {log.out}
        """


rule getcleanreads:
    input:
        rules.removeHuman.output,
    output:
        read1=os.path.join(RESULTS_DIR, "TrimClean/{sample}_1.fq.gz"),
        read2=os.path.join(RESULTS_DIR, "TrimClean/{sample}_2.fq.gz"),
    threads: 10
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    priority: 20
    log:
        out=os.path.join(RESULTS_DIR, "logs/getcleanreads/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/getcleanreads/{sample}.err.log"),
    shell:
        """
        (date
        pixi r -e coverm samtools fastq -@ {threads} -1 {output.read1} -2 {output.read2} {input} &&
        date ) 2> {log.err} > {log.out}
        """


rule getdedupreads:
    input:
        read1=rules.getcleanreads.output.read1,
        read2=rules.getcleanreads.output.read2,
    output:
        read1=os.path.join(RESULTS_DIR, "TrimDedup/{sample}_1.fq.gz"),
        read2=os.path.join(RESULTS_DIR, "TrimDedup/{sample}_2.fq.gz"),
        report_html=os.path.join(RESULTS_DIR, "TrimDedup/{sample}_fastp.html"),
        report_json=os.path.join(RESULTS_DIR, "TrimDedup/{sample}_fastp.json"), 
    threads: 4
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    log:
        out=os.path.join(RESULTS_DIR, "logs/getdedupreads/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/getdedupreads/{sample}.err.log"),
    shell:
        """
        (date
        pixi r -e fastp fastp -i {input.read1} -I {input.read2} -h {output.report_html} -j {output.report_json} -o {output.read1} -O {output.read2} --dedup --thread {threads} &&
        date ) 2> {log.err} > {log.out}
        """
