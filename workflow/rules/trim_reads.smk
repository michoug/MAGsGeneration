rule trim_reads:
    input:
        expand(os.path.join(RESULTS_DIR, "TrimFinal/{sample}_1.fq.gz"), sample=SAMPLES),
    output:
        touch("status/trimming.done"),


rule trimming:
    input:
        read1=os.path.join(READS_DIR, "{sample}_R1_001.fastq.gz"),
        read2=os.path.join(READS_DIR, "{sample}_R2_001.fastq.gz"),
    output:
        trim1=temp(
            os.path.join(RESULTS_DIR, "intermediate/Trim/{sample}_R1_001_val_1.fq.gz")
        ),
        trim2=temp(
            os.path.join(RESULTS_DIR, "intermediate/Trim/{sample}_R2_001_val_2.fq.gz")
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/trimming/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/trimming/{sample}.err.log"),
    priority: 20
    threads: 4
    resources:
        runtime=60 * 3,
        mem_mb=get_mem_mb,
    shell:
        """
        (
            date
            pixi r -e trimgalore trim_galore --paired -j {threads} --output_dir $(dirname {output.trim1}) {input.read1} {input.read2}
            date
        ) 2>{log.err} >{log.out}
        """


rule indexHuman:
    input:
        fa=rules.downloadHuman.output,
    output:
        os.path.join(RESULTS_DIR, "DB/hg38.fa.gz.r{number}.sti"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/indexHuman_{number}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/indexHuman_{number}.err.log"),
    threads: 10
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    params:
        # numb = lambda wildcards, output: output[0][-7:-4]
        numb="{number}",
    shell:
        """
        (
            date
            pixi r -e coverm strobealign --create-index -t {threads} {input.fa} -r {params.numb}
            date
        ) 2>{log.err} >{log.out}
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
        temp(os.path.join(RESULTS_DIR, "intermediate/MappedHuman/{sample}.bam")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/mappingHuman/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/mappingHuman/{sample}.err.log"),
    priority: 20
    threads: 10
    resources:
        runtime=60 * 5,
        mem_mb=get_mem_mb,
    shell:
        """
        (
            date
            pixi r -e coverm strobealign -t {threads} --use-index {input.fasta} {input.read1} {input.read2} \
                | pixi r -e coverm samtools view -@ {threads} -bS - >{output}
            date
        ) 2>{log.err} >{log.out}
        """


rule removeHuman:
    # A pair is human only if a mate aligns over most of its length (see
    # scripts/filter_host_reads.py); `samtools view -f 12` dropped ~50% of
    # pairs on strobealign's short soft-clipped hits. The BAM stays in
    # strobealign's read order so mates remain adjacent for samtools fastq.
    input:
        rules.mappingHuman.output,
    output:
        temp(os.path.join(RESULTS_DIR, "intermediate/TrimHuman/{sample}.noHuman.bam")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/removeHuman/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/removeHuman/{sample}.err.log"),
    priority: 20
    threads: 10
    resources:
        runtime=60 * 3,
        mem_mb=get_mem_mb,
    params:
        script=os.path.join(WORK_DIR, "scripts/filter_host_reads.py"),
        min_aligned_fraction=config.get("human_removal", {}).get(
            "min_aligned_fraction", 0.8
        ),
        max_mismatch_rate=config.get("human_removal", {}).get("max_mismatch_rate", 0.05),
    shell:
        """
        (
            date
            pixi r -e coverm samtools view -h {input} \
                | pixi r -e coverm python {params.script} --min-aligned-fraction {params.min_aligned_fraction} --max-mismatch-rate {params.max_mismatch_rate} \
                | pixi r -e coverm samtools view -@ {threads} -b -o {output} -
            date
        ) 2>{log.err} >{log.out}
        """


rule getcleanreads:
    input:
        rules.removeHuman.output,
    output:
        read1=os.path.join(RESULTS_DIR, "intermediate/TrimClean/{sample}_1.fq.gz"),
        read2=os.path.join(RESULTS_DIR, "intermediate/TrimClean/{sample}_2.fq.gz"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/getcleanreads/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/getcleanreads/{sample}.err.log"),
    priority: 20
    threads: 10
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    shell:
        """
        (
            date
            pixi r -e coverm samtools fastq -@ {threads} -1 {output.read1} -2 {output.read2} {input}
            date
        ) 2>{log.err} >{log.out}
        """


rule getdedupreads:
    input:
        read1=rules.getcleanreads.output.read1,
        read2=rules.getcleanreads.output.read2,
    output:
        read1=os.path.join(RESULTS_DIR, "TrimFinal/{sample}_1.fq.gz"),
        read2=os.path.join(RESULTS_DIR, "TrimFinal/{sample}_2.fq.gz"),
        report_html=os.path.join(RESULTS_DIR, "TrimFinal/Fastp/{sample}_fastp.html"),
        report_json=os.path.join(RESULTS_DIR, "TrimFinal/Fastp/{sample}_fastp.json"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/getdedupreads/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/getdedupreads/{sample}.err.log"),
    priority: 20
    threads: 4
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    shell:
        """
        (
            date
            pixi r -e fastp fastp -i {input.read1} -I {input.read2} -h {output.report_html} -j {output.report_json} -o {output.read1} -O {output.read2} --dedup --thread {threads}
            date
        ) 2>{log.err} >{log.out}
        """
