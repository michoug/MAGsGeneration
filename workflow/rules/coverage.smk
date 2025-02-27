rule coverage_all:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Bins/cat_mags_cov_{sample}.txt"),
            sample=ASSEMBLY,
        ),
        os.path.join(RESULTS_DIR, "cat_mags_cov.txt"),
    output:
        touch("status/coverage.done"),


# CONCATENATING MAGS
rule concatenate_mags:
    input:
        rules.galah.output.final,
    output:
        os.path.join(RESULTS_DIR, "concat/cat_mags.fa"),
    message:
        "Concatenating all dereplicated MAGS"
    localrule: True
    shell:
        "(date && cat {input}/*.fa > {output} && date)"


rule mapping_all:
    input:
        read1=rules.getfinalreads.output.read1,
        read2=rules.getfinalreads.output.read2,
        cont=rules.concatenate_mags.output,
    output:
        temp(os.path.join(RESULTS_DIR, "Bam/cat_mags.fa.{sample}_1.fq.gz.bam")),
    threads: 10
    resources:
        runtime=60 * 12,
        mem_mb=get_mem_mb,
    conda:
        os.path.join(WORK_DIR, "envs/coverm.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/mapping/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/mapping/{sample}.err.log"),
    message:
        "Running coverm to produce sorted bams for: {wildcards.sample}"
    shell:
        "(date && "
        "rm -f {output} && "
        "TMPDIR=. coverm make -r {input.cont} -t {threads} -o $(dirname {output}) --discard-unmapped -c {input.read1} {input.read2} && "
        "date) 2> {log.err} > {log.out}"


# COVERAGE
rule coverage:
    input:
        os.path.join(RESULTS_DIR, "Bam/cat_mags.fa.{sample}_1.fq.gz.bam"),
    output:
        final=os.path.join(RESULTS_DIR, "Bins/cat_mags_cov_{sample}.txt"),
    conda:
        os.path.join(WORK_DIR, "envs/coverm.yaml")
    threads: 10
    resources:
        runtime=60 * 2,
        mem_mb=get_mem_mb,
    priority: 50
    log:
        out=os.path.join(RESULTS_DIR, "logs/coverage/coverage_{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/coverage/coverage_{sample}.err.log"),
    message:
        "Estimating coverage for concatenated MAGs"
    shell:
        "(date && coverm contig -b {input} -m trimmed_mean -t {threads} -o {output.final} && date)"
        " 2> {log.err} > {log.out}"


rule concat_coverage:
    input:
        files=expand(
            os.path.join(RESULTS_DIR, "Bins/cat_mags_cov_{sample}.txt"), sample=SAMPLES
        ),
        script=os.path.join(WORK_DIR, "scripts/mergeTables.pl"),
    output:
        os.path.join(RESULTS_DIR, "cat_mags_cov.txt"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/concat_coverage.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/concat_coverage.err.log"),
    localrule: True
    message:
        "Concatenating the coverage files"
    shell:
        "(date && perl {input.script} {input.files} {output} && date)"
        " 2> {log.err} > {log.out}"
