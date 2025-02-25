# Running CheckM on the bins


rule checkm_all:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Bins/{sample}/checkm2/quality_report.tsv"),
            sample=SAMPLES,
        ),
    output:
        touch("status/checkm2.done"),


# Checking bin quality
rule checkm2:
    input:
        rules.dasTool_run.output,
    output:
        tsv=os.path.join(RESULTS_DIR, "Bins/{sample}/checkm2/quality_report.tsv"),
    conda:
        os.path.join(WORK_DIR, "envs/checkm2.yaml")
    resources:
        runtime=60 * 4,
        mem=lambda wildcards, threads: threads * config["memory_per_thread"],
    log:
        out=os.path.join(RESULTS_DIR, "logs/checkm/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/checkm/{sample}.err.log"),
    threads: 10
    params:
        ext="fa",
        db=config["checkm2"]["db"],
    message:
        "Running CheckM2 for directory renamed bins"
    shell:
        "(date && touch {input} && "
        "export CHECKM2DB={params.db} && "
        "checkm2 predict --threads {threads} -x {params.ext} --input $(dirname {input})/das_DASTool_bins --output-directory $(dirname {output.tsv}) --force && "
        "date) 2> {log.err} > {log.out}"
