# Running CheckM on the bins


rule checkm_all:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/checkm2/quality_report.tsv"
            ),
            sample=ASSEMBLY,
        ),
    output:
        touch("status/checkm2.done"),


# Checking bin quality
# Binette already runs CheckM2 internally to score and select its final
# bins (see rule binette_run), so this reformats its own quality report into
# CheckM2's native schema instead of re-running checkm2 predict from scratch
# for no new information.
rule checkm2:
    input:
        tsv=rules.binette_run.output.quality,
    output:
        tsv=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/checkm2/quality_report.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/checkm/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/checkm/{sample}.err.log"),
    localrule: True
    params:
        script=os.path.join(WORK_DIR, "scripts/binette_quality_to_checkm2.py"),
    message:
        "Reformatting Binette's quality report for {wildcards.sample} into CheckM2's native schema"
    shell:
        """
        (
            date
            pixi r -e binette python {params.script} --binette-tsv {input.tsv} --out {output.tsv}
            date
        ) 2>{log.err} >{log.out}
        """
