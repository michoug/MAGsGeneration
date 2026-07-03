rule mgthermometer:
    input:
        expand(os.path.join(RESULTS_DIR, "MGthermometer/table_OGT.txt")),
    output:
        touch("status/mgthermometer.done"),


def count_mags(wildcards):
    ck_output = checkpoints.countGenome.get(**wildcards).output[0]
    return expand(
        os.path.join(RESULTS_DIR, "MGthermometer/intermediate/proportion_{mags}.txt"),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )


rule MGthermometer_run:
    input:
        fasta=os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.faa"),
        script=os.path.join(WORK_DIR, "scripts/getFrequency.pl"),
    output:
        proportion=os.path.join(
            RESULTS_DIR, "MGthermometer/intermediate/proportion_{mags}.txt"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/MGthermometer/{mags}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/MGthermometer/{mags}.err.log"),
    localrule: True
    params:
        prefix="{mags}",
    message:
        "Calculate proportion of IVYWREL for {wildcards.mags}"
    shell:
        """
        (
            date
            pixi r -e mgthermometer perl {input.script} {input.fasta} {output.proportion}
            date
        ) 2>{log.err} >{log.out}
        """


rule concatResults:
    input:
        count_mags,
    output:
        os.path.join(RESULTS_DIR, "MGthermometer/table_OGT.txt"),
    localrule: True
    shell:
        "(date && cat {input} > {output} && date)"
