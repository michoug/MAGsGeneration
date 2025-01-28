rule mgthermometer:
    input:
        expand(os.path.join(RESULTS_DIR, "mgthermometer/table_OGT.txt")),
    output:
        touch("status/mgthermometer.done"),


def count_mags(wildcards):
    ck_output = checkpoints.countGenome.get(**wildcards).output[0]
    return expand(
        os.path.join(RESULTS_DIR, "mgthermometer/proportion_{mags}.txt"),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )


rule mgthermometer_run:
    input:
        fasta=os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.faa"),
        script=os.path.join(WORK_DIR, "scripts/getFrequency.pl"),
    output:
        proportion=os.path.join(RESULTS_DIR, "mgthermometer/proportion_{mags}.txt"),
    params:
        prefix="{mags}",
    conda:
        os.path.join(WORK_DIR, "envs/mgthermometer.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/mgthermometer/{mags}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/mgthermometer/{mags}.err.log"),
    message:
        "Calculate proportion of IVYWREL for {wildcards.mags}"
    localrule: True
    shell:
        "(date && "
        "perl {input.script} {input.fasta} {output.proportion} && "
        "date) 2> {log.err} > {log.out}"


rule concatResults:
    input:
        count_mags,
    output:
        os.path.join(RESULTS_DIR, "mgthermometer/table_OGT.txt"),
    localrule: True
    shell:
        "(date && cat {input} > {output} && date)"
