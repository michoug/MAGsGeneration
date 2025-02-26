rule annotations:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/Annotation/files.txt")),
    output:
        touch("status/annotations.done"),


checkpoint countGenome:
    input:
        os.path.join(RESULTS_DIR, "Bins/finalBins"),
    output:
        temp(directory(os.path.join(RESULTS_DIR, "Bins/tempBins"))),
    shell:
        "mkdir {output} && cp {input}/*fa {output}/"


def count_mags_prots(wildcards):
    ck_output = checkpoints.countGenome.get(**wildcards).output[0]
    return expand(
        os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.faa"),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )


rule annotation:
    input:
        fasta=os.path.join(RESULTS_DIR, "Bins/tempBins/{mags}.fa"),
    output:
        protein=os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.faa"),
        gff=os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.gff"),
        nucl=os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.ffn"),
    params:
        prefix="{mags}",
    conda:
        os.path.join(WORK_DIR, "envs/annotations.yaml")
    resources:
        runtime=60,
        mem_mb=lambda wildcards, threads: threads * config["memory_per_thread"],
    threads: 4
    log:
        out=os.path.join(RESULTS_DIR, "logs/annotations/{mags}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/annotations/{mags}.err.log"),
    message:
        "Prodigal annotation for {wildcards.mags}"
    shell:
        "(date && "
        "pyrodigal -a {output.protein} -j {threads} -d {output.nucl} -f gff -i {input.fasta} -o {output.gff} && "
        "date) 2> {log.err} > {log.out}"


rule prepare4egg:
    input:
        count_mags_prots,
    output:
        os.path.join(RESULTS_DIR, "Bins/Annotation/files.txt"),
    localrule: True
    shell:
        "ls {input} > {output}"
