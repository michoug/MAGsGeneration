rule annotations:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/Annotation/files.txt")),
    output:
        touch("status/annotations.done"),


checkpoint countGenome:
    input:
        os.path.join(RESULTS_DIR, "Bins/finalBins"),
    output:
        temp(directory(os.path.join(RESULTS_DIR, "Bins/intermediate/tempBins"))),
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
        fasta=os.path.join(RESULTS_DIR, "Bins/intermediate/tempBins/{mags}.fa"),
    output:
        protein=os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.faa"),
        gff=os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.gff"),
        nucl=os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.ffn"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/annotations/{mags}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/annotations/{mags}.err.log"),
    threads: 4
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    params:
        prefix="{mags}",
    message:
        "Prodigal annotation for {wildcards.mags}"
    shell:
        """
        (
            date
            pixi r -e annotations pyrodigal -a {output.protein} -j {threads} -d {output.nucl} -f gff -i {input.fasta} -o {output.gff}
            date
        ) 2>{log.err} >{log.out}
        """


rule prepare4egg:
    input:
        count_mags_prots,
    output:
        os.path.join(RESULTS_DIR, "Bins/Annotation/files.txt"),
    localrule: True
    shell:
        "ls {input} > {output}"
