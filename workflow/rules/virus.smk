rule virus_all:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Virus/{sample}_genomad/summary.log"),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(
                RESULTS_DIR, "Virus/{sample}_checkamg/results/final_results.tsv"
            ),
            sample=SAMPLES,
        ),
    output:
        touch("status/virus.done"),


rule genomad:
    input:
        db=rules.genomad_db.output.done,
        fasta=rules.filter_length.output.out,
    output:
        folder=directory(os.path.join(RESULTS_DIR, "Virus/{sample}_genomad")),
        fna=os.path.join(
            RESULTS_DIR, "Virus/{sample}_genomad/{sample}_summary/{sample}_virus.fna"
        ),
    resources:
        mem_mb=get_mem_mb,
        runtime=4 * 60,
    threads: 10
    log:
        out=os.path.join(RESULTS_DIR, "logs/genomad_{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/genomad_{sample}.err.log"),
    shell:
        """
        (date
        pixi r -e virus genomad end-to-end --cleanup --composition metagenome --threads {threads} {input.fasta} {output.folder} {input.db}
        date) 2> {log.err} > {log.out}
        """


rule check_amg:
    input:
        db=rules.checkamg_db.output.done,
        fna=rules.genomad.output.fna,
    output:
        folder=directory(os.path.join(RESULTS_DIR, "Virus/{sample}_checkamg")),
        done=os.path.join(
            RESULTS_DIR, "Virus/{sample}_checkamg/results/final_results.tsv"
        ),
    resources:
        mem_mb=get_mem_mb,
        runtime=4 * 60,
    threads: 10
    log:
        out=os.path.join(RESULTS_DIR, "logs/checkamg_{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/checkamg_{sample}.err.log"),
    shell:
        """
        (date
        pixi r -e virus checkamg annotate -i {input.fna} -d {input.db} -o {output.folder} -t {threads}
        date) 2> {log.err} > {log.out}
        """
