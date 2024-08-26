rule assembly_all:
    input:
        expand(os.path.join(RESULTS_DIR, "Assembly/{sample_2}.fa"), sample_2 = SAMPLES_2)
    output:
        touch("status/assembly.done")


rule assembly:
    input:
        read1 = rules.getfinalreads.output.read1,
        read2 = rules.getfinalreads.output.read2,
    output:
        fasta = os.path.join(RESULTS_DIR, "Assembly/{sample_2}.fa")
    conda:
        os.path.join(ENV_DIR, "megahit.yaml")
    params:
        directory(os.path.join(RESULTS_DIR, "Assembly/{sample_2}"))
    threads:
        config["mapping"]["threads"]
    log:
        out=os.path.join(RESULTS_DIR, "logs/mapping/{sample_2}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/mapping/{sample_2}.err.log")
    shell:
        """
        (date
        rm -rf {params}
        megahit -1 {input.read1} -2 {input.read2} -o {params} --min-contig-len 1000
        cp {params}/final.contigs.fa {output.fasta}
        date ) 2> {log.err} > {log.out}
        """