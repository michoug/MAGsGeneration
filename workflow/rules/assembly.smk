rule assembly_all:
    input:
        expand(os.path.join(RESULTS_DIR, "Assembly/{sample}_megahit.fa"), sample = SAMPLES)
    output:
        touch("status/assembly.done")


rule assembly_megahit:
    input:
        read1 = rules.getfinalreads.output.read1,
        read2 = rules.getfinalreads.output.read2,
    output:
        fasta = os.path.join(RESULTS_DIR, "Assembly/{sample}_megahit/final.contigs.fa"),
        final = os.path.join(RESULTS_DIR, "Assembly/{sample}_megahit.fa")
    conda:
        os.path.join(ENV_DIR, "megahit.yaml")
    params:
        folder = directory(os.path.join(RESULTS_DIR, "Assembly/{sample}_megahit")), 
        name = "{sample}"
    threads:
        config["assembly"]["threads"]
    log:
        out=os.path.join(RESULTS_DIR, "logs/assembly/{sample}_megahit.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/assembly/{sample}_megahit.err.log")
    shell:
        """
        (date
        rm -rf {params}
        megahit -1 {input.read1} -2 {input.read2} -o {params.folder} --min-contig-len 1000 -t {threads}
        perl -pe "s/>/>{params.name}_/" {output.fasta} > {output.final}
        date ) 2> {log.err} > {log.out}
        """
