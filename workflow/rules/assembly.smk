rule assembly_all:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Assembly/intermediate/{sample}_megahit.fa"),
            sample=ASSEMBLY,
        ),
    output:
        touch("status/assembly.done"),


rule assembly_megahit:
    input:
        read1=rules.getdedupreads.output.read1,
        read2=rules.getdedupreads.output.read2,
    output:
        fasta=os.path.join(
            RESULTS_DIR, "Assembly/intermediate/{sample}_megahit/final.contigs.fa"
        ),
        final=os.path.join(RESULTS_DIR, "Assembly/intermediate/{sample}_megahit.fa"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/assembly/{sample}_megahit.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/assembly/{sample}_megahit.err.log"),
    threads: 20
    resources:
        runtime=60 * 12,
        mem_mb=get_mem_mb,
    params:
        folder=directory(
            os.path.join(RESULTS_DIR, "Assembly/intermediate/{sample}_megahit")
        ),
        name="{sample}",
    shell:
        """
        (
            date
            rm -rf {params}
            pixi r -e megahit megahit -1 {input.read1} -2 {input.read2} -o {params.folder} --min-contig-len 1000 -t {threads}
            perl -pe "s/>/>{params.name}_/" {output.fasta} >{output.final}
            date
        ) 2>{log.err} >{log.out}
        """
