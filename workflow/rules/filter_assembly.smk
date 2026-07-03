rule chekAssembly:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"),
            sample=ASSEMBLY,
        ),
    output:
        touch("status/checkAssembly.done"),


rule filter_length:
    input:
        rules.assembly_megahit.output.final,
    output:
        out=os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"),
        other=temp(
            os.path.join(
                RESULTS_DIR, "Assembly/intermediate/{sample}_filter_name.fasta"
            )
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/filter_length/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/filter_length/{sample}.err.log"),
    threads: 4
    resources:
        mem_mb=get_mem_mb,
    message:
        "Remove contigs less than 1.5kb"
    shell:
        """
        (
            date
            pixi r -e seqkit seqkit seq -j {threads} --remove-gaps -o {output.other} -m 1500 {input}
            perl -pe "s/ .*//g" {output.other} >{output.out}
            date
        ) 2>{log.err} >{log.out}
        """
