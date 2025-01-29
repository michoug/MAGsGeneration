rule chekAssembly:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"), sample=SAMPLES
        ),
    output:
        touch("status/checkAssembly.done"),


rule filter_length:
    input:
        rules.assembly_megahit.output.final,
    output:
        out=os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"),
        other=temp(os.path.join(RESULTS_DIR, "Assembly/{sample}_filter_name.fasta")),
    conda:
        os.path.join(WORK_DIR, "envs/seqkit.yaml")
    threads: 4
    message:
        "Remove contigs less than 1.5kb"
    shell:
        'seqkit seq -j {threads} --remove-gaps -o {output.other} -m 1500 {input} && perl -pe "s/ .*//g" {output.other} > {output.out}'
