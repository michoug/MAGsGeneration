rule chekAssembly:
    input:
        "status/check_inputs.done",
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"),
            sample=ASSEMBLY,
        ),
    output:
        touch("status/checkAssembly.done"),


rule check_inputs_exist:
    input:
        expand(os.path.join(ASSEMBLY_DIR, "{sample}.fna"), sample=ASSEMBLY),
    output:
        touch("status/check_inputs.done"),


rule filter_length:
    input:
        check = rules.check_inputs_exist.output,
        fasta = os.path.join(ASSEMBLY_DIR, "{sample}.fna"),
    output:
        out=os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"),
        other=temp(os.path.join(RESULTS_DIR, "Assembly/{sample}_filter_name.fasta")),
    threads: 4
    message:
        "Remove contigs less than 1.5kb"
    shell:
        'pixi r -e seqkit seqkit seq -j {threads} --remove-gaps -o {output.other} -m 1500 {input.fasta} && perl -pe "s/ .*//g" {output.other} > {output.out}'
