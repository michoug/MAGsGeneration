rule chekAssembly:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"),
            sample=ASSEMBLY,
        ),
        os.path.join(RESULTS_DIR, "qc/assembly/multiqc_report.html"),
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


rule assembly_stats:
    # seqkit stats -a -T (MultiQC seqkit module) on raw and >=1.5 kb assemblies
    input:
        raw=expand(rules.assembly_megahit.output.final, sample=ASSEMBLY),
        filtered=expand(rules.filter_length.output.out, sample=ASSEMBLY),
    output:
        os.path.join(RESULTS_DIR, "qc/assembly/seqkit_stats.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/qc/assembly_stats.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/qc/assembly_stats.err.log"),
    threads: 4
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    message:
        "Computing assembly statistics with seqkit"
    shell:
        """
        (
            date
            pixi r -e seqkit seqkit stats -a -T -b -j {threads} {input.raw} {input.filtered} >{output}
            date
        ) 2>{log.err} >{log.out}
        """


rule multiqc_assembly:
    # MEGAHIT logs and seqkit statistics of the raw and filtered assemblies
    input:
        megahit=expand(rules.assembly_megahit.output.fasta, sample=ASSEMBLY),
        stats=rules.assembly_stats.output,
    output:
        os.path.join(RESULTS_DIR, "qc/assembly/multiqc_report.html"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/qc/multiqc_assembly.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/qc/multiqc_assembly.err.log"),
    threads: 1
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    params:
        logs=expand(
            os.path.join(RESULTS_DIR, "Assembly/intermediate/{sample}_megahit/log"),
            sample=ASSEMBLY,
        ),
    message:
        "Running MultiQC on the assembly reports"
    shell:
        """
        (
            date
            pixi r -e multiqc multiqc {params.logs} {input.stats} \
                -m megahit -m seqkit \
                -o $(dirname {output}) -n $(basename {output}) --force
            date
        ) 2>{log.err} >{log.out}
        """
