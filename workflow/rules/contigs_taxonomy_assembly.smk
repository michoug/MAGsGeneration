rule contigs_taxonomy:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_tax.tsv"),
            sample=ASSEMBLY,
        ),
    output:
        touch("status/taxonomy.done"),


rule taxonomy:
    input:
        fasta=os.path.join(ASSEMBLY_DIR, "{sample}.fna"),
        file=rules.download_tax_db.output.file,
    output:
        tax_report=os.path.join(
            RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_report.tsv"
        ),
        tax_class=os.path.join(
            RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_classifications.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_contig/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_contig/{sample}.err.log"),
    threads: 20
    resources:
        runtime=60 * 5,
        mem_mb=get_mem_mb,
    conda:
        os.path.join(WORK_DIR, "envs/metabuli.yaml")
    shell:
        """
        (date
        name=$(echo $(basename {output.tax_report}) | sed "s/_report.tsv//g")
        metabuli classify {input.fasta} $(dirname {input.file}) $(dirname {output.tax_report}) $name --threads {threads} --max-ram 130 --seq-mode 1
        date) 2> {log.err} > {log.out}
        """


rule taxconverter:
    input:
        report=rules.taxonomy.output.tax_report,
        classif=rules.taxonomy.output.tax_class,
    output:
        os.path.join(RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_tax_temp.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_convert/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_convert/{sample}.err.log"),
    localrule: True
    conda:
        os.path.join(WORK_DIR, "envs/taxconverter.yaml")
    shell:
        "(date && "
        "taxconverter metabuli -c {input.classif} -r {input.report} -o {output} && date) 2> {log.err} > {log.out}"


rule tax_add_missing:
    input:
        fasta=os.path.join(ASSEMBLY_DIR, "{sample}.fna"),
        tax=rules.taxconverter.output,
    output:
        os.path.join(RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_tax_good.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_missing/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_missing/{sample}.err.log"),
    localrule: True
    run:
        import pandas as pd

        contigs = []
        with open(input[0], "r") as f:
            for line in f:
                if line.startswith(">"):
                    contigs.append(line[1:].strip().split()[0])

        tax_df = pd.read_csv(
            input[1], sep="\t", header=None, names=["contigs", "predictions"]
        )

        all_df = pd.DataFrame({"contigs": contigs})

        merged = all_df.merge(tax_df, on="contigs", how="left")
        merged["predictions"].fillna("unknown", inplace=True)

        merged.to_csv(output[0], sep="\t", index=False, header=True)
