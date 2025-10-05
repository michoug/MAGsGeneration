rule contigs_taxonomy:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_tax.tsv"
            ),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_tax.tsv"
            ),
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

use rule taxonomy as taxonomy_euk with:
    input:
        fasta=os.path.join(ASSEMBLY_DIR, "{sample}.fna"),
        file=rules.create_tax_db_euk.output.mapping,
    output:
        taxonomy=os.path.join(
            RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_tax_raw.tsv"
        ),
        tmp=temp(directory(os.path.join(RESULTS_DIR, "Assembly/{sample}_tmp_euk"))),
        dbfasta=temp(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_fasta")
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_contig_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_contig_euk/{sample}.err.log"),
    params: "MicroEuk90"



rule taxconverter:
    input:
        tax=rules.taxonomy.output.taxonomy,
    output:
        temp(
            os.path.join(
                RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_tax_temp.tsv"
            )
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_convert/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_convert/{sample}.err.log"),
    localrule: True
    conda:
        os.path.join(WORK_DIR, "envs/taxconverter.yaml")
    shell:
        "(date && "
        "taxconverter metabuli -c {input.report} -r {input.classif} -o {output} && date) 2> {log.err} > {log.out}"

use rule taxconverter as taxconverter_euk with:
    input:
        tax=rules.taxonomy_euk.output.taxonomy,
    output:
        temp(
            os.path.join(
                RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_tax_temp.tsv"
            )
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_convert_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_convert_euk/{sample}.err.log"),


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

use rule tax_add_missing as tax_add_missing_euk with:
    input:
        fasta=os.path.join(ASSEMBLY_DIR, "{sample}.fna"),
        tax=rules.taxconverter_euk.output,
    output:
        os.path.join(RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_tax.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_missing_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_missing_euk/{sample}.err.log"),