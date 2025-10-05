rule contigs_taxonomy:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_tax_good.tsv"
            ),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(
                RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_tax_good.tsv"
            ),
            sample=SAMPLES,
        ),
    output:
        touch("status/taxonomy.done"),


rule taxonomy:
    input:
        fasta=rules.assembly_megahit.output.final,
        file=os.path.join(RESULTS_DIR, "DB/gtdb/database_report.tsv"),
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
    resources:
        runtime=60 * 5,
        mem_mb=get_mem_mb,
    threads: 20
    conda:
        os.path.join(WORK_DIR, "envs/metabuli.yaml")
    shell:
        """
        (date
        name=$(echo $(basename {output.tax_report}) | sed "s/_report.tsv//g")
        metabuli classify {input.fasta} $(dirname {input.file}) $(dirname {output.tax_report}) $name --threads {threads} --max-ram 130 --seq-mode 1
        date) 2> {log.err} > {log.out}
        """


rule taxonomy_euk:
    input:
        fasta=rules.assembly_megahit.output.final,
        file=rules.create_tax_db_euk.output.mapping,
    output:
        tax_report=os.path.join(
            RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_report.tsv"
        ),
        tmp=temp(directory(os.path.join(RESULTS_DIR, "Assembly/{sample}_tmp_euk"))),
        dbfasta=temp(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_fasta")
        ),
    threads: 55
    resources:
        runtime=60 * 5,
        mem_mb=get_mem_mb,
    conda:
        os.path.join(WORK_DIR, "envs/mmseqs.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_contig_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_contig_euk/{sample}.err.log"),
    params:
        "MicroEuk90",
    shell:
        """
        (date && name=$(echo $(basename {output.tax_report}) | sed "s/_tax.tsv//g")
        rm -rf $(dirname {output.dbfasta})/*
        mmseqs createdb  {input.fasta} {output.dbfasta}
        mmseqs taxonomy {output.dbfasta} $(dirname {input.file})/{params} $(dirname {output.tax_report})/$name {output.tmp} --threads {threads} --tax-lineage 1
        mmseqs createtsv {output.dbfasta} $(dirname {output.tax_report})/$name {output.tax_report}
        date) 2> {log.err} > {log.out}
        """


rule taxconverter:
    input:
        report=rules.taxonomy.output.tax_report,
        classif=rules.taxonomy.output.tax_class,
    output:
        os.path.join(
            RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_tax_temp.tsv"
        )
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_convert/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_convert/{sample}.err.log"),
    localrule: True
    conda:
        os.path.join(WORK_DIR, "envs/taxconverter.yaml")
    shell:
        "(date && "
        "taxconverter metabuli -c {input.classif} -r {input.report} -o {output} && date) 2> {log.err} > {log.out}"


rule taxconverter_euk:
    input:
        tax=rules.taxonomy_euk.output.tax_report,
    output:
            os.path.join(
                RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_tax_conv.tsv"
            )
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_convert_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_convert_euk/{sample}.err.log"),
    localrule: True
    conda:
        os.path.join(WORK_DIR, "envs/taxconverter.yaml")
    shell:
        "(date && "
        "taxconverter mmseqs2 -i {input.tax} -o {output} && date) 2> {log.err} > {log.out}"


rule tax_add_missing:
    input:
        fasta=rules.assembly_megahit.output.final,
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
        fasta=rules.assembly_megahit.output.final,
        tax=rules.taxconverter_euk.output,
    output:
        os.path.join(RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_tax_temp.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_missing_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_missing_euk/{sample}.err.log"),

rule fix_taxonomy_euk:
    input:
        tax= os.path.join(
            RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_tax.tsv"
        ),
        script=os.path.join(WORK_DIR, "scripts/fixTaxEukSimple.pl")
    output:
        os.path.join(
            RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_tax_good.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_fix_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_fix_euk/{sample}.err.log"),
    localrule: True
    shell:
        """
        (date
        perl {input.script} {input.tax} {output}
        date) 2> {log.err} > {log.out}
        """
