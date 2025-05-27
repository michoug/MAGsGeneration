rule contigs_taxonomy:
    input:
        # expand(
        #     os.path.join(
        #         RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_tax_good.tsv"
        #     ),
        #     sample=ASSEMBLY,
        # ),
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
        file=rules.create_tax_db_euk.output.mapping,
    output:
        taxonomy=os.path.join(
            RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_tax.tsv"
        ),
        tmp=temp(directory(os.path.join(RESULTS_DIR, "Assembly/{sample}_tmp"))),
        dbfasta=temp(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_fasta")
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_contig/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_contig/{sample}.err.log"),
    threads: 55
    resources:
        runtime=60 * 5,
        mem_mb=get_mem_mb,
    conda:
        os.path.join(WORK_DIR, "envs/mmseqs.yaml")
    shell:
        '(date && name=$(echo $(basename {output.taxonomy}) | sed "s/_tax.tsv//g") && '
        "rm -rf $(dirname {output.dbfasta})/* && "
        "mmseqs createdb  {input.fasta} {output.dbfasta} && "
        "mmseqs taxonomy {output.dbfasta} $(dirname {input.file})/GTDB_mmseqs $(dirname {output.taxonomy})/$name {output.tmp} --threads {threads} --tax-lineage 1 && "
        "mmseqs createtsv {output.dbfasta} $(dirname {output.taxonomy})/$name {output.taxonomy} && date) 2> {log.err} > {log.out}"


rule taxonomy_euk:
    input:
        fasta=os.path.join(ASSEMBLY_DIR, "{sample}.fna"),
        file=rules.create_tax_db_euk.output.mapping,
    output:
        taxonomy=os.path.join(
            RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_tax.tsv"
        ),
        tmp=temp(directory(os.path.join(RESULTS_DIR, "Assembly/{sample}_tmp_euk"))),
        dbfasta=temp(
            os.path.join(RESULTS_DIR, "Assembly/{sample}_taxonomy_euk/{sample}_fasta")
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_contig_euk/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_contig_euk/{sample}.err.log"),
    threads: 55
    resources:
        runtime=60 * 5,
        mem_mb=get_mem_mb,
    conda:
        os.path.join(WORK_DIR, "envs/mmseqs.yaml")
    shell:
        '(date && name=$(echo $(basename {output.taxonomy}) | sed "s/_tax.tsv//g") && '
        "rm -rf $(dirname {output.dbfasta})/* && "
        "mmseqs createdb  {input.fasta} {output.dbfasta} && "
        "mmseqs taxonomy {output.dbfasta} $(dirname {input.file})/MicroEuk90 $(dirname {output.taxonomy})/$name {output.tmp} --threads {threads} --tax-lineage 1 && "
        "mmseqs createtsv {output.dbfasta} $(dirname {output.taxonomy})/$name {output.taxonomy} && date) 2> {log.err} > {log.out}"


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
    threads: 20
    conda:
        os.path.join(WORK_DIR, "envs/taxconverter.yaml")
    shell:
        "(date && "
        "taxconverter mmseqs2 -i {input.tax} -o {output} && date) 2> {log.err} > {log.out}"


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
