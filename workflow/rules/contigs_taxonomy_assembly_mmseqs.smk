rule contigs_taxonomy:
    input:
        expand(
            os.path.join(
                RESULTS_DIR,
                "Assembly/intermediate/{sample}_taxonomy/{sample}_tax_good.tsv",
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
    threads: 60
    resources:
        runtime=60 * 5,
        mem_mb=get_mem_mb,
    shell:
        """
        (
            date
            name=$(echo $(basename {output.taxonomy}) | sed 's/_tax.tsv//g')
            rm -rf $(dirname {output.dbfasta})/*
            pixi r -e mmseqs mmseqs createdb {input.fasta} {output.dbfasta}
            pixi r -e mmseqs mmseqs taxonomy {output.dbfasta} $(dirname {input.file})/GTDB_mmseqs $(dirname {output.taxonomy})/$name {output.tmp} --threads {threads} --tax-lineage 1
            pixi r -e mmseqs mmseqs createtsv {output.dbfasta} $(dirname {output.taxonomy})/$name {output.taxonomy} && date
        ) 2>{log.err} >{log.out}
        """


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
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        (
            date
            pixi r -e taxconverter taxconverter mmseqs2 -i {input.tax} -o {output} && date
        ) 2>{log.err} >{log.out}
        """


rule tax_add_missing:
    input:
        fasta=os.path.join(ASSEMBLY_DIR, "{sample}.fna"),
        tax=rules.taxconverter.output,
    output:
        os.path.join(
            RESULTS_DIR,
            "Assembly/intermediate/{sample}_taxonomy/{sample}_tax_good.tsv",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_missing/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_missing/{sample}.err.log"),
    localrule: True
    script:
        "../../scripts/tax_add_missing.py"
