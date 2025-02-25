rule contigs_taxonomy:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_tax_good.tsv"
            ),
            sample=SAMPLES,
        ),
    output:
        touch("status/taxonomy.done"),


rule taxonomy:
    input:
        fasta=rules.assembly_megahit.output.final,
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
    resources:
        runtime=60 * 5,
        mem="400G",
    threads: 60
    conda:
        os.path.join(WORK_DIR, "envs/mmseqs.yaml")
    shell:
        '(date && name=$(echo $(basename {output.taxonomy}) | sed "s/_tax.tsv//g") && '
        "rm -rf $(dirname {output.dbfasta})/* && "
        "mmseqs createdb  {input.fasta} {output.dbfasta} && "
        "mmseqs taxonomy {output.dbfasta} $(dirname {input.file})/GTDB_mmseqs $(dirname {output.taxonomy})/$name {output.tmp} --threads {threads} --tax-lineage 1 && "
        "mmseqs createtsv {output.dbfasta} $(dirname {output.taxonomy})/$name {output.taxonomy} && date) 2> {log.err} > {log.out}"


rule taxconverter:
    input:
        tax=rules.taxonomy.output.taxonomy,
        taxconverter=rules.taxconverter_install.output.file,
    output:
        os.path.join(RESULTS_DIR, "Assembly/{sample}_taxonomy/{sample}_tax_good.tsv"),
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
