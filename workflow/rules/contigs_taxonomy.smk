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


rule download_tax_db:
    output:
        file=os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs_mapping"),
        folder=directory(os.path.join(RESULTS_DIR, "DB")),
        tmp=temp(directory(os.path.join(RESULTS_DIR, "tmp"))),
    log:
        out=os.path.join(RESULTS_DIR, "logs/download_tax_db.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/download_tax_db.err.log"),
    threads: 20
    conda:
        os.path.join(ENV_DIR, "mmseqs.yaml")
    shell:
        '(date && name=$(echo $(basename {output.file}) | sed "s/_mapping//g") && '
        "mmseqs databases GTDB {output.folder}/$name {output.tmp} --threads {threads} && date) 2> {log.err} > {log.out}"


rule taxonomy:
    input:
        fasta=rules.assembly_megahit.output.final,
        db=rules.download_tax_db.output.folder,
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
    threads: 20
    conda:
        os.path.join(ENV_DIR, "mmseqs.yaml")
    shell:
        '(date && name=$(echo $(basename {output.taxonomy}) | sed "s/_tax.tsv//g") && '
        "rm -rf $(dirname {output.dbfasta})/* && "
        "mmseqs createdb  {input.fasta} {output.dbfasta} && "
        "mmseqs taxonomy {output.dbfasta} {input.db}/GTDB_mmseqs $(dirname {output.taxonomy})/$name {output.tmp} --threads {threads} --tax-lineage 1 && "
        "mmseqs createtsv {output.dbfasta} $(dirname {output.taxonomy})/$name {output.taxonomy} && date) 2> {log.err} > {log.out}"


rule taxconverter_install:
    output:
        folder=directory(os.path.join(RESULTS_DIR, "DB/taxconverter")),
        file=os.path.join(RESULTS_DIR, "DB/taxconverter/taxconverter.done"),
    conda:
        os.path.join(ENV_DIR, "vamb.yaml")
    localrule: True
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxconverter_install.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxconverter_install.err.log"),
    shell:
        """
        (date
        rm -rf {output.folder}
        rm -rf taxconverter
        git clone https://github.com/RasmussenLab/taxconverter.git
        mv taxconverter {output.folder}
        pip install -e {output.folder}
        unzip {output.folder}/data/lineage.zip -d {output.folder}/data
        touch {output.file}
        date) 2> {log.err} > {log.out}
        """


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
        os.path.join(ENV_DIR, "vamb.yaml")
    shell:
        "(date && "
        "taxconverter mmseqs2 -i {input.tax} -o {output} && date) 2> {log.err} > {log.out}"
