rule install_software:
    input:
        os.path.join(RESULTS_DIR, "DB/taxconverter/taxconverter.done"),
        os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs_mapping"),
        os.path.join(RESULTS_DIR, "DB/MetaBinner/metabinner.done"),
        os.path.join(RESULTS_DIR, "DB/vamb/install.done"),
    output:
        touch("status/install_software.done"),

rule downloadHuman:
    output:
        file = os.path.join(RESULTS_DIR, "DB/hg38.fa.gz"),
    conda:
        os.path.join(ENV_DIR, "datasets.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/human_download.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/human_download.err.log"),
    threads: 10
    shell:
        """
        (date
        rm -f ncbi_dataset.zip
        rm -rf $(dirname {output.file})/ncbi_dataset
        rm -f $(dirname {output.file})/md5sum.txt
        rm -f $(dirname {output.file})/README.md
        datasets download genome accession GCF_000001405.40 --include genome --assembly-source refseq --dehydrated
        unzip ncbi_dataset.zip -d $(dirname {output.file})
        datasets rehydrate --directory $(dirname /scratch/gmichoud/GlacierIce/mags_all/DB/hg38.fa.gz)
        pigz -p 10  $(dirname {output.file})/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna
        mv $(dirname {output.file})/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna.gz {output.file}
        rm -rf $(dirname {output.file})/ncbi_dataset
        rm -f $(dirname {output.file})/md5sum.txt
        rm -f $(dirname {output.file})/README.md
        rm -f ncbi_dataset.zip
        date) 2> {log.err} > {log.out}
        """

rule taxconverter_install:
    # input:
    #     rules.create_db_folder.output.file
    output:
        folder=directory(os.path.join(RESULTS_DIR, "DB/taxconverter")),
        file=os.path.join(RESULTS_DIR, "DB/taxconverter/taxconverter.done"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxconverter_install.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxconverter_install.err.log"),
    conda:
        os.path.join(ENV_DIR, "taxconverter.yaml")
    localrule: True
    shell:
        """
        (date
        rm -rf {output.folder}
        rm -rf taxconverter
        git clone https://github.com/RasmussenLab/taxconverter.git
        mv taxconverter $(dirname {output.folder})
        pip install -e {output.folder}
        unzip {output.folder}/data/lineage.zip -d {output.folder}/data
        touch {output.file}
        date) 2> {log.err} > {log.out}
        """


rule download_tax_db:
    output:
        file=os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs_mapping"),
        tmp=temp(directory(os.path.join(RESULTS_DIR, "tmp"))),
    log:
        out=os.path.join(RESULTS_DIR, "logs/download_tax_db.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/download_tax_db.err.log"),
    threads: 20
    conda:
        os.path.join(ENV_DIR, "mmseqs.yaml")
    shell:
        '(date && name=$(echo $(basename {output.file}) | sed "s/_mapping//g") && '
        "mmseqs databases GTDB $(dirname {output.file})/$name {output.tmp} --threads {threads} && date) 2> {log.err} > {log.out}"


rule metabinner_install:
    input:
        rules.taxconverter_install.output.file
    output:
        dbs=directory(os.path.join(RESULTS_DIR, "DB/MetaBinner")),
        file=os.path.join(RESULTS_DIR, "DB/MetaBinner/metabinner.done"),
    localrule: True
    conda:
        os.path.join(ENV_DIR, "metabinner.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabinner_install.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabinner_install.err.log"),
    shell:
        """
        (date &&  rm -rf MetaBinner
        git clone https://github.com/ziyewang/MetaBinner.git &> /dev/null
        mv MetaBinner $(dirname {output.dbs})
        touch {output.file} ) 2> {log.err} > {log.out}
        """


rule vamb_install:
    input:
        rules.metabinner_install.output.file
    output:
        dbs=directory(os.path.join(RESULTS_DIR, "DB/vamb")),
        testfile=os.path.join(RESULTS_DIR, "DB/vamb/install.done"),
    conda:
        os.path.join(ENV_DIR, "vamb.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/vamb_install.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/vamb_install.err.log"),
    localrule: True
    shell:
        "(date && rm -rf vamb && "
        "rm -rf {output.dbs} && "
        "git clone https://github.com/RasmussenLab/vamb -b master && "
        "mv vamb {output.dbs} && "
        "pip install -e {output.dbs} && "
        "touch {output.testfile} && "
        "date) 2> {log.err} > {log.out}"
