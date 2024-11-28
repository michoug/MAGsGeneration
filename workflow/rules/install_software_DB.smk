rule install_software:
    input:
        os.path.join(RESULTS_DIR, "DB/taxconverter/taxconverter.done"),
        os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs_mapping"),
        os.path.join(RESULTS_DIR, "DB/MetaBinner/metabinner.done"),
        os.path.join(RESULTS_DIR, "DB/vamb/install.done"),
    output:
        touch("status/install_software.done"),


rule taxconverter_install:
    # input:
    #     rules.create_db_folder.output.file
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
        mv taxconverter $(dirname {output.folder})
        pip install -e {output.folder}
        unzip {output.folder}/data/lineage.zip -d {output.folder}/data
        touch {output.file}
        date) 2> {log.err} > {log.out}
        """


rule download_tax_db:
    # input:
    #     db = rules.create_db_folder.output.folder
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
    # input:
    #     rules.create_db_folder.output.file
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
    # input:
    #     rules.create_db_folder.output.file
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
