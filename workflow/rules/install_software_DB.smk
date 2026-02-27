rule install_software:
    input:
        os.path.join(RESULTS_DIR, "DB/gtdb/database_report.tsv"),
        os.path.join(RESULTS_DIR, "DB/MetaBinner/metabinner.done"),
        os.path.join(RESULTS_DIR, "DB/hg38.fa.gz"),
        os.path.join(
            RESULTS_DIR,
            "DB/SingleM/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb/CONTENTS.json",
        ),
        os.path.join(RESULTS_DIR, "DB/Busco_euk/busco_euk.done"),
        os.path.join(RESULTS_DIR, "DB/SynTracker/syntracker.done"),
        os.path.join(RESULTS_DIR, "DB/Genomad/genomad.done"),
        os.path.join(RESULTS_DIR, "DB/AMG/README.txt"),
    output:
        touch("status/install_software.done"),


rule downloadHuman:
    output:
        file=os.path.join(RESULTS_DIR, "DB/hg38.fa.gz"),
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
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
        pixi r -e datasets datasets download genome accession GCF_000001405.40 --include genome --assembly-source refseq --dehydrated
        unzip ncbi_dataset.zip -d $(dirname {output.file})
        pixi r -e datasets datasets rehydrate --directory $(dirname {output.file})
        pigz -p 10  $(dirname {output.file})/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna
        mv $(dirname {output.file})/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna.gz {output.file}
        rm -rf $(dirname {output.file})/ncbi_dataset
        rm -f $(dirname {output.file})/md5sum.txt
        rm -f $(dirname {output.file})/README.md
        rm -f ncbi_dataset.zip
        date) 2> {log.err} > {log.out}
        """


rule download_tax_db:
    output:
        file=os.path.join(RESULTS_DIR, "DB/gtdb/database_report.tsv"),
        folder=directory(os.path.join(RESULTS_DIR, "DB/gtdb")),
        tmp=temp(directory(os.path.join(RESULTS_DIR, "tmp"))),
        done=touch(os.path.join(RESULTS_DIR, "DB/gtdb/gtdb.done")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/download_tax_db.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/download_tax_db.err.log"),
    resources:
        mem_mb=get_mem_mb,
        runtime=6 * 60,
    threads: 20
    shell:
        """
        (date
        pixi r -e metabuli metabuli databases GTDB $(dirname {output.folder}) {output.tmp} --threads {threads}
        touch {output.done}
        date) 2> {log.err} > {log.out}
        """


rule metabinner_install:
    output:
        dbs=directory(os.path.join(RESULTS_DIR, "DB/MetaBinner")),
        file=os.path.join(RESULTS_DIR, "DB/MetaBinner/metabinner.done"),
    localrule: True
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabinner_install.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabinner_install.err.log"),
    priority: 20
    shell:
        """
        (date &&  rm -rf MetaBinner
        git clone https://github.com/ziyewang/MetaBinner.git &> /dev/null
        mv MetaBinner $(dirname {output.dbs})
        touch {output.file} ) 2> {log.err} > {log.out}
        """


rule singlem_data:
    output:
        file=os.path.join(
            RESULTS_DIR,
            "DB/SingleM/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb/CONTENTS.json",
        ),
        folder=directory(
            os.path.join(
                RESULTS_DIR,
                "DB/SingleM/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb",
            )
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_data.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_data.err.log"),
    localrule: True
    priority: 20
    shell:
        """
        (date
        pixi r -e singlem singlem data --output-directory $(dirname {output.folder})
        date) 2> {log.err} > {log.out}
        """


rule busco_db_euk:
    output:
        dir=directory(os.path.join(RESULTS_DIR, "DB/Busco_euk")),
        file=os.path.join(RESULTS_DIR, "DB/Busco_euk/file_versions.tsv"),
        done=touch(os.path.join(RESULTS_DIR, "DB/Busco_euk/busco_euk.done")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/busco_db_euk.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/busco_db_euk.err.log"),
    resources:
        mem_mb=get_mem_mb,
        runtime=2 * 60,
    threads: 2
    shell:
        """
        (date
        rm -rf {output.dir}
        pixi r -e busco busco --download all --download_path {output.dir}
        touch {output.done}
        date) 2> {log.err} > {log.out}
        """


rule syntracker_install:
    output:
        dbs=directory(os.path.join(RESULTS_DIR, "DB/SynTracker")),
        file=os.path.join(RESULTS_DIR, "DB/SynTracker/syntracker.done"),
    localrule: True
    log:
        out=os.path.join(RESULTS_DIR, "logs/syntracker_install.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/syntracker_install.err.log"),
    priority: 20
    shell:
        """
        (date &&  rm -rf SynTracker
        git clone https://github.com/leylabmpi/SynTracker.git &> /dev/null
        mv SynTracker $(dirname {output.dbs})
        touch {output.file} ) 2> {log.err} > {log.out}
        """


rule genomad_db:
    output:
        dir=directory(os.path.join(RESULTS_DIR, "DB/Genomad")),
        # file=os.path.join(RESULTS_DIR, "DB/Genomad/file_versions.tsv"),
        done=os.path.join(RESULTS_DIR, "DB/Genomad/genomad.done"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/genomad_db.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/genomad_db.err.log"),
    resources:
        mem_mb=get_mem_mb,
        runtime=2 * 60,
    threads: 2
    shell:
        """
        (date
        rm -rf {output.dir}
        mkdir -p {output.dir}
        wget https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz -O {output.dir}/genomad_db_v1.9.tar.gz
        tar -xzf {output.dir}/genomad_db_v1.9.tar.gz -C {output.dir} --strip-components=1
        rm {output.dir}/genomad_db_v1.9.tar.gz
        touch {output.done}
        date) 2> {log.err} > {log.out}
        """


rule checkamg_db:
    output:
        file=os.path.join(RESULTS_DIR, "DB/AMG/README.txt"),
        folder=directory(os.path.join(RESULTS_DIR, "DB/AMG")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/check_amg_db.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/check_amg_db.err.log"),
    resources:
        mem_mb=get_mem_mb,
        runtime=2 * 60,
    threads: 2
    shell:
        """
        (date
        rm -rf {output.folder}
        mkdir -p {output.folder}
        wget https://zenodo.org/records/18407279/files/CheckAMG_db_v1_20260128.tar.gz?download=1 -O {output.folder}/CheckAMG_db_v1_20260128.tar.gz
        tar -xzf {output.folder}/CheckAMG_db_v1_20260128.tar.gz -C {output.folder} --strip-components=1
        rm {output.folder}/CheckAMG_db_v1_20260128.tar.gz
        date) 2> {log.err} > {log.out}
        """
