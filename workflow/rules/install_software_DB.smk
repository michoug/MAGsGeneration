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
        # os.path.join(RESULTS_DIR, "DB/MetaCAT.done"),
    output:
        touch("status/install_software.done"),


rule downloadHuman:
    output:
        file=os.path.join(RESULTS_DIR, "DB/hg38.fa.gz"),
    params:
        url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
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
        wget \
            --tries=3 \
            --waitretry=30 \
            --timeout=120 \
            -O {output.file} \
            {params.url}
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


rule metacat_install:
    output:
        file=os.path.join(RESULTS_DIR, "DB/MetaCAT.done"),
    localrule: True
    log:
        out=os.path.join(RESULTS_DIR, "logs/metacat_install.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metacat_install.err.log"),
    priority: 20
    shell:
        """
        (date
        pixi r -e metacat pip install https://github.com/liu-congcong/MetaCAT/releases/download/v1.0.3/metacat-1.0.3-py3-none-any.whl
        touch {output.file}
        date) 2> {log.err} > {log.out}
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
