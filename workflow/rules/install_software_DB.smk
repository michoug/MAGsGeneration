rule install_software:
    input:
        os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs_mapping"),
        # os.path.join(RESULTS_DIR, "DB/gtdb/database_report.tsv"),
        os.path.join(RESULTS_DIR, "DB/hg38.fa.gz"),
        os.path.join(
            RESULTS_DIR,
            "DB/SingleM/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb/CONTENTS.json",
        ),
        os.path.join(RESULTS_DIR, "DB/Busco_euk/busco_euk.done"),
        os.path.join(RESULTS_DIR, "DB/SynTracker/syntracker.done"),
    output:
        touch("status/install_software.done"),


rule downloadHuman:
    output:
        file=os.path.join(RESULTS_DIR, "DB/hg38.fa.gz"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/human_download.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/human_download.err.log"),
    threads: 10
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    params:
        url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
    shell:
        """
        (
            date
            wget \
                --tries=3 \
                --waitretry=30 \
                --timeout=120 \
                -O {output.file} \
                {params.url}
            date
        ) 2>{log.err} >{log.out}
        """


# rule download_tax_db:
#     output:
#         file=os.path.join(RESULTS_DIR, "DB/gtdb/database_report.tsv"),
#         folder=directory(os.path.join(RESULTS_DIR, "DB/gtdb")),
#         tmp=temp(directory(os.path.join(RESULTS_DIR, "tmp"))),
#         done=touch(os.path.join(RESULTS_DIR, "DB/gtdb/gtdb.done")),
#     log:
#         out=os.path.join(RESULTS_DIR, "logs/download_tax_db.out.log"),
#         err=os.path.join(RESULTS_DIR, "logs/download_tax_db.err.log"),
#     threads: 20
#     resources:
#         mem_mb=get_mem_mb,
#         runtime=6 * 60,
#     shell:
#         """
#         (date
#         pixi r -e metabuli metabuli databases GTDB $(dirname {output.folder}) {output.tmp} --threads {threads}
#         touch {output.done}
#         date) 2> {log.err} > {log.out}
#         """


rule download_tax_db:
    output:
        file=os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs_mapping"),
        look=os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs.lookup"),
        tmp=temp(directory(os.path.join(RESULTS_DIR, "tmp"))),
    log:
        out=os.path.join(RESULTS_DIR, "logs/download_tax_db.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/download_tax_db.err.log"),
    threads: 20
    resources:
        mem_mb=get_mem_mb,
        runtime=6 * 60,
    shell:
        """
        (
            date
            name=$(echo $(basename {output.file}) | sed 's/_mapping//g')
            rm -f $(dirname {output.file})/$name*
            pixi r -e mmseqs mmseqs databases GTDB $(dirname {output.file})/$name {output.tmp} --threads {threads}
            date
        ) 2>{log.err} >{log.out}
        """


rule metacat_install:
    output:
        file=os.path.join(RESULTS_DIR, "DB/MetaCAT.done"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/metacat_install.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metacat_install.err.log"),
    priority: 20
    localrule: True
    shell:
        """
        (
            date
            pixi r -e metacat pip install https://github.com/liu-congcong/MetaCAT/releases/download/v1.0.4/metacat-1.0.4-py3-none-any.whl
            touch {output.file}
            date
        ) 2>{log.err} >{log.out}
        """


rule singlem_data:
    output:
        file=os.path.join(
            RESULTS_DIR,
            "DB/SingleM/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb/CONTENTS.json",
        ),
        folder=directory(
            os.path.join(
                RESULTS_DIR,
                "DB/SingleM/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb",
            )
        ),
        filegz=temp(
            os.path.join(
                RESULTS_DIR,
                "DB/SingleM/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb.tar.gz",
            )
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_data.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_data.err.log"),
    priority: 20
    localrule: True
    shell:
        """
        (
            date
            mkdir -p $(dirname {output.filegz})
            wget https://zenodo.org/records/20150069/files/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb.tar.gz -O {output.filegz}
            tar -xzf {output.filegz} -C $(dirname {output.filegz})
            date
        ) 2>{log.err} >{log.out}
        """


rule busco_db_euk:
    output:
        dir=directory(os.path.join(RESULTS_DIR, "DB/Busco_euk")),
        file=os.path.join(RESULTS_DIR, "DB/Busco_euk/file_versions.tsv"),
        done=touch(os.path.join(RESULTS_DIR, "DB/Busco_euk/busco_euk.done")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/busco_db_euk.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/busco_db_euk.err.log"),
    threads: 2
    resources:
        mem_mb=get_mem_mb,
        runtime=10 * 60,
    shell:
        """
        (
            date
            rm -rf {output.dir}
            pixi r -e busco busco --download all --download_path {output.dir}
            touch {output.done}
            date
        ) 2>{log.err} >{log.out}
        """


rule syntracker_install:
    output:
        dbs=directory(os.path.join(RESULTS_DIR, "DB/SynTracker")),
        file=os.path.join(RESULTS_DIR, "DB/SynTracker/syntracker.done"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/syntracker_install.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/syntracker_install.err.log"),
    priority: 20
    localrule: True
    params:
        # Pin to a tagged release instead of tracking the default branch head.
        version="1.4.0",
    shell:
        """
        (
            date
            rm -rf SynTracker
            git clone --branch {params.version} --depth 1 https://github.com/leylabmpi/SynTracker.git
            mv SynTracker $(dirname {output.dbs})
            touch {output.file}
        ) 2>{log.err} >{log.out}
        """
