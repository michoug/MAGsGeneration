rule install_software:
    input:
        os.path.join(RESULTS_DIR, "DB/gtdb/database_report.tsv"),
        os.path.join(RESULTS_DIR, "DB/MetaBinner/metabinner.done"),
        os.path.join(RESULTS_DIR, "DB/vamb.done"),
        os.path.join(RESULTS_DIR, "DB/hg38.fa.gz"),
        os.path.join(
            RESULTS_DIR,
            "DB/SingleM/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb/CONTENTS.json",
        ),
    output:
        touch("status/install_software.done"),


rule downloadHuman:
    output:
        file=os.path.join(RESULTS_DIR, "DB/hg38.fa.gz"),
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    conda:
        os.path.join(WORK_DIR, "envs/datasets.yaml")
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
        datasets rehydrate --directory $(dirname {output.file})
        pigz -p 10  $(dirname {output.file})/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna
        mv $(dirname {output.file})/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna.gz {output.file}
        rm -rf $(dirname {output.file})/ncbi_dataset
        rm -f $(dirname {output.file})/md5sum.txt
        rm -f $(dirname {output.file})/README.md
        rm -f ncbi_dataset.zip
        date) 2> {log.err} > {log.out}
        """


# rule taxconverter_install:
#     # input:
#     #     rules.create_db_folder.output.file
#     output:
#         folder=directory(os.path.join(RESULTS_DIR, "DB/taxconverter")),
#         file=os.path.join(RESULTS_DIR, "DB/taxconverter/taxconverter.done"),
#     log:
#         out=os.path.join(RESULTS_DIR, "logs/taxconverter_install.out.log"),
#         err=os.path.join(RESULTS_DIR, "logs/taxconverter_install.err.log"),
#     conda:
#         os.path.join(WORK_DIR, "envs/taxconverter.yaml")
#     priority: 20
#     localrule: True
#     shell:
#         """
#         (date
#         rm -rf {output.folder}
#         rm -rf taxconverter
#         git lfs install
#         git clone https://github.com/RasmussenLab/taxconverter.git
#         mv taxconverter $(dirname {output.folder})
#         pip install -e {output.folder}
#         unzip {output.folder}/data/lineage.zip -d {output.folder}/data
#         touch {output.file}
#         date) 2> {log.err} > {log.out}
#         """


# rule download_tax_db:
#     output:
#         file=os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs_mapping"),
#         look=os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs.lookup"),
#         tmp=temp(directory(os.path.join(RESULTS_DIR, "tmp"))),
#     log:
#         out=os.path.join(RESULTS_DIR, "logs/download_tax_db.out.log"),
#         err=os.path.join(RESULTS_DIR, "logs/download_tax_db.err.log"),
#     resources:
#         mem_mb=get_mem_mb,
#         runtime=6 * 60,
#     threads: 20
#     conda:
#         os.path.join(WORK_DIR, "envs/mmseqs.yaml")
#     shell:
#         '(date && name=$(echo $(basename {output.file}) | sed "s/_mapping//g") && '
#         "rm -f $(dirname {output.file})/$name* && "
#         "mmseqs databases GTDB $(dirname {output.file})/$name {output.tmp} --threads {threads} && date) 2> {log.err} > {log.out}"


rule download_tax_db:
    output:
        file=os.path.join(RESULTS_DIR, "DB/gtdb/database_report.tsv"),
        folder=directory(os.path.join(RESULTS_DIR, "DB/gtdb")),
        tmp=temp(directory(os.path.join(RESULTS_DIR, "tmp"))),
    log:
        out=os.path.join(RESULTS_DIR, "logs/download_tax_db.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/download_tax_db.err.log"),
    resources:
        mem_mb=get_mem_mb,
        runtime=6 * 60,
    threads: 20
    conda:
        os.path.join(WORK_DIR, "envs/metabuli.yaml")
    shell:
        """
        (date
        metabuli databases GTDB $(dirname {output.folder}) {output.tmp} --threads {threads}
        date) 2> {log.err} > {log.out}
        """


rule metabinner_install:
    output:
        dbs=directory(os.path.join(RESULTS_DIR, "DB/MetaBinner")),
        file=os.path.join(RESULTS_DIR, "DB/MetaBinner/metabinner.done"),
    localrule: True
    conda:
        os.path.join(WORK_DIR, "envs/metabinner.yaml")
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


rule vamb_install:
    input:
        rules.metabinner_install.output.file,
    output:
        os.path.join(RESULTS_DIR, "DB/vamb.done"),
    conda:
        os.path.join(WORK_DIR, "envs/vamb.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/vamb_install.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/vamb_install.err.log"),
    priority: 20
    localrule: True
    shell:
        """
        (date
        rm -rf $CONDA_PREFIX/compiler_compat/
        pip install vamb
        touch {output}
        date) 2> {log.err} > {log.out}
        """


rule singlem_data:
    input:
        rules.vamb_install.output,
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
    conda:
        os.path.join(WORK_DIR, "envs/singlem.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_data.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_data.err.log"),
    localrule: True
    priority: 20
    shell:
        """
        (date
        singlem data --output-directory $(dirname {output.folder})
        date) 2> {log.err} > {log.out}
        """
