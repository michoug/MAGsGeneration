rule install_software:
    input:
        # os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs_mapping"),
        os.path.join(RESULTS_DIR, "DB/MetaBinner/metabinner.done"),
        os.path.join(RESULTS_DIR, "DB/MicroEuk_v3/MicroEuk90"),
        os.path.join(RESULTS_DIR, "DB/db_euk.done"),
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


rule download_tax_db:
    output:
        file=os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs_mapping"),
        look=os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs.lookup"),
        tmp=temp(directory(os.path.join(RESULTS_DIR, "GTDB_tmp"))),
    log:
        out=os.path.join(RESULTS_DIR, "logs/download_tax_db.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/download_tax_db.err.log"),
    resources:
        mem_mb=get_mem_mb,
        runtime=6 * 60,
    threads: 20
    conda:
        os.path.join(WORK_DIR, "envs/mmseqs.yaml")
    shell:
        '(date && name=$(echo $(basename {output.file}) | sed "s/_mapping//g") && '
        "rm -f $(dirname {output.file})/$name* && "
        "mmseqs databases GTDB $(dirname {output.file})/$name {output.tmp} --threads {threads} && date) 2> {log.err} > {log.out}"


rule download_tax_db_euk:
    output:
        gzip=temp(os.path.join(RESULTS_DIR, "DB/MicroEuk_v3.tar.gz")),
        db=directory(os.path.join(RESULTS_DIR, "DB/MicroEuk_v3")),
        final=os.path.join(RESULTS_DIR, "DB/MicroEuk_v3/MicroEuk90"),
        prot_id=os.path.join(
            RESULTS_DIR,
            "DB/MicroEuk_v3/MicroEuk_v3/identifier_mapping.proteins.tsv.gz",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/download_tax_db_euk.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/download_tax_db_euk.err.log"),
    resources:
        mem_mb=get_mem_mb,
        runtime=6 * 60,
    threads: 20
    conda:
        os.path.join(WORK_DIR, "envs/mmseqs.yaml")
    shell:
        """
        (date 
        wget -v -O {output.gzip} https://zenodo.org/records/10139451/files/MicroEuk_v3.tar.gz?download=1 
        tar xvzf {output.gzip} -C {output.db}
        gzip -v -d {output.db}/MicroEuk_v3/MicroEuk100.faa.gz
        gzip -v -d -c {output.db}/MicroEuk_v3/MicroEuk90_clusters.tsv.gz | cut -f1  > {output.db}/MicroEuk90.list
        seqkit faidx -j 4 -l {output.db}/MicroEuk90.list {output.db}/MicroEuk_v3/MicroEuk100.faa > {output.db}/MicroEuk90.faa
        mmseqs createdb --compressed 1 {output.db}/MicroEuk90.faa {output.final}
        date) 2> {log.err} > {log.out}
        """


rule create_tax_db_euk:
    input:
        seq=rules.download_tax_db_euk.output.final,
        file=rules.download_tax_db_euk.output.prot_id,
        tax=rules.download_tax_db_euk.output.db,
    output:
        file=os.path.join(RESULTS_DIR, "DB/db_euk.done"),
        ids_micro=temp(os.path.join(RESULTS_DIR, "DB/MicroEuk_v3/ids_micro.txt")),
        ids=os.path.join(RESULTS_DIR, "DB/MicroEuk_v3/ids.txt"),
        tmp=temp(directory(os.path.join(RESULTS_DIR, "tmp"))),
        mapping=os.path.join(RESULTS_DIR, "DB/MicroEuk_v3/MicroEuk90_mapping"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/create_tax_db_euk.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/create_tax_db_euk.err.log"),
    resources:
        mem_mb=get_mem_mb,
        runtime=60 * 12,
    conda:
        os.path.join(WORK_DIR, "envs/mmseqs.yaml")
    threads: 20
    shell:
        """
        (date &&
        zcat {input.file} | cut -f3,4 > {output.ids_micro}
        perl -ne 'print "$2\t$1\n" if /jgi\\|.*\\|(\\d+)\\|[^\t]+\t([a-f0-9]{{32}})/' {output.ids_micro} > {output.ids}
        rm -rf {input.tax}/taxdump.tar.gz
        wget -O {input.tax}/taxdump.tar.gz ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        tar -xzvf {input.tax}/taxdump.tar.gz -C {input.tax}
        mmseqs createtaxdb {input.seq} {output.tmp} --ncbi-tax-dump {input.tax} --tax-mapping-file {output.ids}
        mv {input.tax}/nodes.dmp {input.tax}/MicroEuk90_nodes.dmp
        mv {input.tax}/names.dmp {input.tax}/MicroEuk90_names.dmp
        mv {input.tax}/merged.dmp {input.tax}/MicroEuk90_merged.dmp
        touch {output.file}
        ) 2> {log.err} > {log.out}
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
