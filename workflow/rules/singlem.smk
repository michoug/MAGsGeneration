rule singlem_all:
    input:
        expand(os.path.join(RESULTS_DIR, "SingleM/all/{sample}_otu.tsv"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "SingleM/genomes_otu.tsv"))
    output:
        touch("status/singlem.done")

rule singlem_data:
    output:
        os.path.join(DB_DIR, "SingleM/S3.0.5.metapackage20220806.smpkg.zb")
    conda:
        os.path.join(ENV_DIR, "singlem.yaml")
    localrule: True
    shell:
        """
        singlem data --output-directory ($dirname {output})
        """

rule singlem_reads:
    input:
        read1=os.path.join(READS_DIR, "{sample}_R1_001.fastq.gz"),
        read2=os.path.join(READS_DIR, "{sample}_R2_001.fastq.gz"),
        db=rules.singlem_data.output
    output:
        table=os.path.join(RESULTS_DIR, "SingleM/all/{sample}_otu.tsv"),
        profile=os.path.join(RESULTS_DIR, "SingleM/all/{sample}_profile.tsv")
    conda:
        os.path.join(ENV_DIR, "singlem.yaml")
    threads:
        config["singlem"]["threads"]
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_raw/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_raw/{sample}.err.log")
    shell:
        """
        (date
        export SINGLEM_METAPACKAGE_PATH={input.db}
        singlem pipe --threads {threads} -1 {input.read1} -2 {input.read2} --otu-table {output.table} -p {output.profile}
        date ) 2> {log.err} > {log.out}
        """

rule singlem_genomes:
    input:
        genomes=os.path.join(RESULTS_DIR, "Bins/cluster_galah.txt")
    output:
        table=os.path.join(RESULTS_DIR, "SingleM/genomes_otu.tsv")
    conda:
        os.path.join(ENV_DIR, "singlem.yaml")
    threads:
        config["singlem"]["threads"]
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_genomes.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_genomes.err.log")
    shell:
        """
        (date
        singlem pipe --threads {threads} --otu-table {output.table} --genome-fasta-files $(dirname {input.genomes})/finalBins/*fa &&
        date ) 2> {log.err} > {log.out}
        """