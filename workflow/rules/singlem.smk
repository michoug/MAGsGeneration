rule singlem_all:
    input:
        expand(
            os.path.join(RESULTS_DIR, "SingleM/all/{sample}_otu.tsv"), sample=SAMPLES
        ),
        expand(os.path.join(RESULTS_DIR, "SingleM/genomes_otu.tsv")),
        expand(os.path.join(RESULTS_DIR, "SingleM/SingleM_MF.tsv")),
        expand(os.path.join(RESULTS_DIR, "SingleM/SingleM_appraise.tsv")),
    output:
        touch("status/singlem.done"),


rule singlem_data:
    output:
        file=os.path.join(
            DB_DIR,
            "SingleM/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb/CONTENTS.json",
        ),
        folder=os.path.join(
            DB_DIR, "SingleM/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb"
        ),
    conda:
        os.path.join(ENV_DIR, "singlem.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_data.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_data.err.log"),
    localrule: True
    shell:
        """
        singlem data --output-directory $(dirname {output.folder})
        """


rule singlem_reads:
    input:
        read1=os.path.join(READS_DIR, "{sample}_R1_001.fastq.gz"),
        read2=os.path.join(READS_DIR, "{sample}_R2_001.fastq.gz"),
        db=rules.singlem_data.output.folder,
    output:
        table=os.path.join(RESULTS_DIR, "SingleM/all/{sample}_otu.tsv"),
        profile=os.path.join(RESULTS_DIR, "SingleM/all/{sample}_profile.tsv"),
    conda:
        os.path.join(ENV_DIR, "singlem.yaml")
    threads: config["singlem"]["threads"]
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_raw/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_raw/{sample}.err.log"),
    shell:
        """
        (date
        export SINGLEM_METAPACKAGE_PATH={input.db}
        singlem pipe --threads {threads} -1 {input.read1} -2 {input.read2} --otu-table {output.table} -p {output.profile}
        date ) 2> {log.err} > {log.out}
        """


rule singlem_microbial:
    input:
        read1=expand(
            os.path.join(READS_DIR, "{sample}_R1_001.fastq.gz"), sample=SAMPLES
        ),
        read2=expand(
            os.path.join(READS_DIR, "{sample}_R2_001.fastq.gz"), sample=SAMPLES
        ),
        db=rules.singlem_data.output.folder,
        profile=expand(rules.singlem_reads.output.profile, sample=SAMPLES),
    output:
        profile=os.path.join(RESULTS_DIR, "SingleM/allProfiles.tsv"),
        mf=os.path.join(RESULTS_DIR, "SingleM/SingleM_MF.tsv"),
    conda:
        os.path.join(ENV_DIR, "singlem.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_MF.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_MF.err.log"),
    shell:
        """
        (date
        export SINGLEM_METAPACKAGE_PATH={input.db}
        singlem summarise --input-taxonomic-profiles {input.profile} --output-taxonomic-profile {output.profile}
        singlem microbial_fraction -1 {input.read1} -2 {input.read2} -p {output.profile} > {output.mf}
        date ) 2> {log.err} > {log.out}
        """


rule singlem_genomes:
    input:
        genomes=os.path.join(RESULTS_DIR, "Bins/cluster_galah.txt"),
        db=rules.singlem_data.output.folder,
    output:
        table=os.path.join(RESULTS_DIR, "SingleM/genomes_otu.tsv"),
    conda:
        os.path.join(ENV_DIR, "singlem.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_genomes.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_genomes.err.log"),
    shell:
        """
        (date
        export SINGLEM_METAPACKAGE_PATH={input.db}
        singlem pipe --otu-table {output.table} --genome-fasta-files $(dirname {input.genomes})/finalBins/*fa &&
        date ) 2> {log.err} > {log.out}
        """


rule singlem_appraise:
    input:
        metagenomes=expand(
            os.path.join(RESULTS_DIR, "SingleM/all/{sample}_otu.tsv"), sample=SAMPLES
        ),
        genomes=os.path.join(RESULTS_DIR, "SingleM/genomes_otu.tsv"),
    output:
        os.path.join(RESULTS_DIR, "SingleM/SingleM_appraise.tsv"),
    conda:
        os.path.join(ENV_DIR, "singlem.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_appraise.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_appraise.err.log"),
    localrule: True
    shell:
        """
        (date
        singlem appraise --metagenome-otu-tables {input.metagenomes} --genome-otu-tables {input.genomes} --imperfect > {output}
        date ) 2> {log.err} > {log.out}
        """
