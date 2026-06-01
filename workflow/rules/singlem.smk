rule singlem_all:
    input:
        expand(
            os.path.join(RESULTS_DIR, "SingleM/intermediate/{sample}_otu.tsv"),
            sample=SAMPLES,
        ),
        expand(os.path.join(RESULTS_DIR, "SingleM/genomes_otu.tsv")),
        expand(os.path.join(RESULTS_DIR, "SingleM/SingleM_MF.tsv")),
        expand(os.path.join(RESULTS_DIR, "SingleM/SingleM_appraise.tsv")),
        expand(os.path.join(RESULTS_DIR, "SingleM/clustered_all.otu_table.tsv")),
        expand(
            os.path.join(
                RESULTS_DIR,
                "SingleM/intermediate/Unifrac/otu_table.S3.1.ribosomal_protein_L2_rplB.unifrac",
            )
        ),
        expand(
            os.path.join(
                RESULTS_DIR,
                "SingleM/Betadiv/otu_table.S3.1.ribosomal_protein_L2_rplB.ebd",
            )
        ),
    output:
        touch("status/singlem.done"),


rule singlem_reads:
    input:
        read1=os.path.join(READS_DIR, "{sample}_R1_001.fastq.gz"),
        read2=os.path.join(READS_DIR, "{sample}_R2_001.fastq.gz"),
        db=rules.singlem_data.output.folder,
    output:
        table=os.path.join(RESULTS_DIR, "SingleM/intermediate/{sample}_otu.tsv"),
        profile=os.path.join(RESULTS_DIR, "SingleM/intermediate/{sample}_profile.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_raw/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_raw/{sample}.err.log"),
    threads: 12
    resources:
        runtime=60 * 8,
        mem_mb=get_mem_mb,
    shell:
        """
        (date
        export SINGLEM_METAPACKAGE_PATH={input.db}
        pixi r -e singlem singlem pipe --threads {threads} -1 {input.read1} -2 {input.read2} --otu-table {output.table} -p {output.profile} --metapackage {input.db}
        date ) 2> {log.err} > {log.out}
        """


rule singlem_cluster:
    input:
        expand(
            os.path.join(RESULTS_DIR, "SingleM/intermediate/{sample}_otu.tsv"),
            sample=SAMPLES,
        ),
    output:
        os.path.join(RESULTS_DIR, "SingleM/clustered_all.otu_table.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem/singlem_cluster.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem/singlem_cluster.err.log"),
    threads: 2
    resources:
        runtime=60 * 6,
        mem_mb=get_mem_mb,
    shell:
        """
        (date
        pixi r -e singlem singlem summarise --input-otu-tables {input} --cluster --clustered-output-otu-table {output}
        date ) 2> {log.err} > {log.out}
        """


rule singlem_unifrac:
    input:
        rules.singlem_cluster.output,
    output:
        os.path.join(
            RESULTS_DIR,
            "SingleM/intermediate/Unifrac/otu_table.S3.1.ribosomal_protein_L2_rplB.unifrac",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem/singlem_unifrac.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem/singlem_unifrac.err.log"),
    threads: 2
    resources:
        runtime=60 * 1,
        mem_mb=get_mem_mb,
    params:
        "otu_table",
    shell:
        """
        (date
        pixi r -e singlem singlem summarise --input-otu-table {input} --unifrac-by-otu $(dirname {output})/{params}
        date ) 2> {log.err} > {log.out}
        """


rule singlem_betadiv:
    input:
        rules.singlem_unifrac.output,
    output:
        ebd=os.path.join(
            RESULTS_DIR, "SingleM/Betadiv/otu_table.S3.1.ribosomal_protein_L2_rplB.ebd"
        ),
        tree=os.path.join(
            RESULTS_DIR, "SingleM/Betadiv/otu_table.S3.1.ribosomal_protein_L2_rplB.tre"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem/singlem_betadiv.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem/singlem_betadiv.err.log"),
    threads: 2
    resources:
        runtime=60 * 1,
        mem_mb=get_mem_mb,
    params:
        "otu_table",
    shell:
        """
        (date
        for file in $(dirname {input})/*.unifrac; do
            pixi r -e singlem convertToEBD.py $file $(dirname {output.ebd})/$(basename $file .unifrac).ebd
            pixi r -e singlem ExpressBetaDiversity -s $(dirname {output.tree})/$(basename $file .unifrac).ebd -c Bray-Curtis -p $(dirname {output.tree})/$(basename $file .unifrac)
        done
        ) 2> {log.err} > {log.out}
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
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_MF.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_MF.err.log"),
    resources:
        runtime=60 * 5,
        mem_mb=get_mem_mb,
    shell:
        """
        (date
        export SINGLEM_METAPACKAGE_PATH={input.db}
        pixi r -e singlem singlem summarise --input-taxonomic-profiles {input.profile} --output-taxonomic-profile {output.profile} --metapackage {input.db}
        pixi r -e singlem singlem microbial_fraction -1 {input.read1} -2 {input.read2} -p {output.profile} > {output.mf} --metapackage {input.db}
        date ) 2> {log.err} > {log.out}
        """


rule singlem_genomes:
    input:
        genomes=os.path.join(RESULTS_DIR, "Bins/intermediate/cluster_galah.txt"),
        folder=os.path.join(RESULTS_DIR, "Bins/finalBins"),
        db=rules.singlem_data.output.folder,
    output:
        table=os.path.join(RESULTS_DIR, "SingleM/genomes_otu.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_genomes.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_genomes.err.log"),
    threads: 16
    resources:
        runtime=60 * 2,
        mem_mb=get_mem_mb,
    shell:
        """
        (date
        export SINGLEM_METAPACKAGE_PATH={input.db}
        pixi r -e singlem singlem pipe --threads {threads} --otu-table {output.table} --genome-fasta-files {input.folder}/*fa --metapackage {input.db} &&
        date ) 2> {log.err} > {log.out}
        """


rule singlem_appraise:
    input:
        metagenomes=expand(
            os.path.join(RESULTS_DIR, "SingleM/intermediate/{sample}_otu.tsv"),
            sample=SAMPLES,
        ),
        genomes=os.path.join(RESULTS_DIR, "SingleM/genomes_otu.tsv"),
        db=rules.singlem_data.output.folder,
    output:
        os.path.join(RESULTS_DIR, "SingleM/SingleM_appraise.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/singlem_appraise.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/singlem_appraise.err.log"),
    localrule: True
    shell:
        """
        (date
        export SINGLEM_METAPACKAGE_PATH={input.db}
        pixi r -e singlem singlem appraise --metagenome-otu-tables {input.metagenomes} --genome-otu-tables {input.genomes} --imperfect --metapackage {input.db} > {output}
        date ) 2> {log.err} > {log.out}
        """
