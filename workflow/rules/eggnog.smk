# Running EggNog-Mapper on the concatenated bins


rule eggnog_all:
    input:
        expand(os.path.join(RESULTS_DIR, "EggNOG/allMAGs_egg.txt.emapper.annotations")),
    output:
        touch("status/eggnog.done"),


# Preparing the MAG files
rule multi_to_single:
    input:
        rules.prepare4egg.output,
    output:
        fixed_fa=os.path.join(RESULTS_DIR, "EggNOG/intermediate/cat_proteins.faa"),
    message:
        "Converting multi-line FASTA to single line FASTA files"
    shell:
        "(date && pixi r -e seqkit seqkit seq $(dirname {input})/*faa > {output.fixed_fa} && date)"


# Preparing EGGNOG run
checkpoint prep_eggnog:
    input:
        fa=rules.multi_to_single.output.fixed_fa,
    output:
        dir=directory(os.path.join(RESULTS_DIR, "EggNOG/intermediate/split")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/split/split.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/split/split.err.log"),
    localrule: True
    message:
        "Splitting the large concactenated protein file into chunks of 1000000"
    shell:
        "(date && rm -f {output.dir}/* && pixi r -e seqkit seqkit split -s 1000000 -O {output.dir} {input.fa} "
        "&& date) 2> {log.err} > {log.out}"


def aggregate_emapper(wildcards):
    checkpoint_output = checkpoints.prep_eggnog.get(**wildcards).output[0]
    return expand(
        os.path.join(
            RESULTS_DIR,
            "EggNOG/intermediate/split/temp/cat_proteins.part_{i}.emapper.seed_orthologs",
        ),
        i=glob_wildcards(os.path.join(checkpoint_output, "cat_proteins.part_{i}")).i,
    )


# EGGNOG mapping to annotations
rule emapper:
    input:
        os.path.join(RESULTS_DIR, "EggNOG/intermediate/split/cat_proteins.part_{i}"),
    output:
        hits=os.path.join(
            RESULTS_DIR, "EggNOG/intermediate/split/temp/cat_proteins.part_{i}.emapper.hits"
        ),
        seed=os.path.join(
            RESULTS_DIR,
            "EggNOG/intermediate/split/temp/cat_proteins.part_{i}.emapper.seed_orthologs",
        ),
    threads: 10
    resources:
        mem_mb=lambda wildcards, threads: threads * config["memory_per_thread"],
        runtime=60 * 12,
    log:
        out=os.path.join(RESULTS_DIR, "logs/eggnog/emapper_{i}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/eggnog/emapper_{i}.err.log"),
    params:
        db=config["eggnog"]["db"],
    message:
        "Running EggNog-mapper on {wildcards.i}"
    shell:
        "(date && "
        "rm -f {input}.emapper.hits && "
        "pixi r -e eggnog emapper.py -m diamond --data_dir {params.db} --no_annot --no_file_comments --cpu {threads} -i {input} -o {input} && "
        "mv {input}.emapper.hits {output.hits} && mv {input}.emapper.seed_orthologs {output.seed} && "
        "date) 2> {log.err} > {log.out}"


rule cat_emapper_out:
    input:
        aggregate_emapper,
    output:
        os.path.join(RESULTS_DIR, "EggNOG/allMAGs.emapper_orthologs"),
    message:
        "Cat EggNog-mapper output"
    localrule: True
    shell:
        "(date && cat {input} > {output} && date)"


rule emapper_final:
    input:
        os.path.join(RESULTS_DIR, "EggNOG/allMAGs.emapper_orthologs"),
    output:
        os.path.join(RESULTS_DIR, "EggNOG/allMAGs_egg.txt.emapper.annotations"),
    threads: 10
    resources:
        runtime=60 * 12,
        mem_mb=get_mem_mb,
    log:
        out=os.path.join(RESULTS_DIR, "logs/eggnog/emapper.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/eggnog/emapper.err.log"),
    params:
        db=config["eggnog"]["db"],
    message:
        "Running EggNog-mapper on All MAGs"
    shell:
        "(date && "
        "pixi r -e eggnog emapper.py --data_dir {params.db} --annotate_hits_table {input} --no_file_comments -o $(echo {output} | sed 's/.emapper.annotations//g' ) --cpu {threads} --dbmem && "
        "date) 2> {log.err} > {log.out}"
