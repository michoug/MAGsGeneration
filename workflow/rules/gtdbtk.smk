# GTDBTK taxonomy
rule gtdbtk_all:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/gtdbtk_final")),
    output:
        touch("status/gtdbtk.done"),


rule gtdbtk:
    input:
        rules.galah.output.final,
    output:
        directory(os.path.join(RESULTS_DIR, "Bins/gtdbtk_final")),
    log:
        os.path.join(RESULTS_DIR, "logs/gtdbtk.log"),
    conda:
        os.path.join(WORK_DIR, "envs/gtdbtk.yaml")
    params:
        db=config["gtdbtk"]["path"],
    threads: 30
    resources:
        runtime=60 * 12,
        mem="250G",
    message:
        "Running GTDB on MAGs"
    shell:
        "(date && export GTDBTK_DATA_PATH={params.db} && gtdbtk classify_wf --cpus {threads} -x fa --genome_dir {input} --out_dir {output} --mash_db $(dirname {params.db})/gtdb_ref_sketch.msh && date) &> >(tee {log})"
