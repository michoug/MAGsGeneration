# GTDBTK taxonomy
rule quality:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/gtdbtk_final")),
        expand(os.path.join(RESULTS_DIR, "Bins/checkm2/quality_report.tsv")),
    output:
        touch("status/quality.done"),


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
        mem_mb=get_mem_mb,
    message:
        "Running GTDB on MAGs"
    shell:
        "(date && export GTDBTK_DATA_PATH={params.db} && gtdbtk classify_wf --cpus {threads} -x fa --genome_dir {input} --out_dir {output} --mash_db $(dirname {params.db})/gtdb_ref_sketch.msh && date) &> >(tee {log})"


use rule checkm2_refine as checkm2_final with:
    input:
        rules.galah.output.final,
    output:
        tsv=os.path.join(RESULTS_DIR, "Bins/checkm2/quality_report.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/checkm/checkm2_final.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/checkm/checkm2_final.err.log"),
    threads: 20
