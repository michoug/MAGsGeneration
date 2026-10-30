# GTDBTK taxonomy
rule quality:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/GDTBTk")),
        expand(os.path.join(RESULTS_DIR, "Bins/checkm2/quality_report.tsv")),
        os.path.join(RESULTS_DIR, "qc/mags/multiqc_report.html"),
    output:
        touch("status/quality.done"),


rule gtdbtk:
    input:
        rules.galah.output.final,
    output:
        directory(os.path.join(RESULTS_DIR, "Bins/GDTBTk")),
    log:
        os.path.join(RESULTS_DIR, "logs/gtdbtk.log"),
    threads: 30
    resources:
        runtime=60 * 12,
        mem_mb=get_mem_mb,
    params:
        db=config["gtdbtk"]["path"],
    message:
        "Running GTDB on MAGs"
    shell:
        """
        (
            date
            export GTDBTK_DATA_PATH={params.db}
            pixi r -e gtdbtk gtdbtk classify_wf --cpus {threads} -x fa --genome_dir {input} --out_dir {output}
            date
        ) &> >(tee {log})
        """


use rule checkm2_refine as checkm2_final with:
    input:
        rules.galah.output.final,
    output:
        tsv=os.path.join(RESULTS_DIR, "Bins/checkm2/quality_report.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/checkm/checkm2_final.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/checkm/checkm2_final.err.log"),
    threads: 20


_MAG_QC_BINNERS = (
    "binner_selection" in STEPS and config["binners"]["selection"]["enabled"]
)


rule mag_counts_table:
    # Refined bins per sample by CheckM2 tier, and the binner selection ranking
    input:
        reports=expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/checkm2/quality_report.tsv"
            ),
            sample=ASSEMBLY,
        ),
        binners=(
            os.path.join(RESULTS_DIR, "Bins/intermediate/binner_selection_report.tsv")
            if _MAG_QC_BINNERS
            else []
        ),
    output:
        mags=os.path.join(RESULTS_DIR, "qc/mags/mag_counts_mqc.tsv"),
        binners=(
            os.path.join(RESULTS_DIR, "qc/mags/binner_selection_mqc.tsv")
            if _MAG_QC_BINNERS
            else []
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/qc/mag_counts_table.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/qc/mag_counts_table.err.log"),
    localrule: True
    params:
        script=os.path.join(WORK_DIR, "scripts/mag_counts_mqc.py"),
        binners=lambda wc, input, output: (
            f"--binner-report {input.binners} --binner-out {output.binners}"
            if _MAG_QC_BINNERS
            else ""
        ),
    message:
        "Counting refined bins per sample and quality tier"
    shell:
        """
        (
            date
            pixi r -e multiqc python {params.script} --reports {input.reports} --out {output.mags} {params.binners}
            date
        ) 2>{log.err} >{log.out}
        """


rule multiqc_mags:
    # CheckM2 and GTDB-Tk on the dereplicated MAGs, plus the per-sample tables
    input:
        checkm2=rules.checkm2_final.output.tsv,
        gtdbtk=rules.gtdbtk.output,
        tables=rules.mag_counts_table.output,
    output:
        os.path.join(RESULTS_DIR, "qc/mags/multiqc_report.html"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/qc/multiqc_mags.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/qc/multiqc_mags.err.log"),
    threads: 1
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    message:
        "Running MultiQC on the MAG quality and taxonomy reports"
    shell:
        """
        (
            date
            pixi r -e multiqc multiqc {input.checkm2} {input.gtdbtk} {input.tables} \
                -m checkm2 -m gtdbtk -m custom_content \
                -o $(dirname {output}) -n $(basename {output}) --force
            date
        ) 2>{log.err} >{log.out}
        """
