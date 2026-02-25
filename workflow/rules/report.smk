rule report:
    input:
        qmd=os.path.join(WORK_DIR, "report/MAGs_summary.qmd"),
    output:
        html=os.path.join(RESULTS_DIR, "report/MAGs_summary.html"),
    params:
        results_dir=RESULTS_DIR,
        outdir=os.path.join(RESULTS_DIR, "report"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/report/report.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/report/report.err.log"),
    threads: 1
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    message:
        "Rendering Quarto summary report"
    shell:
        """
        (date
        mkdir -p {params.outdir}
        pixi r -e report quarto render {input.qmd} \
            -P results_dir:{params.results_dir} \
            --output-dir {params.outdir}
        date) 2> {log.err} > {log.out}
        """
