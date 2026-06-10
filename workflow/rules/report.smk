_STEP_DONE_FILES = {
    "trim_reads": "status/trimming.done",
    "assembly_all": "status/assembly.done",
    "contigs_taxonomy": "status/taxonomy.done",
    "check_assembly": "status/checkAssembly.done",
    "mags_generation": "status/mags_generation.done",
    "binning_refinement": "status/binning_refinement.done",
    "eukaryotes_filter": "status/eukaryotes_filter.done",
    "checkm": "status/checkm2.done",
    "refine_all": "status/refine.done",
    "filter_checkm": "status/filter_checkm2.done",
    "dereplicate": "status/dereplicate.done",
    "syntracker": "status/syntracker.done",
    "annotations": "status/annotations.done",
    "coverage": "status/coverage.done",
    "quality": "status/quality.done",
    "mgthermometer": "status/mgthermometer.done",
    "eggnog": "status/eggnog.done",
    "singlem_all": "status/singlem.done",
}


rule report:
    input:
        qmd=os.path.join(WORK_DIR, "report/MAGs_summary.qmd"),
        done=[v for k, v in _STEP_DONE_FILES.items() if k in STEPS],
    output:
        html=os.path.join(RESULTS_DIR, "report/MAGs_summary.html"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/report/report.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/report/report.err.log"),
    threads: 1
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    params:
        results_dir=RESULTS_DIR,
        outdir=os.path.join(RESULTS_DIR, "report"),
        sample_list=os.path.join(WORK_DIR, "config/sample_list.txt"),
        steps=",".join(STEPS),
    message:
        "Rendering Quarto summary report"
    shell:
        """
        (
            date
            mkdir -p {params.outdir}
            pixi r -e report quarto render {input.qmd} \
                -P results_dir:{params.results_dir} \
                -P sample_list:{params.sample_list} \
                -P steps:"{params.steps}" \
                --output-dir {params.outdir}
            date
        ) 2>{log.err} >{log.out}
        """
