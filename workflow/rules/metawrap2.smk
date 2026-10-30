if "binning_refinement" in STEPS and "binner_selection" not in STEPS:
    raise Exception(
        "'binning_refinement' requires 'binner_selection' earlier in config "
        "'steps' - it provides Bins/intermediate/selected_binners.txt, even "
        "when binners.selection.enabled is false."
    )


METAWRAP2_CFG = config["metawrap2"]
METAWRAP2_REFINED = "metawrap_%d_%d_bins" % (
    int(METAWRAP2_CFG["completeness"]),
    int(METAWRAP2_CFG["contamination"]),
)


rule binning_refinement:
    input:
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Metawrap2/BIN_REFINEMENT",
                METAWRAP2_REFINED + ".stats.tsv",
            ),
            sample=ASSEMBLY,
        ),
    output:
        touch("status/binning_refinement.done"),


rule metawrap2_run:
    # metawrap2 bin_refinement on the selected binners (-A/-B/-C..., max 6)
    input:
        bins=expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{{sample}}/BinnerEval/{binner}/bins"
            ),
            binner=BINNERS,
        ),
        sel=os.path.join(RESULTS_DIR, "Bins/intermediate/selected_binners.txt"),
        gunc_db=GUNC_DB if GUNC_ENABLED else [],
    output:
        bins=directory(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Metawrap2/BIN_REFINEMENT",
                METAWRAP2_REFINED,
            )
        ),
        stats=os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/{sample}/Metawrap2/BIN_REFINEMENT",
            METAWRAP2_REFINED + ".stats.tsv",
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/metawrap2/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metawrap2/{sample}.err.log"),
    threads: 24
    resources:
        runtime=60 * 12,
        mem_mb=get_mem_mb,
    params:
        db=config["checkm2"]["db"],
        bins_dir=os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/BinnerEval"),
        completeness=METAWRAP2_CFG["completeness"],
        contamination=METAWRAP2_CFG["contamination"],
        gunc_db_line='GUNC_DB = "%s"\\n' % GUNC_DB if GUNC_ENABLED else "",
        gunc_flag="" if GUNC_ENABLED else "--skip-gunc",
        max_disk_fraction=METAWRAP2_CFG.get("max_disk_fraction", 0.95),
    message:
        "Running Metawrap2 bin_refinement on the binners selected by binner_selection"
    shell:
        """
        (
            date
            mapfile -t B <{input.sel}
            if [[ ${{#B[@]}} -eq 0 ]]; then
                echo "ERROR: {input.sel} is empty - no binners selected for Metawrap2" >&2
                exit 1
            fi
            if [[ ${{#B[@]}} -gt 6 ]]; then
                echo "ERROR: ${{#B[@]}} binners selected, metawrap2 bin_refinement accepts at most 6" >&2
                exit 1
            fi
            outdir=$(dirname {output.bins})
            mkdir -p $outdir
            cfg=$(mktemp --suffix=.metawrap2.toml)
            trap 'rm -f $cfg' EXIT
            printf '[settings]\nuse_conda_envs = false\n\n' >$cfg
            printf '[databases]\nCHECKM2DB = "%s"\n' {params.db} >>$cfg
            printf '{params.gunc_db_line}' >>$cfg
            printf '\n[limits]\nmax_threads = %s\nmax_memory_gb = %s\nmax_disk_fraction = %s\n' \
                {threads} $(({resources.mem_mb} / 1024)) {params.max_disk_fraction} >>$cfg
            cat $cfg
            ARGS=()
            letters=(A B C D E F)
            for i in "${{!B[@]}}"; do
                ARGS+=(-${{letters[$i]}} {params.bins_dir}/${{B[$i]}}/bins)
            done
            status=0
            METAWRAP2_CONFIG=$cfg pixi run -e metawrap2 metawrap2 bin_refinement --force --config $cfg {params.gunc_flag} -o $outdir -t {threads} -c {params.completeness} -x {params.contamination} "${{ARGS[@]}}" || status=$?
            if [[ $status -ne 0 ]]; then
                if grep -q "Consolidation produced no bins" $outdir/metawrap2.log 2>/dev/null; then
                    echo "WARNING: no bin of {wildcards.sample} passed -c {params.completeness} -x {params.contamination}, writing an empty bin set" >&2
                    mkdir -p {output.bins}
                    printf 'bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n' >{output.stats}
                else
                    exit $status
                fi
            fi
            date
        ) 2>{log.err} >{log.out}
        """


rule metawrap2_bins:
    # Renames bins to {sample}_bin_NNN.fa and drops GUNC failures
    input:
        bins=rules.metawrap2_run.output.bins,
        stats=rules.metawrap2_run.output.stats,
    output:
        bins=directory(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/Metawrap2/final_bins"
            )
        ),
        table=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Metawrap2/final_bins.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/metawrap2_bins/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metawrap2_bins/{sample}.err.log"),
    localrule: True
    params:
        script=os.path.join(WORK_DIR, "scripts/collect_metawrap2_bins.py"),
        keep_chimeric="" if GUNC_CFG.get("remove_chimeric", True) else "--keep-chimeric",
    message:
        "Collecting Metawrap2 refined bins for {wildcards.sample}"
    shell:
        """
        (
            date
            pixi r -e metawrap2 python {params.script} --bins {input.bins} --stats {input.stats} \
                --outdir {output.bins} --table {output.table} --sample {wildcards.sample} {params.keep_chimeric}
            date
        ) 2>{log.err} >{log.out}
        """
