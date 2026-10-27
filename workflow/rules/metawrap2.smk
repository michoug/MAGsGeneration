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
    # Feeds only the binners listed in selected_binners.txt to
    # `metawrap2 bin_refinement`, as -A/-B/-C... (at most 6 bin sets).
    # CheckM2's database is read from the [databases] table of the
    # metawrap2 TOML, not from $CHECKM2DB, so a per-run copy of
    # config/metawrap.toml is written with it, plus [limits] max_threads so
    # metawrap2 uses every allocated core (its default budget is 50%).
    # When no bin passes -c/-x, metawrap2 exits non-zero ("Consolidation
    # produced no bins"); that is a valid result for a sample, so it is
    # turned into an empty bin folder and a header-only stats table.
    input:
        bins=expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{{sample}}/BinnerEval/{binner}/bins"
            ),
            binner=BINNERS,
        ),
        sel=os.path.join(RESULTS_DIR, "Bins/intermediate/selected_binners.txt"),
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
        config=os.path.join(WORK_DIR, "config/metawrap.toml"),
        bins_dir=os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/BinnerEval"),
        completeness=METAWRAP2_CFG["completeness"],
        contamination=METAWRAP2_CFG["contamination"],
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
            cfg=$outdir/metawrap2.toml
            cat {params.config} >$cfg
            printf '\n\n[databases]\nCHECKM2DB = "%s"\n' {params.db} >>$cfg
            printf '\n[limits]\nmax_threads = %s\n' {threads} >>$cfg
            export METAWRAP2_CONFIG=$cfg
            ARGS=()
            letters=(A B C D E F)
            for i in "${{!B[@]}}"; do
                ARGS+=(-${{letters[$i]}} {params.bins_dir}/${{B[$i]}}/bins)
            done
            status=0
            pixi run -e metawrap2 metawrap2 bin_refinement --force --config $cfg --skip-gunc -o $outdir -t {threads} -c {params.completeness} -x {params.contamination} "${{ARGS[@]}}" || status=$?
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
    # metawrap2 writes bin_NNN.fasta; downstream CheckM2 and Rosella refine
    # read `-x fa`, so copy them under a stable path with that extension.
    input:
        rules.metawrap2_run.output.bins,
    output:
        directory(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/Metawrap2/final_bins"
            )
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/metawrap2_bins/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metawrap2_bins/{sample}.err.log"),
    localrule: True
    message:
        "Collecting Metawrap2 refined bins for {wildcards.sample}"
    shell:
        """
        (
            date
            mkdir -p {output}
            n=0
            for f in {input}/*.fasta; do
                [[ -e "$f" ]] || continue
                cp "$f" {output}/$(basename "$f" .fasta).fa
                n=$((n + 1))
            done
            if [[ $n -eq 0 ]]; then
                echo "WARNING: no refined bins found in {input}" >&2
            fi
            date
        ) 2>{log.err} >{log.out}
        """
