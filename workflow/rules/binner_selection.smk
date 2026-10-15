# Score the 7 binners from mags_generation.smk with CheckM2 and keep only
# the best `binners.selection.n` (default 3) to feed DAS Tool.
#
# BINNER_C2B mirrors the outputs of `prepare_dasTool` in dastool.smk as
# literal paths rather than `rules.prepare_dasTool.output.*`: this step is
# included before dastool.smk in the pipeline order, so the `rules.`
# attribute is not yet available when this file is parsed. Keep the two in
# sync if binner names or the DasTool output layout change.

BINNER_CATALOGUE = config["binners"]["catalogue"]
BINNERS = list(BINNER_CATALOGUE)
SELECTION_CFG = config["binners"]["selection"]
SELECTION_ENABLED = SELECTION_CFG["enabled"]

BINNER_C2B = {
    b: os.path.join(RESULTS_DIR, f"Bins/intermediate/{{sample}}/DasTool/{b}_das.tsv")
    for b in BINNERS
}

_binner_selection_targets = [
    os.path.join(RESULTS_DIR, "Bins/intermediate/selected_binners.txt")
]
if SELECTION_ENABLED:
    _binner_selection_targets.append(
        os.path.join(RESULTS_DIR, "Bins/intermediate/binner_selection_report.tsv")
    )


rule binner_selection:
    input:
        _binner_selection_targets,
    output:
        touch("status/binner_selection.done"),


if SELECTION_ENABLED:

    rule explode_binner_bins:
        input:
            c2b=lambda wc: BINNER_C2B[wc.binner].format(sample=wc.sample),
            assembly=rules.filter_length.output.out,
        output:
            directory(
                os.path.join(
                    RESULTS_DIR, "Bins/intermediate/{sample}/BinnerEval/{binner}/bins"
                )
            ),
        log:
            out=os.path.join(
                RESULTS_DIR, "logs/explode_binner_bins/{sample}_{binner}.out.log"
            ),
            err=os.path.join(
                RESULTS_DIR, "logs/explode_binner_bins/{sample}_{binner}.err.log"
            ),
        wildcard_constraints:
            binner="|".join(BINNERS),
        threads: 2
        resources:
            runtime=60,
            mem_mb=get_mem_mb,
        params:
            script=os.path.join(WORK_DIR, "scripts/contig2bin_to_fasta.py"),
            min_size=SELECTION_CFG["min_bin_size"],
        message:
            "Rebuilding {wildcards.binner} bins for {wildcards.sample} for independent scoring"
        shell:
            """
            (
                date
                python3 {params.script} --contig2bin {input.c2b} --assembly {input.assembly} --outdir {output} --min-size {params.min_size}
                status=$?
                date
                exit $status
            ) 2>{log.err} >{log.out}
            """

    rule checkm2_binner:
        input:
            bins=rules.explode_binner_bins.output,
        output:
            tsv=os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/BinnerEval/{binner}/checkm2/quality_report.tsv",
            ),
        log:
            out=os.path.join(
                RESULTS_DIR, "logs/checkm2_binner/{sample}_{binner}.out.log"
            ),
            err=os.path.join(
                RESULTS_DIR, "logs/checkm2_binner/{sample}_{binner}.err.log"
            ),
        wildcard_constraints:
            binner="|".join(BINNERS),
        threads: 10
        resources:
            runtime=60 * 4,
            mem_mb=get_mem_mb,
        params:
            ext="fa",
            db=config["checkm2"]["db"],
        message:
            "Running CheckM2 on {wildcards.binner} bins for {wildcards.sample}"
        shell:
            """
            (
                date
                if compgen -G "{input.bins}/*.{params.ext}" >/dev/null 2>&1; then
                    export CHECKM2DB={params.db}
                    pixi r -e checkm2 checkm2 predict --threads {threads} -x {params.ext} --input {input.bins} --output-directory $(dirname {output.tsv}) --force
                else
                    mkdir -p $(dirname {output.tsv})
                    echo -e "Name\tCompleteness\tContamination" >{output.tsv}
                fi
                date
            ) 2>{log.err} >{log.out}
            """

    rule rank_binners:
        input:
            reports=expand(
                os.path.join(
                    RESULTS_DIR,
                    "Bins/intermediate/{sample}/BinnerEval/{binner}/checkm2/quality_report.tsv",
                ),
                sample=ASSEMBLY,
                binner=BINNERS,
            ),
        output:
            os.path.join(RESULTS_DIR, "Bins/intermediate/binner_ranking.tsv"),
        localrule: True
        message:
            "Scoring binners on CheckM2 MQ/NC/NC-strict counts across all samples"
        run:
            import sys

            sys.path.insert(0, os.path.join(WORK_DIR, "scripts"))
            from rank_binners import load_counts

            counts_df = load_counts(
                list(input.reports), tiers=SELECTION_CFG["tiers"]
            )
            counts_df.to_csv(output[0], sep="\t", index=False)

    rule select_binners:
        input:
            ranking=rules.rank_binners.output[0],
        output:
            selected=os.path.join(RESULTS_DIR, "Bins/intermediate/selected_binners.txt"),
            report=os.path.join(
                RESULTS_DIR, "Bins/intermediate/binner_selection_report.tsv"
            ),
        localrule: True
        message:
            "Selecting the top binners for DAS Tool"
        run:
            import sys

            sys.path.insert(0, os.path.join(WORK_DIR, "scripts"))
            from rank_binners import select
            import pandas as pd

            counts_df = pd.read_csv(input.ranking, sep="\t")
            result = select(
                counts_df,
                BINNER_CATALOGUE,
                n=SELECTION_CFG["n"],
                one_per_family=SELECTION_CFG["one_per_family"],
                tier_names=tuple(SELECTION_CFG["tiers"].keys()),
            )
            result["report"].to_csv(output.report, sep="\t", index=False)
            with open(output.selected, "w") as fh:
                for binner in result["selected"]:
                    fh.write(binner + "\n")
            for binner in result["selected"]:
                row = result["report"][result["report"]["binner"] == binner].iloc[0]
                print(
                    f"selected {binner}: family={row['family']} score={row['score']:.3f}"
                )

else:

    rule select_binners_noop:
        output:
            selected=os.path.join(RESULTS_DIR, "Bins/intermediate/selected_binners.txt"),
        localrule: True
        params:
            binners=BINNERS,
        message:
            "binners.selection.enabled is false - passing all binners through to DAS Tool"
        run:
            with open(output.selected, "w") as fh:
                for binner in params.binners:
                    fh.write(binner + "\n")
