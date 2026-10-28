# Turn the 7 binners' outputs from mags_generation.smk into per-binner
# contig2bin tables and FASTA folders under BinnerEval/{binner}/. Those
# folders are scored with CheckM2 to keep only the best
# `binners.selection.n` (default 3), and are the inputs of metawrap2_run
# (metawrap2.smk) whether or not selection is enabled.

BINNER_CATALOGUE = config["binners"]["catalogue"]
BINNERS = list(BINNER_CATALOGUE)
SELECTION_CFG = config["binners"]["selection"]
SELECTION_ENABLED = SELECTION_CFG["enabled"]

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


rule prepare_contig2bin:
    # One contig<TAB>bin table per binner, bin names prefixed with the sample
    input:
        metabat=rules.metabat2.output,
        rosella=rules.rosella.output,
        taxvamb=rules.taxvamb.output,
        concoct=rules.concoct.output,
        comebin=rules.comebin.output.tsv,
        semibin=rules.semibin.output,
        metacat=rules.metacat_cluster.output.folder,
        metacat_done=rules.metacat_cluster.output.done,
    output:
        metabat=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/BinnerEval/metabat/contig2bin.tsv"
        ),
        rosella=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/BinnerEval/rosella/contig2bin.tsv"
        ),
        taxvamb=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/BinnerEval/taxvamb/contig2bin.tsv"
        ),
        concoct=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/BinnerEval/concoct/contig2bin.tsv"
        ),
        comebin=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/BinnerEval/comebin/contig2bin.tsv"
        ),
        semibin=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/BinnerEval/semibin/contig2bin.tsv"
        ),
        metacat=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/BinnerEval/metacat/contig2bin.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/prepare_contig2bin/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prepare_contig2bin/{sample}.err.log"),
    localrule: True
    params:
        script=os.path.join(WORK_DIR, "scripts/Fasta_to_Scaffolds2Bin.sh"),
    message:
        "Preparing per-binner contig2bin tables for {wildcards.sample}"
    shell:
        r"""
        (
            date
            bash {params.script} -e fa -i {input.metabat} \
                | perl -pe 's/^(\S+)[^\t]*(\t.*)$/$1$2/' \
                | perl -pe 's/metabat./{wildcards.sample}_metabat_/g' >{output.metabat}
            bash {params.script} -e fna -i $(dirname {input.rosella}) \
                | perl -pe 's/rosella_(bin|refined)_/{wildcards.sample}_rosella_/g' >{output.rosella}
            bash {params.script} -e fna -i $(dirname {input.taxvamb})/bins \
                | perl -pe 's/\t/\t{wildcards.sample}_taxvamb_/g' >{output.taxvamb}
            perl -pe 's/,/\t{wildcards.sample}_concoct_/g' {input.concoct} | tail -n +2 >{output.concoct}
            perl -pe 's/\tgroup/\t{wildcards.sample}_comebin_/g' {input.comebin} >{output.comebin}
            perl -pe 's/\t/\t{wildcards.sample}_semibin_/g' {input.semibin} | tail -n +2 >{output.semibin}
            bash {params.script} -e fasta -i {input.metacat} \
                | perl -pe 's/_metacat\./_metacat_/g' >{output.metacat}
            date
        ) 2>{log.err} >{log.out}
        """


# Always defined: these per-binner FASTA folders are scored here when
# selection is enabled, and are the -A/-B/-C... inputs of metawrap2_run.
rule explode_binner_bins:
    input:
        c2b=lambda wc: getattr(rules.prepare_contig2bin.output, wc.binner).format(
            sample=wc.sample
        ),
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
    message:
        "Rebuilding {wildcards.binner} bins for {wildcards.sample}"
    shell:
        """
        (
            date
            pixi r -e coverm python {params.script} --contig2bin {input.c2b} --assembly {input.assembly} --outdir {output}
            status=$?
            date
            exit $status
        ) 2>{log.err} >{log.out}
        """


if SELECTION_ENABLED:

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
