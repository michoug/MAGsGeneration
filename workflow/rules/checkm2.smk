# Running CheckM on the bins


rule checkm_all:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/checkm2/quality_report.tsv"
            ),
            sample=ASSEMBLY,
        ),
    output:
        touch("status/checkm2.done"),


# Checking bin quality
rule checkm2:
    input:
        bins=rules.metawrap2_bins.output,
    output:
        tsv=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/checkm2/quality_report.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/checkm/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/checkm/{sample}.err.log"),
    threads: 10
    resources:
        runtime=60 * 4,
        mem_mb=get_mem_mb,
    params:
        ext="fa",
        db=config["checkm2"]["db"],
    message:
        "Running CheckM2 for directory renamed bins"
    shell:
        """
        (
            date
            if compgen -G "{input.bins}/*.{params.ext}" >/dev/null 2>&1; then
                export CHECKM2DB={params.db}
                pixi r -e checkm2 checkm2 predict --threads {threads} -x {params.ext} --input {input.bins} --output-directory $(dirname {output.tsv}) --force
            else
                echo "WARNING: no bins in {input.bins}, writing an empty CheckM2 report" >&2
                mkdir -p $(dirname {output.tsv})
                echo -e "Name\\tCompleteness\\tContamination" >{output.tsv}
            fi
            date
        ) 2>{log.err} >{log.out}
        """
