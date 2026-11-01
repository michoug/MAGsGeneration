# Filter Binette bins with GUNC (chimerism) and contamination


rule refine_all:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/Refine/gunc.done"),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Refine_rename/genomes_list.txt",
            ),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/checkm2_refine/quality_report.tsv"
            )
        ),
    output:
        touch("status/refine.done"),


rule gunc_run:
    # GUNC chimerism check on the Binette bins
    input:
        bins=rules.binette_run.output.bins,
        db=GUNC_DB,
    output:
        os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Refine/gunc/GUNC.maxCSS_level.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/gunc/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/gunc/{sample}.err.log"),
    threads: 10
    resources:
        runtime=60 * 4,
        mem_mb=get_mem_mb,
    params:
        ext=".fa",
    message:
        "Running GUNC on the Binette bins of {wildcards.sample}"
    shell:
        """
        (
            date
            outdir=$(dirname {output})
            rm -rf $outdir
            mkdir -p $outdir/tmp
            if compgen -G "{input.bins}/*{params.ext}" >/dev/null; then
                pixi r -e gunc gunc run --input_dir {input.bins} --file_suffix {params.ext} --db_file {input.db} --threads {threads} --out_dir $outdir --temp_dir $outdir/tmp
                cp $outdir/GUNC.*.maxCSS_level.tsv {output}
            else
                echo "WARNING: no bins in {input.bins}, writing an empty GUNC table" >&2
                printf 'genome\tpass.GUNC\n' >{output}
            fi
            rm -rf $outdir/tmp
            date
        ) 2>{log.err} >{log.out}
        """


rule gunc_filter:
    # Keeps bins with pass.GUNC not False and contamination <= gunc.max_contamination
    input:
        bins=rules.binette_run.output.bins,
        quality=rules.checkm2.output.tsv,
        gunc=rules.gunc_run.output,
    output:
        bins=directory(
            os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/Refine/gunc_bins")
        ),
        table=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Refine/gunc_filter.tsv"
        ),
        done=os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/Refine/gunc.done"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/gunc_filter/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/gunc_filter/{sample}.err.log"),
    localrule: True
    params:
        script=os.path.join(WORK_DIR, "scripts/gunc_filter_bins.py"),
        max_contamination=GUNC_CFG.get("max_contamination", 10),
    message:
        "Filtering the Binette bins of {wildcards.sample} on GUNC and contamination"
    shell:
        """
        (
            date
            refine=$(dirname {output.done})
            # leftovers of the former Rosella refine step would be picked up by rename_refine
            rm -rf $refine/refined_bins $refine/unchanged_bins $refine/coverage.tsv $refine/kmer_frequencies.tsv $refine/rosella.done
            pixi r -e gunc python {params.script} --bins {input.bins} --quality {input.quality} --gunc {input.gunc} \
                --outdir {output.bins} --table {output.table} --max-contamination {params.max_contamination}
            touch {output.done}
            date
        ) 2>{log.err} >{log.out}
        """


rule rename_refine:
    input:
        ref=rules.gunc_filter.output.done,
    output:
        list=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Refine_rename/genomes_list.txt"
        ),
        list_temp=temp(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Refine_rename/genomes_list_temp.txt",
            )
        ),
        bash_temp=temp(
            os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/Refine_rename/t.sh")
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/refine_rename/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/refine_rename/{sample}.err.log"),
    localrule: True
    params:
        sample="{sample}",
    shell:
        """
        (
            date
            rm -f $(dirname {output.list})/*
            : >{output.list_temp}

            shopt -s nullglob
            for f in $(dirname {input.ref})/*bins/*fa; do
                name=$(basename $f)
                echo cp $f $(dirname {output.list})/{params.sample}_$name >>{output.list_temp}
            done
            shopt -u nullglob

            if [ -s {output.list_temp} ]; then
                perl -pe 's/{params.sample}_{params.sample}/{params.sample}/g' {output.list_temp} \
                    | perl -pe 's/.*rosella_unbinned.*//g' >{output.bash_temp}
                bash {output.bash_temp}
            else
                echo "WARNING: no bins found for sample {params.sample}" >&2
                : >{output.bash_temp}
            fi

            ls $(dirname {output.list})/*fa >>{output.list} 2>/dev/null \
                || echo "WARNING: sample {params.sample} produced 0 bins, writing empty list" >&2

            touch {output.list}
            date
        ) 2>{log.err} >{log.out}
        """


rule cp_mags:
    input:
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Refine_rename/genomes_list.txt",
            ),
            sample=ASSEMBLY,
        ),
    output:
        directory(os.path.join(RESULTS_DIR, "Bins/intermediate/renamed_mags")),
    localrule: True
    message:
        "Copy MAGs to another folder"
    shell:
        """
        (
            date
            shopt -s nullglob
            mkdir -p {output}
            for dir in $(dirname {input}); do
                for fa in "$dir"/*fa; do
                    cp "$fa" {output}
                done
            done
            date
        )
        """


# Checking bin quality
rule checkm2_refine:
    input:
        os.path.join(RESULTS_DIR, "Bins/intermediate/renamed_mags"),
    output:
        tsv=os.path.join(
            RESULTS_DIR, "Bins/intermediate/checkm2_refine/quality_report.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/checkm_refine.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/checkm_refine.err.log"),
    threads: 20
    resources:
        runtime=60 * 24,
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
            export CHECKM2DB={params.db}
            pixi r -e checkm2 checkm2 predict --threads {threads} -x {params.ext} --input {input} --output-directory $(dirname {output.tsv}) --force
            date
        ) 2>{log.err} >{log.out}
        """
