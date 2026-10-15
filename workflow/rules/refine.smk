# Refine DasTool Bins with rosella


rule refine_all:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/Refine/rosella.done"
            ),
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


rule rosella_ref:
    input:
        cont=rules.filter_length.output.out,
        cov=expand(
            os.path.join(
                RESULTS_DIR,
                "intermediate/BamFiles/{{sample}}/{{sample}}_{sample_2}.bam",
            ),
            sample_2=SAMPLES,
        ),
        check=rules.checkm2.output,
        bins=rules.filter_dastool_bins.output,
    output:
        os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/Refine/rosella.done"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/refine/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/refine/{sample}.err.log"),
    threads: 30
    resources:
        runtime=60 * 48,
        mem_mb=get_mem_mb,
    params:
        ext="fa",
    message:
        "Refine bins with Rosella"
    shell:
        """
        (
            date
            rm -rf $(dirname {output})/*
            pixi r -e rosella rosella refine --min-bin-size 200000 --checkm-results {input.check} -d {input.bins} --max-contamination 15 -x {params.ext} -b {input.cov} -r {input.cont} --threads {threads} -o $(dirname {output})
            touch {output}
            date
        ) 2>{log.err} >{log.out}
        """


rule rename_refine:
    input:
        ref=rules.rosella_ref.output,
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
            mkdir -p {output} && for file in $(dirname {input}); do cp $file/*fa {output}; done
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
