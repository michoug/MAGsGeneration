if "binning_refinement" in STEPS and "binner_selection" not in STEPS:
    raise Exception(
        "'binning_refinement' requires 'binner_selection' earlier in config "
        "'steps' - it provides Bins/intermediate/selected_binners.txt, even "
        "when binners.selection.enabled is false."
    )


rule binning_refinement:
    input:
        expand(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/Binette/final_bins_quality_reports.tsv",
            ),
            sample=ASSEMBLY,
        ),
    output:
        touch("status/binning_refinement.done"),


rule prepare_Binette:
    # One contig<TAB>bin table per binner
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
        metabatout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Binette/metabat_binette.tsv"
        ),
        rosellaout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Binette/rosella_binette.tsv"
        ),
        taxvambout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Binette/taxvamb_binette.tsv"
        ),
        concoctout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Binette/concoct_binette.tsv"
        ),
        comebinout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Binette/comebin_binette.tsv"
        ),
        semibinout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Binette/semibin_binette.tsv"
        ),
        metacatout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/Binette/metacat_binette.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/prepare_Binette/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prepare_Binette/{sample}.err.log"),
    localrule: True
    params:
        value="{sample}",
        script=os.path.join(WORK_DIR, "scripts/Fasta_to_Scaffolds2Bin.sh"),
    message:
        "Preparing per-binner contig2bin tables for {wildcards.sample}"
    shell:
        r"""
        (
            date
            bash {params.script} -e fa -i {input.metabat} >{output.metabatout}
            perl -pe 's/^(\S+)[^\t]*(\t.*)$/$1$2/' {output.metabatout} >$(dirname {output.metabatout})/t.txt
            mv $(dirname {output.metabatout})/t.txt {output.metabatout}
            perl -pe 's/metabat./{params.value}_metabat_/g' {output.metabatout} >$(dirname {output.metabatout})/t.txt
            mv $(dirname {output.metabatout})/t.txt {output.metabatout}
            bash {params.script} -e fna -i $(dirname {input.rosella}) >{output.rosellaout}
            perl -pe 's/rosella_(bin|refined)_/{params.value}_rosella_/g' {output.rosellaout} >$(dirname {output.metabatout})/t2.txt
            mv $(dirname {output.metabatout})/t2.txt {output.rosellaout}
            bash {params.script} -e fna -i $(dirname {input.taxvamb})/bins >{output.taxvambout}
            perl -pe 's/\t/\t{params.value}_taxvamb_/g' {output.taxvambout} >$(dirname {output.metabatout})/t3.txt
            mv $(dirname {output.metabatout})/t3.txt {output.taxvambout}
            perl -pe 's/,/\t{params.value}_concoct_/g' {input.concoct} | tail -n +2 >{output.concoctout}
            perl -pe 's/\tgroup/\t{params.value}_comebin_/g' {input.comebin} >{output.comebinout}
            perl -pe 's/\t/\t{params.value}_semibin_/g' {input.semibin} | tail -n +2 >{output.semibinout}
            bash {params.script} -e fasta -i {input.metacat} >{output.metacatout}
            perl -pe 's/_metacat\./_metacat_/g' {output.metacatout} >$(dirname {output.metabatout})/t4.txt
            mv $(dirname {output.metabatout})/t4.txt {output.metacatout}
            date
        ) 2>{log.err} >{log.out}
        """


rule binette_run:
    input:
        metabat=rules.prepare_Binette.output.metabatout,
        taxvamb=rules.prepare_Binette.output.taxvambout,
        rosella=rules.prepare_Binette.output.rosellaout,
        concoct=rules.prepare_Binette.output.concoctout,
        comebin=rules.prepare_Binette.output.comebinout,
        semibin=rules.prepare_Binette.output.semibinout,
        metacat=rules.prepare_Binette.output.metacatout,
        cont=rules.filter_length.output.out,
        sel=os.path.join(RESULTS_DIR, "Bins/intermediate/selected_binners.txt"),
    output:
        quality=os.path.join(
            RESULTS_DIR,
            "Bins/intermediate/{sample}/Binette/final_bins_quality_reports.tsv",
        ),
        bins=directory(
            os.path.join(RESULTS_DIR, "Bins/intermediate/{sample}/Binette/final_bins")
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/binette/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/binette/{sample}.err.log"),
    threads: 10
    resources:
        runtime=60 * 12,
        mem_mb=get_mem_mb,
    params:
        db=config["checkm2"]["db"],
    message:
        "Running Binette on the binners selected by binner_selection"
    shell:
        """
        (
            date
            declare -A C2B=(
                [metabat]={input.metabat}
                [taxvamb]={input.taxvamb}
                [rosella]={input.rosella}
                [concoct]={input.concoct}
                [comebin]={input.comebin}
                [semibin]={input.semibin}
                [metacat]={input.metacat}
            )
            mapfile -t B <{input.sel}
            if [[ ${{#B[@]}} -eq 0 ]]; then
                echo "ERROR: {input.sel} is empty - no binners selected for Binette" >&2
                exit 1
            fi
            PATHS=()
            for binner in "${{B[@]}}"; do
                if [[ -z "${{C2B[$binner]+x}}" ]]; then
                    echo "ERROR: unknown binner '$binner' listed in {input.sel}" >&2
                    exit 1
                fi
                PATHS+=("${{C2B[$binner]}}")
            done
            rm -rf $(dirname {output.quality})
            pixi r -e binette binette -b "${{PATHS[@]}}" -c {input.cont} -o $(dirname {output.quality}) -t {threads} --checkm2_db {params.db} --min_completeness 70 --max_length 50000000 --verbose
            date
        ) 2>{log.err} >{log.out}
        """
