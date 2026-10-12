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
                "Bins/intermediate/{sample}/DasTool/das_DASTool_summary.tsv",
            ),
            sample=ASSEMBLY,
        ),
    output:
        touch("status/binning_refinement.done"),


rule prepare_dasTool:
    input:
        metabat=rules.prok_metabat2.output,
        rosella=rules.prok_rosella.output,
        taxvamb=rules.prok_taxvamb.output,
        concoct=rules.prok_concoct.output,
        comebin=rules.comebin.output,
        semibin=rules.prok_semibin.output,
        metacat=rules.metacat_cluster.output.folder,
        metacat_done=rules.metacat_cluster.output.done,
    output:
        metabatout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/DasTool/metabat_das.tsv"
        ),
        rosellaout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/DasTool/rosella_das.tsv"
        ),
        taxvambout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/DasTool/taxvamb_das.tsv"
        ),
        concoctout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/DasTool/concoct_das.tsv"
        ),
        comebinout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/DasTool/comebin_das.tsv"
        ),
        semibinout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/DasTool/semibin_das.tsv"
        ),
        metacatout=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/DasTool/metacat_das.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/prepare_dasTool/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prepare_dasTool/{sample}.err.log"),
    localrule: True
    params:
        value="{sample}",
    message:
        "Preparing files for DasTool"
    shell:
        r"""
        (
            date
            pixi r -e dastool Fasta_to_Contig2Bin.sh -e fa -i {input.metabat} >{output.metabatout}
            perl -pe 's/^(\S+)[^\t]*(\t.*)$/$1$2/' {output.metabatout} >$(dirname {output.metabatout})/t.txt
            mv $(dirname {output.metabatout})/t.txt {output.metabatout}
            perl -pe 's/metabat./{params.value}_metabat_/g' {output.metabatout} >$(dirname {output.metabatout})/t.txt
            mv $(dirname {output.metabatout})/t.txt {output.metabatout}
            pixi r -e dastool Fasta_to_Contig2Bin.sh -e fna -i $(dirname {input.rosella}) >{output.rosellaout}
            perl -pe 's/rosella_(bin|refined)_/{params.value}_rosella_/g' {output.rosellaout} >$(dirname {output.metabatout})/t2.txt
            mv $(dirname {output.metabatout})/t2.txt {output.rosellaout}
            pixi r -e dastool Fasta_to_Contig2Bin.sh -e fna -i $(dirname {input.taxvamb})/bins >{output.taxvambout}
            perl -pe 's/\t/\t{params.value}_taxvamb_/g' {output.taxvambout} >$(dirname {output.metabatout})/t3.txt
            mv $(dirname {output.metabatout})/t3.txt {output.taxvambout}
            perl -pe 's/,/\t{params.value}_concoct_/g' {input.concoct} | tail -n +2 >{output.concoctout}
            perl -pe 's/\tgroup/\t{params.value}_comebin_/g' {input.comebin} >{output.comebinout}
            perl -pe 's/\t/\t{params.value}_semibin_/g' {input.semibin} | tail -n +2 >{output.semibinout}
            pixi r -e dastool Fasta_to_Contig2Bin.sh -e fasta -i {input.metacat} >{output.metacatout}
            perl -pe 's/_metacat\./_metacat_/g' {output.metacatout} >$(dirname {output.metabatout})/t4.txt
            mv $(dirname {output.metabatout})/t4.txt {output.metacatout}
            date
        ) 2>{log.err} >{log.out}
        """


rule dasTool_run:
    input:
        metabat=rules.prepare_dasTool.output.metabatout,
        taxvamb=rules.prepare_dasTool.output.taxvambout,
        rosella=rules.prepare_dasTool.output.rosellaout,
        concoct=rules.prepare_dasTool.output.concoctout,
        comebin=rules.prepare_dasTool.output.comebinout,
        semibin=rules.prepare_dasTool.output.semibinout,
        metacat=rules.prepare_dasTool.output.metacatout,
        cont=rules.filter_length.output.out,
        sel=os.path.join(RESULTS_DIR, "Bins/intermediate/selected_binners.txt"),
    output:
        os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/DasTool/das_DASTool_summary.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/dasTool/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/dasTool/{sample}.err.log"),
    threads: 10
    resources:
        runtime=60 * 12,
        mem_mb=get_mem_mb,
    message:
        "Running DasTool on the binners selected by binner_selection"
    shell:
        """
        (
            date
            mapfile -t B <{input.sel}
            if [[ ${{#B[@]}} -eq 0 ]]; then
                echo "ERROR: {input.sel} is empty - no binners selected for DAS Tool" >&2
                exit 1
            fi
            IFS=,
            LABELS="${{B[*]}}"
            unset IFS
            PATHS=$(printf "$(dirname {output})/%s_das.tsv," "${{B[@]}}")
            PATHS=${{PATHS%,}}
            pixi r -e dastool DAS_Tool -i "$PATHS" --score_threshold -42 -l "$LABELS" -c {input.cont} -o $(dirname {output})/das --write_bins --search_engine diamond --threads {threads}
            date
        ) 2>{log.err} >{log.out}
        """


rule filter_dastool_bins:
    # DAS Tool's --write_bins writes a FASTA for every candidate bin it
    # evaluated (including SCG-split "_sub" candidates and original bins
    # that lost every contig to a better-scoring competitor), not just its
    # final non-redundant selection. Only the bin names listed in
    # das_DASTool_summary.tsv are DAS Tool's actual answer - feeding the
    # raw das_DASTool_bins/ directory to CheckM2/Rosella refine (as the
    # unfiltered directory does) lets rejected bins that still score well
    # on their own (e.g. a near-complete bin that lost its contigs to a
    # different winning bin) pass every quality gate and reach Galah as a
    # bogus "distinct" genome that still shares raw contigs with the bin
    # that actually won them.
    input:
        summary=rules.dasTool_run.output,
    output:
        directory(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/DasTool/das_DASTool_bins_filtered",
            )
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/filter_dastool_bins/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/filter_dastool_bins/{sample}.err.log"),
    localrule: True
    params:
        ext="fa",
    message:
        "Keeping only DAS Tool's selected bins for {wildcards.sample}"
    shell:
        """
        (
            date
            mkdir -p {output}
            src=$(dirname {input.summary})/das_DASTool_bins
            missing=0
            while read -r bin; do
                if [[ -f "$src/$bin.{params.ext}" ]]; then
                    cp "$src/$bin.{params.ext}" {output}/
                else
                    echo "ERROR: expected winning bin file missing: $src/$bin.{params.ext}" >&2
                    missing=$((missing + 1))
                fi
            done < <(tail -n +2 {input.summary} | cut -f1)
            if [[ $missing -gt 0 ]]; then
                echo "ERROR: $missing winning bin(s) listed in das_DASTool_summary.tsv were missing from $src" >&2
                exit 1
            fi
            date
        ) 2>{log.err} >{log.out}
        """
