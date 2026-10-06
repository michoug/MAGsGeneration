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
        metabat=rules.metabat2.output,
        rosella=rules.rosella.output,
        taxvamb=rules.taxvamb.output,
        concoct=rules.concoct.output,
        comebin=rules.comebin.output,
        semibin=rules.semibin.output,
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
        "Running DasTool"
    shell:
        """
        (
            date
            pixi r -e dastool DAS_Tool -i {input.metabat},{input.taxvamb},{input.rosella},{input.concoct},{input.comebin},{input.semibin},{input.metacat} --score_threshold -42 -l metabat,taxvamb,rosella,concoct,comebin,semibin,metacat -c {input.cont} -o $(dirname {output})/das --write_bins --search_engine diamond --threads {threads}
            date
        ) 2>{log.err} >{log.out}
        """
