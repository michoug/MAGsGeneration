rule binning_refinement:
    input:
        expand(
            os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/das_DASTool_summary.tsv"),
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
    output:
        metabatout=os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/metabat_das.tsv"),
        rosellaout=os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/rosella_das.tsv"),
        taxvambout=os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/taxvamb_das.tsv"),
        concoctout=os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/concoct_das.tsv"),
        comebinout=os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/comebin_das.tsv"),
        semibinout=os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/semibin_das.tsv"),
    params:
        value="{sample}",
    conda:
        os.path.join(WORK_DIR, "envs/dastool.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/prepare_dasTool/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prepare_dasTool/{sample}.err.log"),
    message:
        "Preparing files for DasTool"
    localrule: True
    shell:
        """
        (date && Fasta_to_Contig2Bin.sh -e fa -i {input.metabat} > {output.metabatout}
        perl -pe 's/metabat./{params.value}_metabat_/g' {output.metabatout} > $(dirname {output.metabatout})/t.txt && mv $(dirname {output.metabatout})/t.txt {output.metabatout} 
        Fasta_to_Contig2Bin.sh -e fna -i $(dirname {input.rosella}) >  {output.rosellaout}
        perl -pe 's/rosella_(bin|refined)_/{params.value}_rosella_/g' {output.rosellaout} > $(dirname {output.metabatout})/t2.txt && mv $(dirname {output.metabatout})/t2.txt {output.rosellaout}
        Fasta_to_Contig2Bin.sh -e fna -i $(dirname {input.taxvamb})/bins > {output.taxvambout}
        perl -pe 's/\t/\t{params.value}_taxvamb_/g' {output.taxvambout} > $(dirname {output.metabatout})/t3.txt && mv $(dirname {output.metabatout})/t3.txt {output.taxvambout} 
        perl -pe 's/,/\t{params.value}_concoct_/g' {input.concoct} | tail -n +2 > {output.concoctout}
        perl -pe 's/\tgroup/\t{params.value}_comebin_/g' {input.comebin} > {output.comebinout}
        perl -pe 's/\t/\t{params.value}_semibin_/g' {input.semibin} | tail -n +2 > {output.semibinout}
        date) 2> {log.err} > {log.out}
        """


rule dasTool_run:
    input:
        metabat=rules.prepare_dasTool.output.metabatout,
        # metabinner=rules.prepare_dasTool.output.metabinnerout,
        taxvamb=rules.prepare_dasTool.output.taxvambout,
        rosella=rules.prepare_dasTool.output.rosellaout,
        concoct=rules.prepare_dasTool.output.concoctout,
        comebin=rules.prepare_dasTool.output.comebinout,
        semibin=rules.prepare_dasTool.output.semibinout,
        cont=rules.filter_length.output.out,
    output:
        os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/das_DASTool_summary.tsv"),
    conda:
        os.path.join(WORK_DIR, "envs/dastool.yaml")
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
        (date && 
        DAS_Tool -i {input.metabat},{input.taxvamb},{input.rosella},{input.concoct},{input.comebin},{input.semibin} --score_threshold -42 -l metabat,taxvamb,rosella,concoct,comebin,semibin -c {input.cont} -o $(dirname {output})/das --write_bins --search_engine diamond --threads {threads} && 
        date) 2> {log.err} > {log.out}        
        """
