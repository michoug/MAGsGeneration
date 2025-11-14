rule binning_refinement:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/{sample}/BinRef/final_bins_quality_reports.tsv"
            ),
            sample=ASSEMBLY,
        ),
    output:
        touch("status/binning_refinement.done"),


rule prepare_binning_refinement:
    input:
        metabat=rules.metabat2.output,
        # metabinner=rules.metabinner.output,
        rosella=rules.rosella.output,
        taxvamb=rules.taxvamb.output,
        concoct=rules.concoct.output,
        comebin=rules.comebin.output,
        semibin=rules.semibin.output,
    output:
        metabatout=os.path.join(RESULTS_DIR, "Bins/{sample}/BinRef/metabat_das.tsv"),
        # metabinnerout=os.path.join(
        #     RESULTS_DIR, "Bins/{sample}/BinRef/metabinner_das.tsv"
        # ),
        rosellaout=os.path.join(RESULTS_DIR, "Bins/{sample}/BinRef/rosella_das.tsv"),
        taxvambout=os.path.join(RESULTS_DIR, "Bins/{sample}/BinRef/taxvamb_das.tsv"),
        concoctout=os.path.join(RESULTS_DIR, "Bins/{sample}/BinRef/concoct_das.tsv"),
        comebinout=os.path.join(RESULTS_DIR, "Bins/{sample}/BinRef/comebin_das.tsv"),
        semibinout=os.path.join(RESULTS_DIR, "Bins/{sample}/BinRef/semibin_das.tsv"),
    params:
        value="{sample}",
        script=os.path.join(WORK_DIR, "scripts/Fasta_to_Scaffolds2Bin.sh"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/prepare_BinRef/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prepare_BinRef/{sample}.err.log"),
    message:
        "Preparing files for BinRef"
    localrule: True
    shell:
        """
        (date && bash {params.script} -e fa -i {input.metabat} > {output.metabatout}
        perl -pe 's/metabat./{params.value}_metabat_/g' {output.metabatout} > $(dirname {output.metabatout})/t.txt && mv $(dirname {output.metabatout})/t.txt {output.metabatout} 
        bash {params.script} -e fna -i $(dirname {input.rosella}) >  {output.rosellaout}
        perl -pe 's/rosella_(bin|refined)_/{params.value}_rosella_/g' {output.rosellaout} > $(dirname {output.metabatout})/t2.txt && mv $(dirname {output.metabatout})/t2.txt {output.rosellaout}
        bash {params.script} -e fna -i $(dirname {input.taxvamb})/bins > {output.taxvambout}
        perl -pe 's/\t/\t{params.value}_taxvamb_/g' {output.taxvambout} > $(dirname {output.metabatout})/t3.txt && mv $(dirname {output.metabatout})/t3.txt {output.taxvambout} 
        perl -pe 's/,/\t{params.value}_concoct_/g' {input.concoct} | tail -n +2 > {output.concoctout}
        perl -pe 's/\tgroup/\t{params.value}_comebin_/g' {input.comebin} > {output.comebinout}
        perl -pe 's/\t/\t{params.value}_semibin_/g' {input.semibin} | tail -n +2 > {output.semibinout}
        date) 2> {log.err} > {log.out}
        """
        # perl -pe 's/\t/\t{params.value}_metabinner_/g' {input.metabinner} > {output.metabinnerout}


rule binref:
    input:
        metabat=rules.prepare_binning_refinement.output.metabatout,
        # metabinner=rules.prepare_binning_refinement.output.metabinnerout,
        taxvamb=rules.prepare_binning_refinement.output.taxvambout,
        rosella=rules.prepare_binning_refinement.output.rosellaout,
        concoct=rules.prepare_binning_refinement.output.concoctout,
        comebin=rules.prepare_binning_refinement.output.comebinout,
        semibin=rules.prepare_binning_refinement.output.semibinout,
        cont=rules.filter_length.output.out,
    output:
        os.path.join(RESULTS_DIR, "Bins/{sample}/BinRef/final_bins_quality_reports.tsv"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/binref/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/binref/{sample}.err.log"),
    threads: 20
    resources:
        runtime=60 * 12,
        mem_mb=get_mem_mb,
    message:
        "Running BinRef"
    params:
        db=config["checkm2"]["db"],
    shell:
        """
        (date && 
        pixi r -e binette binette -b {input.metabat} {input.taxvamb} {input.rosella} {input.concoct} {input.comebin} {input.semibin} -c {input.cont} -o $(dirname {output}) -t {threads} --checkm2_db {params.db} --min_completeness 50 --max_length 20000000 --verbose && 
        date) 2> {log.err} > {log.out}        
        """


#
