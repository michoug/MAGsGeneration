rule dastool:
    input:
       expand(os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/das_DASTool_summary.tsv"), sample = SAMPLES)
    output:
       touch("status/dastool.done")

rule prepare_dasTool:
    input:
        metabat=rules.metabat2.output,
        metabinner=rules.metabinner.output,
        rosella=rules.rosella.output
    output:
        metabatout=os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/metabat_das.tsv"),
        metabinnerout=os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/metabinner_das.tsv"),
        rosellaout=os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/rosella_das.tsv")
    params:
        value="{sample}"
    conda: 
        os.path.join(ENV_DIR, "dastool.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/prepare_dasTool/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prepare_dasTool/{sample}.err.log")
    message:
        "Preparing files for DasTool"
    shell:
        """
        (date && Fasta_to_Contig2Bin.sh -e fa -i {input.metabat} > {output.metabatout} && 
        perl -pe 's/metabat./{params.value}_metabat_/g' {output.metabatout} > t.txt && mv t.txt {output.metabatout} && 
        perl -pe 's/\t/\t{params.value}_metabinner_/g' {input.metabinner} > {output.metabinnerout} && 
        Fasta_to_Contig2Bin.sh -e fna -i $(dirname {input.rosella}) >  {output.rosellaout} && 
        perl -pe 's/rosella_bin./{params.value}_rosella_/g' {output.rosellaout} > t2.txt && mv t2.txt {output.rosellaout} && 
        date) 2> {log.err} > {log.out}
        """


rule dasTool_run:
    input:
        metabat=rules.prepare_dasTool.output.metabatout,
        metabinner=rules.prepare_dasTool.output.metabinnerout,
        rosella=rules.prepare_dasTool.output.rosellaout,
        cont=rules.filter_length.output
    output:
        os.path.join(RESULTS_DIR, "Bins/{sample}/DasTool/das_DASTool_summary.tsv")
    conda:
        os.path.join(ENV_DIR, "dastool.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/dasTool/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/dasTool/{sample}.err.log")
    threads:
        config["dasTool"]["threads"]
    message:
        "Running DasTool"
    shell:
        """
        (date && 
        DAS_Tool -i {input.metabat},{input.metabinner},{input.rosella} --score_threshold -42 -l metabat,metabinner,rosella -c {input.cont} -o $(dirname {output})/das --write_bins --search_engine diamond --threads {threads} && 
        date) 2> {log.err} > {log.out}        
        """
