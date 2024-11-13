# Refine DasTool Bins with rosella

rule refine_all:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Refine/rosella.done"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Refine_rename/genomes_list.txt"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/checkm2_refine/quality_report.tsv"))
    output:
        touch("status/refine.done")


rule rosella_ref:
    input:
        cont=rules.filter_length.output.out,
        cov=expand(os.path.join(RESULTS_DIR,"BamFiles/{{sample}}/{{sample}}_{sample_2}.bam"), sample_2 = SAMPLES),
        check=rules.checkm2.output,
        bins=rules.dasTool_run.output
    output:
        os.path.join(RESULTS_DIR, "Bins/{sample}/Refine/rosella.done")
    conda:
        os.path.join(ENV_DIR, "rosella.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/refine/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/refine/{sample}.err.log")
    threads:
        config["rosella"]["threads"]
    message:
        "Refine bins with Rosella"
    params:
        ext=config["checkm2"]["extension"]
    shell:
        """
        (date && 
        rm -rf $(dirname {output})/* && rosella refine --min-bin-size 200000 --checkm-results {input.check} -d $(dirname {input.bins})/das_DASTool_bins --max-contamination 15 -x {params.ext} -b {input.cov} -r {input.cont} --threads {threads} -o $(dirname {output}) && 
        touch {output} &&
        date) 2> {log.err} > {log.out}
        """


rule rename_refine:
    input:
        ref=rules.rosella_ref.output
    params:
        sample="{sample}"
    output:
        list=os.path.join(RESULTS_DIR, "Bins/{sample}/Refine_rename/genomes_list.txt"),
        list_temp=temp(os.path.join(RESULTS_DIR, "Bins/{sample}/Refine_rename/genomes_list_temp.txt")),
        bash_temp=temp(os.path.join(RESULTS_DIR, "Bins/{sample}/Refine_rename/t.sh"))
    log:
        out=os.path.join(RESULTS_DIR, "logs/refine_rename/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/refine_rename/{sample}.err.log")
    localrule: True
    shell:
        """
        (date 
        rm -f $(dirname {output.list})/*
        for f in $(dirname {input.ref})/*bins/*fa; 
        do name=$(basename $f);
        echo cp $f $(dirname {output.list})/{params.sample}_$name >> {output.list_temp};
        done && 
        perl -pe 's/{params.sample}_{params.sample}/{params.sample}/g' {output.list_temp} | perl -pe 's/.*rosella_unbinned.*//g'> {output.bash_temp}
        bash {output.bash_temp}
        ls $(dirname {output.list})/*fa >> {output.list}; 
        date) 2> {log.err} > {log.out}
        """

rule cp_mags:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Refine_rename/genomes_list.txt"), sample = SAMPLES)
    output:
        directory(os.path.join(RESULTS_DIR, "Bins/renamed_mags"))
    message:
        "Copy MAGs to another folder"
    localrule: True
    shell:
        "(date && mkdir {output} && for file in $(dirname {input}); do cp $file/*fa {output}; done && date)"

# Checking bin quality
rule checkm2_refine:
    input:
        os.path.join(RESULTS_DIR, "Bins/renamed_mags")
    output:
        tsv=os.path.join(RESULTS_DIR, "Bins/checkm2_refine/quality_report.tsv")
    conda:
        os.path.join(ENV_DIR, "checkm2.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/checkm_refine.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/checkm_refine.err.log")
    threads:
        config["checkm2"]["threads"]
    params:
        ext=config["checkm2"]["extension"],
        db=config["checkm2"]["db"]
    message:
        "Running CheckM2 for directory renamed bins"
    shell:
        "(date && "
        "export CHECKM2DB={params.db} && "
        "checkm2 predict --threads {threads} -x fa --input {input} --output-directory $(dirname {output.tsv}) --force && "
        "date) 2> {log.err} > {log.out}"
