# Refine DasTool Bins with rosella

localrules: 

rule refine_all:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Refine/rosella_bins.json"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Refine_rename/genomes_list.txt"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/checkm2_refine/quality_report.tsv"))
    output:
        touch("status/refine.done")


rule rosella_ref:
    input:
        cont=rules.filter_length.output,
        cov=rules.metabat2_mapping.output,
        check=rules.checkm2.output,
        bins=rules.dasTool_run.output
    output:
        os.path.join(RESULTS_DIR, "Bins/{sample}/Refine/rosella_bins.json")
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
        rm -f $(dirname {output})/* && rosella refine --min-bin-size 200000 --checkm-file {input.check} -d $(dirname {input.bins})/das_DASTool_bins --max-contamination 15 -x {params.ext} -i {input.cov} -a {input.cont} --threads {threads} -o $(dirname {output}) && 
        date) 2> {log.err} > {log.out}
        """


rule rename_refine:
    input:
        ref=rules.rosella_ref.output
    params:
        sample="{sample}"
    output:
        list=os.path.join(RESULTS_DIR, "Bins/{sample}/Refine_rename/genomes_list.txt")
    log:
        out=os.path.join(RESULTS_DIR, "logs/refine_rename/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/refine_rename/{sample}.err.log")
    shell:
        "(date && "
        "rm -f {output.list} && "
        "for f in $(dirname {input.ref})/*fna; "
        "do "
        "name=$(basename $f); "
        "cp $f $(dirname {output})/{params.sample}_$name; "
        "echo $(dirname {output})/{params.sample}_$name >> {output.list}; " 
        "done && date) 2> {log.err} > {log.out}"

rule cp_mags:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Refine_rename/genomes_list.txt"), sample = SAMPLES)
    output:
        directory(os.path.join(RESULTS_DIR, "Bins/renamed_mags"))
    message:
        "Copy MAGs to another folder"
    shell:
        "(date && mkdir {output} && for file in $(dirname {input}); do cp $file/*fna {output}; done && date)"

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
        db=config["checkm2"]["db"],
        checkm2=os.path.join(SUBMODULES, "bin/checkm2"),
    message:
        "Running CheckM2 for directory renamed bins"
    shell:
        "(date && "
        "export CHECKM2DB={params.db} && "
        "checkm2 predict --threads {threads} -x fna --input {input} --output-directory $(dirname {output.tsv}) --force && "
        "date) 2> {log.err} > {log.out}"
