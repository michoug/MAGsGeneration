rule dereplicate:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/finalBins")),
    output:
        touch("status/dereplicate.done"),


rule galah_prepare:
    input:
        check=rules.checkm2_refine.output.tsv,
    output:
        final=os.path.join(RESULTS_DIR, "Bins/intermediate/checkm2BeforeGalah.tsv"),
        temp=temp(os.path.join(RESULTS_DIR, "Bins/intermediate/checkm2_temp.tsv")),
    localrule: True
    message:
        "Adjusting CheckM output for input to Galah"
    shell:
        """
        cat {input.check} | sed '/^Name/d' | awk '{{print $1","$2","$3}}' > {output.temp} &&
        (echo "genome,completeness,contamination" && cat {output.temp}) > {output.final}
        """


## galah cluster --genome-fasta-list genomes_list.txt --genome-info t2.txt --min-completeness 70 --max-contamination 10 --output-representative-fasta-directory-copy final -t 10


# Dereplicating the genomes
checkpoint galah:
    input:
        bins=rules.cp_mags.output,
        check=rules.galah_prepare.output.final,
    output:
        final=directory(os.path.join(RESULTS_DIR, "Bins/finalBins")),
        cluster=os.path.join(RESULTS_DIR, "Bins/intermediate/cluster_galah.txt"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/galah/galah.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/galah/galah.err.log"),
    threads: 10
    resources:
        runtime=60 * 2,
        mem_mb=get_mem_mb,
    params:
        comp=70,
        cont=10,
    message:
        "Running Galah on all NOMIS MAGs"
    shell:
        "(date && "
        "pixi r -e galah galah cluster --genome-fasta-directory {input.bins} -x fa --genome-info {input.check} --min-completeness {params.comp} --max-contamination {params.cont} --output-representative-fasta-directory-copy {output.final} --output-cluster-definition {output.cluster} -t {threads} && "
        "date) 2> {log.err} > {log.out}"
