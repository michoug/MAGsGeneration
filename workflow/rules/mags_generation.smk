
rule mags_generation:
    input:
        expand(os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "TempReads/{sample_2}.paired.bcsp"), sample_2 = SAMPLES_2),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Metabat"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Metabinner/metabinner_result.tsv"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Rosella/UMAP_projection_with_clusters_final.png"), sample = SAMPLES)
    output:
        touch("status/mags_generation.done")

rule filter_length:
    input:
        os.path.join(RESULTS_DIR, "Assembly/{sample}.fa")
    output:
        out = os.path.join(RESULTS_DIR,"Assembly/{sample}_filter.fasta"),
        other = temp(os.path.join(RESULTS_DIR,"Assembly/{sample}_filter_name.fasta"))
    conda:
        os.path.join(ENV_DIR, "seqkit.yaml")
    threads:
        config["filter_length"]["threads"]
    message:
        "Remove contigs less than 1.5kb"
    shell:
        "seqkit seq -j {threads} --remove-gaps -o {output.other} -m 1500 {input} && perl -pe \"s/ .*//g\" {output.other} > {output.out}"

rule mapping_index:
    input:
        readsF=expand(os.path.join(RESULTS_DIR, "TrimFinal/{sample_2}_1.fq.gz"),sample_2 = SAMPLES_2),
        readsR=expand(os.path.join(RESULTS_DIR, "TrimFinal/{sample_2}_2.fq.gz"),sample_2 = SAMPLES_2),
        files=os.path.join(WORK_DIR,"config/sample_list.txt")
    output:
        folder=directory(os.path.join(RESULTS_DIR, "TempReads/")),
        bcsp=expand(os.path.join(RESULTS_DIR, "TempReads/{sample_2}.paired.bcsp"),sample_2 = SAMPLES_2)
    log:
        out=os.path.join(RESULTS_DIR, "logs/fairy.index.err.log"),
        err=os.path.join(RESULTS_DIR, "logs/fairy.index.out.log")
    threads:
        config["fairy"]["threads"]
    conda:
        os.path.join(ENV_DIR, "fairy.yaml")
    message:
        "Mapping: Fairy index for assembly mapping"
    shell:
        """
        (date 
        fairy sketch -1 {input.readsF} -2 {input.readsR} -d {output.folder} -t {threads} --lS {input.files}
        date) 2> {log.err} > {log.out}
        """

rule fairy_mapping:
    input:
        bcsp=rules.mapping_index.output.bcsp,
        cont=rules.filter_length.output.out
    output:
        os.path.join(RESULTS_DIR,"Bins/{sample}/{sample}_cov.txt")
    conda:
        os.path.join(ENV_DIR, "fairy.yaml")
    threads:
        config["fairy"]["threads"]
    log:
        out=os.path.join(RESULTS_DIR, "logs/fairy_mapping/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/fairy_mapping/{sample}.err.log")
    message:
        "Run Fairy to estimate the coverage"
    shell:
        "(date && fairy coverage {input.bcsp} {input.cont} -t {threads} -o {output} && "
        "date) 2> {log.err} > {log.out}"


rule metabat2:
    priority: 20
    input:
        contig=rules.filter_length.output.out,
        cov=rules.fairy_mapping.output
    output:
        directory(os.path.join(RESULTS_DIR,"Bins/{sample}/Metabat/"))
    conda:
        os.path.join(ENV_DIR, "metabat.yaml")
    threads:
        config["metabat2"]["threads"]
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabat2/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabat2/{sample}.err.log")
    message:
        "Running metabat2 to obtain bins"
    shell:
        "(date && metabat2 -i {input.contig} -a {input.cov} -o {output}/metabat -t {threads} && "
        "date) 2> {log.err} > {log.out}"


rule metabinner_install:
    output:
        dbs=directory(os.path.join(DB_DIR, "MetaBinner"))
    conda:
        os.path.join(ENV_DIR, "metabinner.yaml")
    shell:
        "git clone https://github.com/ziyewang/MetaBinner.git &> /dev/null && "
        "mv MetaBinner $(dirname {output.dbs})"


rule metabinner_coverage:
    input:
        rules.fairy_mapping.output
    output:
        os.path.join(RESULTS_DIR,"Bins/{sample}/{sample}_metabinner_cov.txt")
    threads:
        config["metabinner"]["threads"]
    params:
        config["metabinner"]["length"]
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabinner_mapping/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabinner_mapping/{sample}.err.log")
    message:
        "Running the first part of metabinner to estimate the coverage"
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep = '\t')
        df_filtered = df.loc[:, ~df.columns.str.endswith('var')]
        df_filtered = df_filtered[df_filtered['contigLen'] > params[0]]
        df_filtered = df_filtered.drop(columns=['contigLen', 'totalAvgDepth'])
        df_filtered.to_csv(output[0], sep='\t', index=False)

rule metabinner_prepare:
    input:
        contig=rules.filter_length.output.out,
        bin=rules.metabinner_install.output.dbs
    output:
        cont_t=temp(RESULTS_DIR + "/Assembly/{sample}_filter_" + str(config["metabinner"]["length"]) + ".fa"),       
        cont=RESULTS_DIR + "/Bins/{sample}/{sample}_metabinner_" + str(config["metabinner"]["length"]) +".fa",        
        kmer_t=temp(RESULTS_DIR + "/Assembly/{sample}_filter_kmer_4_f" + str(config["metabinner"]["length"]) + ".csv"),
        kmer=RESULTS_DIR + "/Bins/{sample}/{sample}_kmer_4_f" + str(config["metabinner"]["length"]) + ".csv"
    conda:
        os.path.join(ENV_DIR, "metabinner.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabinner_prepare/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabinner_prepare/{sample}.err.log")
    message:
        "Prepare fasta and kmer file for Metabinner"
    shell:
        "(date && python {input.bin}/scripts/Filter_tooshort.py {input.contig} {config[metabinner][length]} && "
        "cp {output.cont_t} {output.cont} && " 
        "python {input.bin}/scripts/gen_kmer.py {input.contig} {config[metabinner][length]} 4 && "
        "cp {output.kmer_t} {output.kmer} && date) 2> {log.err} > {log.out}"


rule metabinner:
    priority: 20
    input:
        cont=rules.metabinner_prepare.output.cont,
        cov=rules.metabinner_coverage.output,
        kmer=rules.metabinner_prepare.output.kmer,    
        path=rules.metabinner_install.output.dbs
    output:
        os.path.join(RESULTS_DIR, "Bins/{sample}/Metabinner/metabinner_result.tsv")
    conda:
        os.path.join(ENV_DIR, "metabinner.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabinner/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabinner/{sample}.err.log")
    threads:
        config["metabinner"]["threads"]
    message:
        "Running Metabinner"
    shell:
        "(date && run_metabinner.sh -a {input.cont} -d {input.cov} -k {input.kmer} -p {input.path} -o $(dirname {output}) -t {threads} -s huge && "
        "mv $(dirname {output})/metabinner_res/* $(dirname {output}) && "
        "rm -r $(dirname {output})/metabinner_res/ $(dirname {output})/intermediate_result $(dirname {output})/unitem_profile $(dirname {output})/ensemble_res "
        "&& date) 2> {log.err} > {log.out}"

rule rosella:
    input:
        cont=rules.filter_length.output.out,
        cov=rules.fairy_mapping.output
    output:
        os.path.join(RESULTS_DIR, "Bins/{sample}/Rosella/UMAP_projection_with_clusters_final.png")
    conda:
        os.path.join(ENV_DIR, "rosella.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/rosella/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/rosella/{sample}.err.log")
    threads:
        config["rosella"]["threads"]
    message:
        "Running Rosella"
    shell:
        "(date && rm -r $(dirname {output}) && "
        "rosella recover -C {input.cov} -r {input.cont} --threads {threads} -o $(dirname {output}) "
        "&& date) 2> {log.err} > {log.out}"
