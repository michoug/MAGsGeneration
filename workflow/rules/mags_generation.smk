
rule mags_generation:
    input:
        expand(os.path.join(RESULTS_DIR, "Assembly/{sample}_filter.fasta"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "TempReads/{sample}.paired.bcsp"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Metabat"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Metabinner/metabinner_result.tsv"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Bins/{sample}/Rosella/UMAP_projection_with_clusters_final.png"), sample = SAMPLES),
        expand(os.path.join(RESULTS_DIR,"Bins/{sample}/Vamb/Vamb_Results/vaevae_clusters_unsplit.tsv"), sample = SAMPLES)
    output:
        touch("status/mags_generation.done")


rule filter_length:
    input:
        rules.assembly_megahit.output.final
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
        readsF=expand(os.path.join(RESULTS_DIR, "TrimFinal/{sample}_1.fq.gz"),sample = SAMPLES),
        readsR=expand(os.path.join(RESULTS_DIR, "TrimFinal/{sample}_2.fq.gz"),sample = SAMPLES),
        files=os.path.join(WORK_DIR,"config/sample_list.txt")
    output:
        folder=directory(os.path.join(RESULTS_DIR, "TempReads/")),
        bcsp=expand(os.path.join(RESULTS_DIR, "TempReads/{sample}.paired.bcsp"),sample = SAMPLES)
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
        dbs=directory(os.path.join(RESULTS_DIR, "DB/MetaBinner")),
        file = os.path.join(RESULTS_DIR, "DB/MetaBinner/metabinner.done")
    localrule: True
    conda:
        os.path.join(ENV_DIR, "metabinner.yaml")
    shell:
        """
        rm -rf MetaBinner
        git clone https://github.com/ziyewang/MetaBinner.git &> /dev/null
        mv MetaBinner $(dirname {output.dbs})
        touch {output.file}
        """


rule metabinner_coverage:
    input:
        rules.fairy_mapping.output
    output:
        os.path.join(RESULTS_DIR,"Bins/{sample}/{sample}_metabinner_cov.txt")
    localrule: True
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
        bin=rules.metabinner_install.output.file
    output:
        cont_t=temp(RESULTS_DIR + "/Assembly/{sample}_filter_" + str(config["metabinner"]["length"]) + ".fa"),       
        cont=RESULTS_DIR + "/Bins/{sample}/{sample}_metabinner_" + str(config["metabinner"]["length"]) +".fa",        
        kmer_t=temp(RESULTS_DIR + "/Assembly/{sample}_filter_kmer_4_f" + str(config["metabinner"]["length"]) + ".csv"),
        kmer=RESULTS_DIR + "/Bins/{sample}/{sample}_kmer_4_f" + str(config["metabinner"]["length"]) + ".csv"
    localrule: True
    conda:
        os.path.join(ENV_DIR, "metabinner.yaml")
    params:
        config["metabinner"]["length"]
    log:
        out=os.path.join(RESULTS_DIR, "logs/metabinner_prepare/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/metabinner_prepare/{sample}.err.log")
    message:
        "Prepare fasta and kmer file for Metabinner"
    shell:
        "(date && python $(dirname {input.bin})/scripts/Filter_tooshort.py {input.contig} {params} && "
        "cp {output.cont_t} {output.cont} && " 
        "python $(dirname {input.bin})/scripts/gen_kmer.py {input.contig} {params} 4 && "
        "cp {output.kmer_t} {output.kmer} && date) 2> {log.err} > {log.out}"


rule metabinner:
    priority: 20
    input:
        cont=rules.metabinner_prepare.output.cont,
        cov=rules.metabinner_coverage.output,
        kmer=rules.metabinner_prepare.output.kmer,    
        path=rules.metabinner_install.output.file
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
        "(date && run_metabinner.sh -a {input.cont} -d {input.cov} -k {input.kmer} -p $(dirname {input.path}) -o $(dirname {output}) -t {threads} -s huge && "
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

rule vamb_install:
    output:
        dbs=directory(os.path.join(RESULTS_DIR, "vamb")),
        testfile=os.path.join(RESULTS_DIR, "vamb/install.done")
    conda:
        os.path.join(ENV_DIR, "vamb.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/vamb_install.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/vamb_install.err.log")
    localrule: True
    shell:
        "(date && rm -rf vamb && "
        "rm -rf {output.dbs} && "
        "git clone https://github.com/RasmussenLab/vamb -b master && "
        "mv vamb $(dirname {output.dbs}) && "
        "pip install -e {output.dbs} && "
        "touch {output.testfile} && "
        "date) 2> {log.err} > {log.out}"

rule vamb_coverage:
    input:
        rules.metabinner_coverage.output
    output:
        os.path.join(RESULTS_DIR,"Bins/{sample}/Vamb/{sample}_aemb_cov.txt")
    params:
        name = "{sample}"
    localrule: True
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep = '\t')
        df.rename(columns={"contigName": "contigname"}, inplace=True)
        df.to_csv(output[0], sep='\t', index=False, header=True)

rule vamb_length:
    input:
        rules.filter_length.output.out
    output:
        temp(os.path.join(RESULTS_DIR,"Assembly/{sample}_filter_vamb.fasta"))
    conda:
        os.path.join(ENV_DIR, "seqkit.yaml")
    threads:
        config["filter_length"]["threads"]
    localrule: True
    message:
        "Remove contigs less than 1501 bp for Vamb"
    shell:
        "seqkit seq -j {threads} --remove-gaps -o {output} -m 1501 {input}"

rule vamb_filter:
    input:
        rules.taxconverter.output,
        rules.vamb_length.output
    output:
        temp(os.path.join(RESULTS_DIR,"Assembly/{sample}_tax.tsv"))
    log:
        out=os.path.join(RESULTS_DIR, "logs/vamb_filter/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/vamb_filter/{sample}.err.log")
    localrule: True
    run:
        import pandas as pd
        ids = []
        with open(input[1], 'r') as file:
            for line in file:
                if line.startswith('>'):
                    ids.append(line[1:].strip())
        df = pd.read_csv(input[0], sep = '\t')
        filtered_df = df[df.iloc[:, 0].isin(ids)]
        filtered_df.to_csv(output[0], sep='\t', index=False, header=True)

rule taxvamb:
    input:
        file=rules.vamb_install.output.testfile,
        contig=rules.vamb_length.output,
        tax=rules.vamb_filter.output,
        cov = rules.vamb_coverage.output
    output:
        os.path.join(RESULTS_DIR,"Bins/{sample}/Vamb/Vamb_Results/vaevae_clusters_unsplit.tsv")
    conda:
        os.path.join(ENV_DIR, "vamb.yaml")
    threads:
        config["vamb"]["threads"]
    log:
        out=os.path.join(RESULTS_DIR, "logs/vamb/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/vamb/{sample}.err.log")
    shell:
        "(date && "
        "rm -r $(dirname {output}) && "
        "vamb bin taxvamb --outdir $(dirname {output}) -m 1501 --fasta {input.contig} --taxonomy {input.tax} --minfasta 500000 --abundance_tsv {input.cov} -p {threads} -o && "
        "date) 2> {log.err} > {log.out}"
