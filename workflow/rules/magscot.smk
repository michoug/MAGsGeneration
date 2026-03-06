rule binning_refinement:
    input:
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/MAGScoT/magscot_input.tsv"
            ),
            sample=ASSEMBLY,
        ),
        expand(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/MAGScoT/magscot_hmm.tsv"
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
        metacat=rules.metacat_cluster.output.done,
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
    params:
        value="{sample}",
    log:
        out=os.path.join(RESULTS_DIR, "logs/prepare_dasTool/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prepare_dasTool/{sample}.err.log"),
    message:
        "Preparing files for DasTool"
    localrule: True
    shell:
        """
        (date &&
        pixi r -e dastool Fasta_to_Contig2Bin.sh -e fa -i {input.metabat} > {output.metabatout}
        perl -pe 's/metabat./{params.value}_metabat_/g' {output.metabatout} > $(dirname {output.metabatout})/t.txt && mv $(dirname {output.metabatout})/t.txt {output.metabatout} 
        pixi r -e dastool Fasta_to_Contig2Bin.sh -e fna -i $(dirname {input.rosella}) >  {output.rosellaout}
        perl -pe 's/rosella_(bin|refined)_/{params.value}_rosella_/g' {output.rosellaout} > $(dirname {output.metabatout})/t2.txt && mv $(dirname {output.metabatout})/t2.txt {output.rosellaout}
        pixi r -e dastool Fasta_to_Contig2Bin.sh -e fna -i $(dirname {input.taxvamb})/bins > {output.taxvambout}
        perl -pe 's/\t/\t{params.value}_taxvamb_/g' {output.taxvambout} > $(dirname {output.metabatout})/t3.txt && mv $(dirname {output.metabatout})/t3.txt {output.taxvambout} 
        perl -pe 's/,/\t{params.value}_concoct_/g' {input.concoct} | tail -n +2 > {output.concoctout}
        perl -pe 's/\tgroup/\t{params.value}_comebin_/g' {input.comebin} > {output.comebinout}
        perl -pe 's/\t/\t{params.value}_semibin_/g' {input.semibin} | tail -n +2 > {output.semibinout}
        perl -pe 's/metacat\./metacat_/g' {input.metacat} | tail -n +2 > {output.metacatout}
        date) 2> {log.err} > {log.out}
        """


rule prepare_magscot_bin:
    input:
        metabat=rules.prepare_dasTool.output.metabatout,
        taxvamb=rules.prepare_dasTool.output.taxvambout,
        rosella=rules.prepare_dasTool.output.rosellaout,
        concoct=rules.prepare_dasTool.output.concoctout,
        comebin=rules.prepare_dasTool.output.comebinout,
        semibin=rules.prepare_dasTool.output.semibinout,
        metacat=rules.prepare_dasTool.output.metacatout,
    output:
        os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/MAGScoT/magscot_input.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/prepare_magscot/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prepare_magscot/{sample}.err.log"),
    message:
        "Preparing files for MAGScoT"
    localrule: True
    shell:
        """
        (date &&
        awk '{{print $2"\t"$1"\tmetabat"}}' {input.metabat} > {output} &&
        awk '{{print $2"\t"$1"\ttaxvamb"}}' {input.taxvamb} >> {output} &&
        awk '{{print $2"\t"$1"\trosella"}}' {input.rosella} >> {output} &&
        awk '{{print $2"\t"$1"\tconcoct"}}' {input.concoct} >> {output} &&
        awk '{{print $2"\t"$1"\tcomebin"}}' {input.comebin} >> {output} &&
        awk '{{print $2"\t"$1"\tsemibin"}}' {input.semibin} >> {output} &&
        awk '{{print $2"\t"$1"\tmetacat"}}' {input.metacat} >> {output} &&
        date) 1> {log.out} 2> {log.err}
        """


rule prepare_magscot_hmm:
    input:
        cont=rules.filter_length.output.out,
        db=rules.magscot_install.output.file,
    output:
        protein=temp(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/MAGScoT/magscot_proteins.faa"
            )
        ),
        nucl=temp(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/MAGScoT/magscot_nucl.ffn"
            )
        ),
        gff=temp(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/MAGScoT/magscot_genes.gff"
            )
        ),
        tigr=temp(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/MAGScoT/magscot_hmm.tigr.out"
            )
        ),
        tigr_hit=temp(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/MAGScoT/magscot_hmm.tigr.hit.out",
            )
        ),
        pfam=temp(
            os.path.join(
                RESULTS_DIR, "Bins/intermediate/{sample}/MAGScoT/magscot_hmm.pfam.out"
            )
        ),
        pfam_hit=temp(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/MAGScoT/magscot_hmm.pfam.hit.out",
            )
        ),
        tigr_clean=temp(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/MAGScoT/magscot_hmm.tigr.clean.tsv",
            )
        ),
        pfam_clean=temp(
            os.path.join(
                RESULTS_DIR,
                "Bins/intermediate/{sample}/MAGScoT/magscot_hmm.pfam.clean.tsv",
            )
        ),
        hmm=os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/MAGScoT/magscot_hmm.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/prepare_magscot_hmm/{sample}.hmm.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prepare_magscot_hmm/{sample}.hmm.err.log"),
    resources:
        mem_mb=get_mem_mb,
        runtime=60 * 4,
    message:
        "Preparing HMM files for MAGScoT"
    threads: 10
    shell:
        """
        (date
        pixi r -e magscot pyrodigal -a {output.protein} -j {threads} -d {output.nucl} -f gff -i {input.cont} -o {output.gff}
        pixi r -e magscot hmmsearch -o {output.tigr} --tblout {output.tigr_hit} --noali --notextw --cut_nc --cpu {threads} $(dirname {input.db})/hmm/gtdbtk_rel207_tigrfam.hmm {output.protein}
        pixi r -e magscot hmmsearch -o {output.pfam} --tblout {output.pfam_hit} --noali --notextw --cut_nc --cpu {threads} $(dirname {input.db})/hmm/gtdbtk_rel207_Pfam-A.hmm {output.protein}
        cat {output.tigr_hit} | grep -v "^#" | awk '{{print $1"\t"$3"\t"$5}}' > {output.tigr_clean}
        cat {output.pfam_hit} | grep -v "^#" | awk '{{print $1"\t"$4"\t"$5}}' > {output.pfam_clean}
        cat {output.pfam_clean} {output.tigr_clean} > {output.hmm}
        date) 2> {log.err} > {log.out}
        """


rule magscot_run:
    input:
        bins=rules.prepare_magscot_bin.output,
        hmm=rules.prepare_magscot_hmm.output.hmm,
        soft=rules.magscot_install.output.file,
    output:
        os.path.join(
            RESULTS_DIR, "Bins/intermediate/{sample}/MAGScoT/{sample}_magscot_bins.tsv"
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/magscot/{sample}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/magscot/{sample}.err.log"),
    resources:
        mem_mb=get_mem_mb,
        runtime=60 * 4,
    message:
        "Running MAGScoT"
    threads: 10
    shell:
        """
        (date &&
        pixi r -e magscot Rscript {input.soft} -i {input.bins} --hmm {input.hmm} -o {output} -t {threads}
        date) 2> {log.err} > {log.out}
        """
