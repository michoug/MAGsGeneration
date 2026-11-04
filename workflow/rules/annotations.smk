rule annotations:
    input:
        expand(os.path.join(RESULTS_DIR, "Bins/Annotation/files.txt")),
        os.path.join(RESULTS_DIR, "qc/mags/mag_rna_counts.tsv"),
    output:
        touch("status/annotations.done"),


checkpoint countGenome:
    input:
        os.path.join(RESULTS_DIR, "Bins/finalBins"),
    output:
        temp(directory(os.path.join(RESULTS_DIR, "Bins/intermediate/tempBins"))),
    shell:
        "mkdir {output} && cp {input}/*fa {output}/"


def count_mags_prots(wildcards):
    ck_output = checkpoints.countGenome.get(**wildcards).output[0]
    return expand(
        os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.faa"),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )


rule annotation:
    input:
        fasta=os.path.join(RESULTS_DIR, "Bins/intermediate/tempBins/{mags}.fa"),
    output:
        protein=os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.faa"),
        gff=os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.gff"),
        nucl=os.path.join(RESULTS_DIR, "Bins/Annotation/{mags}.ffn"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/annotations/{mags}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/annotations/{mags}.err.log"),
    threads: 4
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    params:
        prefix="{mags}",
    message:
        "Prodigal annotation for {wildcards.mags}"
    shell:
        """
        (
            date
            pixi r -e annotations pyrodigal -a {output.protein} -j {threads} -d {output.nucl} -f gff -i {input.fasta} -o {output.gff}
            date
        ) 2>{log.err} >{log.out}
        """


rule prepare4egg:
    input:
        count_mags_prots,
    output:
        os.path.join(RESULTS_DIR, "Bins/Annotation/files.txt"),
    localrule: True
    shell:
        "ls {input} > {output}"


def count_mags_rna(wildcards):
    ck_output = checkpoints.countGenome.get(**wildcards).output[0]
    return expand(
        os.path.join(RESULTS_DIR, "Bins/Annotation/RNA/{mags}_trna.txt"),
        mags=glob_wildcards(os.path.join(ck_output, "{mags}.fa")).mags,
    )


rule mag_rna:
    # rRNAs with barrnap (bacterial and archaeal models, the one with more hits
    # is kept when counting), tRNAs with tRNAscan-SE (general model, so no
    # domain is needed) and CRISPR arrays with MinCED
    input:
        fasta=os.path.join(RESULTS_DIR, "Bins/intermediate/tempBins/{mags}.fa"),
    output:
        rrna_bac=os.path.join(RESULTS_DIR, "Bins/Annotation/RNA/{mags}_rrna_bac.gff"),
        rrna_arc=os.path.join(RESULTS_DIR, "Bins/Annotation/RNA/{mags}_rrna_arc.gff"),
        trna=os.path.join(RESULTS_DIR, "Bins/Annotation/RNA/{mags}_trna.txt"),
        crispr=os.path.join(RESULTS_DIR, "Bins/Annotation/RNA/{mags}_crispr.gff"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/annotations/{mags}_rna.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/annotations/{mags}_rna.err.log"),
    threads: 4
    resources:
        runtime=60,
        mem_mb=get_mem_mb,
    message:
        "Annotating rRNAs, tRNAs and CRISPR arrays of {wildcards.mags}"
    shell:
        """
        (
            date
            pixi r -e annotations barrnap --kingdom bac --threads {threads} --quiet {input.fasta} >{output.rrna_bac}
            pixi r -e annotations barrnap --kingdom arc --threads {threads} --quiet {input.fasta} >{output.rrna_arc}
            rm -f {output.trna}
            pixi r -e annotations tRNAscan-SE -G --thread {threads} -q -o {output.trna} {input.fasta}
            touch {output.trna}
            pixi r -e annotations minced {input.fasta} $(dirname {output.crispr})/{wildcards.mags}_crispr.txt {output.crispr}
            touch {output.crispr}
            date
        ) 2>{log.err} >{log.out}
        """


rule mag_rna_counts:
    input:
        count_mags_rna,
    output:
        os.path.join(RESULTS_DIR, "qc/mags/mag_rna_counts.tsv"),
    localrule: True
    message:
        "Counting rRNAs, tRNAs and CRISPR arrays per MAG"
    run:
        import pandas as pd

        def rrna_counts(gff):
            counts = {"16S_rRNA": 0, "23S_rRNA": 0, "5S_rRNA": 0}
            with open(gff) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    for gene in counts:
                        if f"Name={gene}" in line:
                            counts[gene] += 1
            return counts

        rows = []
        for trna in input:
            prefix = trna.removesuffix("_trna.txt")
            genome = os.path.basename(prefix)
            rrna = max(
                (rrna_counts(f"{prefix}_rrna_{k}.gff") for k in ("bac", "arc")),
                key=lambda c: sum(c.values()),
            )
            # tRNAscan-SE table: 3 header lines, then one tRNA per line;
            # pseudogenes are flagged in the last (Note) column
            with open(trna) as fh:
                lines = fh.read().splitlines()[3:]
            n_trna = sum(1 for l in lines if l.strip() and "pseudo" not in l)
            with open(f"{prefix}_crispr.gff") as fh:
                n_crispr = sum(1 for l in fh if not l.startswith("#") and l.strip())
            rows.append(
                {
                    "genome": genome,
                    "rRNAs": sum(rrna.values()),
                    **rrna,
                    "tRNAs": n_trna,
                    "CRISPR_arrays": n_crispr,
                }
            )
        columns = [
            "genome",
            "rRNAs",
            "16S_rRNA",
            "23S_rRNA",
            "5S_rRNA",
            "tRNAs",
            "CRISPR_arrays",
        ]
        pd.DataFrame(rows, columns=columns).sort_values("genome").to_csv(
            output[0], sep="\t", index=False
        )
