EUK_TAX_RANKS = {
    "d": "domain",
    "k": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
EUK_TAX_DOMAINS = {"Eukaryota", "Bacteria", "Archaea", "Viruses"}
EUK_TAX_COLUMNS = [
    "genome",
    "lineage",
    "deepest_rank",
    "support",
    "classified_fraction",
    "eukaryota_fraction",
]


rule eukaryotes_taxonomy:
    input:
        os.path.join(RESULTS_DIR, "Bins/Eukaryotic_mags_taxonomy/euk_taxonomy.tsv"),
    output:
        touch("status/eukaryotes_taxonomy.done"),


checkpoint euk_contig_lengths:
    input:
        bins=rules.galah_euk.output.final,
        done=rules.galah_euk.output.done,
    output:
        os.path.join(RESULTS_DIR, "Bins/Eukaryotic_mags_taxonomy/contig_lengths.tsv"),
    localrule: True
    message:
        "Listing the contigs and lengths of the dereplicated eukaryotic MAGs"
    run:
        from pathlib import Path

        # One row per contig; also the MAG list for count_euk_mags, so no
        # temporary copy of the MAGs is needed
        with open(output[0], "w") as out:
            out.write("genome\tcontig\tlength\n")
            for fasta in sorted(Path(input.bins).glob("*.fa")):
                lengths = {}
                with open(fasta) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            contig = line[1:].split()[0]
                            lengths[contig] = 0
                        else:
                            lengths[contig] += len(line.strip())
                for contig, length in lengths.items():
                    out.write(f"{fasta.stem}\t{contig}\t{length}\n")


def count_euk_mags(wildcards):
    import pandas as pd

    ck_output = checkpoints.euk_contig_lengths.get(**wildcards).output[0]
    return expand(
        os.path.join(RESULTS_DIR, "Bins/Eukaryotic_mags_taxonomy/{mags}_tax.tsv"),
        mags=pd.read_csv(ck_output, sep="\t")["genome"].unique(),
    )


rule euk_taxonomy_query:
    input:
        bins=rules.galah_euk.output.final,
        done=rules.galah_euk.output.done,
    output:
        temp(os.path.join(RESULTS_DIR, "Bins/intermediate/euk_taxonomy_query.fa")),
    localrule: True
    message:
        "Pooling the contigs of all dereplicated eukaryotic MAGs for one MMseqs2 search"
    run:
        from pathlib import Path

        # Prefix contigs with "<mag>|" so finishEukTaxonomy can split the hits back
        with open(output[0], "w") as out:
            for fasta in sorted(Path(input.bins).glob("*.fa")):
                with open(fasta) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            line = f">{fasta.stem}|{line[1:]}"
                        out.write(line)


rule euk_taxonomy_mmseqs:
    input:
        fasta=rules.euk_taxonomy_query.output[0],
        # Path rather than rules.download_tax_euk_db, as install_software_DB.smk
        # is only included with the install_software step
        db=os.path.join(RESULTS_DIR, "DB/UniRef90_mmseqs_mapping"),
    output:
        taxonomy=os.path.join(
            RESULTS_DIR, "Bins/Eukaryotic_mags_taxonomy/all_contigs_tax_mmseqs.tsv"
        ),
        tmp=temp(
            directory(os.path.join(RESULTS_DIR, "Bins/intermediate/euk_taxonomy_tmp"))
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/taxonomy_euk/mmseqs.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/taxonomy_euk/mmseqs.err.log"),
    threads: 40
    resources:
        runtime=60 * 24,
        mem_mb=get_mem_mb,
    params:
        # MMseqs sizes its index from the node's RAM, not the SLURM allocation:
        # without a limit the UniRef90 index (~650G) exceeds mem_mb and is killed
        split_memory=lambda wildcards, resources: f"{int(resources.mem_mb*0.7/1024)}G",
    message:
        "Classifying the contigs of all eukaryotic MAGs with MMseqs2 against UniRef90"
    shell:
        """
        (
            date
            rm -rf {output.tmp}
            mkdir -p {output.tmp}
            if [ -s {input.fasta} ]; then
                pixi r -e mmseqs mmseqs createdb {input.fasta} {output.tmp}/query
                pixi r -e mmseqs mmseqs taxonomy {output.tmp}/query $(dirname {input.db})/UniRef90_mmseqs {output.tmp}/result {output.tmp}/work --threads {threads} --tax-lineage 1 --split-memory-limit {params.split_memory}
                pixi r -e mmseqs mmseqs createtsv {output.tmp}/query {output.tmp}/result {output.taxonomy}
            else
                echo "No eukaryotic MAGs to classify; writing an empty taxonomy table."
                : >{output.taxonomy}
            fi
            date
        ) 2>{log.err} >{log.out}
        """


rule finishEukTaxonomy:
    input:
        taxonomy=rules.euk_taxonomy_mmseqs.output.taxonomy,
        lengths=rules.euk_contig_lengths.output[0],
    output:
        os.path.join(RESULTS_DIR, "Bins/Eukaryotic_mags_taxonomy/{mags}_tax.tsv"),
    localrule: True
    message:
        "Building the consensus taxonomy of {wildcards.mags}"
    run:
        import pandas as pd

        contigs = pd.read_csv(input.lengths, sep="\t", dtype={"contig": str})
        contigs = contigs[contigs["genome"] == wildcards.mags]
        lengths = dict(zip(contigs["contig"], contigs["length"]))
        mag_length = sum(lengths.values())
        # createtsv columns: contig, taxid, rank, name, fragment counts and
        # support, then the lineage last (--tax-lineage 1)
        if os.path.getsize(input.taxonomy) > 0:
            raw = pd.read_csv(input.taxonomy, sep="\t", header=None, dtype=str)
            tax = pd.DataFrame(
                {"contig": raw[0], "taxid": raw[1], "lineage": raw.iloc[:, -1]}
            )
        else:
            tax = pd.DataFrame(columns=["contig", "taxid", "lineage"], dtype=str)
        # Keep this MAG's contigs from the pooled search and drop the "<mag>|" prefix
        prefix = f"{wildcards.mags}|"
        tax = tax[tax["contig"].str.startswith(prefix)].copy()
        tax["contig"] = tax["contig"].str[len(prefix) :]
        tax = tax[(tax["taxid"] != "0") & tax["lineage"].notna()]
        # Keep the ranked levels only (d_, k_, ..., s_); "-_" entries are unranked,
        # except the domain: recent NCBI taxonomy ranks it "domain", which MMseqs
        # writes as "-_Eukaryota" rather than "d_Eukaryota"
        parsed = []
        for lineage in tax["lineage"]:
            levels = {}
            for taxon in lineage.split(";"):
                prefix, _, name = taxon.partition("_")
                if prefix in EUK_TAX_RANKS:
                    levels[prefix] = taxon
                elif name in EUK_TAX_DOMAINS:
                    levels["d"] = f"d_{name}"
            parsed.append(levels)
        weights = [lengths.get(c, 0) for c in tax["contig"]]
        classified = sum(weights)
        # Walk down the ranks, keeping a taxon while it holds >50% of the
        # classified length among contigs agreeing with the consensus so far.
        # A rank that most of that length lacks (often kingdom for protists)
        # is skipped rather than ending the walk
        consensus, deepest_rank, support = [], "none", 0.0
        consistent = [True] * len(parsed)
        for rank in EUK_TAX_RANKS:
            if not classified:
                break
            votes, missing = {}, 0
            for levels, weight, ok in zip(parsed, weights, consistent):
                if not ok:
                    continue
                if rank in levels:
                    votes[levels[rank]] = votes.get(levels[rank], 0) + weight
                else:
                    missing += weight
            if not votes or missing / classified > 0.5:
                continue
            taxon, weight = max(votes.items(), key=lambda kv: kv[1])
            if weight / classified <= 0.5:
                break
            consensus.append(taxon)
            deepest_rank, support = EUK_TAX_RANKS[rank], weight / classified
            consistent = [
                ok and levels.get(rank) == taxon
                for levels, ok in zip(parsed, consistent)
            ]
        eukaryota = sum(
            weight
            for levels, weight in zip(parsed, weights)
            if levels.get("d") == "d_Eukaryota"
        )
        pd.DataFrame(
            [
                {
                    "genome": wildcards.mags,
                    "lineage": ";".join(consensus) if consensus else "Unclassified",
                    "deepest_rank": deepest_rank,
                    "support": round(support, 4),
                    "classified_fraction": (
                        round(classified / mag_length, 4) if mag_length else 0.0
                    ),
                    "eukaryota_fraction": (
                        round(eukaryota / classified, 4) if classified else 0.0
                    ),
                }
            ],
            columns=EUK_TAX_COLUMNS,
        ).to_csv(output[0], sep="\t", index=False)


rule get_euk_taxonomy:
    input:
        count_euk_mags,
    output:
        os.path.join(RESULTS_DIR, "Bins/Eukaryotic_mags_taxonomy/euk_taxonomy.tsv"),
    localrule: True
    message:
        "Combining the taxonomy of all dereplicated eukaryotic MAGs"
    run:
        import pandas as pd

        frames = [pd.read_csv(f, sep="\t") for f in input]
        if frames:
            df = pd.concat(frames, ignore_index=True).sort_values("genome")
        else:
            df = pd.DataFrame(columns=EUK_TAX_COLUMNS)
        df.to_csv(output[0], sep="\t", index=False)
