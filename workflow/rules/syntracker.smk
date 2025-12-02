rule syntracker:
    input:
        expand(os.path.join(RESULTS_DIR, "SynTracker/blastDB")),
        expand(os.path.join(RESULTS_DIR, "SynTracker/syntracker_finalized.done")),
    output:
        touch("status/syntracker.done"),


rule cp_assemblies:
    input:
        expand(rules.filter_length.output.out, sample=SAMPLES),
    output:
        dir=temp(directory(os.path.join(RESULTS_DIR, "SynTracker/assemblies"))),
    localrule: True
    log:
        out=os.path.join(RESULTS_DIR, "logs/cp_assemblies.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/cp_assemblies.err.log"),
    shell:
        """
        (date
        mkdir -p {output.dir}
        cp {input} {output.dir}/
        )
        2> {log.err} > {log.out}
        """


rule syntracker_blast:
    input:
        syntracker=rules.syntracker_install.output.dbs,
        assembly=rules.cp_assemblies.output.dir,
    output:
        dir=directory(os.path.join(RESULTS_DIR, "SynTracker/blastDB")),
        done=touch(os.path.join(RESULTS_DIR, "SynTracker/blastDB.done")),
    log:
        out=os.path.join(RESULTS_DIR, "logs/syntracker_blast.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/syntracker_blast.err.log"),
    resources:
        mem_mb=get_mem_mb,
        runtime=4 * 60,
    threads: 4
    shell:
        """
        (date
        rm -rf {output.dir}
        pixi r -e syntracker python {input.syntracker}/syntracker_makeDB.py -target {input.assembly} -out {output.dir}
        touch {output.done}
        date ) 2> {log.err} > {log.out}
        """


def get_syn_mags(wildcards):
    checkpoint_output = checkpoints.galah.get(**wildcards).output[0]
    return expand(
        os.path.join(
            RESULTS_DIR,
            "SynTracker/temp/syntracker_output_{mags}/summary_output/synteny_scores_per_region.csv",
        ),
        mags=glob_wildcards(os.path.join(checkpoint_output, "{mags}.fa")).mags,
    )


rule syntracker_run:
    input:
        syntracker=rules.syntracker_install.output.dbs,
        blastdb=rules.syntracker_blast.output.dir,
        blastdb_done=rules.syntracker_blast.output.done,
        mags=os.path.join(RESULTS_DIR, "Bins/finalBins/{mags}.fa"),
    output:
        done=os.path.join(
            RESULTS_DIR,
            "SynTracker/temp/syntracker_output_{mags}/summary_output/synteny_scores_per_region.csv",
        ),
        folder=
            directory(
                os.path.join(RESULTS_DIR, "SynTracker/temp/syntracker_output_{mags}")
        ),
    log:
        out=os.path.join(RESULTS_DIR, "logs/syntracker/syntracker_run_{mags}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/syntracker/syntracker_run_{mags}.err.log"),
    resources:
        mem_mb=get_mem_mb,
        runtime=60,
    # group:
    #     "syntracker"
    threads: 5
    shell:
        """
        (date
        mkdir -p {output.folder}_mag
        cp {input.mags} {output.folder}_mag/
        pixi r -e syntracker python {input.syntracker}/syntracker.py -ref {output.folder}_mag -mode new -blastDB {input.blastdb} -out {output.folder} -cores {threads} -length 2500
        rm -r {output.folder}_mag
        date ) 2> {log.err} > {log.out}
        """


rule syntracker_finalize:
    input:
        get_syn_mags,
    output:
        score_100=os.path.join(
            RESULTS_DIR, "SynTracker/avg_synteny_scores_100_regions.csv"
        ),
        score_200=os.path.join(
            RESULTS_DIR, "SynTracker/avg_synteny_scores_200_regions.csv"
        ),
        score_40=os.path.join(
            RESULTS_DIR, "SynTracker/avg_synteny_scores_40_regions.csv"
        ),
        score_60=os.path.join(
            RESULTS_DIR, "SynTracker/avg_synteny_scores_60_regions.csv"
        ),
        score_80=os.path.join(
            RESULTS_DIR, "SynTracker/avg_synteny_scores_80_regions.csv"
        ),
        score_all=os.path.join(
            RESULTS_DIR, "SynTracker/avg_synteny_scores_all_regions.csv"
        ),
        score_per=os.path.join(RESULTS_DIR, "SynTracker/synteny_scores_per_region.csv"),
        done=os.path.join(RESULTS_DIR, "SynTracker/syntracker_finalized.done"),
    log:
        out=os.path.join(RESULTS_DIR, "logs/syntracker/syntracker_finalize.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/syntracker/syntracker_finalize.err.log"),
    message:
        "Finalizing SynTracker output"
    localrule: True
    shell:
        """
        (
        date

        out_dir=$(dirname "{output.score_all}")

        for csv in avg_synteny_scores_100_regions.csv \
                   avg_synteny_scores_40_regions.csv \
                   avg_synteny_scores_60_regions.csv \
                   avg_synteny_scores_80_regions.csv \
                   avg_synteny_scores_200_regions.csv \
                   avg_synteny_scores_all_regions.csv \
                   synteny_scores_per_region.csv;
        do
            first=1
            dest_file="${{out_dir}}/${{csv}}"

            # Ensure destination directory exists
            mkdir -p "${{out_dir}}"

            for file in {input}; do
                src_file="$(dirname "$file")/${{csv}}"
                if [ -f "$src_file" ]; then
                    if [ $first -eq 1 ]; then
                        cat "$src_file" > "$dest_file"
                        first=0
                    else
                        tail -n +2 "$src_file" >> "$dest_file"
                    fi
                fi
            done
        done

        touch {output.done}
        date
        ) > {log.out} 2> {log.err}
        """
