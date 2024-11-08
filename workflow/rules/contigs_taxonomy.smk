rule contigs_taxonomy:
    input:
       expand(os.path.join(RESULTS_DIR, "Assembly/{sample}_taxonomy/taxonomy.txt"),sample = SAMPLES)
    output:
       touch("status/taxonomy.done")


rule download_tax_db:
	output:
		file = os.path.join(RESULTS_DIR, "DB/GTDB_mmseqs_mapping"),
		folder = os.path.join(RESULTS_DIR, "DB"),
		tmp = temp(directory(os.path.join(RESULTS_DIR, "tmp")))
	log:
		out=os.path.join(RESULTS_DIR, "logs/download_tax_db.out.log"),
		err=os.path.join(RESULTS_DIR, "logs/download_tax_db.err.log")
	threads:
		20
	conda:
		os.path.join(ENV_DIR, "mmseqs.yaml")
	shell:
		"(date && name=$(echo $(basename {output.file}) | sed \"s/_mapping//g\") && "
        "mmseqs databases GTDB {output.folder}/$name {output.tmp} --threads {threads} && date) 2> {log.err} > {log.out}"

rule taxonomy:
	input:
		fasta = rules.assembly_megahit.output.final,
		db = rules.download_tax_db.output.file
	output:
		taxonomy = os.path.join(RESULTS_DIR, "Assembly/{sample}_taxonomy/taxonomy.txt"),
		tmp = temp(os.path.join(RESULTS_DIR, "Assembly/{sample}_tmp"))
	log:
		out=os.path.join(RESULTS_DIR, "logs/taxonomy_contig/{sample}.out.log"),
		err=os.path.join(RESULTS_DIR, "logs/taxonomy_contig/{sample}.err.log")
	threads:
		20
	conda:
		os.path.join(ENV_DIR, "mmseqs.yaml")
	shell:
		"mmseqs taxonomy {input.fasta} $(dirname {input.db}) $(dirname {output.taxonomy}) {output.tmp} -threads {threads}"
