ref = snakemake.wildcards.ref
sample = snakemake.wildcards.sample
junction_overhang = snakemake.config["junction_overhang"]
max_insertions = snakemake.config["max_insertions"]
read_support = snakemake.config["read_support"]

input_fasta_prefix = f"results/retrogenes/{ref}/{sample}/long/consensus/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.retrogene."
input_fasta_suffix = ".consensus.fasta"

with open(snakemake.input[0], "r") as input_file, open(snakemake.output[0], "w") as output_file:
	for line in input_file:
		geneid, _, retrogene_genome_location = line.strip().split()[:3]
		retrogene_chrom, retrogene_span = retrogene_genome_location.split(":")
		if "chroms" in snakemake.config:
			chrom_to_nicechrom = snakemake.config["chroms"]
			nicechrom_to_chrom = {v: k for k, v in chrom_to_nicechrom.items()}
			original_chrom = nicechrom_to_chrom[retrogene_chrom]
			#retrogene_genome_location = f"{original_chrom}:{retrogene_span}"
		else:
			original_chrom = retrogene_chrom
		input_fasta = input_fasta_prefix + geneid + "." + original_chrom + "." + retrogene_span + input_fasta_suffix
		with open(input_fasta, "r") as input_fasta_file:
			for line2 in input_fasta_file:
				if line2.startswith(">"):
					output_file.write(f">{geneid}_{retrogene_genome_location}\n")
				else:
					output_file.write(line2)
