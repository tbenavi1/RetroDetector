chrom_to_nicechrom = snakemake.config["chroms"]

with open(snakemake.input[0], "r") as input_results_file, open(snakemake.output[0], "w") as output_results_file:
	for line in input_results_file:
		geneid, transcript_genome_location, retrogene_genome_location, retrogene_transcriptome_location, strand, retrogene_length, coverage_pct, pct_identity = line.strip().split()
		chrom, span = transcript_genome_location.split(":")
		nicechrom = chrom_to_nicechrom[chrom]
		transcript_genome_location = f"{nicechrom}:{span}"
		chrom, span = retrogene_genome_location.split(":")
		nicechrom = chrom_to_nicechrom[chrom]
		retrogene_genome_location = f"{nicechrom}:{span}"
		output_results_file.write(f"{geneid}\t{transcript_genome_location}\t{retrogene_genome_location}\t{retrogene_transcriptome_location}\t{strand}\t{retrogene_length}\t{coverage_pct}\t{pct_identity}\n")
