with open(snakemake.input[0], "r") as input_results_file, open(snakemake.output[0], "w") as output_file:
	for line in input_results_file:
		parent_location, retrogene_location = line.strip().split()[1:3]
		parent_chrom = parent_location.split(":")[0]
		retro_chrom = retrogene_location.split(":")[0]
		if "chroms" not in snakemake.config or ("chroms" in snakemake.config and parent_chrom in snakemake.config["chroms"].values() and retro_chrom in snakemake.config["chroms"].values()):
			output_file.write(line)
