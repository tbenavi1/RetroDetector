with open(snakemake.input[0], "r") as input_file, open(snakemake.output[0], "w") as output_file:
	for line in input_file:
		genome_location = line.strip().split()[2]
		genome_range = genome_location.split(":")[1]
		start, stop = genome_range.split("-")
		start, stop = int(start), int(stop)
		length = stop - start + 1
		output_file.write(f"{line.strip()}\t{length}\n")
