import gzip

output_folder = snakemake.output[0][:-26]

with gzip.open(snakemake.input[0], "rt") as input_file:
	for line in input_file:
		if line.startswith(">"):
			transcript = line.strip().split()[0][1:]
			output_file = open(f"{output_folder}/{transcript}.fa", "w")
			output_file.write(line)
		else:
			output_file.write(line)

with open(snakemake.output[0], "w") as output_file:
	output_file.write(f"All transcripts processed.")
