print(snakemake.output[0])
with open(snakemake.input[0], "r") as input_file:
	for line in input_file:
		GeneID = line.strip().split()[0]
		with open(snakemake.output[0]+f"/{GeneID}.txt", "w") as output_file:
			output_file.write(f"{GeneID}\n")
