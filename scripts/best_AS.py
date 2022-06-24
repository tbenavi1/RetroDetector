geneid_to_best = {}

with open(snakemake.input[0], "r") as input_file:
	for line in input_file:
		geneid, location, AS = line.strip().split()
		AS = int(AS)
		if geneid in geneid_to_best:
			previous_location, previous_AS = geneid_to_best[geneid]
			if AS > previous_AS:
				geneid_to_best[geneid] = (location, AS)
		else:
			geneid_to_best[geneid] = (location, AS)

with open(snakemake.output[0], "w") as output_file:
	for geneid in geneid_to_best:
		location, AS = geneid_to_best[geneid]
		output_file.write(f"{geneid}\t{location}\t{AS}\n")
