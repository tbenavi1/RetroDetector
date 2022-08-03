geneid_to_best = {}

with open(snakemake.input[0], "r") as input_file:
	for line in input_file:
		geneid, transcript, location, AS = line.strip().split()
		AS = int(AS)
		if geneid in geneid_to_best:
			previous_transcript, previous_location, previous_AS = geneid_to_best[geneid]
			if AS > previous_AS:
				geneid_to_best[geneid] = (transcript, location, AS)
		else:
			geneid_to_best[geneid] = (transcript, location, AS)

with open(snakemake.output[0], "w") as output_file:
	for geneid in geneid_to_best:
		transcript, location, AS = geneid_to_best[geneid]
		output_file.write(f"{geneid}\t{transcript}\t{location}\t{AS}\n")
