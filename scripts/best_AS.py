geneid_to_best = {}

with open(snakemake.input[0], "r") as input_file:
	for line in input_file:
		geneid, transcript_location, location, AS, direction, readnames = line.strip().split()
		AS = int(AS)
		if geneid in geneid_to_best:
			previous_transcript_location, previous_location, previous_AS, previous_direction, previous_readnames = geneid_to_best[geneid]
			if AS > previous_AS:
				geneid_to_best[geneid] = (transcript_location, location, AS, direction, readnames)
		else:
			geneid_to_best[geneid] = (transcript_location, location, AS, direction, readnames)

with open(snakemake.output[0], "w") as output_file:
	for geneid in geneid_to_best:
		transcript_location, location, AS, direction, readnames = geneid_to_best[geneid]
		output_file.write(f"{geneid}\t{transcript_location}\t{location}\t{AS}\t{direction}\t{readnames}\n")
