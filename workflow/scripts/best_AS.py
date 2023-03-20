read_support = snakemake.params.read_support

previous_geneid = ""
with open(snakemake.input[0], "r") as input_file, open(snakemake.output[0], "w") as output_file:
	for line in input_file:
		geneid, transcript_location, location, AS, direction, readnames = line.strip().split()
		#if we are starting a new cluster
		if geneid != previous_geneid:
			geneid_readnames = set()
		readnames = set(readnames.split(","))
		#remove reads from previous clusters
		readnames = readnames - geneid_readnames
		num_reads = len(readnames)
		#if this cluster has sufficient support
		if num_reads >= read_support:
			geneid_readnames.update(readnames)
			readnames = ",".join(readnames)
			output_file.write(f"{geneid}\t{transcript_location}\t{location}\t{AS}\t{direction}\t{readnames}\n")
		previous_geneid = geneid
