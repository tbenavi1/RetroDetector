supporting_reads = set()
with open(snakemake.input[0], "r") as input_junctions_file:
	for line in input_junctions_file:
		transcript, junction, read = line.strip().split()
		supporting_reads.add(read)

with open(snakemake.output[0], "w") as output_reads_file:
	for read in supporting_reads:
		output_reads_file.write(f"{read}\n")
