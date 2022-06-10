from collections import defaultdict

transcript_to_gene = {}
with open(snakemake.input[0], "r") as input_coords_file:
	for line in input_coords_file:
		transcript, gene = line.strip().split()[:2]
		#for flybase_dyak
		if "-" in gene:
			gene = gene.split("-")[0]
		transcript_to_gene[transcript] = gene

transcript_to_junctions = {}
with open(snakemake.input[1], "r") as input_junctions_file:
	for line in input_junctions_file:
		transcript, *junctions = line.strip().split()
		transcript_to_junctions[transcript] = junctions

transcript_junction_to_reads = defaultdict(set)
with open(snakemake.input[2], "r") as input_spannedjunctions_file:
	for line in input_spannedjunctions_file:
		transcript, junction, readname = line.strip().split()
		transcript_junction_to_reads[(transcript, junction)].add(readname)

read_threshold = snakemake.config["read_threshold"]

with open(snakemake.output[0], "w") as output_coverage_file, open(snakemake.output[1], "w") as output_numspanned_file:
	for transcript in sorted(transcript_to_junctions):
		gene = transcript_to_gene[transcript]
		junctions = transcript_to_junctions[transcript]
		num_junctions = len(junctions)
		num_junctions_spanned_by_threshold = 0
		output_coverage_line = f"{transcript}\t{gene}"
		for junction in junctions:
			reads = transcript_junction_to_reads[(transcript, junction)]
			num_reads = len(reads)
			output_coverage_line += f"\t{junction}\t{num_reads}"
			if num_reads >= read_threshold:
				num_junctions_spanned_by_threshold += 1
		output_coverage_line += "\n"
		if num_junctions_spanned_by_threshold > 0:
			output_coverage_file.write(output_coverage_line)
			output_numspanned_file.write(f"{transcript}\t{gene}\t{num_junctions}\t{num_junctions_spanned_by_threshold}\n")
