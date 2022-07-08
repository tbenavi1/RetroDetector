from collections import defaultdict

transcript_to_junctions = {}
transcript_to_geneid = {}
with open(snakemake.input[0], "r") as input_junctions_file:
	for line in input_junctions_file:
		geneid, transcript, *junctions = line.strip().split()
		#for flybase_dyak
		if "-" in geneid:
			geneid = geneid.split("-")[0]
		transcript_to_junctions[transcript] = junctions
		transcript_to_geneid[transcript] = geneid

transcript_junction_to_reads = defaultdict(set)
with open(snakemake.input[1], "r") as input_spannedjunctions_file:
	for line in input_spannedjunctions_file:
		geneid, transcript, junction, readname = line.strip().split()
		transcript_junction_to_reads[(transcript, junction)].add(readname)

read_threshold = snakemake.config["read_threshold"]

with open(snakemake.output[0], "w") as output_coverage_file, open(snakemake.output[1], "w") as output_numspanned_file:
	for transcript in sorted(transcript_to_junctions):
		geneid = transcript_to_geneid[transcript]
		junctions = transcript_to_junctions[transcript]
		num_junctions = len(junctions)
		num_junctions_spanned_by_threshold = 0
		output_coverage_line = f"{geneid}\t{transcript}"
		for junction in junctions:
			reads = transcript_junction_to_reads[(transcript, junction)]
			num_reads = len(reads)
			output_coverage_line += f"\t{junction}\t{num_reads}"
			if num_reads >= read_threshold:
				num_junctions_spanned_by_threshold += 1
		output_coverage_line += "\n"
		if num_junctions_spanned_by_threshold > 0:
			output_coverage_file.write(output_coverage_line)
			output_numspanned_file.write(f"{geneid}\t{transcript}\t{num_junctions}\t{num_junctions_spanned_by_threshold}\n")
