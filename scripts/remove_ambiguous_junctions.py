from collections import defaultdict

transcript_to_spannedjunctions = defaultdict(set)
with open(snakemake.input[0], "r") as input_spannedjunctions_file:
	for line in input_spannedjunctions_file:
		transcript, spanned_junction, readname = line.strip().split()
		transcript_to_spannedjunctions[transcript].add(spanned_junction)

with open(snakemake.input[1], "r") as input_junctions_file, open(snakemake.output[0], "w") as output_junctions_file:
	for line in input_junctions_file:
		transcript, *junctions = line.strip().split()
		spanned_junctions = transcript_to_spannedjunctions[transcript]
		output_line = transcript
		keep_transcript = False
		for junction in junctions:
			if junction not in spanned_junctions:
				keep_transcript = True
				output_line += f"\t{junction}"
		output_line += "\n"
		if keep_transcript:
			output_junctions_file.write(output_line)
