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
transcript_junction_to_has_overlapping = defaultdict(bool)
transcript_junction_to_has_nonoverlapping = defaultdict(bool)
transcript_junction_to_has_large_tlen = defaultdict(bool)
with open(snakemake.input[1], "r") as input_spannedjunctions_file:
	for line in input_spannedjunctions_file:
		geneid, transcript, junction, readname, overlapping_var, expected_genome_tlen = line.strip().split()
		if expected_genome_tlen != ".":
			expected_genome_tlen = int(expected_genome_tlen)
		else:
			expected_genome_tlen = 0
		transcript_junction_to_reads[(transcript, junction)].add(readname)
		if overlapping_var == "overlapping":
			transcript_junction_to_has_overlapping[(transcript, junction)] = True
		else:
			assert overlapping_var == "non-overlapping", line
			transcript_junction_to_has_nonoverlapping[(transcript, junction)] = True
		if expected_genome_tlen > 1000:
			transcript_junction_to_has_large_tlen[(transcript, junction)] = True

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
			has_overlapping = transcript_junction_to_has_overlapping[(transcript, junction)]
			has_nonoverlapping = transcript_junction_to_has_nonoverlapping[(transcript, junction)]
			has_large_tlen = transcript_junction_to_has_large_tlen[(transcript, junction)]
			num_reads = len(reads)
			output_coverage_line += f"\t{junction}\t{num_reads}"
			if num_reads >= read_threshold and (has_overlapping or has_large_tlen) and has_nonoverlapping:
				num_junctions_spanned_by_threshold += 1
		output_coverage_line += "\n"
		if num_junctions_spanned_by_threshold > 0:
			output_coverage_file.write(output_coverage_line)
			output_numspanned_file.write(f"{geneid}\t{transcript}\t{num_junctions}\t{num_junctions_spanned_by_threshold}\n")
