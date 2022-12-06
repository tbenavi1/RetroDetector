from collections import defaultdict

short_read_strong_spanning_alignment_expected_genomic_insertion_size_threshold = snakemake.params.short_read_strong_spanning_alignment_expected_genomic_insertion_size_threshold
junction_total_read_support_threshold = snakemake.params.junction_total_read_support_threshold
junction_strong_short_read_support_threshold = snakemake.params.junction_strong_short_read_support_threshold

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
transcript_junction_to_strong_reads = defaultdict(set)
no_non_overlapping_transcript_junction_to_reads = defaultdict(set)
no_non_overlapping_transcript_junction_to_strong_reads = defaultdict(set)
no_alternate_transcript_junction_to_reads = defaultdict(set)
no_alternate_transcript_junction_to_strong_reads = defaultdict(set)
with open(snakemake.input[1], "r") as input_spannedjunctions_file:
	for line in input_spannedjunctions_file:
		geneid, transcript, junction, readname, overlapping_var, expected_genome_tlen = line.strip().split()
		#if this is an alternate alignment
		if "lcl|" in readname:
			is_alternate = True
		else:
			is_alternate = False
		#if this is a non overlapping alignment
		if overlapping_var == "non-overlapping":
			is_nonoverlapping = True
		else:
			is_nonoverlapping = False
		if expected_genome_tlen != ".":
			expected_genome_tlen = int(expected_genome_tlen)
		else:
			expected_genome_tlen = 0
		#if the alignment is overlapping or the expected genome tlen is at least the strong threshold
		if not is_nonoverlapping or expected_genome_tlen >= short_read_strong_spanning_alignment_expected_genomic_insertion_size_threshold:
			is_strong = True
		else:
			is_strong = False
		transcript_junction_to_reads[(transcript, junction)].add(readname)
		if is_strong:
			transcript_junction_to_strong_reads[(transcript, junction)].add(readname)
		if not is_nonoverlapping:
			no_non_overlapping_transcript_junction_to_reads[(transcript, junction)].add(readname)
			if is_strong:
				no_non_overlapping_transcript_junction_to_strong_reads[(transcript, junction)].add(readname)
		if not is_alternate:
			no_alternate_transcript_junction_to_reads[(transcript, junction)].add(readname)
			if is_strong:
				no_alternate_transcript_junction_to_strong_reads[(transcript, junction)].add(readname)

with open(snakemake.output[0], "w") as output_coverage_file, open(snakemake.output[1], "w") as output_numspanned_file, open(snakemake.output[2], "w") as output_no_non_overlapping_coverage_file, open(snakemake.output[3], "w") as output_no_non_overlapping_numspanned_file, open(snakemake.output[4], "w") as output_no_alternate_coverage_file, open(snakemake.output[5], "w") as output_no_alternate_numspanned_file:
	for transcript in sorted(transcript_to_junctions):
		geneid = transcript_to_geneid[transcript]
		junctions = transcript_to_junctions[transcript]
		num_junctions = len(junctions)
		num_junctions_spanned_by_thresholds = 0
		no_non_overlapping_num_junctions_spanned_by_thresholds = 0
		no_alternate_num_junctions_spanned_by_thresholds = 0
		output_coverage_line = f"{geneid}\t{transcript}"
		no_non_overlapping_output_coverage_line = f"{geneid}\t{transcript}"
		no_alternate_output_coverage_line = f"{geneid}\t{transcript}"
		for junction in junctions:
			reads = transcript_junction_to_reads[(transcript, junction)]
			strong_reads = transcript_junction_to_strong_reads[(transcript, junction)]
			no_non_overlapping_reads = no_non_overlapping_transcript_junction_to_reads[(transcript, junction)]
			no_non_overlapping_strong_reads = no_non_overlapping_transcript_junction_to_strong_reads[(transcript, junction)]
			no_alternate_reads = no_alternate_transcript_junction_to_reads[(transcript, junction)]
			no_alternate_strong_reads = no_alternate_transcript_junction_to_strong_reads[(transcript, junction)]
			num_reads = len(reads)
			num_strong_reads = len(strong_reads)
			num_no_non_overlapping_reads = len(no_non_overlapping_reads)
			num_no_non_overlapping_strong_reads = len(no_non_overlapping_strong_reads)
			num_no_alternate_reads = len(no_alternate_reads)
			num_no_alternate_strong_reads = len(no_alternate_strong_reads)
			output_coverage_line += f"\t{junction}\t{num_reads}"
			no_non_overlapping_output_coverage_line += f"\t{junction}\t{num_no_non_overlapping_reads}"
			no_alternate_output_coverage_line += f"\t{junction}\t{num_no_alternate_reads}"
			#if it has enough total read support and enough strong read support
			if num_reads >= junction_total_read_support_threshold and num_strong_reads >= junction_strong_short_read_support_threshold:
				num_junctions_spanned_by_thresholds += 1
			if num_no_non_overlapping_reads >= junction_total_read_support_threshold and num_no_non_overlapping_strong_reads >= junction_strong_short_read_support_threshold:
				no_non_overlapping_num_junctions_spanned_by_thresholds += 1
			if num_no_alternate_reads >= junction_total_read_support_threshold and num_no_alternate_strong_reads >= junction_strong_short_read_support_threshold:
				no_alternate_num_junctions_spanned_by_thresholds += 1
		output_coverage_line += "\n"
		no_non_overlapping_output_coverage_line += "\n"
		no_alternate_output_coverage_line += "\n"
		if num_junctions_spanned_by_thresholds > 0:
			output_coverage_file.write(output_coverage_line)
			output_numspanned_file.write(f"{geneid}\t{transcript}\t{num_junctions}\t{num_junctions_spanned_by_thresholds}\n")
		if no_non_overlapping_num_junctions_spanned_by_thresholds > 0:
			output_no_non_overlapping_coverage_file.write(no_non_overlapping_output_coverage_line)
			output_no_non_overlapping_numspanned_file.write(f"{geneid}\t{transcript}\t{num_junctions}\t{no_non_overlapping_num_junctions_spanned_by_thresholds}\n")
		if no_alternate_num_junctions_spanned_by_thresholds > 0:
			output_no_alternate_coverage_file.write(no_alternate_output_coverage_line)
			output_no_alternate_numspanned_file.write(f"{geneid}\t{transcript}\t{num_junctions}\t{no_alternate_num_junctions_spanned_by_thresholds}\n")
