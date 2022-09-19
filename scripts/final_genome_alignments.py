def is_supporting_alignment(ref_transcript, ref_chrom, ref_start, ref_stop, retrogene_locations):
	for retro_transcript, retro_chrom, retro_start, retro_stop in retrogene_locations:
		if ref_transcript == retro_transcript and ref_chrom == retro_chrom and retro_start <= ref_start <= ref_stop <= retro_stop:
			return True
	return False

retrogene_locations = []

with open(snakemake.input[0], "r") as input_AS_file:
	for line in input_AS_file:
		geneID, transcript_location, anchor_range = line.strip().split()
		transcript = transcript_location.split(":")[0]
		chrom, anchor_span = anchor_range.split(":")
		start, stop = anchor_span.split("-")
		start, stop = int(start), int(stop)
		retrogene_locations.append((transcript, chrom, start, stop))

with open(snakemake.input[1], "r") as input_headers_file, open(snakemake.output[0], "w") as output_sam_file:
	for line in input_headers_file:
		if line.startswith("@"):
			output_sam_file.write(line)
		else:
			break

with open(snakemake.input[2], "r") as input_as_sam_file, open(snakemake.output[0], "a") as output_sam_file:
	for line in input_as_sam_file:
		geneid, transcript, readname, read_start, read_stop, _, transcript_start, transcript_stop, chrom, ref_start, ref_stop = line.strip().split()[:11]
		ref_start, ref_stop = int(ref_start), int(ref_stop)
		if is_supporting_alignment(transcript, chrom, ref_start, ref_stop, retrogene_locations):
			alignment = "\t".join(line.strip().split()[13:])
			output_sam_file.write(f"{alignment}\n")
