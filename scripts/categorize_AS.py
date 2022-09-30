dist_threshold = 10000

GeneID_to_location = {}
with open(snakemake.input[1], "r") as input_coords_file:
  for line in input_coords_file:
    GeneID, transcript, chrom, start, stop = line.strip().split()
    #for flybase_dyak
    if "-" in GeneID:
      GeneID = GeneID.split("-")[0]
    start, stop = int(start), int(stop)
    if GeneID in GeneID_to_location:
      previous_chrom, previous_start, previous_stop = GeneID_to_location[GeneID]
      GeneID_to_location[GeneID] = [chrom, min(previous_start, start), max(previous_stop, stop)]
    else:
      GeneID_to_location[GeneID] = [chrom, start, stop]

with open(snakemake.input[0], "r") as input_AS_file, open(snakemake.output[0], "w") as output_diff_file, open(snakemake.output[1], "w") as output_same_file:
	for line in input_AS_file:
		geneid, transcript_location, anchor_range, anchor_AS, direction = line.strip().split()
		gene_chrom, gene_start, gene_stop = GeneID_to_location[geneid]
		anchor_chrom, anchor_startstop = anchor_range.split(":")
		anchor_start, anchor_stop = anchor_startstop.split("-")
		anchor_start, anchor_stop = int(anchor_start), int(anchor_stop)
		if gene_chrom != anchor_chrom or anchor_start > gene_stop + dist_threshold or gene_start > anchor_stop + dist_threshold:
			output_diff_file.write(f"{geneid}\t{transcript_location}\t{anchor_range}\t{direction}\n")
		else:
			output_same_file.write(f"{geneid}\t{transcript_location}\t{anchor_range}\t{direction}\n")
