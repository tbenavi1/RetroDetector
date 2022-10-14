from collections import defaultdict

num_samples = len(snakemake.input)-1
dist_threshold = 10000

GeneID_to_location = {}
with open(snakemake.input[-1], "r") as input_coords_file:
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

geneid_to_diffs = defaultdict(list)
geneid_to_sames = defaultdict(list)
for i in range(num_samples):
	geneid_to_cluster_found = {}
	geneid_to_anchor_range = {}
	geneid_to_best_AS = {}
	with open(snakemake.input[i], "r") as input_AS_file:
		for line in input_AS_file:
			geneid, transcript_range, anchor_range, anchor_AS, strand = line.strip().split()
			if geneid not in geneid_to_cluster_found:
				geneid_to_cluster_found[geneid] = False
			if geneid not in geneid_to_best_AS:
				geneid_to_best_AS[geneid] = 0
			anchor_AS = int(anchor_AS)
			best_AS = geneid_to_best_AS[geneid]
			if anchor_AS > best_AS:
				geneid_to_cluster_found[geneid] = True
				geneid_to_anchor_range[geneid] = anchor_range
				geneid_to_best_AS[geneid] = anchor_AS
	for geneid in geneid_to_best_AS:
		cluster_found = geneid_to_cluster_found[geneid]
		if cluster_found:
			gene_chrom, gene_start, gene_stop = GeneID_to_location[geneid]
			anchor_range = geneid_to_anchor_range[geneid]
			anchor_chrom, anchor_startstop = anchor_range.split(":")
			anchor_start, anchor_stop = anchor_startstop.split("-")
			anchor_start, anchor_stop = int(anchor_start), int(anchor_stop)
			if gene_chrom != anchor_chrom or anchor_start > gene_stop + dist_threshold or gene_start > anchor_stop + dist_threshold:
				diff = True
				same = False
			else:
				same = True
				diff = False
			geneid_to_diffs[geneid].append(diff)
			geneid_to_sames[geneid].append(same)

with open(snakemake.output[0], "w") as output_alldiff_file, open(snakemake.output[1], "w") as output_allsame_file, open(snakemake.output[2], "w") as output_somesome_file:
	for geneid in geneid_to_diffs:
		sames = geneid_to_sames[geneid]
		diffs = geneid_to_diffs[geneid]
		if all(diffs):
			output_alldiff_file.write(f"{geneid}\n")
		elif all(sames):
			output_allsame_file.write(f"{geneid}\n")
		else:
			output_somesome_file.write(f"{geneid}\n")
