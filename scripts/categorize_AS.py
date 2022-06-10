geneids = snakemake.config["geneids"]
num_geneids = len(geneids)
#assert num_geneids == 213
num_samples = 10

GeneID_to_location = {}
with open(snakemake.input[-1], "r") as input_coords_file:
	for line in input_coords_file:
		transcript, GeneID, chrom, start, stop = line.strip().split()
		#for flybase_dyak
		if "-" in GeneID:
			GeneID = GeneID.split("-")[0]
		start, stop = int(start), int(stop)
		if GeneID in GeneID_to_location:
			previous_chrom, previous_start, previous_stop = GeneID_to_location[GeneID]
			GeneID_to_location[GeneID] = [chrom, min(previous_start, start), max(previous_stop, stop)]
		else:
			GeneID_to_location[GeneID] = [chrom, start, stop]

with open(snakemake.output[0], "w") as output_alldiff_file, open(snakemake.output[1], "w") as output_allsame_file, open(snakemake.output[2], "w") as output_somesome_file:
	for i, geneid in enumerate(geneids):
		gene_chrom, gene_start, gene_stop = GeneID_to_location[geneid]
		sames = []
		diffs = []
		for j in range(num_samples):
			best_AS = 0
			cluster_found = False
			with open(snakemake.input[i*num_samples + j], "r") as input_AS_file:
				for line in input_AS_file:
					anchor_range, anchor_AS = line.strip().split()
					anchor_AS = int(anchor_AS)
					if anchor_AS > best_AS:
						cluster_found = True
						anchor_chrom, anchor_startstop = anchor_range.split(":")
						anchor_start, anchor_stop = anchor_startstop.split("-")
						anchor_start, anchor_stop = int(anchor_start), int(anchor_stop)
						best_AS = anchor_AS
				if cluster_found:
					if gene_chrom != anchor_chrom or anchor_start > gene_stop + 10000 or gene_start > anchor_stop + 10000:
						diff = True
						same = False
					else:
						same = True
						diff = False
					diffs.append(diff)
					sames.append(same)
		if all(diffs):
			output_alldiff_file.write(f"{geneid}\n")
		elif all(sames):
			output_allsame_file.write(f"{geneid}\n")
		else:
			output_somesome_file.write(f"{geneid}\n")
