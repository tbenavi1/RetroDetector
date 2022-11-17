from collections import defaultdict

transcript_to_GeneID = {}

with open(snakemake.input[0], "r") as input_coords:
	for line in input_coords:
		GeneID, transcript, chrom, start, stop = line.strip().split()
		transcript_to_GeneID[transcript] = GeneID

read_to_GeneIDs = defaultdict(set)
with open(snakemake.input[1], "r") as input_transcriptome_sam, open(snakemake.output[0], "w") as output_file:
	for line in input_transcriptome_sam:
		if not line.startswith("@"):
			read, _, transcript = line.strip().split()[:3]
			GeneID = transcript_to_GeneID[transcript]
			read_to_GeneIDs[read].add(GeneID)
			output_file.write(f"{GeneID}\t{line}")

with open(snakemake.input[2], "r") as input_genome_sam, open(snakemake.output[1], "w") as output_file:
	for line in input_genome_sam:
		if not line.startswith("@"):
			read = line.strip().split()[0]
			GeneIDs = read_to_GeneIDs[read]
			for GeneID in GeneIDs:
				output_file.write(f"{GeneID}\t{line}")
