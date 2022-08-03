from collections import defaultdict

long_samples = snakemake.params["long_samples"]
sample = snakemake.wildcards["sample"]
sample_index = long_samples.index(sample)

GeneIDs_of_interest = set()

with open(snakemake.input[0], "r") as input_SFS:
	for line in input_SFS:
		GeneID, *genotypes = line.strip().split()
		genotype = genotypes[sample_index]
		if genotype == "1":
			GeneIDs_of_interest.add(GeneID)

transcript_to_GeneID = {}

with open(snakemake.input[1], "r") as input_coords:
	for line in input_coords:
		GeneID, transcript, chrom, start, stop = line.strip().split()
		if GeneID in GeneIDs_of_interest:
			transcript_to_GeneID[transcript] = GeneID

read_to_GeneIDs = defaultdict(set)
with open(snakemake.input[2], "r") as input_transcript_sam, open(snakemake.output[0], "w") as output_file:
	for line in input_transcript_sam:
		if not line.startswith("@"):
			read, _, transcript = line.strip().split()[:3]
			if transcript in transcript_to_GeneID:
				GeneID = transcript_to_GeneID[transcript]
				read_to_GeneIDs[read].add(GeneID)
				output_file.write(f"{GeneID}\t{line}")

with open(snakemake.input[3], "r") as input_genome_sam, open(snakemake.output[1], "w") as output_file:
	for line in input_genome_sam:
		if not line.startswith("@"):
			read = line.strip().split()[0]
			if read in read_to_GeneIDs:
				GeneIDs = read_to_GeneIDs[read]
				for GeneID in GeneIDs:
					output_file.write(f"{GeneID}\t{line}")
