from collections import defaultdict

#GeneID_of_interest = snakemake.params[0]

transcript_to_indices = defaultdict(set)
transcript_to_geneid = {}
with open(snakemake.input[0], "r") as input_coverage_file:
	for line in input_coverage_file:
		transcript, GeneID, *junctions = line.strip().split()
		transcript_to_geneid[transcript] = GeneID
		#if GeneID == GeneID_of_interest:
		num_junctions = len(junctions)//2
		for i in range(num_junctions):
			junction = int(junctions[2*i])
			coverage = int(junctions[2*i+1])
			if coverage >= 3:
				transcript_to_indices[transcript].add((i, junction))

spanned_introns = set()
with open(snakemake.input[1], "r") as input_introns_file:
	for line in input_introns_file:
		transcript, *introns = line.strip().split()
		if transcript in transcript_to_indices:
			geneid = transcript_to_geneid[transcript]
			indices = transcript_to_indices[transcript]
			for i, junction in indices:
				intron = introns[i]
				spanned_introns.add((geneid, transcript, junction, intron))

with open(snakemake.output[0], "w") as output_introns_file:
	for intron_info in spanned_introns:
		geneid, transcript, junction, intron = intron_info
		output_introns_file.write(f"{geneid}\t{transcript}\t{junction}\t{intron}\n")
