from collections import defaultdict

long_samples = snakemake.params["long_samples"]

geneids = set()
geneid_to_samples = defaultdict(set)

for i, sample in enumerate(long_samples):
	with open(snakemake.input[i], "r") as input_AS_file:
		for line in input_AS_file:
			geneid, transcript, anchor_range = line.strip().split()
			geneids.add(geneid)
			geneid_to_samples[geneid].add(sample)

with open(snakemake.output[0], "w") as output_file:
	samples_text = "\t".join(long_samples)
	output_file.write(f"GeneID\t{samples_text}\n")
	for geneid in sorted(geneids):
		geneid_samples = geneid_to_samples[geneid]
		genotypes = "\t".join("1" if sample in geneid_samples else "0" for sample in long_samples)
		output_file.write(f"{geneid}\t{genotypes}\n")
