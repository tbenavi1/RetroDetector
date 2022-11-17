from collections import defaultdict

samples = []
for sample in snakemake.config["fastqs"]:
	techs = snakemake.config["fastqs"][sample]
	if "pacbio_hifi" in techs or "pacbio_clr" in techs or "ont" in techs:
		samples.append(sample)

num_samples = len(samples)

GeneID_to_samples = defaultdict(set)
multiply_spanned_in_at_least_one_sample = set()

for i, sample in enumerate(samples):
	with open(snakemake.input[i], "r") as input_junctions_file:
		for line in input_junctions_file:
			GeneID, transcript, num_junctions, num_junctions_spanned_by_threshold = line.strip().split()
			num_junctions_spanned_by_threshold = int(num_junctions_spanned_by_threshold)
			if num_junctions_spanned_by_threshold >= 2:
				multiply_spanned_in_at_least_one_sample.add(GeneID)
			if num_junctions_spanned_by_threshold >= 1:
				GeneID_to_samples[GeneID].add(sample)

with open(snakemake.output[0], "w") as output_freqtable_file, open(snakemake.output[1], "w") as output_freqtable_singly_file, open(snakemake.output[2], "w") as output_freqtable_multiply_file:
	for GeneID in sorted(GeneID_to_samples):
		samples_at_threshold = GeneID_to_samples[GeneID]
		num_samples_at_threshold = 0
		for sample in samples:
			if sample in samples_at_threshold:
				num_samples_at_threshold += 1
		freq_vector = "\t".join("1" if sample in samples_at_threshold else "0" for sample in samples)
		output_freqtable_file.write(f"{GeneID}\t{freq_vector}\t{num_samples_at_threshold}\n")
		if GeneID not in multiply_spanned_in_at_least_one_sample:
			output_freqtable_singly_file.write(f"{GeneID}\t{freq_vector}\t{num_samples_at_threshold}\n")
		else:
			output_freqtable_multiply_file.write(f"{GeneID}\t{freq_vector}\t{num_samples_at_threshold}\n")
