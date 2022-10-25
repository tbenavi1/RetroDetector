from statistics import mean
from collections import defaultdict

geneids = []
with open(snakemake.input[0], "r") as input_geneids_file:
	for line in input_geneids_file:
		geneid = line.strip()
		geneids.append(geneid)

geneid_intron_to_percents = defaultdict(list)
with open(snakemake.input[1], "r") as input_intron_file:
	for line in input_intron_file:
		if not line.startswith("Read"):
			geneid, transcript, intron, readname, missing_percent, cigar = line.strip().split()
			if geneid in geneids:
				missing_percent = float(missing_percent)
				geneid_intron_to_percents[(geneid, intron)].append(missing_percent)

geneid_to_averages = defaultdict(list)
for geneid, intron in geneid_intron_to_percents:
	percents = geneid_intron_to_percents[(geneid, intron)]
	average_missing = mean(percents)
	geneid_to_averages[geneid].append(average_missing)

with open(snakemake.output[0], "w") as output_file, open(snakemake.output[1], "w") as output_weirdos_file:
	for geneid in geneids:
		#if none of the read alignments span across the entire intron
		if geneid not in geneid_to_averages:
			output_file.write(f"{geneid}\tinconclusive\n")
			continue
		averages = geneid_to_averages[geneid]
		if all(average >= 0.5 for average in averages):
			output_file.write(f"{geneid}\tmissing\n")
		elif all(average < 0.5 for average in averages):
			output_file.write(f"{geneid}\tpresent\n")
		else:
			output_weirdos_file.write(f"{geneid} is inconsistent.\n")
