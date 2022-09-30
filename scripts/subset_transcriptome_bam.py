from collections import defaultdict

read_threshold = snakemake.config["read_threshold"]

transcript_junction_to_reads = defaultdict(set)
with open(snakemake.input[0], "r") as input_spannedjunctions_file:
	for line in input_spannedjunctions_file:
		geneid, transcript, junction, readname = line.strip().split()
		transcript_junction_to_reads[(transcript, junction)].add(readname)

supporting_reads = set()
supported_junctions = []

for (transcript, junction) in transcript_junction_to_reads:
	reads = transcript_junction_to_reads[(transcript, junction)]
	num_reads = len(reads)
	if num_reads >= read_threshold:
		supporting_reads.update(reads)
		supported_junctions.append((transcript, junction))

with open(snakemake.output[0], "w") as output_readnames_file:
	for read in supporting_reads:
		output_readnames_file.write(f"{read}\n")

with open(snakemake.output[1], "w") as output_subset_file:
	#write headers
	with open(snakemake.input[1], "r") as input_spanningalignments_sam:
		for line in input_spanningalignments_sam:
			if line.startswith("@"):
				output_subset_file.write(line)
			else:
				break
	#write alignments
	with open(snakemake.input[2], "r") as input_spanningalignments_txt:
		for line in input_spanningalignments_txt:
			geneid, transcript, junctionlist, *alignment = line.strip().split()
			alignment_supported = False
			for junction in junctionlist.split(","):
				if (transcript, junction) in supported_junctions:
					alignment_supported = True
			if alignment_supported:
				alignment = "\t".join(alignment)
				output_subset_file.write(f"{alignment}\n")
