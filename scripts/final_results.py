import gzip
from collections import defaultdict

ref = snakemake.wildcards.ref
sample = snakemake.wildcards.sample

output_prefix = f"RESULTS/{ref}/{sample}/long/consensus/{ref}.{sample}.retrogene."
output_needle_suffix = ".needle"
output_consensus_suffix = ".consensus.fasta"

transcript_to_genome_location = {}
geneid_to_chrom = {}
geneid_to_start = defaultdict(int)
geneid_to_stop = defaultdict(int)
with open(snakemake.input[0], "r") as input_coords_file:
	for line in input_coords_file:
			geneid, transcript, chrom, start, stop = line.strip().split()
			geneid_to_chrom[geneid] = chrom
			start, stop = int(start), int(stop)
			previous_start = geneid_to_start[geneid]
			if previous_start > 0:
				start = min(previous_start, start)
			previous_stop = geneid_to_stop[geneid]
			stop = max(previous_stop, stop)
			geneid_to_start[geneid] = start
			geneid_to_stop[geneid] = stop
			transcript_to_genome_location[transcript] = f"{chrom}:{start}-{stop}"

transcript_to_length = {}
with gzip.open(snakemake.input[1], "rt") as input_transcriptome_file:
	transcript_length = 0
	for line in input_transcriptome_file:
		if line.startswith(">"):
			if transcript_length > 0:
				transcript_to_length[transcript] = transcript_length
			transcript = line.strip().split()[0][1:]
			transcript_length = 0
		else:
			line_length = len(line.strip())
			transcript_length += line_length

def get_overlapping_geneids(retrogene_genome_location):
	overlapping_geneids = []
	retrogene_chrom, retrogene_span = retrogene_genome_location.split(":")
	retrogene_start, retrogene_stop = retrogene_span.split("-")
	retrogene_start, retrogene_stop = int(retrogene_start), int(retrogene_stop)
	for geneid in geneid_to_chrom:
		geneid_chrom = geneid_to_chrom[geneid]
		geneid_start = geneid_to_start[geneid]
		geneid_stop = geneid_to_stop[geneid]
		if geneid_chrom == retrogene_chrom and geneid_start <= retrogene_stop and geneid_stop >= retrogene_start:
			overlapping_geneids.append(geneid)
	return overlapping_geneids

with open(snakemake.input[2], "r") as input_clustered_file, open(snakemake.output[0], "w") as output_results_file:
	for line in input_clustered_file:
		geneid, retrogene_transcriptome_location, retrogene_genome_location, strand = line.strip().split()
		transcript = retrogene_transcriptome_location.split(":")[0]
		transcript_genome_location = transcript_to_genome_location[transcript]
		transcript_length = transcript_to_length[transcript]
		output_consensus = output_prefix + geneid + output_consensus_suffix
		with open(output_consensus, "r") as input_consensus_file:
			retrogene_length = 0
			for line in input_consensus_file:
				if not line.startswith(">"):
					line_length = len(line.strip())
					retrogene_length += line_length
		coverage_pct = retrogene_length/transcript_length
		output_needle = output_prefix + geneid + output_needle_suffix
		try:
			with open(output_needle, "r") as input_needle_file:
				pct_identity = input_needle_file.readlines()[23].split("(")[1].split(")")[0][:-1]
		except:
			pct_identity = "."
		overlapping_geneids = get_overlapping_geneids(retrogene_genome_location)
		overlapping_geneids = ",".join(overlapping_geneids)
		output_results_file.write(f"{geneid}\t{transcript_genome_location}\t{retrogene_genome_location}\t{retrogene_transcriptome_location}\t{strand}\t{retrogene_length}\t{coverage_pct:.2f}\t{pct_identity}\t{overlapping_geneids}\n")
