from collections import defaultdict
import re
from itertools import groupby

def ref_len(cigar_string):
	"""
	Given a CIGAR string, return the number of
	bases consumed in the reference sequence.
	"""
	ref_consuming_ops = {"M", "D", "N", "=", "X"}
	result = 0
	cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
	for _, length_digits in cig_iter:
		length = int(''.join(length_digits))
		op = next(next(cig_iter)[1])
		if op in ref_consuming_ops:
			result += length
	return result

def subset_cigar_string(full_cigar, alignment_pos, ref_start, ref_stop):
	"""
	This function will return the cigar string that overlaps
	the reference coordinates from ref_start to ref_stop inclusive.
	For now we assume that the alignment fully covers this region.
	This function is using the 1-based coordinate system.
	So, we assume that alignment_pos is 1-based (as it should be from SAM), 
	and that ref_start and ref_stop are 1-based. 
	"""
	consumes_ref = {"M", "D", "N", "=", "X"}
	assert alignment_pos <= ref_start
	X = re.finditer("(\d+)([MIDNSHPX=])", full_cigar)
	previous_block_stop = alignment_pos - 1
	subset_cigar = ""
	overlapping_region = False
	for block in X:
		cigar_digits = int(block.group(1))
		cigar_operation = block.group(2)
		if cigar_operation in consumes_ref:
			block_start = previous_block_stop + 1
			block_stop = block_start + cigar_digits - 1
			#if we overlap the region of interest
			if block_stop >= ref_start:
				overlapping_region = True
				overlap_start = max(ref_start, block_start)
				overlap_stop = min(ref_stop, block_stop)
				overlap_length = overlap_stop - overlap_start + 1
				subset_cigar += f"{overlap_length}{cigar_operation}"
			#if we have finished the region of interest
			if block_stop >= ref_stop:
				break
		else:
			if overlapping_region:
				subset_cigar += f"{cigar_digits}{cigar_operation}"
			block_stop = previous_block_stop
		previous_block_stop = block_stop
	return subset_cigar

transcript_junction_to_intron = {}
with open(snakemake.input[0], "r") as input_intron_file:
	for line in input_intron_file:
		geneid, transcript, junction, intron = line.strip().split()
		transcript_junction_to_intron[(transcript, junction)] = intron

readname_to_introns = defaultdict(set)
intron_to_transcripts = defaultdict(set)
with open(snakemake.input[1], "r") as input_spannedjunctions_file:
	for line in input_spannedjunctions_file:
		geneie, transcript, junction, readname = line.strip().split()
		if (transcript, junction) in transcript_junction_to_intron:
			chrom = "_".join(transcript.split("|")[1].split("_")[:2])
			intron = transcript_junction_to_intron[(transcript, junction)]
			#readname_to_introns[readname].add(f"{chrom}:{intron}")
			readname_to_introns[readname].add(intron)
			#intron_to_transcripts[f"{chrom}:{intron}"].add(transcript)
			intron_to_transcripts[intron].add(transcript)

def cigar_missing_percent(cigar):
	consumes_ref = {"M", "D", "N", "=", "X"}
	X = re.finditer("(\d+)([MIDNSHPX=])", cigar)
	numerator = 0
	denominator = 0
	for block in X:
		cigar_digits = int(block.group(1))
		cigar_operation = block.group(2)
		if cigar_operation in consumes_ref:
			denominator += cigar_digits
			if cigar_operation in {"N", "D"}:
				numerator += cigar_digits
	missing_percent = numerator/denominator
	return missing_percent

with open(snakemake.input[2], "r") as input_subsetsam_file, open(snakemake.output[0], "w") as output_file:
	for line in input_subsetsam_file:
		if not line.startswith("@"):
			geneid, readname, _, chrom, pos, _, cigar = line.strip().split()[:7]
			pos = int(pos)
			introns = readname_to_introns[readname]
			cigar_ref_len = ref_len(cigar)
			for intron in introns:
				intron_chrom = intron.split(":")[0]
				if intron_chrom == chrom:
					intron_start, intron_stop = intron.split(":")[1].split("-")
					intron_start, intron_stop = int(intron_start), int(intron_stop)
					#if the read spans across the entire intron
					if pos <= intron_start and intron_stop <= pos + cigar_ref_len:
						transcripts = intron_to_transcripts[intron]
						region_cigar = subset_cigar_string(cigar, pos, intron_start, intron_stop)
						missing_percent = cigar_missing_percent(region_cigar)
						for transcript in transcripts:
							output_file.write(f"{geneid}\t{transcript}\t{intron}\t{readname}\t{missing_percent}\t{region_cigar}\n")
					else:
						output_file.write("Read doesn't span across entire intron.\n")
