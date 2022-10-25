from itertools import groupby
import re
import sys

junction_overhang = int(sys.argv[1])
insertions_threshold = int(sys.argv[2])

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

def is_intronless(cigar, allowed_insertions = 10):
	"""
	This function will return true if there are <=
	allowed_insertions insertions in the cigar string
	and false otherwise.
	"""
	number_insertions = 0
	X = re.finditer("(\d+)([MIDNSHPX=])", cigar)
	for block in X:
		cigar_digits = int(block.group(1))
		cigar_operation = block.group(2)
		if cigar_operation == "I":
			number_insertions += cigar_digits
			if number_insertions > allowed_insertions:
				return False
	return True

def get_flank_regions(cigar):
	flank_ops = {"S", "H"}
	query_ops = {"M", "I", "=", "X"}
	query_pos = 0
	X = re.finditer("(\d+)([MIDNSHPX=])", cigar)
	finished_left_flank = False
	finished_right_flank = False
	regions = []
	for block in X:
		cigar_digits = int(block.group(1))
		cigar_operation = block.group(2)
		if cigar_operation in flank_ops:
			if finished_left_flank and not finished_right_flank:
				right_flank_pos = query_pos + 1
				finished_right_flank = True
			query_pos += cigar_digits
		else:
			if not finished_left_flank:
				finished_left_flank = True
				left_flank_pos = query_pos
			if cigar_operation in query_ops:
				query_pos += cigar_digits
	if not finished_right_flank:
		right_flank_pos = query_pos + 1
	if left_flank_pos >= 1:
		region = f"1-{left_flank_pos}"
		regions.append(region)
	if right_flank_pos <= query_pos:
		region = f"{right_flank_pos}-{query_pos}"
		regions.append(region)
	return regions

transcript_to_junctions = {}
transcript_to_geneid = {}
with open(sys.argv[3], "r") as input_junctions_file:
	for line in input_junctions_file:
		geneid, transcript, *junctions = line.strip().split()
		transcript_to_junctions[transcript] = junctions
		transcript_to_geneid[transcript] = geneid

with open(sys.argv[4], "w") as output_spanningalignments_file, open(sys.argv[5], "w") as output_spanningalignments_txtfile, open(sys.argv[6], "w") as output_spannedjunctions_file, open(sys.argv[7], "w") as output_flankregions_file:
	for line in sys.stdin:
		if line.startswith("@"):
			output_spanningalignments_file.write(line)
		else:
			readname, _, transcript, pos, _, cigar = line.strip().split()[:6]
			pos = int(pos)
			if transcript in transcript_to_junctions:
				junctions = transcript_to_junctions[transcript]
				geneid = transcript_to_geneid[transcript]
				cigar_ref_len = ref_len(cigar)
				is_supporting_alignment = False
				spanned_junctions = []
				for junction in junctions:
					junction = int(junction)
					#if there is an alignment that spans junction_overhang base pairs on both sides of junction
					if pos + junction_overhang - 1 <= junction <= pos + cigar_ref_len - junction_overhang:
						region_cigar = subset_cigar_string(cigar, pos, junction - junction_overhang + 1, junction + junction_overhang)
						#if this alignment does not have an intron
						#aka, if the cigar has at most insertions_threshold insertions
						if is_intronless(region_cigar, insertions_threshold):
							output_spannedjunctions_file.write(f"{geneid}\t{transcript}\t{junction}\t{readname}\n")
							is_supporting_alignment = True
							spanned_junctions.append(str(junction))
				if is_supporting_alignment:
					output_spanningalignments_file.write(line)
					spanned_junctions = ",".join(spanned_junctions)
					output_spanningalignments_txtfile.write(f"{geneid}\t{transcript}\t{spanned_junctions}\t{line}")
					flank_regions = get_flank_regions(cigar)
					for flank_region in flank_regions:
						read_flank = f"{readname}:{flank_region}"
						output_flankregions_file.write(f"{geneid}\t{transcript}\t{read_flank}\n")
