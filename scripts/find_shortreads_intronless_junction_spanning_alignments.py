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
		#if cigar_operation in ["I", "D", "N", "S"]: #changed to this to match long read, also removed "X"
			number_insertions += cigar_digits
			if number_insertions > allowed_insertions:
				return False
	return True

transcript_to_junctions = {}
transcript_to_geneid = {}
with open(sys.argv[3], "r") as input_junctions_file:
	for line in input_junctions_file:
		geneid, transcript, *junctions = line.strip().split()
		transcript_to_junctions[transcript] = junctions
		transcript_to_geneid[transcript] = geneid

with open(sys.argv[4], "w") as output_spanningalignments_file, open(sys.argv[5], "w") as output_spannedjunctions_file:
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
				if is_supporting_alignment:
					output_spanningalignments_file.write(line)
#			#also need to parse alternate hits
#			optional_fields = line.strip().split()[11:]
#			for optional_field in optional_fields:
#				#if there are optional hits for this alignment
#				#optional_field may look like XA:Z:FBtr0257228,+1385,63M88S,5;FBtr0402045,+1383,63M88S,5;
#				if optional_field.startswith("XA:Z:"):
#					alternate_hits = optional_field.split(":")[2][:-1].split(";")
#					#alternate_hit may look like FBtr0257228,+1385,63M88S,5
#					for alternate_hit in alternate_hits:
#						transcript, pos, cigar, NM = alternate_hit.split(",")
#						#strip off the plus or minus, convert to int
#						pos = int(pos[1:])
#						if transcript in transcript_to_junctions:
#							junctions = transcript_to_junctions[transcript]
#							cigar_ref_len = ref_len(cigar)
#							for junction in junctions:
#								junction = int(junction)
#								#if there is an alignment that spans 10bp on either side of junction
#								if pos + 9 <= junction <= pos + cigar_ref_len - 10:
#									region_cigar = subset_cigar_string(cigar, pos, junction - 9, junction + 10)
#									#if this alignment does not have an intron
#									#aka, if the cigar has at most 10 insertions
#									if is_intronless(region_cigar, 10):
#									#if is_intronless(region_cigar, 2): #changed from 10 to 2 to match long read
#										output_spanningalignments_file.write(line)
#										output_spannedjunctions_file.write(f"{transcript}\t{junction}\t{readname}\n")
