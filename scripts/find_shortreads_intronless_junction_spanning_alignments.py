from itertools import groupby
from collections import defaultdict
import re
import sys

junction_overhang = int(sys.argv[1])
insertions_threshold = int(sys.argv[2])
insert_size = 350

def asbin(n):
	return str(bin(n))[2:].zfill(16)

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

def get_flank_lengths(cigar, side):
	flank_ops = {"S", "H"}
	query_ops = {"M", "I", "=", "X"}
	query_pos = 0
	X = re.finditer("(\d+)([MIDNSHPX=])", cigar)
	finished_left_flank = False
	finished_right_flank = False
	#regions = []
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
	if side == "left":
		if left_flank_pos >= 1:
			length = left_flank_pos
			#region = f"1-{left_flank_pos}"
			#regions.append(region)
		else:
			length = 0
	elif side == "right":
		if right_flank_pos <= query_pos:
			length = query_pos - right_flank_pos + 1
			#region = f"{right_flank_pos}-{query_pos}"
			#regions.append(region)
		else:
			length = 0
	return length

transcript_to_junctions = {}
transcript_to_geneid = {}
with open(sys.argv[3], "r") as input_junctions_file:
	for line in input_junctions_file:
		geneid, transcript, *junctions = line.strip().split()
		transcript_to_junctions[transcript] = junctions
		transcript_to_geneid[transcript] = geneid

transcript_to_longjunctions = defaultdict(list)
with open(sys.argv[4], "r") as input_longjunctions_file:
	for line in input_longjunctions_file:
		geneid, transcript, *junctions = line.strip().split()
		transcript_to_longjunctions[transcript] = junctions

transcript_junction_to_intron_length = {}
with open(sys.argv[5], "r") as input_intronlengths_file:
	for line in input_intronlengths_file:
		geneid, transcript, junction, intron_length = line.strip().split()
		junction = int(junction)
		intron_length = int(intron_length)
		transcript_junction_to_intron_length[(transcript, junction)] = intron_length

read_to_alignment = defaultdict(list)

with open(sys.argv[6], "w") as output_spanningalignments_file, open(sys.argv[7], "w") as output_spannedjunctions_file, open("true_postiives.txt", "w") as output_truepositives_file:
	for line in sys.stdin:
		if line.startswith("@"):
			output_spanningalignments_file.write(line)
		else:
			readname, flag, transcript, pos, _, cigar, _, pnext, tlen = line.strip().split()[:9]
			pos = int(pos)
			pnext = int(pnext)
			flag = int(flag)
			tlen = abs(int(tlen))
			if readname == "K00282:232:HN2CLBBXX:6:2123:28148:13763":
				print("hi")
				print(transcript)
			if transcript in transcript_to_junctions:
				if readname == "K00282:232:HN2CLBBXX:6:2123:28148:13763":
					print("hello")
				junctions = transcript_to_junctions[transcript]
				geneid = transcript_to_geneid[transcript]
				is_supporting_alignment = False
				#if this read is part of a proper pair
				if asbin(flag)[-2] == '1':
					optional_fields = line.strip().split()[11:]
					has_MC = False
					for optional_field in optional_fields:
						if optional_field.startswith("MC:Z:"):
							has_MC = True
							cigar_next = optional_field[5:]
					assert has_MC, line
					if pos <= pnext:
						start = pos
						cigar_ref_len = ref_len(cigar_next)
						stop = pnext + cigar_ref_len - 1
						left_overhang = get_flank_lengths(cigar, "left")
						right_overhang = get_flank_lengths(cigar_next, "right")
					else:
						start = pnext
						cigar_ref_len = ref_len(cigar)
						stop = pos + cigar_ref_len - 1
						left_overhang = get_flank_lengths(cigar_next, "left")
						right_overhang = get_flank_lengths(cigar, "right")
					spanned_junctions = []
					total_spanned_intron_length = 0
					long_junctions = transcript_to_longjunctions[transcript]
					is_true_positive = False
					for junction in junctions:
						junction = int(junction)
						#if there is an paired alignment that spans junction_overhang base pairs on both sides of junction
						if start + junction_overhang - 1 <= junction <= stop - junction_overhang:
							spanned_junctions.append(junction)
							if str(junction) in long_junctions:
								is_true_positive = True
							spanned_intron_length = transcript_junction_to_intron_length[(transcript, junction)]
							total_spanned_intron_length += spanned_intron_length
							#output_spannedjunctions_file.write(f"{geneid}\t{transcript}\t{junction}\t{readname}\n")
							#is_supporting_alignment = True
					expected_genome_tlen = tlen + left_overhang + total_spanned_intron_length + right_overhang
					if is_true_positive:
						output_truepositives_file.write(f"{geneid}\t{transcript}\t{tlen}\t{expected_genome_tlen}\t{left_overhang}\t{total_spanned_intron_length}\t{right_overhang}\t{line}")
					#if readname == "K00282:232:HN2CLBBXX:6:1125:30837:12234":
					#if readname == "K00282:232:HN2CLBBXX:6:1102:18609:37361":
					#if readname == "K00282:232:HN2CLBBXX:6:1105:11129:8084":
					if readname == "K00282:232:HN2CLBBXX:6:2123:28148:13763":
						print(f"tlen is {tlen}")
						print(f"left overhang is {left_overhang}")
						print(f"total spanned intron length is {total_spanned_intron_length}")
						print(f"right overhang is {right_overhang}")
						print(f"expected genome tlen is {expected_genome_tlen}")
					#if there was a spanned junction, and the transcript tlen is closer to insert size than expected genome tlen
					if total_spanned_intron_length > 0 and abs(tlen - insert_size) < abs(expected_genome_tlen - insert_size):
						is_supporting_alignment = True
						for spanned_junction in spanned_junctions:
							output_spannedjunctions_file.write(f"{geneid}\t{transcript}\t{spanned_junction}\t{readname}\n")
				else:
					cigar_ref_len = ref_len(cigar)
					for junction in junctions:
						junction = int(junction)
						#if there is an alignment that spans junction_overhang base pairs on both sides of junction
						if pos + junction_overhang - 1 <= junction <= pos + cigar_ref_len - junction_overhang:
							region_cigar = subset_cigar_string(cigar, pos, junction - junction_overhang + 1, junction + junction_overhang)
							#assert region_cigar[-1] == "M", line
							#if this alignment does not have an intron
							#aka, if the cigar has at most insertions_threshold insertions
							if is_intronless(region_cigar, insertions_threshold):
								output_spannedjunctions_file.write(f"{geneid}\t{transcript}\t{junction}\t{readname}\n")
								is_supporting_alignment = True
				if is_supporting_alignment:
					output_spanningalignments_file.write(line)
			#also need to parse alternate hits
			#if this is the first of the pair
			if readname not in read_to_alignment:
				read_to_alignment[readname] = line
			#else this is the second of the pair
			else:
				line2 = read_to_alignment[readname]
				optional_fields = line.strip().split()[11:]
				optional_fields2 = line2.strip().split()[11:]
				alternate_hits = []
				alternate_hits2 = []
				for optional_field in optional_fields:
					#if there are optional hits for this alignment
					#optional_field may look like XA:Z:FBtr0257228,+1385,63M88S,5;FBtr0402045,+1383,63M88S,5;
					if optional_field.startswith("XA:Z:"):
						alternate_hits = optional_field.split(":")[2][:-1].split(";")
						#alternate_hit may look like FBtr0257228,+1385,63M88S,5
						break
				#alternate_hits.sort()
				for optional_field2 in optional_fields2:
					if optional_field2.startswith("XA:Z:"):
						alternate_hits2 = optional_field2.split(":")[2][:-1].split(";")
						break
				#alternate_hits2.sort()
				#assert len(alternate_hits) == len(alternate_hits2), line+line2
				for alternate_hit in alternate_hits:
					transcript, pos, cigar, NM = alternate_hit.split(",")
					found_match = False
					for alternate_hit2 in alternate_hits2:
						transcript2, pnext, cigar_next, NM2 = alternate_hit2.split(",")
						if transcript == transcript2:
							found_match = True
							break
					if not found_match:
						break
					assert transcript == transcript2, line+line2
					#strip off the plus or minus, convert to int
					pos = int(pos[1:])
					pnext = int(pnext[1:])
					if transcript in transcript_to_junctions:
						junctions = transcript_to_junctions[transcript]
						geneid = transcript_to_geneid[transcript]
						is_supporting_alignment = False
						if pos <= pnext:
							start = pos
							cigar_ref_len = ref_len(cigar_next)
							stop = pnext + cigar_ref_len - 1
							left_overhang = get_flank_lengths(cigar, "left")
							right_overhang = get_flank_lengths(cigar_next, "right")
						else:
							start = pnext
							cigar_ref_len = ref_len(cigar)
							stop = pos + cigar_ref_len - 1
							left_overhang = get_flank_lengths(cigar_next, "left")
							right_overhang = get_flank_lengths(cigar, "right")
						tlen = stop - start + 1
						spanned_junctions = []
						total_spanned_intron_length = 0
						long_junctions = transcript_to_longjunctions[transcript]
						is_true_positive = False
						for junction in junctions:
							junction = int(junction)
							#if there is an alignment that spans 10bp on either side of junction
							if start + junction_overhang - 1 <= junction <= stop - junction_overhang:
								spanned_junctions.append(junction)
								if str(junction) in long_junctions:
									is_true_positive = True
								spanned_intron_length = transcript_junction_to_intron_length[(transcript, junction)]
								total_spanned_intron_length += spanned_intron_length
						expected_genome_tlen = tlen + left_overhang + total_spanned_intron_length + right_overhang
						if is_true_positive:
							output_truepositives_file.write(f"{geneid}\t{transcript}\t{tlen}\t{expected_genome_tlen}\t{left_overhang}\t{total_spanned_intron_length}\t{right_overhang}\t{line2}")
							output_truepositives_file.write(f"{geneid}\t{transcript}\t{tlen}\t{expected_genome_tlen}\t{left_overhang}\t{total_spanned_intron_length}\t{right_overhang}\t{line}")
						if total_spanned_intron_length > 0 and abs(tlen - insert_size) < abs(expected_genome_tlen - insert_size):
							is_supporting_alignment = True
							for spanned_junction in spanned_junctions:
								output_spannedjunctions_file.write(f"{geneid}\t{transcript}\t{spanned_junction}\t{readname}\n") #I didn't double this
						if is_supporting_alignment:
							output_spanningalignments_file.write(line2)
							output_spanningalignments_file.write(line)
#								region_cigar = subset_cigar_string(cigar, pos, junction - 9, junction + 10)
#								#if this alignment does not have an intron
#								#aka, if the cigar has at most 10 insertions
#								if is_intronless(region_cigar, 10):
#								#if is_intronless(region_cigar, 2): #changed from 10 to 2 to match long read
#									output_spanningalignments_file.write(line)
#									output_spannedjunctions_file.write(f"{transcript}\t{junction}\t{readname}\n")
				del read_to_alignment[readname]
