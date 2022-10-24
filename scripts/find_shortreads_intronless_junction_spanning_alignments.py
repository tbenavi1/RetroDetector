from collections import defaultdict
from itertools import groupby
import re
import sys

junction_overhang = int(sys.argv[1])
insertions_threshold = int(sys.argv[2])
short_read_spanning_alignment_minimum_matching_bp = int(sys.argv[3])
short_read_spanning_alignment_minimum_transcriptomic_insertion_size = int(sys.argv[4])
short_read_spanning_alignment_maximum_transcriptomic_insertion_size = int(sys.argv[5])

def asbin(n):
	"""
	Given an integer corresponding to a sam flag,
	return the corresponding string of bits.
	"""
	return str(bin(n))[2:].zfill(16)

def get_cigar_next(optional_fields):
	"""
	Given the optional fields of a sam alignment,
	return the cigar of the mate, if there is one.
	"""
	has_MC = False
	for optional_field in optional_fields:
		if optional_field.startswith("MC:Z:"):
			has_MC = True
			cigar_next = optional_field[5:]
			break
	if not has_MC:
		cigar_next = "*"
	return cigar_next

def get_num_matches(cigar):
	"""
	Given a cigar string,
	return the number of matches.
	"""
	if cigar == "*":
		return 0
	num_matches = 0
	X = re.finditer("(\d+)([MIDNSHPX=])", cigar)
	for block in X:
		cigar_digits = int(block.group(1))
		cigar_operation = block.group(2)
		if cigar_operation == "M":
			num_matches += cigar_digits
	return num_matches

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

def get_flank_lengths(cigar, side):
	"""
	Given a cigar string and a side "left" or "right",
	return the combined length of any S or H operations on that side.
	"""
	flank_ops = {"S", "H"}
	query_ops = {"M", "I", "=", "X"}
	query_pos = 0
	X = re.finditer("(\d+)([MIDNSHPX=])", cigar)
	finished_left_flank = False
	finished_right_flank = False
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
		else:
			length = 0
	elif side == "right":
		if right_flank_pos <= query_pos:
			length = query_pos - right_flank_pos + 1
		else:
			length = 0
	return length

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

def is_intronless(cigar, allowed_insertions):
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

def process_overlapping_alignment(readname, transcript, pos, cigar, junction):
	junction_is_supported = False
	geneid = transcript_to_geneid[transcript]
	assert pos <= junction - junction_overhang + 1, f"{readname}, {transcript}"
	region_cigar = subset_cigar_string(cigar, pos, junction - junction_overhang + 1, junction + junction_overhang)
	#if this alignment does not have an intron
	#aka, if the cigar has at most insertions_threshold insertions
	if is_intronless(region_cigar, insertions_threshold):
		junction_is_supported = True
		output_spannedjunctions_file.write(f"{geneid}\t{transcript}\t{junction}\t{readname}\toverlapping\t.\n")
	return junction_is_supported

def process_spanning_alignment(readname, transcript, tlen, left_overhang, right_overhang, spanned_overlapped_junctions, spanned_nonoverlapped_junctions):
	geneid = transcript_to_geneid[transcript]
	total_spanned_intron_length = 0
	for junction in spanned_overlapped_junctions:
		spanned_intron_length = transcript_junction_to_intron_length[(transcript, junction)]
		total_spanned_intron_length += spanned_intron_length
	for junction in spanned_nonoverlapped_junctions:
		spanned_intron_length = transcript_junction_to_intron_length[(transcript, junction)]
		total_spanned_intron_length += spanned_intron_length
	expected_genome_tlen = tlen + left_overhang + total_spanned_intron_length + right_overhang
	for junction in spanned_nonoverlapped_junctions:
		output_spannedjunctions_file.write(f"{geneid}\t{transcript}\t{junction}\t{readname}\tnon-overlapping\t{expected_genome_tlen}\n")
	return

processed_paired_alignment_reads = set()
supporting_paired_alignment_reads = set()
def process_paired_alignment(readname, flag, transcript, pos, cigar, pnext, tlen, cigar_next, line):
	#if we have already processed the mate
	if readname in processed_paired_alignment_reads:
		processed_paired_alignment_reads.remove(readname)
		#if the mate was a supporting read
		if readname in supporting_paired_alignment_reads:
			supporting_paired_alignment_reads.remove(readname)
			output_spanningalignments_file.write(line)
	#else we have not processed this read before
	else:
		if pos <= pnext:
			start = pos
			start_inner = pnext
			left_cigar = cigar
			right_cigar = cigar_next
		else:
			start = pnext
			start_inner = pos
			left_cigar = cigar_next
			right_cigar = cigar
		stop_inner = start + ref_len(left_cigar) - 1
		stop = start_inner + ref_len(right_cigar) - 1
		left_overhang = get_flank_lengths(left_cigar, "left")
		right_overhang = get_flank_lengths(right_cigar, "right")
		junctions = transcript_to_junctions[transcript]
		is_supporting_alignment = False
		spanned_overlapped_junctions = []
		spanned_nonoverlapped_junctions = []
		for junction in junctions:
			junction_is_supported = False
			#if left read overlaps junction
			if start + junction_overhang - 1 <= junction <= stop_inner - junction_overhang:
				junction_is_supported = process_overlapping_alignment(readname, transcript, start, left_cigar, junction)
				spanned_overlapped_junctions.append(junction)
			#else if right read overlaps junction
			elif start_inner + junction_overhang - 1 <= junction <= stop - junction_overhang:
				junction_is_supported = process_overlapping_alignment(readname, transcript, start_inner, right_cigar, junction)
				spanned_overlapped_junctions.append(junction)
			#else if paired alignment spans junction
			elif start + junction_overhang - 1 <= junction <= stop - junction_overhang:
				junction_is_supported = True
				spanned_nonoverlapped_junctions.append(junction)
			if junction_is_supported:
				is_supporting_alignment = True
		if spanned_nonoverlapped_junctions:
			process_spanning_alignment(readname, transcript, tlen, left_overhang, right_overhang, spanned_overlapped_junctions, spanned_nonoverlapped_junctions)
		processed_paired_alignment_reads.add(readname)
		if is_supporting_alignment:
			supporting_paired_alignment_reads.add(readname)
			output_spanningalignments_file.write(line)
	return

def process_unpaired_alignment(readname, flag, transcript, pos, cigar, pnext, tlen, cigar_next, line):
	start = pos
	stop = start + ref_len(cigar)  - 1
	junctions = transcript_to_junctions[transcript]
	is_supporting_alignment = False
	for junction in junctions:
		junction_is_supported = False
		#We are interested in the base pairs from junction - junction_overhang + 1 to junction
		#and the base pairs from junction + 1 to junction + junction_overhang
		if start + junction_overhang - 1 <= junction <= stop - junction_overhang:
			junction_is_supported = process_overlapping_alignment(readname, transcript, pos, cigar, junction)
		if junction_is_supported:
			is_supporting_alignment = True
	if is_supporting_alignment:
		output_spanningalignments_file.write(line)
	return

def process_alignment(readname, flag, transcript, pos, cigar, pnext, tlen, cigar_next, line):
	reads_mapped_to_parent_gene = transcript_to_reads_mapped_to_parent_gene[transcript]
	if readname in reads_mapped_to_parent_gene:
		read_is_mapped_to_parent_gene = True
	else:
		read_is_mapped_to_parent_gene = False
	if transcript in transcript_to_junctions:
		transcript_has_exon_exon_junctions = True
	else:
		transcript_has_exon_exon_junctions = False
	#if this alignment could potentially be evidence of a retrogene
	if not read_is_mapped_to_parent_gene and transcript_has_exon_exon_junctions:
		#if alignment has expected tlen
		if short_read_spanning_alignment_minimum_transcriptomic_insertion_size <= abs(tlen) <= short_read_spanning_alignment_maximum_transcriptomic_insertion_size:
			alignment_has_expected_tlen = True
		else:
			alignment_has_expected_tlen = False
		#if alignment is a proper pair
		if asbin(flag)[-2] == '1':
			alignment_is_proper_pair = True
		else:
			alignment_is_proper_pair = False
		num_matches = get_num_matches(cigar)
		num_matches_next = get_num_matches(cigar_next)
		if num_matches >= short_read_spanning_alignment_minimum_matching_bp and num_matches_next >= short_read_spanning_alignment_minimum_matching_bp:
			reads_have_sufficient_matches = True
		else:
			reads_have_sufficient_matches = False
		#if we are treating this as a paired alignment
		if alignment_has_expected_tlen and alignment_is_proper_pair and reads_have_sufficient_matches:
			process_paired_alignment(readname, flag, transcript, pos, cigar, pnext, tlen, cigar_next, line)
		else:
			process_unpaired_alignment(readname, flag, transcript, pos, cigar, pnext, tlen, cigar_next, line)
	return

def calculate_tlen(pos, cigar, pnext, cigar_next):
	if pos <= pnext:
		start = pos
		cigar_ref_len = ref_len(cigar_next)
		stop = pnext + cigar_ref_len - 1
		tlen = stop - start + 1
	else:
		start = pnext
		cigar_ref_len = ref_len(cigar)
		stop = pos + cigar_ref_len - 1
		tlen = start - stop - 1
	return tlen

read_to_alignment = {}
def get_alternate_alignments(line):
	alternate_alignments = []
	readname, flag = line.strip().split()[:2]
	flag = int(flag)
	#if the mate is unmapped
	if asbin(flag)[-4] == "1":
		alternate_hits = []
		optional_fields = line.strip().split()[11:]
		for optional_field in optional_fields:
			if optional_field.startswith("XA:Z:"):
				alternate_hits = optional_field.split(":")[2][:-1].split(";")
				break
		alternate_transcripts = set()
		transcript_to_alternate_hits = defaultdict(list)
		for alternate_hit in alternate_hits:
			transcript = alternate_hit.split(",")[0]
			alternate_transcripts.add(transcript)
			transcript_to_alternate_hits[transcript].append(alternate_hit)
		for transcript in alternate_transcripts:
			alternate_hits = transcript_to_alternate_hits[transcript]
			if len(alternate_hits) == 1:
				alternate_hit = alternate_hits[0]
				#set an "alternate" readname that includes the transcript
				alternate_readname = f"{readname}_{transcript}"
				pos, cigar = alternate_hit.split(",")[1:3]
				strand = pos[0]
				if strand == "+":
					flag = 0
				else:
					#set flag to read reverse strand
					flag = 16
				pos = int(pos[1:])
				pnext = 0
				tlen = 0
				cigar_next = "*"
				alternate_line = f"{alternate_readname}\t{flag}\t{transcript}\t{pos}\t255\t{cigar}\t*\t{pnext}\t{tlen}\t*\t*\n"
				alternate_alignment = (alternate_readname, flag, transcript, pos, cigar, pnext, tlen, cigar_next, alternate_line)
				alternate_alignments.append(alternate_alignment)
	#else if the mate is mapped and this is the first of the pair that we've encountered
	elif readname not in read_to_alignment:
		read_to_alignment[readname] = line
	#else this is the second of the pair that we've encountered
	else:
		line2 = line
		line = read_to_alignment[readname]
		alternate_hits = []
		alternate_hits2 = []
		optional_fields = line.strip().split()[11:]
		optional_fields2 = line2.strip().split()[11:]
		for optional_field in optional_fields:
			if optional_field.startswith("XA:Z:"):
				alternate_hits = optional_field.split(":")[2][:-1].split(";")
				break
		for optional_field2 in optional_fields2:
			if optional_field2.startswith("XA:Z:"):
				alternate_hits2 = optional_field2.split(":")[2][:-1].split(";")
				break
		alternate_transcripts = set()
		transcript_to_alternate_hits = defaultdict(list)
		transcript_to_alternate_hits2 = defaultdict(list)
		for alternate_hit in alternate_hits:
			transcript = alternate_hit.split(",")[0]
			alternate_transcripts.add(transcript)
			transcript_to_alternate_hits[transcript].append(alternate_hit)
		for alternate_hit2 in alternate_hits2:
			transcript2 = alternate_hit2.split(",")[0]
			alternate_transcripts.add(transcript2)
			transcript_to_alternate_hits2[transcript2].append(alternate_hit2)
		for transcript in alternate_transcripts:
			#set an "alternate" readname that has the transcript
			alternate_readname = f"{readname}_{transcript}"
			alternate_hits = transcript_to_alternate_hits[transcript]
			alternate_hits2 = transcript_to_alternate_hits2[transcript]
			if len(alternate_hits) == 1:
				first_is_singly_mapped = True
			else:
				first_is_singly_mapped = False
			if len(alternate_hits2) == 1:
				second_is_singly_mapped = True
			else:
				second_is_singly_mapped = False
			#if both are singly mapped to this transcript
			if first_is_singly_mapped and second_is_singly_mapped:
				alternate_hit = alternate_hits[0]
				alternate_hit2 = alternate_hits2[0]
				pos, cigar = alternate_hit.split(",")[1:3]
				second_pos, second_cigar = alternate_hit2.split(",")[1:3]
				first_strand = pos[0]
				second_strand = second_pos[0]
				#strip off the plus or minus, convert to int
				pos = int(pos[1:])
				second_pos = int(second_pos[1:])
				#if these map to different strands, so still possible a pair
				if first_strand != second_strand:
					if first_strand == "+":
						#set first flag to "read paired", "read mapped in proper pair", "mate reverse strand", "first in pair"
						first_flag = 99
						#set second flag to "read paired", "read mapped in proper pair", "read reverse strand", "second in pair"
						second_flag = 147
					else:
						#set first flag to "read paired", "read mapped in proper pair", "read reverse strand", "first in pair"
						first_flag = 83
						#set second flag to "read paired", "read mapped in proper pair", "mate reverse strand", "second in pair"
						second_flag = 163
					first_tlen = calculate_tlen(pos, cigar, second_pos, second_cigar)
					second_tlen = -first_tlen
					first_alternate_line = f"{alternate_readname}\t{first_flag}\t{transcript}\t{pos}\t255\t{cigar}\t=\t{second_pos}\t{first_tlen}\t*\t*\n"
					second_alternate_line = f"{alternate_readname}\t{second_flag}\t{transcript}\t{second_pos}\t255\t{second_cigar}\t=\t{pos}\t{second_tlen}\t*\t*\n"
					alternate_alignment = (alternate_readname, first_flag, transcript, pos, cigar, second_pos, first_tlen, second_cigar, first_alternate_line)
					alternate_alignments.append(alternate_alignment)
					alternate_alignment = (alternate_readname, second_flag, transcript, second_pos, second_cigar, pos, second_tlen, cigar, second_alternate_line)
					alternate_alignments.append(alternate_alignment)
				#else these are not possibly a pair
				else:
					if first_strand == "+":
						first_flag = 0
					else:
						first_flag = 16
					if second_strand == "+":
						second_flag = 0
					else:
						second_flag = 16
					pnext = 0
					tlen = 0
					cigar_next = "*"
					first_alternate_line = f"{alternate_readname}\t{first_flag}\t{transcript}\t{pos}\t255\t{cigar}\t*\t{pnext}\t{tlen}\t*\t*\n"
					second_alternate_line = f"{alternate_readname}\t{second_flag}\t{transcript}\t{second_pos}\t255\t{second_cigar}\t*\t{pnext}\t{tlen}\t*\t*\n"
					alternate_alignment = (alternate_readname, first_flag, transcript, pos, cigar, pnext, tlen, cigar_next, first_alternate_line)
					alternate_alignments.append(alternate_alignment)
					alternate_alignment = (alternate_readname, second_flag, transcript, second_pos, second_cigar, pnext, tlen, cigar_next, second_alternate_line)
					alternate_alignments.append(alternate_alignment)
			#else if only the first is singly mapped to this transcript
			elif first_is_singly_mapped and not second_is_singly_mapped:
				alternate_hit = alternate_hits[0]
				pos, cigar = alternate_hit.split(",")[1:3]
				strand = pos[0]
				if strand == "+":
					flag = 0
				else:
					#set flag to read reverse strand
					flag = 16
				#strip off the plus or minus, convert to int
				pos = int(pos[1:])
				pnext = 0
				tlen = 0
				cigar_next = "*"
				alternate_line = f"{alternate_readname}\t{flag}\t{transcript}\t{pos}\t255\t{cigar}\t*\t{pnext}\t{tlen}\t*\t*\n"
				alternate_alignment = (alternate_readname, flag, transcript, pos, cigar, pnext, tlen, cigar_next, alternate_line)
				alternate_alignments.append(alternate_alignment)
			#else if only the second is singly mapped to this transcript
			elif not first_is_singly_mapped and second_is_singly_mapped:
				alternate_hit2 = alternate_hits2[0]
				pos, cigar = alternate_hit2.split(",")[1:3]
				strand = pos[0]
				if strand == "+":
					flag = 0
				else:
					#set flag to read reverse strand
					flag = 16
				#strip off the plus or minus, convert to int
				pos = int(pos[1:])
				pnext = 0
				tlen = 0
				cigar_next = "*"
				alternate_line = f"{alternate_readname}\t{flag}\t{transcript}\t{pos}\t255\t{cigar}\t*\t{pnext}\t{tlen}\t*\t*\n"
				alternate_alignment = (alternate_readname, flag, transcript, pos, cigar, pnext, tlen, cigar_next, alternate_line)
				alternate_alignments.append(alternate_alignment)
			#else neither are singly mapped to this transcript
			else:
				pass
		del read_to_alignment[readname]
	return alternate_alignments

transcript_to_junctions = {}
transcript_to_geneid = {}
with open(sys.argv[6], "r") as input_junctions_file:
	for line in input_junctions_file:
		geneid, transcript, *junctions = line.strip().split()
		transcript_to_junctions[transcript] = [int(junction) for junction in junctions]
		transcript_to_geneid[transcript] = geneid

transcript_junction_to_intron_length = {}
with open(sys.argv[7], "r") as input_intronlengths_file:
	for line in input_intronlengths_file:
		geneid, transcript, junction, intron_length = line.strip().split()
		junction, intron_length = int(junction), int(intron_length)
		transcript_junction_to_intron_length[(transcript, junction)] = intron_length

transcript_to_reads_mapped_to_parent_gene = defaultdict(set)
with open(sys.argv[8], "r") as input_overlapping_file:
	for line in input_overlapping_file:
		geneid, transcript, *reads = line.strip().split()
		transcript_to_reads_mapped_to_parent_gene[transcript] = set(reads)

with open(sys.argv[9], "w") as output_spanningalignments_file, open(sys.argv[10], "w") as output_spannedjunctions_file:
	for line in sys.stdin:
		if line.startswith("@"):
			output_spanningalignments_file.write(line)
		else:
			flag = int(line.strip().split()[1])
			if asbin(flag)[-12] == '1':
				alignment_is_supplementary = True
			else:
				alignment_is_supplementary = False
			if not alignment_is_supplementary:
				#process primary alignment
				readname, flag, transcript, pos, _, cigar, _, pnext, tlen, _, _, *optional_fields = line.strip().split()
				flag, pos, pnext, tlen = int(flag), int(pos), int(pnext), int(tlen)
				cigar_next = get_cigar_next(optional_fields)
				process_alignment(readname, flag, transcript, pos, cigar, pnext, tlen, cigar_next, line)
				#process any alternate alignments
				alternate_alignments = get_alternate_alignments(line)
				for alternate_alignment in alternate_alignments:
					readname, flag, transcript, pos, cigar, pnext, tlen, cigar_next, line = alternate_alignment
					process_alignment(readname, flag, transcript, pos, cigar, pnext, tlen, cigar_next, line)
