from collections import defaultdict
import re
from itertools import groupby

def asbin(n):
	return str(bin(n))[2:].zfill(16)

def query_len(cigar_string):
  """
  Given a CIGAR string, return the number of
  bases consumed in the query sequence.
  """
  query_consuming_ops = {"M", "I", "S", "H", "=", "X"}
  result = 0
  cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
  for _, length_digits in cig_iter:
    length = int(''.join(length_digits))
    op = next(next(cig_iter)[1])
    if op in query_consuming_ops:
      result += length
  return result

def get_read_alignment_positions(cigar):
  flank_ops = {"S", "H"}
  query_ops = {"M", "I", "=", "X"}
  query_pos = 1
  X = re.finditer("(\d+)([MIDNSHPX=])", cigar)
  finished_left_flank = False
  finished_right_flank = False
  regions = []
  for block in X:
    cigar_digits = int(block.group(1))
    cigar_operation = block.group(2)
    if cigar_operation in flank_ops:
      if finished_left_flank and not finished_right_flank:
        right_pos = query_pos - 1
        finished_right_flank = True
      query_pos += cigar_digits
    else:
      if not finished_left_flank:
        finished_left_flank = True
        left_pos = query_pos
      if cigar_operation in query_ops:
        query_pos += cigar_digits
  if not finished_right_flank:
    right_pos = query_pos - 1
  return (left_pos, right_pos)

def convert_query_to_ref(full_cigar, read_reversed, alignment_pos, query_start, query_stop):
  if read_reversed:
    read_length = query_len(full_cigar)
    query_start, query_stop = read_length + 1 - query_start, read_length + 1 - query_stop
  if query_start <= query_stop:
    ref_reversed = False
  else:
    ref_reversed = True
    query_start, query_stop = query_stop, query_start
  consumes_ref = {"M", "D", "N", "=", "X"}
  consumes_query = {"M", "I", "S", "H", "=", "X"}
  X = re.finditer("(\d+)([MIDNSHPX=])", full_cigar)
  previous_ref_block_stop = alignment_pos - 1
  previous_query_block_stop = 1 - 1
  overlapping_region = False
  for block in X:
    cigar_digits = int(block.group(1))
    cigar_operation = block.group(2)
    if cigar_operation in consumes_ref:
      ref_block_start = previous_ref_block_stop + 1
      ref_block_stop = ref_block_start + cigar_digits - 1
    else:
      ref_block_stop = previous_ref_block_stop
    if cigar_operation in consumes_query:
      query_block_start = previous_query_block_stop + 1
      query_block_stop = query_block_start + cigar_digits - 1
      #if we overlap the region of interest for the first time
      if query_block_stop >= query_start and not overlapping_region:
        overlapping_region = True
        if cigar_operation in consumes_ref:
          ref_start = (query_start - query_block_start) + ref_block_start
        else:
          ref_start = previous_ref_block_stop
      #if we have finished the region of interest
      if query_block_stop >= query_stop:
        if cigar_operation in consumes_ref:
          ref_stop = (query_stop - query_block_start) + ref_block_start
        else:
          ref_stop = previous_ref_block_stop
        if ref_reversed:
          ref_start, ref_stop = ref_stop, ref_start
        return (ref_start, ref_stop)
    else:
      query_block_stop = previous_query_block_stop
    previous_ref_block_stop = ref_block_stop
    previous_query_block_stop = query_block_stop

read_geneid_to_positions = defaultdict(list)
with open(snakemake.input[0], "r") as input_transcript_sam_file:
	for line in input_transcript_sam_file:
		if not line.startswith("@"):
			geneid, qname, flag, chrom, pos, mapq, cigar = line.strip().split()[:7]
			flag = int(flag)
			read_start, read_stop = get_read_alignment_positions(cigar)
			if asbin(flag)[-5] == '1':
				read_length = query_len(cigar)
				read_start, read_stop = read_length + 1 - read_start, read_length + 1 - read_stop
			read_geneid_to_positions[(qname, geneid)].append((read_start, read_stop))

with open(snakemake.input[1], "r") as input_genome_sam_file, open(snakemake.output[0], "w") as output_file:
	for line in input_genome_sam_file:
		if not line.startswith("@"):
			geneid, qname, flag, chrom, pos, mapq, cigar = line.strip().split()[:7]
			pos = int(pos)
			flag = int(flag)
			try:
				AS = int(line.strip().split()[13][5:])
			except:
				continue
			read_positions = read_geneid_to_positions[(qname, geneid)]
			for read_position in read_positions:
				read_start, read_stop = read_position
				if asbin(flag)[-5] == '1':
					read_reversed = True
				else:
					read_reversed = False
				ref_start, ref_stop = convert_query_to_ref(cigar, read_reversed, pos, read_start, read_stop)
				output_file.write(f"{geneid}\t{qname}\t{read_start}\t{read_stop}\t{chrom}\t{ref_start}\t{ref_stop}\t{AS}\n")
