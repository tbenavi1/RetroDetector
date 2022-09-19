from collections import defaultdict
import sys

distance_threshold = 1000 #changed to 100 from 10000 since are locations are now more specific

cluster_chrom = ""
previous_chrom = ""
previous_transcript = ""
previous_pos = 0
qname_to_AS = defaultdict(int)
#with open(sys.argv[1], "r") as input_sam_file, open(sys.argv[2], "w") as output_AS_file:
with open(snakemake.input[0], "r") as input_AS_file, open(snakemake.output[0], "w") as output_AS_file:
	for line in input_AS_file:
		geneid, transcript, qname, read_start, read_stop, transcript, transcript_start, transcript_stop, chrom, ref_start, ref_stop, AS = line.strip().split()
		ref_start, ref_stop, AS = int(ref_start), int(ref_stop), int(AS)
		transcript_start, transcript_stop = int(transcript_start), int(transcript_stop)
		if ref_start > ref_stop:
			ref_start, ref_stop = ref_stop, ref_start
		#if this is the same as previous cluster
		if transcript == previous_transcript and chrom == previous_chrom and ref_stop <= previous_pos + distance_threshold:
			#print("hi")
			cluster_stop = ref_stop #pos
			cluster_transcript_start = min(cluster_transcript_start, transcript_start)
			cluster_transcript_stop = max(cluster_transcript_stop, transcript_stop)
			previous_AS = qname_to_AS[qname]
			qname_to_AS[qname] = max(previous_AS, AS)
		#else this is a new cluster
		else:
			#write previous cluster
			cluster_AS = 0
			#print(previous_chrom, chrom, previous_pos, pos)
			for qname in qname_to_AS:
				read_AS = qname_to_AS[qname]
				cluster_AS += read_AS
			if cluster_chrom:
				output_AS_file.write(f"{cluster_geneid}\t{cluster_transcript}:{cluster_transcript_start}-{cluster_transcript_stop}\t{cluster_chrom}:{cluster_start}-{cluster_stop}\t{cluster_AS}\n")
			#initialize new cluster
			qname_to_AS = defaultdict(int)
			cluster_chrom = chrom
			cluster_geneid = geneid
			cluster_transcript = transcript
			cluster_start = ref_start #pos
			cluster_stop = ref_stop #pos
			cluster_transcript_start = transcript_start
			cluster_transcript_stop = transcript_stop
			qname_to_AS[qname] = AS
		previous_chrom = chrom
		previous_transcript = transcript
		previous_pos = ref_stop #pos
	#write last cluster
	cluster_AS = 0
	for qname in qname_to_AS:
		read_AS = qname_to_AS[qname]
		cluster_AS += read_AS
	if cluster_chrom:
		output_AS_file.write(f"{cluster_geneid}\t{cluster_transcript}:{cluster_transcript_start}-{cluster_transcript_stop}\t{cluster_chrom}:{cluster_start}-{cluster_stop}\t{cluster_AS}\n")
