from collections import defaultdict
import sys

distance_threshold = 10000 #changed to 100 from 10000 since are locations are now more specific

cluster_chrom = ""
previous_chrom = ""
previous_geneid = ""
previous_pos = 0
qname_to_AS = defaultdict(int)
#with open(sys.argv[1], "r") as input_sam_file, open(sys.argv[2], "w") as output_AS_file:
with open(snakemake.input[0], "r") as input_AS_file, open(snakemake.output[0], "w") as output_AS_file:
	for line in input_AS_file:
		geneid, qname, read_start, read_stop, chrom, ref_start, ref_stop, AS = line.strip().split()
		ref_start, ref_stop, AS = int(ref_start), int(ref_stop), int(AS)
		if ref_start > ref_stop:
			ref_start, ref_stop = ref_stop, ref_start
		#if this is the same as previous cluster
		if geneid == previous_geneid and chrom == previous_chrom and ref_stop <= previous_pos + distance_threshold:
			#print("hi")
			cluster_stop = ref_stop #pos
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
				output_AS_file.write(f"{cluster_geneid}\t{cluster_chrom}:{cluster_start}-{cluster_stop}\t{cluster_AS}\n")
			#initialize new cluster
			qname_to_AS = defaultdict(int)
			cluster_chrom = chrom
			cluster_geneid = geneid
			cluster_start = ref_start #pos
			cluster_stop = ref_stop #pos
			qname_to_AS[qname] = AS
		previous_chrom = chrom
		previous_geneid = geneid
		previous_pos = ref_stop #pos
	#write last cluster
	cluster_AS = 0
	for qname in qname_to_AS:
		read_AS = qname_to_AS[qname]
		cluster_AS += read_AS
	if cluster_chrom:
		output_AS_file.write(f"{cluster_geneid}\t{cluster_chrom}:{cluster_start}-{cluster_stop}\t{cluster_AS}\n")
