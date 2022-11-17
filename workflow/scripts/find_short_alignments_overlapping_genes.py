import os
import subprocess

bam_file = snakemake.input[1] 

with open(snakemake.input[0], "r") as input_coords_file, open(snakemake.output[0], "w") as output_overlapping_file:
	for line in input_coords_file:
		geneid, transcript, chrom, start, stop = line.strip().split()
		region = f"{chrom}:{start}-{stop}"
		subprocess.run(f"samtools view {bam_file} {region} > gene.sam", shell=True)
		overlapping_reads = set()
		with open("gene.sam", "r") as input_sam_file:
			for line in input_sam_file:
				readname = line.strip().split()[0]
				overlapping_reads.add(readname)
		overlapping_reads = "\t".join(overlapping_reads)
		output_overlapping_file.write(f"{geneid}\t{transcript}\t{overlapping_reads}\n")
os.remove("gene.sam")
