import subprocess

ref = snakemake.wildcards.ref
sample = snakemake.wildcards.sample
junction_overhang = snakemake.wildcards.junction_overhang
insertions_threshold = snakemake.wildcards.insertions_threshold
junction_total_read_support_threshold = snakemake.wildcards.junction_total_read_support_threshold

input_diff = snakemake.input[0]
input_bam = snakemake.input[1]
output = snakemake.output[0]
output_fasta_prefix = f"results/retrogenes/{ref}/{sample}/long/consensus/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogene."
output_fasta_suffix = f".consensus.fasta"

num_retrogenes = 0
with open(input_diff, "r") as input_file:
	for line in input_file:
		num_retrogenes += 1
		geneid, transcript, region, direction, readnames = line.strip().split()
		output_fasta = output_fasta_prefix + geneid + "." + region.split(":")[0] + "." + region.split(":")[1] + output_fasta_suffix
		subprocess.run(f"samtools consensus -r {region} -o '{output_fasta}' {input_bam}", shell=True)

with open(output, "w") as output_file:
	output_file.write(f"All {num_retrogenes} retrogene(s) processed successfully.")
