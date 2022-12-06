import subprocess

ref = snakemake.wildcards.ref
sample = snakemake.wildcards.sample
junction_overhang = snakemake.wildcards.junction_overhang
insertions_threshold = snakemake.wildcards.insertions_threshold
junction_total_read_support_threshold = snakemake.wildcards.junction_total_read_support_threshold

transcriptome = f"results/Transcriptome/{ref}/{ref}.transcriptome.fna.gz"

output_transcript_prefix = f"results/Consensus/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}."
output_transcript_suffix = ".transcript.fasta"

subprocess.run(f"mkdir -p results/Consensus/{ref}/{sample}/long", shell=True)

consensus_prefix = f"results/retrogenes/{ref}/{sample}/long/consensus/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogene."
consensus_suffix = ".consensus.fasta"

output_needle_prefix = f"results/retrogenes/{ref}/{sample}/long/consensus/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogene."
output_needle_suffix = ".needle"

output = snakemake.output[0]

num_retrogenes = 0
with open(snakemake.input[0], "r") as input_file:
	for line in input_file:
		num_retrogenes += 1
		geneid, transcript_location, region, direction, readnames = line.strip().split()
		region = region.split(":")[0] + "." + region.split(":")[1]
		if direction == "forward":
			option = ""
		else:
			assert direction == "reverse", line
			option = " --reverse-complement"
		output_transcript_fasta = output_transcript_prefix + geneid + output_transcript_suffix
		consensus_fasta = consensus_prefix + geneid + "." + region + consensus_suffix
		output_needle = output_needle_prefix + geneid + "." + region + output_needle_suffix
		subprocess.run(f"samtools faidx {transcriptome} '{transcript_location}' -o '{output_transcript_fasta}'{option}", shell=True)
		subprocess.run(f"needle -asequence '{output_transcript_fasta}' -bsequence '{consensus_fasta}' -gapopen 10.0 -gapextend 0.5 -outfile '{output_needle}'", shell=True)

with open(output, "w") as output_file:
	output_file.write(f"All {num_retrogenes} retrogene(s) processed successfully.")
