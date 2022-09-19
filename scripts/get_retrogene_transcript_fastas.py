import subprocess

ref = snakemake.wildcards.ref
sample = snakemake.wildcards.sample

transcriptome = f"Reference/{ref}.rna_from_genomic.fna.gz"

#output_name_prefix = f"Consensus/{ref}/{sample}/long/{ref}.{sample}."
#output_name_suffix = ".transcript.name"

output_transcript_prefix = f"Consensus/{ref}/{sample}/long/{ref}.{sample}."
output_transcript_suffix = ".transcript.fasta"

subprocess.run(f"mkdir -p Consensus/{ref}/{sample}/long", shell=True)

consensus_prefix = f"RESULTS/{ref}/{sample}/long/consensus/{ref}.{sample}.retrogene."
consensus_suffix = ".consensus.fasta"

output_needle_prefix = f"RESULTS/{ref}/{sample}/long/consensus/{ref}.{sample}.retrogene."
output_needle_suffix = ".needle"

output = snakemake.output[0]

num_retrogenes = 0
with open(snakemake.input[0], "r") as input_file:
	for line in input_file:
		num_retrogenes += 1
		geneid, transcript_location, region = line.strip().split()
		#output_name = output_name_prefix + geneid + output_name_suffix
		output_transcript_fasta = output_transcript_prefix + geneid + output_transcript_suffix
		consensus_fasta = consensus_prefix + geneid + consensus_suffix
		output_needle = output_needle_prefix + geneid + output_needle_suffix
		#with open(output_name, "w") as output_name_file:
		#	output_name_file.write(f"{transcript}\n")
		#subprocess.run(f"seqtk subseq {transcriptome} {output_name} > {output_transcript_fasta}", shell=True)
		subprocess.run(f"samtools faidx {transcriptome} '{transcript_location}' -o {output_transcript_fasta}", shell=True)
		subprocess.run(f"needle -asequence {output_transcript_fasta} -bsequence {consensus_fasta} -gapopen 10.0 -gapextend 0.5 -outfile {output_needle}", shell=True)

with open(output, "w") as output_file:
	output_file.write(f"All {num_retrogenes} retrogene(s) processed successfully.")
