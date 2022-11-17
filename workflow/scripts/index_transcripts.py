import gzip
import subprocess

output_folder = snakemake.output[0][:-24]

with gzip.open(snakemake.input[0], "rt") as input_file:
	for line in input_file:
		if line.startswith(">"):
			transcript = line.strip().split()[0][1:]
			transcript_file = f"{output_folder}/{transcript}.fa"
			subprocess.run(f"bwa index '{transcript_file}'", shell=True)

with open(snakemake.output[0], "w") as output_file:
	output_file.write(f"All transcripts indexed.")
