import gzip
import subprocess

output_folder = snakemake.output[0][:-29]

with gzip.open(snakemake.input[0], "rt") as input_file:
	for line in input_file:
		if line.startswith(">"):
			transcript = line.strip().split()[0][1:]
			transcript_file = f"{output_folder}/{transcript}.fa"
			mmi_file = f"{output_folder}/{transcript}-hifi.mmi"
			subprocess.run(f"minimap2 -x map-hifi -d '{mmi_file}' '{transcript_file}'", shell=True)

with open(snakemake.output[0], "w") as output_file:
	output_file.write(f"All transcripts indexed.")
