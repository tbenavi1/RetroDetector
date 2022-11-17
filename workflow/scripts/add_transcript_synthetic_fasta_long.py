with open(snakemake.input[0], "r") as input_transcripts_file, open(snakemake.input[1], "r") as input_fasta_file, open(snakemake.output[0], "w") as output_fasta_file:
	for line in input_fasta_file:
		if line.startswith(">"):
			transcript = input_transcripts_file.readline().strip()
			region = line.strip()[1:]
			new_line = f">{transcript}:{region}\n"
			output_fasta_file.write(new_line)
		else:
			output_fasta_file.write(line)
