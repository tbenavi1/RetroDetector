import gzip

with gzip.open(snakemake.input[0], "rt") as input_transcript_file, open(snakemake.output[0], "w") as output_multiexon_file, open(snakemake.output[1], "w") as output_coords_file, open(snakemake.output[2], "w") as output_junctions_file, open(snakemake.output[3], "w") as output_introns_file, open(snakemake.output[4], "w") as output_intronlengths_file:
	for line in input_transcript_file:
		if line.startswith(">"):
			transcript = line.strip().split(" ")[0][1:]
			if transcript.startswith("FB"):
				ncbi_transcript = False
			else:
				assert transcript.startswith("lcl"), "We currently only support NCBI or FlyBase transcriptome files. Please submit a GitHub issue."
				ncbi_transcript = True
			qualifiers = line.strip().split(" ")
			#don't include any partial transcripts, any transcripts with exceptions, or any ncRNA
			#these may occur in the NCBI transcript files
			if ncbi_transcript and any(qualifier.startswith("[partial=") or qualifier.startswith("[exception=") or qualifier.startswith("[ncRNA_class=") for qualifier in qualifiers):
				write_line = False
				continue
			GeneID = ""
			locus_tag = ""
			if not ncbi_transcript:
				loc = line.strip().split(";")[1]
				chrom = loc.split(":")[0].split("=")[1]
				ranges = loc.split(":")[1]
				for qualifier in qualifiers:
					if qualifier.startswith("dbxref="):
						items = qualifier[7:].split(",")
						for item in items:
							if item.startswith("FlyBase_Annotation_IDs:GE"):
								GeneID = item.split(":")[1]
			else:
				chrom = transcript.split("|")[1].split("_")[0] + "_" + transcript.split("|")[1].split("_")[1]
				for qualifier in qualifiers:
					if qualifier.startswith("[location="):
						ranges = qualifier.split("=")[1][:-1]
					if qualifier.startswith("[db_xref="):
						if "," in qualifier:
							GeneID = qualifier.split(",")[1].split(":")[1][:-1]
						else:
							GeneID = qualifier.split(":")[1][:-1]
					#if qualifier.startswith("[locus_tag="):
					#	assert qualifier.startswith("[locus_tag=Dmel_"), "We currently only process Drosophila melanogaster locus tags. Please submit a GitHub issue."
					if qualifier.startswith("[locus_tag=Dmel_"):
						locus_tag = qualifier.split("_")[2][:-1]
			if locus_tag:
				GeneID = locus_tag
			#There is a glitch on the Drosophila melanogaster transcriptome file where gene 49228 doesn't have a locus tag
			if GeneID == "49228":
				GeneID = "CG32491"
			#if there are at least two exons for this transcript
			if "join" in ranges:
				#if this transcript occurs on the forward strand
				if "complement" not in ranges:
					positions = ranges[5:-1]
				#else this transcript occurs on the reverse strand
				else:
					positions = ranges[16:-2]
			#else there is only one exon for this transcript and we do not include
			else:
				write_line = False
				continue
			write_line = True
			output_multiexon_file.write(line)
			exons = positions.split(",")
			num_exons = len(exons)
			transcript_start = exons[0].split("..")[0]
			transcript_stop = exons[-1].split("..")[-1]
			output_coords_file.write(f"{GeneID}\t{transcript}\t{chrom}\t{transcript_start}\t{transcript_stop}\n")
			startstops = []
			skip_line = False
			for exon in exons:
				#some of the lines seem malformatted where the exon doesn't have a start and stop position, we skip these
				try:
					start, stop = exon.split("..")
					start, stop = int(start), int(stop)
					startstops.append((start, stop))
				except:
					skip_line = True
			if skip_line:
				continue
			junctions = []
			ref_junctions = []
			ref_introns = []
			ref_intronlengths = []
			current_length = 0
			for i, (start, stop) in enumerate(startstops):
				if i > 0:
					ref_junctions.append(start - 1)
					previous_stop = startstops[i-1][1]
					ref_introns.append(f"{previous_stop+1}-{start-1}")
					intron_length = start - previous_stop - 1
					ref_intronlengths.append(intron_length)
				if i < num_exons - 1:
					ref_junctions.append(stop)
				dist = stop - start + 1
				current_length += dist
				#we append the 1-indexed position of the last base pair before the junction
				junctions.append(current_length)
			#we remove the last item, since this will be the 1-indexed position of the
			#last base in the transcript but this is not a junction between two exons
			junctions = junctions[:-1]
			if ranges.startswith("complement"):
				junctions = [current_length - junction for junction in reversed(junctions)]
				ref_introns = reversed(ref_introns)
				ref_intronlengths.reverse()
			for i, junction in enumerate(junctions):
				intron_length = ref_intronlengths[i]
				output_intronlengths_file.write(f"{GeneID}\t{transcript}\t{junction}\t{intron_length}\n")
			junctions = "\t".join(str(junction) for junction in junctions)
			ref_introns = "\t".join(f"{chrom}:{ref_intron}" for ref_intron in ref_introns)
			output_junctions_file.write(f"{GeneID}\t{transcript}\t{junctions}\n")
			output_introns_file.write(f"{GeneID}\t{transcript}\t{ref_introns}\n")
		else:
			if write_line:
				output_multiexon_file.write(line)
