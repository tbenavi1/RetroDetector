import gzip

with gzip.open(snakemake.input[0], "rt") as input_transcript_file, open(snakemake.output[0], "w") as output_coords_file, open(snakemake.output[1], "w") as output_junctions_file, open(snakemake.output[2], "w") as output_introns_file, open(snakemake.output[3], "w") as output_intronlengths_file, open(snakemake.output[4], "w") as output_longregions_file, open(snakemake.output[5], "w") as output_longregiontranscripts_file:
	for line in input_transcript_file:
		if line.startswith(">"):
			transcript = line.strip().split(" ")[0][1:]
			if transcript.startswith("FB"):
				ncbi_transcript = False
			else:
				assert transcript.startswith("lcl")
				ncbi_transcript = True
			qualifiers = line.strip().split(" ")
			#don't include any partial transcripts, any transcripts with exceptions, or any ncRNA
			#these may occur in the NCBI transcript files
			if ncbi_transcript and any(qualifier.startswith("[partial=") or qualifier.startswith("[exception=") or qualifier.startswith("[ncRNA_class=") for qualifier in qualifiers):
				continue
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
				continue
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
				#ref_intronlengths = reversed(ref_intronlengths)
				ref_intronlengths.reverse()
			for i, junction in enumerate(junctions):
				intron_length = ref_intronlengths[i]
				output_intronlengths_file.write(f"{GeneID}\t{transcript}\t{junction}\t{intron_length}\n")
			junctions = "\t".join(str(junction) for junction in junctions)
			ref_introns = "\t".join(f"{chrom}:{ref_intron}" for ref_intron in ref_introns)
			output_junctions_file.write(f"{GeneID}\t{transcript}\t{junctions}\n")
			output_introns_file.write(f"{GeneID}\t{transcript}\t{ref_introns}\n")
			for ref_junction in ref_junctions:
				#short_start = ref_junction - 74
				#short_stop = ref_junction + 75
				#output_shortregions_file.write(f"{chrom}:{short_start}-{short_stop}\n")
				#short_start = ref_junction - 129
				#short_stop = ref_junction + 20
				#output_shortregions_file.write(f"{chrom}:{short_start}-{short_stop}\n")
				#short_start = ref_junction - 19
				#short_stop = ref_junction + 130
				#output_shortregions_file.write(f"{chrom}:{short_start}-{short_stop}\n")
				#short_start = ref_junction - 139
				#short_stop = ref_junction + 10
				#output_shortregions_file.write(f"{chrom}:{short_start}-{short_stop}\n")
				#short_start = ref_junction - 9
				#short_stop = ref_junction + 140
				#output_shortregions_file.write(f"{chrom}:{short_start}-{short_stop}\n")
				long_start = max(ref_junction - 4999, 1)
				long_stop = ref_junction + 5000
				output_longregions_file.write(f"{chrom}:{long_start}-{long_stop}\n")
				output_longregiontranscripts_file.write(f"{transcript}\n")
				long_start = max(ref_junction - 9899, 1)
				long_stop = ref_junction + 100
				output_longregions_file.write(f"{chrom}:{long_start}-{long_stop}\n")
				output_longregiontranscripts_file.write(f"{transcript}\n")
				long_start = max(ref_junction - 99, 1)
				long_stop = ref_junction + 9900
				output_longregions_file.write(f"{chrom}:{long_start}-{long_stop}\n")
				output_longregiontranscripts_file.write(f"{transcript}\n")
