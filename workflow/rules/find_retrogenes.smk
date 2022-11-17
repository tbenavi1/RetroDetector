#Find transcriptome alignments which span exon-exon junctions but don't have intronic sequence

rule find_retrogene_alignments_long:
  input:
    junctions="results/Transcriptome/{ref}/{ref}.transcriptome.junctions.tsv",
    bam="results/BAMS/{ref}/{sample}/transcriptome/long/{ref}.{sample}.transcriptome.long.sorted.bam"
  output:
    spanningalignments="results/Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spanningalignments.sam",
    spanningalignmentstxt="results/Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spanningalignments.txt",
    spannedjunctions="results/Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.tsv",
    flankregions="results/Flanks/{ref}/{sample}/long/{ref}.{sample}.longreads.flankregions.tsv"
  params:
    junction_overhang=config["junction_overhang"],
    insertions_threshold=config["insertions_threshold"]
  shell:
    "samtools view -h {input.bam} | python workflow/scripts/find_retrogene_alignments_long.py {params.junction_overhang} {params.insertions_threshold} {input.junctions} {output.spanningalignments} {output.spanningalignmentstxt} {output.spannedjunctions} {output.flankregions}"

#Only include transcriptome alignments that support a spanned exon-exon junction with sufficient read coverage

rule subset_transcriptome_bam:
  input:
    "results/Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.tsv",
    "results/Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spanningalignments.sam",
    "results/Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spanningalignments.txt"
  output:
    "results/Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.readnames.txt",
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.sam"
  script:
    "../scripts/subset_transcriptome_bam.py"

#Only include genome alignments that have read names corresponding to the transcriptome alignments above

rule subset_genome_bam:
  input:
    reads="results/Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.readnames.txt",
    bam="results/BAMS/{ref}/{sample}/genome/long/{ref}.{sample}.genome.long.sorted.bam"
  output:
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.sam"
  shell:
    "samtools view -h -N {input.reads} {input.bam} > {output}"

#Add parent GeneIDs to the alignments so that we can sort and process everything one parent gene at a time

rule add_geneid_to_subset_sams:
  input:
    "results/Transcriptome/{ref}/{ref}.transcriptome.coords.tsv",
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.sam",
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.sam"
  output:
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.geneid.sam",
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sam"
  script:
    "../scripts/add_geneid_to_subset_sam.py"

#sort by parent GeneID

rule sort_genome_geneid_sam:
  input:
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sam"
  output:
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sorted.sam"
  shell:
    "sort -k1,1 {input} > {output}"

rule sort_transcriptome_geneid_sam:
  input:
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.geneid.sam"
  output:
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.geneid.sorted.sam"
  shell:
    "sort -k1,1 {input} > {output}"
