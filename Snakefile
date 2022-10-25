configfile: "config.yaml"

def get_all_short_fastq(wildcards):
  files = []
  for sample in config["fastqs"]:
    if "short_paired_end" in config["fastqs"][sample]:
      i_s = config["fastqs"][sample]["short_paired_end"]
      for i in i_s:
        file = f"FASTQS/{sample}/short_paired_end/{sample}.short.{i}.R1.fastq.gz"
        files.append(file)
        file = f"FASTQS/{sample}/short_paired_end/{sample}.short.{i}.R2.fastq.gz"
        files.append(file)
  return files

def get_all_long_fastq(wildcards):
  files = []
  for sample in config["fastqs"]:
    if "pacbio_hifi" in config["fastqs"][sample]:
      i_s = config["fastqs"][sample]["pacbio_hifi"]
      for i in i_s:
        file = f"FASTQS/{sample}/pacbio_hifi/{sample}.pacbio_hifi.{i}.fastq.gz"
        files.append(file)
    if "pacbio_clr" in config["fastqs"][sample]:
      i_s = config["fastqs"][sample]["pacbio_clr"]
      for i in i_s:
        file = f"FASTQS/{sample}/pacbio_clr/{sample}.pacbio_clr.{i}.fastq.gz"
        files.append(file)
    if "ont" in config["fastqs"][sample]:
      i_s = config["fastqs"][sammple]["ont"]
      for i in i_s:
        file = f"FASTQS/{sample}/ont/{sample}.ont.{i}.fastq.gz"
        files.append(file)
  return files

short_samples = []
for sample in config["fastqs"]:
  if "short_paired_end" in config["fastqs"][sample]:
    short_samples.append(sample)

long_samples = []
for sample in config["fastqs"]:
  if "pacbio_hifi" in config["fastqs"][sample] or "pacbio_clr" in config["fastqs"][sample] or "ont" in config["fastqs"][sample]:
    long_samples.append(sample)

rule all:
  input:
    #expand("BAMS/{ref}/{sample}/transcriptome/short_paired_end/{ref}.{sample}.transcriptome.short.sorted.bam", ref=config["ref"], sample=short_samples),
    #expand("BAMS/{ref}/{sample}/transcriptome/long/{ref}.{sample}.transcriptome.long.sorted.bam", ref=config["ref"], sample=long_samples),
    #expand("BAMS/{ref}/{sample}/genome/short_paired_end/{ref}.{sample}.genome.short.sorted.bam", ref=config["ref"], sample=short_samples),
    #expand("BAMS/{ref}/{sample}/genome/long/{ref}.{sample}.genome.long.sorted.bam", ref=config["ref"], sample=long_samples)
    #expand("Spanned/{ref}/{sample}/short_paired_end/{ref}.{sample}.shortreads.overlapping.genes.tsv", ref=config["ref"], sample=["1300_13"])
    #expand("Spanned/{ref}/{sample}/short_paired_end/{ref}.{sample}.shortreads.spannedjunctions.tsv", ref=config["ref"], sample=["1300_13"])
    #expand("Summary/{ref}/{sample}/short_paired_end/{ref}.{sample}.{strongthreshold}.shortreads.numspannedjunctions.tsv", ref=config["ref"], sample=["1300_13"], strongthreshold=config["short_read_strong_spanning_alignment_expected_genomic_insertion_size_threshold"])
    #expand("Summary/{ref}/{sample}/{ref}.{sample}.{strongthreshold}.results.txt", ref=config["ref"], sample=["1300_13"], strongthreshold=config["short_read_strong_spanning_alignment_expected_genomic_insertion_size_threshold"])
    #expand("Junctions/{ref}.synthetic.longreads.fasta", ref=config["ref"])
    #expand("Junctions/{ref}.junctions.unambiguous.long.tsv", ref=config["ref"])
    #expand("Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.tsv", ref=config["ref"], sample=["1300_13"])
    expand("AS/{ref}/{sample}/{ref}.{sample}.genome.AS", ref=config["ref"], sample=["1300_13"])
    #expand("AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.best.AS.diff", ref=config["ref"], sample=["1300_13"])
    #expand("RESULTS/{ref}/{sample}/long/consensus/{ref}.{sample}.transcript.fastas.txt", ref=config["ref"], sample=["1300_13"])
    #expand("Analysis/{ref}/{sample}/{ref}.{sample}.introns.missing_percent", ref=config["ref"], sample=["1300_13"])
    #expand("Analysis/{ref}/{sample}/{ref}.{sample}.intronsummary", ref=config["ref"], sample=["1300_13"])
    #expand("SFS/{ref}/summary/long/{ref}.longreads.freqtable.multiplyspanned.tsv", ref=config["ref"])
    #expand("SFS/{ref}/summary/short_paired_end/{ref}.shortreads.freqtable.multiplyspanned.tsv", ref=config["ref"])

#bgzip reference genome

def get_genomic_fna(wildcards):
  return config["ref"][wildcards.ref]["genomic_fna"]

rule unzip_genomic_fna:
  input:
    get_genomic_fna
  output:
    "Reference/{ref}.genomic.fna"
  shell:
    "gunzip -c {input} > {output}"

rule bgzip_genomic_fna:
  input:
    "Reference/{ref}.genomic.fna"
  output:
    "Reference/{ref}.genomic.fna.gz"
  shell:
    "bgzip {input}"

#bgzip reference gtf

def get_genomic_gtf(wildcards):
  return config["ref"][wildcards.ref]["genomic_gtf"]

rule unzip_genomic_gtf:
  input:
    get_genomic_gtf
  output:
    "Reference/{ref}.genomic.gtf"
  shell:
    "gunzip -c {input} > {output}"

rule bgzip_genomic_gtf:
  input:
    "Reference/{ref}.genomic.gtf"
  output:
    "Reference/{ref}.genomic.gtf.gz"
  shell:
    "bgzip {input}"

#bgzip reference transcriptome

def get_rna_from_genomic(wildcards):
  return config["ref"][wildcards.ref]["rna_from_genomic"]

rule unzip_rna_from_genomic:
  input:
    get_rna_from_genomic
  output:
    "Reference/{ref}.rna_from_genomic.fna"
  shell:
    "gunzip -c {input} > {output}"

rule bgzip_rna_from_genomic:
  input:
    "Reference/{ref}.rna_from_genomic.fna"
  output:
    "Reference/{ref}.rna_from_genomic.fna.gz"
  shell:
    "bgzip {input}"

#index reference genome

rule index_genome:
  input:
    "Reference/{ref}.genomic.fna.gz"
  output:
    "Reference/{ref}.genomic.fna.gz.fai",
    "Reference/{ref}.genomic.fna.gz.gzi"
  shell:
    "samtools faidx {input}"

#link FASTQ files

def get_short_fastq(wildcards):
  return config["fastqs"][wildcards.sample]["short_paired_end"][int(wildcards.i)][int(wildcards.j)]

rule link_short_fastq:
  input:
    get_short_fastq
  output:
    "FASTQS/{sample}/short_paired_end/{sample}.short.{i}.R{j}.fastq.gz"
  shell:
    "ln -sr {input} {output}"

def get_long_fastq(wildcards):
	return config["fastqs"][wildcards.sample][wildcards.tech][int(wildcards.i)]

rule link_long_fastq:
	input:
		get_long_fastq
	output:
		"FASTQS/{sample}/{tech}/{sample}.{tech}.{i}.fastq.gz"
	shell:
		"ln -sr {input} {output}"

#Index transcriptome according to short-read/long-read sequencing technology

rule bwa_index_transcriptome:
  input:
    "Reference/{ref}.rna_from_genomic.fna.gz"
  output:
    "Reference/{ref}.rna_from_genomic.fna.gz.amb",
    "Reference/{ref}.rna_from_genomic.fna.gz.ann",
    "Reference/{ref}.rna_from_genomic.fna.gz.bwt",
    "Reference/{ref}.rna_from_genomic.fna.gz.pac",
    "Reference/{ref}.rna_from_genomic.fna.gz.sa"
  shell:
    "bwa index {input}"

rule minimap2_hifi_index_transcriptome:
  input:
    "Reference/{ref}.rna_from_genomic.fna.gz"
  output:
    "Reference/{ref}.rna_from_genomic.pacbio_hifi.mmi"
  threads: 3
  shell:
    "minimap2 -x map-hifi -d {output} {input}"

rule minimap2_pb_index_transcriptome:
  input:
    "Reference/{ref}.rna_from_genomic.fna.gz"
  output:
    "Reference/{ref}.rna_from_genomic.pacbio_clr.mmi"
  threads: 3
  shell:
    "minimap2 -x map-pb -d {output} {input}"

rule minimap2_ont_index_transcriptome:
  input:
    "Reference/{ref}.rna_from_genomic.fna.gz"
  output:
    "Reference/{ref}.rna_from_genomic.ont.mmi"
  threads: 3
  shell:
    "minimap2 -x map-ont -d {output} {intput}"

#Index genome according to short-read/long-read sequencing technology

rule bwa_index_genome:
  input:
    "Reference/{ref}.genomic.fna.gz"
  output:
    "Reference/{ref}.genomic.fna.gz.amb",
    "Reference/{ref}.genomic.fna.gz.ann",
    "Reference/{ref}.genomic.fna.gz.bwt",
    "Reference/{ref}.genomic.fna.gz.pac",
    "Reference/{ref}.genomic.fna.gz.sa"
  shell:
    "bwa index {input}"

rule minimap2_hifi_index_genome:
  input:
    "Reference/{ref}.genomic.fna.gz"
  output:
    "Reference/{ref}.genomic.pacbio_hifi.mmi"
  threads: 3
  shell:
    "minimap2 -x map-hifi -d {output} {input}"

rule minimap2_pb_index_genome:
  input:
    "Reference/{ref}.genomic.fna.gz"
  output:
    "Reference/{ref}.genomic.pacbio_clr.mmi"
  threads: 3
  shell:
    "minimap2 -x map-pb -d {output} {input}"

rule minimap2_ont_index_genome:
  input:
    "Reference/{ref}.genomic.fna.gz"
  output:
    "Reference/{ref}.genomic.ont.mmi"
  threads: 3
  shell:
    "minimap2 -x map-ont -d {output} {intput}"

#Analyze transcriptome to get coords, junctions, introns, intron lengths, and synthetic long-read regions

def get_scripts(wildcards):
  return config["scripts_directory"]

rule get_transcriptome_junctions:
  input:
    "Reference/{ref}.rna_from_genomic.fna.gz"
  output:
    "Transcriptome/{ref}.transcriptome.coords.tsv",
    "Junctions/{ref}.junctions.tsv",
    "Junctions/{ref}.introns.tsv",
    "Junctions/{ref}.intron_lengths.tsv",
    "Junctions/{ref}.synthetic.longregions.txt",
    "Junctions/{ref}.synthetic.longregiontranscripts.txt"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/get_transcriptome_junctions.py"

#Make synethetic long reads

rule get_synthetic_fasta_long:
  input:
    "Reference/{ref}.genomic.fna.gz.fai",
    "Reference/{ref}.genomic.fna.gz.gzi",
    ref="Reference/{ref}.genomic.fna.gz",
    regions="Junctions/{ref}.synthetic.longregions.txt"
  output:
    "Junctions/{ref}.synthetic.temp.longreads.fasta"
  shell:
    "samtools faidx {input.ref} -r {input.regions} > {output}"

rule add_transcript_synthetic_fasta_long:
  input:
    "Junctions/{ref}.synthetic.longregiontranscripts.txt",
    "Junctions/{ref}.synthetic.temp.longreads.fasta"
  output:
    "Junctions/{ref}.synthetic.longreads.fasta"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/add_transcript_synthetic_fasta_long.py"

#Map synthetic long reads to transcriptome

rule minimap2_synthetic_long_to_transcriptome:
	input:
		mmi="Reference/{ref}.rna_from_genomic.pacbio_hifi.mmi",
		reads="Junctions/{ref}.synthetic.longreads.fasta"
	output:
		temp("BAMS/{ref}/synthetic/transcriptome/pacbio_hifi/{ref}.synthetic.longreads.transcriptome.sam")
	threads: 15
	shell:
		"minimap2 -t {threads} -a {input.mmi} {input.reads} > {output}"

#Sort transcriptome bams for synthetic long reads

rule sam_to_bam_transcriptome_synthetic_long:
  input:
    "BAMS/{ref}/synthetic/transcriptome/pacbio_hifi/{ref}.synthetic.longreads.transcriptome.sam"
  output:
    temp("BAMS/{ref}/synthetic/transcriptome/pacbio_hifi/{ref}.synthetic.longreads.transcriptome.bam")
  threads: 16
  shell:
    "samtools view -F 4 -@ {threads} -b -o {output} {input}"

rule sort_bam_transcriptome_synthetic_long:
  input:
    "BAMS/{ref}/synthetic/transcriptome/pacbio_hifi/{ref}.synthetic.longreads.transcriptome.bam"
  output:
    "BAMS/{ref}/synthetic/transcriptome/pacbio_hifi/{ref}.synthetic.longreads.transcriptome.sorted.bam"
  threads: 16
  shell:
    "samtools sort -@ {threads} -o {output} {input}"

#Find alignments of synthetic long reads that span exon/exon junctions without corresponding intronic sequence

rule find_synthetic_long_alignments:  
  input:
    junctions="Junctions/{ref}.junctions.tsv",
    bam="BAMS/{ref}/synthetic/transcriptome/pacbio_hifi/{ref}.synthetic.longreads.transcriptome.sorted.bam"
  output:
    spanningalignments="Ambiguous/{ref}.ambiguous.long.spanningalignments.sam",
    spanningalignmentstxt="Ambiguous/{ref}.ambiguous.long.spanningalignments.txt",
    spannedjunctions="Ambiguous/{ref}.ambiguous.long.spannedjunctions.tsv",
    flankregions="Flanks/{ref}/synthetic/pacbio_hifi/{ref}.synthetic.flankregions.tsv"
  params:
    scripts=get_scripts,
    junction_overhang=config["junction_overhang"],
    insertions_threshold=config["insertions_threshold"]
  shell:
    "samtools view -h {input.bam} | python {params.scripts}/find_synthetic_longreads_intronless_junction_spanning_alignments.py {params.junction_overhang} {params.insertions_threshold} {input.junctions} {output.spanningalignments} {output.spanningalignmentstxt} {output.spannedjunctions} {output.flankregions}"

#Remove ambiguous junctions (based on alignments of synthetic long reads)

rule remove_ambiguous_junctions_long:
  input:
    "Ambiguous/{ref}.ambiguous.long.spannedjunctions.tsv",
    "Junctions/{ref}.junctions.tsv"
  output:
    "Junctions/{ref}.junctions.unambiguous.long.tsv"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/remove_ambiguous_junctions.py"

#Map real short and long reads to transcriptome

rule bwa_mem_transcriptome_short:
  input:
    "Reference/{ref}.rna_from_genomic.fna.gz.amb",
    "Reference/{ref}.rna_from_genomic.fna.gz.ann",
    "Reference/{ref}.rna_from_genomic.fna.gz.bwt",
    "Reference/{ref}.rna_from_genomic.fna.gz.pac",
    "Reference/{ref}.rna_from_genomic.fna.gz.sa",
    ref="Reference/{ref}.rna_from_genomic.fna.gz",
    read1="FASTQS/{sample}/short_paired_end/{sample}.short.{i}.R1.fastq.gz",
    read2="FASTQS/{sample}/short_paired_end/{sample}.short.{i}.R2.fastq.gz"
  output:
    temp("BAMS/{ref}/{sample}/transcriptome/short_paired_end/{ref}.{sample}.transcriptome.illumina.{i}.sam")
  threads: 16
  shell:
    "bwa mem -t {threads} {input.ref} {input.read1} {input.read2} > {output}"

rule minimap_transcriptome_long:
	input:
		mmi="Reference/{ref}.rna_from_genomic.{tech}.mmi",
		reads="FASTQS/{sample}/{tech}/{sample}.{tech}.{i}.fastq.gz"
	output:
		temp("BAMS/{ref}/{sample}/transcriptome/{tech}/{ref}.{sample}.transcriptome.{tech}.{i}.sam")
	threads: 15
	shell:
		"minimap2 -t {threads} -a {input.mmi} {input.reads} > {output}"

#Sort and merge transcriptome bams for real short and long reads

rule sam_to_bam_transcriptome_short:
  input:
    "BAMS/{ref}/{sample}/transcriptome/short_paired_end/{ref}.{sample}.transcriptome.illumina.{i}.sam"
  output:
    temp("BAMS/{ref}/{sample}/transcriptome/short_paired_end/{ref}.{sample}.transcriptome.illumina.{i}.bam")
  threads: 16
  shell:
    "samtools view -F 4 -@ {threads} -b -o {output} {input}"

rule sort_bam_transcriptome_short:
  input:
    "BAMS/{ref}/{sample}/transcriptome/short_paired_end/{ref}.{sample}.transcriptome.illumina.{i}.bam"
  output:
    temp("BAMS/{ref}/{sample}/transcriptome/short_paired_end/{ref}.{sample}.transcriptome.illumina.{i}.sorted.bam")
  threads: 16
  shell:
    "samtools sort -@ {threads} -o {output} {input}"

def get_sorted_bam_transcriptome_short(wildcards):
  files = []
  if "short_paired_end" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["short_paired_end"]
    for i in i_s:
      file = f"BAMS/{wildcards.ref}/{wildcards.sample}/transcriptome/short_paired_end/{wildcards.ref}.{wildcards.sample}.transcriptome.illumina.{i}.sorted.bam"
      files.append(file)
  return files

rule merge_bam_transcriptome_short:
  input:
    get_sorted_bam_transcriptome_short
  output:
    "BAMS/{ref}/{sample}/transcriptome/short_paired_end/{ref}.{sample}.transcriptome.short.sorted.bam"
  threads: 16
  shell:
    "samtools merge -@ {threads} -o {output} {input}"

rule sam_to_bam_transcriptome_long:
	input:
		"BAMS/{ref}/{sample}/transcriptome/{tech}/{ref}.{sample}.transcriptome.{tech}.{i}.sam"
	output:
		temp("BAMS/{ref}/{sample}/transcriptome/{tech}/{ref}.{sample}.transcriptome.{tech}.{i}.bam")
	threads: 16
	shell:
		"samtools view -F 4 -@ {threads} -b -o {output} {input}"

rule sort_bam_transcriptome_long:
	input:
		"BAMS/{ref}/{sample}/transcriptome/{tech}/{ref}.{sample}.transcriptome.{tech}.{i}.bam"
	output:
		temp("BAMS/{ref}/{sample}/transcriptome/{tech}/{ref}.{sample}.transcriptome.{tech}.{i}.sorted.bam")
	threads: 16
	shell:
		"samtools sort -@ {threads} -o {output} {input}"

def get_sorted_bam_transcriptome_long(wildcards):
  files = []
  if "pacbio_hifi" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["pacbio_hifi"]
    for i in i_s:
      file = f"BAMS/{wildcards.ref}/{wildcards.sample}/transcriptome/pacbio_hifi/{wildcards.ref}.{wildcards.sample}.transcriptome.pacbio_hifi.{i}.sorted.bam"
      files.append(file)
  if "pacbio_clr" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["pacbio_clr"]
    for i in i_s:
      file = f"BAMS/{wildcards.ref}/{wildcards.sample}/transcriptome/pacbio_clr/{wildcards.ref}.{wildcards.sample}.transcriptome.pacbio_clr.{i}.sorted.bam"
      files.append(file)
  if "ont" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["ont"]
    for i in i_s:
      file = f"BAMS/{wildcards.ref}/{wildcards.sample}/transcriptome/ont/{wildcards.ref}.{wildcards.sample}.transcriptome.ont.{i}.sorted.bam"
      files.append(file)
  return files

rule merge_bam_transcriptome_long:
	input:
		get_sorted_bam_transcriptome_long
	output:
		"BAMS/{ref}/{sample}/transcriptome/long/{ref}.{sample}.transcriptome.long.sorted.bam"
	threads: 16
	shell:
		"samtools merge -@ {threads} -o {output} {input}"

#Map real short and long reads to genome

rule bwa_mem_genome_short:
  input:
    "Reference/{ref}.genomic.fna.gz.amb",
    "Reference/{ref}.genomic.fna.gz.ann",
    "Reference/{ref}.genomic.fna.gz.bwt",
    "Reference/{ref}.genomic.fna.gz.pac",
    "Reference/{ref}.genomic.fna.gz.sa",
    ref="Reference/{ref}.genomic.fna.gz",
    read1="FASTQS/{sample}/short_paired_end/{sample}.short.{i}.R1.fastq.gz",
    read2="FASTQS/{sample}/short_paired_end/{sample}.short.{i}.R2.fastq.gz"
  output:
    temp("BAMS/{ref}/{sample}/genome/short_paired_end/{ref}.{sample}.genome.illumina.{i}.sam")
  threads: 16
  shell:
    "bwa mem -t {threads} {input.ref} {input.read1} {input.read2} > {output}"

rule minimap_genome_long:
	input:
		mmi="Reference/{ref}.genomic.{tech}.mmi",
		reads="FASTQS/{sample}/{tech}/{sample}.{tech}.{i}.fastq.gz"
	output:
		temp("BAMS/{ref}/{sample}/genome/{tech}/{ref}.{sample}.genome.{tech}.{i}.sam")
	threads: 15
	shell:
		"minimap2 -t {threads} -a {input.mmi} {input.reads} > {output}"

#Sort and merge genome bams for real short and long reads

rule sam_to_bam_genome_short:
  input:
    "BAMS/{ref}/{sample}/genome/short_paired_end/{ref}.{sample}.genome.illumina.{i}.sam"
  output:
    temp("BAMS/{ref}/{sample}/genome/short_paired_end/{ref}.{sample}.genome.illumina.{i}.bam")
  threads: 16
  shell:
    "samtools view -F 4 -@ {threads} -b -o {output} {input}"

rule sort_bam_genome_short:
  input:
    "BAMS/{ref}/{sample}/genome/short_paired_end/{ref}.{sample}.genome.illumina.{i}.bam"
  output:
    temp("BAMS/{ref}/{sample}/genome/short_paired_end/{ref}.{sample}.genome.illumina.{i}.sorted.bam")
  threads: 16
  shell:
    "samtools sort -@ {threads} -o {output} {input}"

def get_sorted_bam_genome_short(wildcards):
  files = []
  i_s = config["fastqs"][wildcards.sample]["short_paired_end"]
  for i in i_s:
    file = f"BAMS/{wildcards.ref}/{wildcards.sample}/genome/short_paired_end/{wildcards.ref}.{wildcards.sample}.genome.illumina.{i}.sorted.bam"
    files.append(file)
  return files

rule merge_bam_genome_short:
  input:
    get_sorted_bam_genome_short
  output:
    "BAMS/{ref}/{sample}/genome/short_paired_end/{ref}.{sample}.genome.short.sorted.bam"
  threads: 16
  shell:
    "samtools merge -@ {threads} -o {output} {input}"

rule index_bam_genome_short:
  input:
    "BAMS/{ref}/{sample}/genome/short_paired_end/{ref}.{sample}.genome.short.sorted.bam"
  output:
    "BAMS/{ref}/{sample}/genome/short_paired_end/{ref}.{sample}.genome.short.sorted.bam.bai"
  threads: 16
  shell:
    "samtools index -@ {threads} {input}"

rule sam_to_bam_genome_long:
	input:
		"BAMS/{ref}/{sample}/genome/{tech}/{ref}.{sample}.genome.{tech}.{i}.sam"
	output:
		temp("BAMS/{ref}/{sample}/genome/{tech}/{ref}.{sample}.genome.{tech}.{i}.bam")
	threads: 16
	shell:
		"samtools view -F 4 -@ {threads} -b -o {output} {input}"

rule sort_bam_genome_long:
	input:
		"BAMS/{ref}/{sample}/genome/{tech}/{ref}.{sample}.genome.{tech}.{i}.bam"
	output:
		temp("BAMS/{ref}/{sample}/genome/{tech}/{ref}.{sample}.genome.{tech}.{i}.sorted.bam")
	threads: 16
	shell:
		"samtools sort -@ {threads} -o {output} {input}"

def get_sorted_bam_genome_long(wildcards):
  files = []
  if "pacbio_hifi" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["pacbio_hifi"]
    for i in i_s:
      file = f"BAMS/{wildcards.ref}/{wildcards.sample}/genome/pacbio_hifi/{wildcards.ref}.{wildcards.sample}.genome.pacbio_hifi.{i}.sorted.bam"
      files.append(file)
  if "pacbio_clr" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["pacbio_clr"]
    for i in i_s:
      file = f"BAMS/{wildcards.ref}/{wildcards.sample}/genome/pacbio_clr/{wildcards.ref}.{wildcards.sample}.genome.pacbio_clr.{i}.sorted.bam"
      files.append(file)
  if "ont" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["ont"]
    for i in i_s:
      file = f"BAMS/{wildcards.ref}/{wildcards.sample}/genome/ont/{wildcards.ref}.{wildcards.sample}.genome.ont.{i}.sorted.bam"
      files.append(file)
  return files

rule merge_bam_genome_long:
	input:
		get_sorted_bam_genome_long
	output:
		"BAMS/{ref}/{sample}/genome/long/{ref}.{sample}.genome.long.sorted.bam"
	threads: 16
	shell:
		"samtools merge -@ {threads} -o {output} {input}"

#Find alignments of real short reads that overlap genes

rule find_short_alignments_overlapping_genes:
  input:
    "Transcriptome/{ref}.transcriptome.coords.tsv",
    "BAMS/{ref}/{sample}/genome/short_paired_end/{ref}.{sample}.genome.short.sorted.bam",
    "BAMS/{ref}/{sample}/genome/short_paired_end/{ref}.{sample}.genome.short.sorted.bam.bai"
  output:
    "Spanned/{ref}/{sample}/short_paired_end/{ref}.{sample}.shortreads.overlapping.genes.tsv"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/find_short_alignments_overlapping_genes.py"

#Find alignments of real short and long reads that span exon/exon junctions without corresponding intronic sequence

rule find_shortreads_alignments:
  input:
    bam="BAMS/{ref}/{sample}/transcriptome/short_paired_end/{ref}.{sample}.transcriptome.short.sorted.bam",
    junctions="Junctions/{ref}.junctions.tsv",
    intronlengths="Junctions/{ref}.intron_lengths.tsv",
    overlapping="Spanned/{ref}/{sample}/short_paired_end/{ref}.{sample}.shortreads.overlapping.genes.tsv"
  output:
    spanningalignments="Spanned/{ref}/{sample}/short_paired_end/{ref}.{sample}.shortreads.spanningalignments.sam",
    spannedjunctions="Spanned/{ref}/{sample}/short_paired_end/{ref}.{sample}.shortreads.spannedjunctions.tsv"
  params:
    scripts=get_scripts,
    junction_overhang=config["junction_overhang"],
    insertions_threshold=config["insertions_threshold"],
    short_read_spanning_alignment_minimum_matching_bp=config["short_read_spanning_alignment_minimum_matching_bp"],
    short_read_spanning_alignment_minimum_transcriptomic_insertion_size=config["short_read_spanning_alignment_minimum_transcriptomic_insertion_size"],
    short_read_spanning_alignment_maximum_transcriptomic_insertion_size=config["short_read_spanning_alignment_maximum_transcriptomic_insertion_size"]
  shell:
    "samtools view -h {input.bam} | python {params.scripts}/find_shortreads_intronless_junction_spanning_alignments.py {params.junction_overhang} {params.insertions_threshold} {params.short_read_spanning_alignment_minimum_matching_bp} {params.short_read_spanning_alignment_minimum_transcriptomic_insertion_size} {params.short_read_spanning_alignment_maximum_transcriptomic_insertion_size} {input.junctions} {input.intronlengths} {input.overlapping} {output.spanningalignments} {output.spannedjunctions}"

rule find_longreads_alignments:
  input:
    junctions="Junctions/{ref}.junctions.unambiguous.long.tsv",
    bam="BAMS/{ref}/{sample}/transcriptome/long/{ref}.{sample}.transcriptome.long.sorted.bam"
  output:
    spanningalignments="Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spanningalignments.sam",
    spanningalignmentstxt="Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spanningalignments.txt",
    spannedjunctions="Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.tsv",
    flankregions="Flanks/{ref}/{sample}/long/{ref}.{sample}.longreads.flankregions.tsv"
  params:
    scripts=get_scripts,
    junction_overhang=config["junction_overhang"],
    insertions_threshold=config["insertions_threshold"]
  shell:
    "samtools view -h {input.bam} | python {params.scripts}/find_longreads_intronless_junction_spanning_alignments.py {params.junction_overhang} {params.insertions_threshold} {input.junctions} {output.spanningalignments} {output.spanningalignmentstxt} {output.spannedjunctions} {output.flankregions}"

#Summarize alignments of real short and long reads

rule summarize_alignments_short:
  input:
    "Junctions/{ref}.junctions.tsv",
    "Spanned/{ref}/{sample}/short_paired_end/{ref}.{sample}.shortreads.spannedjunctions.tsv"
  output:
    "Summary/{ref}/{sample}/short_paired_end/{ref}.{sample}.{strongthreshold}.shortreads.coverage.across.junctions.tsv",
    "Summary/{ref}/{sample}/short_paired_end/{ref}.{sample}.{strongthreshold}.shortreads.numspannedjunctions.tsv",
    "Summary/{ref}/{sample}/short_paired_end/{ref}.{sample}.{strongthreshold}.shortreads.no_non_overlapping.coverage.across.junctions.tsv",
    "Summary/{ref}/{sample}/short_paired_end/{ref}.{sample}.{strongthreshold}.shortreads.no_non_overlapping.numspannedjunctions.tsv",
    "Summary/{ref}/{sample}/short_paired_end/{ref}.{sample}.{strongthreshold}.shortreads.no_alternate.coverage.across.junctions.tsv",
    "Summary/{ref}/{sample}/short_paired_end/{ref}.{sample}.{strongthreshold}.shortreads.no_alternate.numspannedjunctions.tsv"
  params:
    scripts=get_scripts,
  script:
    "{params.scripts}/summarize_spanned_junctions_short.py"

rule paper_results:
  input:
    "long_truepositives.txt",
    "Summary/{ref}/{sample}/short_paired_end/{ref}.{sample}.{strongthreshold}.shortreads.numspannedjunctions.tsv",
    "Summary/{ref}/{sample}/short_paired_end/{ref}.{sample}.{strongthreshold}.shortreads.no_non_overlapping.numspannedjunctions.tsv",
    "Summary/{ref}/{sample}/short_paired_end/{ref}.{sample}.{strongthreshold}.shortreads.no_alternate.numspannedjunctions.tsv"
  output:
    "Summary/{ref}/{sample}/{ref}.{sample}.{strongthreshold}.results.txt"
  params:
    scripts=get_scripts,
  script:
    "{params.scripts}/paper_results.py"

rule summarize_alignments_long:
  input:
    "Junctions/{ref}.junctions.tsv",
    "Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.tsv"
  output:
    "Summary/{ref}/{sample}/long/{ref}.{sample}.longreads.coverage.across.junctions.tsv",
    "Summary/{ref}/{sample}/long/{ref}.{sample}.longreads.numspannedjunctions.tsv"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/summarize_spanned_junctions.py"

#Make frequency table for putative retrogenes

rule make_freqtable_short:
  input:
    expand("Summary/{{ref}}/{sample}/short_paired_end/{{ref}}.{sample}.shortreads.numspannedjunctions.tsv", sample=short_samples)
  output:
    "SFS/{ref}/summary/short_paired_end/{ref}.shortreads.freqtable.tsv",
    "SFS/{ref}/summary/short_paired_end/{ref}.shortreads.freqtable.singlyspanned.tsv",
    "SFS/{ref}/summary/short_paired_end/{ref}.shortreads.freqtable.multiplyspanned.tsv"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/make_freqtable_shortreads.py"

rule make_freqtable_long:
  input:
    expand("Summary/{{ref}}/{sample}/long/{{ref}}.{sample}.longreads.numspannedjunctions.tsv", sample=long_samples)
  output:
    "SFS/{ref}/summary/long/{ref}.longreads.freqtable.tsv",
    "SFS/{ref}/summary/long/{ref}.longreads.freqtable.singlyspanned.tsv",
    "SFS/{ref}/summary/long/{ref}.longreads.freqtable.multiplyspanned.tsv"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/make_freqtable_longreads.py"

#For every spanned junction that is supported by enough reads, get the supporting readnames

#rule get_supporting_readnames:
#  input:
#    "Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.tsv"
#  output:
#    "Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.readnames.txt"
#  params:
#    scripts=get_scripts
#  script:
#    "{params.scripts}/get_supporting_readnames.py"

#rule subset_transcriptome_bam:
#  input:
#    reads="Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.readnames.txt",
#    #bam="BAMS/{ref}/{sample}/transcriptome/long/{ref}.{sample}.transcriptome.long.sorted.bam"
#    bam="Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spanningalignments.sam"
#  output:
#    "Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.sam"
#  shell:
#    "samtools view -h -N {input.reads} {input.bam} > {output}"

rule subset_transcriptome_bam:
  input:
    "Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.tsv",
    "Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spanningalignments.sam",
    "Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spanningalignments.txt"
  output:
    "Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.readnames.txt",
    "Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.sam"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/subset_transcriptome_bam.py"

rule subset_genome_bam:
  input:
    reads="Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.readnames.txt",
    bam="BAMS/{ref}/{sample}/genome/long/{ref}.{sample}.genome.long.sorted.bam"
  output:
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.sam"
  shell:
    "samtools view -h -N {input.reads} {input.bam} > {output}"

rule add_geneid_to_subset_sams:
  input:
    "SFS/{ref}/summary/long/{ref}.longreads.freqtable.tsv",
    "Transcriptome/{ref}.transcriptome.coords.tsv",
    "Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.sam",
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.sam"
  output:
    "Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.geneid.sam",
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sam"
  params:
    scripts=get_scripts,
    long_samples=long_samples
  script:
    "{params.scripts}/add_geneid_to_subset_sam.py"

rule sort_genome_geneid_sam:
  input:
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sam"
  output:
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sorted.sam"
  shell:
    "sort -k1,1 {input} > {output}"

rule sort_transcriptome_geneid_sam:
  input:
    "Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.geneid.sam"
  output:
    "Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.geneid.sorted.sam"
  shell:
    "sort -k1,1 {input} > {output}"

rule analyze_anchors:
  input:
    "Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.geneid.sorted.sam",
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sorted.sam",
  output:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.AS"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/analyze_anchors.py"

rule sort_anchors:
  input:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.AS"
  output:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.sorted.AS"
  shell:
    "cut -f1-13 {input} | sort -k1,1 -k2,2 -k9,9 -k10,10n -k11,11n > {output}"

rule cluster_anchors:
  input:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.sorted.AS"
  output:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.AS"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/cluster_anchors.py"

rule best_AS:
  input:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.AS"
  output:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.best.AS"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/best_AS.py"

rule categorize_AS:
  input:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.best.AS",
    "Transcriptome/{ref}.transcriptome.coords.tsv"
  output:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.best.AS.diff",
    "AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.best.AS.same"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/categorize_AS.py"

#rule summary_short_retrogene_genotypes:
#  input:

rule summary_long_retrogene_genotypes:
  input:
    expand("AS/{{ref}}/{sample}/{{ref}}.{sample}.genome.clustered.best.AS.diff", sample=long_samples)
  output:
    "RESULTS/{ref}/summary/long/{ref}.summary.long.retrogene.genotypes.tsv"
  params:
    scripts=get_scripts,
    long_samples=long_samples
  script:
    "{params.scripts}/summary_retrogene_genotypes.py"

rule final_genome_alignments:
  input:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.best.AS.diff",
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.sam",
    "AS/{ref}/{sample}/{ref}.{sample}.genome.AS"
  output:
    temp("RESULTS/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supporting_alignments.genome.sam")
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/final_genome_alignments.py"

rule bam_final_alignments:
  input:
    "RESULTS/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supporting_alignments.genome.sam"
  output:
    temp("RESULTS/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supporting_alignments.genome.bam")
  threads: 16
  shell:
    "samtools view -@ {threads} -b -o {output} {input}"

rule sort_final_alignments:
  input:
    "RESULTS/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supporting_alignments.genome.bam"
  output:
    "RESULTS/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supporting_alignments.genome.sorted.bam"
  threads: 16
  shell:
    "samtools sort -@ {threads} -o {output} {input}"

rule index_final_alignments:
  input:
    "RESULTS/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supporting_alignments.genome.sorted.bam"
  output:
    "RESULTS/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supporting_alignments.genome.sorted.bam.bai"
  threads: 16
  shell:
    "samtools index -@ {threads} {input}"

rule get_retrogene_consensus_sequences:
  input:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.best.AS.diff",
    "RESULTS/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supporting_alignments.genome.sorted.bam",
    "RESULTS/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supporting_alignments.genome.sorted.bam.bai"
  output:
    "RESULTS/{ref}/{sample}/long/consensus/{ref}.{sample}.retrogene.consensus.txt"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/get_retrogene_consensus_sequences.py"

rule get_retrogene_transcript_fastas:
  input:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.best.AS.diff",
    "Reference/{ref}.rna_from_genomic.fna.gz",
    "RESULTS/{ref}/{sample}/long/consensus/{ref}.{sample}.retrogene.consensus.txt"
  output:
    "RESULTS/{ref}/{sample}/long/consensus/{ref}.{sample}.transcript.fastas.txt"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/get_retrogene_transcript_fastas.py"

rule sample_long_retrogene_results:
  input:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.best.AS.diff"
  output:
    "RESULTS/{ref}/{sample}/long/{ref}.{sample}.long.retrogenes.tsv"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/sample_long_retrogene_results.py"

rule categorize_AS_summary:
  input:
    expand("AS/{{ref}}/{sample}/{{ref}}.{sample}.genome.clustered.AS", sample=long_samples),
    "Transcriptome/{ref}.transcriptome.coords.tsv"
  output:
    "AS/{ref}/SUMMARY/{ref}.SUMMARY.alldiff.txt",
    "AS/{ref}/SUMMARY/{ref}.SUMMARY.allsame.txt",
    "AS/{ref}/SUMMARY/{ref}.SUMMARY.somediffsomesame.txt"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/categorize_AS_summary.py"

rule get_GeneID_introns:
  input:
    "Summary/{ref}/{sample}/long/{ref}.{sample}.longreads.coverage.across.junctions.tsv",
    "Junctions/{ref}.introns.tsv"
  output:
    "Analysis/{ref}/{sample}/long/{ref}.{sample}.introns"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/get_GeneID_introns.py"

rule missing_intron_percent:
  input:
    "Analysis/{ref}/{sample}/long/{ref}.{sample}.introns",
    "Spanned/{ref}/{sample}/long/{ref}.{sample}.longreads.spannedjunctions.tsv",
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sorted.sam"
  output:
    "Analysis/{ref}/{sample}/{ref}.{sample}.introns.missing_percent"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/missing_intron_percent.py"

rule summarize_missing_intron:
  input:
    "AS/{ref}/SUMMARY/{ref}.SUMMARY.allsame.txt",
    "Analysis/{ref}/{sample}/{ref}.{sample}.introns.missing_percent"
  output:
    "Analysis/{ref}/{sample}/{ref}.{sample}.intronsummary",
    "Analysis/{ref}/{sample}/{ref}.{sample}.intronweirdos"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/summarize_missing_intron.py"

#checkpoint get_GeneIDs:
#  input:
#    "SFS/{ref}/long/{ref}.longreads.freqtable.tsv"
#  output:
#    directory(GeneIDs)
#  params:
#    scripts=get_scripts
#  script:
#    "{params.scripts}/get_GeneIDs.py"

#rule get_GeneID_readnames:
#  input:
#    "Transcriptome/{ref}.transcriptome.coords.tsv",
#    "Spanned/{ref}/long/{ref}.{sample}.longreads.spannedjunctions.tsv"
#  output:
#    directory(GeneIDs)
#  params:
#    scripts=get_scripts
#  script:
#    "{params.scripts}/get_GeneID_readnames.py"
  
#def aggregate_input(wildcards):
#  checkpoint_output = checkpoints.get_GeneID_readnames.get(**wildcards).output[0]
#  return expand("", sample=, geneID=glob_wildcards(os.path.join(checkpoint_output, "")).geneID)

#rule get_GeneID_readnames:
#  input:
#    "../RetroGenes_short/Transcript/{ref}.transcript.coords.tsv",
#    "Spanned/{ref}/long/{ref}.{sample}.longreads.spannedjunctions.tsv",
#  output:
#    "Analysis/{ref}/{sample}/{ref}.{sample}.{geneid}.names"
#  params:
#    "{geneid}"
#  script:
#    "scripts/get_GeneID_readnames.py"

#def aggregate_input(wildcards):
#	checkpoint_output = checkpoints.get_GeneIDS.get(**wildcards).output[0]
#	return expand("AS/{ref}/{sample}/{ref}.{sample}.{geneid}.genome.clustered.AS", ref=wildcards.ref, sample=wildcards.sample, geneid=glob_wildcards(os.path.join(checkpoint_output, "{geneid}.txt")).geneid)

#rule categorize_AS2:
#  input:
#    expand("AS/{{ref}}/{sample}/{{ref}}.{sample}.{geneid}.genome.clustered.AS", geneid=config["geneids"], sample=long_samples),
#    "../RetroGenes_short/Transcript/{transcript}.transcript.coords.tsv"
#  output:
#    "AS/{transcript}/SUMMARY/{transcript}.SUMMARY.alldiff.txt",
#    "AS/{transcript}/SUMMARY/{transcript}.SUMMARY.allsame.txt",
#    "AS/{transcript}/SUMMARY/{transcript}.SUMMARY.somediffsomesame.txt"
#  script:
#    "scripts/categorize_AS.py"
