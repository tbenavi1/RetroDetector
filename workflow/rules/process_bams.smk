#Sort and merge transcriptome bams for real short reads

rule sam_to_bam_transcriptome_short:
  input:
    "results/BAMS/{ref}/{sample}/transcriptome/short/{ref}.{sample}.transcriptome.illumina.{i}.sam"
  output:
    temp("results/BAMS/{ref}/{sample}/transcriptome/short/{ref}.{sample}.transcriptome.illumina.{i,[0-9]+}.bam")
  threads: 32
  shell:
    "samtools view -F 4 -@ {threads} -b -o {output} {input}"

rule sort_bam_transcriptome_short:
  input:
    "results/BAMS/{ref}/{sample}/transcriptome/short/{ref}.{sample}.transcriptome.illumina.{i}.bam"
  output:
    temp("results/BAMS/{ref}/{sample}/transcriptome/short/{ref}.{sample}.transcriptome.illumina.{i,[0-9]+}.sorted.bam")
  threads: 32
  shell:
    "samtools sort -@ {threads} -o {output} {input}"

def get_sorted_bam_transcriptome_short(wildcards):
  files = []
  if "short" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["short"]
    for i in i_s:
      file = f"results/BAMS/{wildcards.ref}/{wildcards.sample}/transcriptome/short/{wildcards.ref}.{wildcards.sample}.transcriptome.illumina.{i}.sorted.bam"
      files.append(file)
  return files

rule merge_bam_transcriptome_short:
  input:
    get_sorted_bam_transcriptome_short
  output:
    "results/BAMS/{ref}/{sample}/transcriptome/short/{ref}.{sample}.transcriptome.short.sorted.bam"
  threads: 32
  shell:
    "samtools merge -@ {threads} -o {output} {input}"

#sort and merge transcriptome bams for real long reads

rule sam_to_bam_transcriptome_long:
	input:
		"results/BAMS/{ref}/{sample}/transcriptome/{tech}/{ref}.{sample}.transcriptome.{tech}.{i}.sam"
	output:
		temp("results/BAMS/{ref}/{sample}/transcriptome/{tech}/{ref}.{sample}.transcriptome.temp.{tech}.{i,[0-9]+}.bam")
	threads: 32
	shell:
		"samtools view -F 4 -@ {threads} -b -o {output} {input}"

rule sort_bam_transcriptome_long:
	input:
		"results/BAMS/{ref}/{sample}/transcriptome/{tech}/{ref}.{sample}.transcriptome.temp.{tech}.{i}.bam"
	output:
		temp("results/BAMS/{ref}/{sample}/transcriptome/{tech}/{ref}.{sample}.transcriptome.temp.{tech}.{i,[0-9]+}.sorted.bam")
	threads: 32
	shell:
		"samtools sort -@ {threads} -o {output} {input}"

def get_sorted_bam_transcriptome_long(wildcards):
  files = []
  if "pacbio_hifi" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["pacbio_hifi"]
    for i in i_s:
      file = f"results/BAMS/{wildcards.ref}/{wildcards.sample}/transcriptome/pacbio_hifi/{wildcards.ref}.{wildcards.sample}.transcriptome.temp.pacbio_hifi.{i}.sorted.bam"
      files.append(file)
  if "pacbio_clr" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["pacbio_clr"]
    for i in i_s:
      file = f"results/BAMS/{wildcards.ref}/{wildcards.sample}/transcriptome/pacbio_clr/{wildcards.ref}.{wildcards.sample}.transcriptome.temp.pacbio_clr.{i}.sorted.bam"
      files.append(file)
  if "ont" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["ont"]
    for i in i_s:
      file = f"results/BAMS/{wildcards.ref}/{wildcards.sample}/transcriptome/ont/{wildcards.ref}.{wildcards.sample}.transcriptome.temp.ont.{i}.sorted.bam"
      files.append(file)
  return files

rule merge_bam_transcriptome_long:
	input:
		get_sorted_bam_transcriptome_long
	output:
		"results/BAMS/{ref}/{sample}/transcriptome/long/{ref}.{sample}.transcriptome.long.sorted.bam"
	threads: 32
	shell:
		"samtools merge -@ {threads} -o {output} {input}"

#Sort, merge, and index genome bams for real short reads

rule sam_to_bam_genome_short:
  input:
    "results/BAMS/{ref}/{sample}/genome/short/{ref}.{sample}.genome.illumina.{i}.sam"
  output:
    temp("results/BAMS/{ref}/{sample}/genome/short/{ref}.{sample}.genome.illumina.{i,[0-9]+}.bam")
  threads: 32
  shell:
    "samtools view -F 4 -@ {threads} -b -o {output} {input}"

rule sort_bam_genome_short:
  input:
    "results/BAMS/{ref}/{sample}/genome/short/{ref}.{sample}.genome.illumina.{i}.bam"
  output:
    temp("results/BAMS/{ref}/{sample}/genome/short/{ref}.{sample}.genome.illumina.{i,[0-9]+}.sorted.bam")
  threads: 32
  shell:
    "samtools sort -@ {threads} -o {output} {input}"

def get_sorted_bam_genome_short(wildcards):
  files = []
  i_s = config["fastqs"][wildcards.sample]["short"]
  for i in i_s:
    file = f"results/BAMS/{wildcards.ref}/{wildcards.sample}/genome/short/{wildcards.ref}.{wildcards.sample}.genome.illumina.{i}.sorted.bam"
    files.append(file)
  return files

rule merge_bam_genome_short:
  input:
    get_sorted_bam_genome_short
  output:
    "results/BAMS/{ref}/{sample}/genome/short/{ref}.{sample}.genome.short.sorted.bam"
  threads: 32
  shell:
    "samtools merge -@ {threads} -o {output} {input}"

rule index_bam_genome_short:
  input:
    "results/BAMS/{ref}/{sample}/genome/short/{ref}.{sample}.genome.short.sorted.bam"
  output:
    "results/BAMS/{ref}/{sample}/genome/short/{ref}.{sample}.genome.short.sorted.bam.bai"
  threads: 32
  shell:
    "samtools index -@ {threads} {input}"

#Sort, merge, and index genome bams for real long reads

rule sam_to_bam_genome_long:
	input:
		"results/BAMS/{ref}/{sample}/genome/{tech}/{ref}.{sample}.genome.{tech}.{i}.sam"
	output:
		temp("results/BAMS/{ref}/{sample}/genome/{tech}/{ref}.{sample}.genome.{tech}.{i,[0-9]+}.bam")
	threads: 32
	shell:
		"samtools view -F 4 -@ {threads} -b -o {output} {input}"

rule sort_bam_genome_long:
	input:
		"results/BAMS/{ref}/{sample}/genome/{tech}/{ref}.{sample}.genome.{tech}.{i}.bam"
	output:
		temp("results/BAMS/{ref}/{sample}/genome/{tech}/{ref}.{sample}.genome.{tech}.{i,[0-9]+}.sorted.bam")
	threads: 32
	shell:
		"samtools sort -@ {threads} -o {output} {input}"

def get_sorted_bam_genome_long(wildcards):
  files = []
  if "pacbio_hifi" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["pacbio_hifi"]
    for i in i_s:
      file = f"results/BAMS/{wildcards.ref}/{wildcards.sample}/genome/pacbio_hifi/{wildcards.ref}.{wildcards.sample}.genome.pacbio_hifi.{i}.sorted.bam"
      files.append(file)
  if "pacbio_clr" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["pacbio_clr"]
    for i in i_s:
      file = f"results/BAMS/{wildcards.ref}/{wildcards.sample}/genome/pacbio_clr/{wildcards.ref}.{wildcards.sample}.genome.pacbio_clr.{i}.sorted.bam"
      files.append(file)
  if "ont" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["ont"]
    for i in i_s:
      file = f"results/BAMS/{wildcards.ref}/{wildcards.sample}/genome/ont/{wildcards.ref}.{wildcards.sample}.genome.ont.{i}.sorted.bam"
      files.append(file)
  return files

rule merge_bam_genome_long:
	input:
		get_sorted_bam_genome_long
	output:
		"results/BAMS/{ref}/{sample}/genome/long/{ref}.{sample}.genome.long.sorted.bam"
	threads: 32
	shell:
		"samtools merge -@ {threads} -o {output} {input}"

rule index_bam_genome_long:
  input:
    "results/BAMS/{ref}/{sample}/genome/long/{ref}.{sample}.genome.long.sorted.bam"
  output:
    "results/BAMS/{ref}/{sample}/genome/long/{ref}.{sample}.genome.long.sorted.bam.bai"
  threads: 32
  shell:
    "samtools index -@ {threads} {input}"
