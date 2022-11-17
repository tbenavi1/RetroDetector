#bgzip reference genome

def get_genome(wildcards):
  return config["ref"][wildcards.ref]["genome"]

rule unzip_genome:
  input:
    get_genome
  output:
    "results/Genome/{ref}/{ref}.genome.fna"
  shell:
    "gunzip -c {input} > {output}"

rule bgzip_genome:
  input:
    "results/Genome/{ref}/{ref}.genome.fna"
  output:
    "results/Genome/{ref}/{ref}.genome.fna.gz"
  shell:
    "bgzip {input}"

#bgzip reference transcriptome

def get_transcriptome(wildcards):
  return config["ref"][wildcards.ref]["transcriptome"]

rule unzip_transcriptome:
  input:
    get_transcriptome
  output:
    "results/Transcriptome/{ref}/{ref}.transcriptome.fna"
  shell:
    "gunzip -c {input} > {output}"

rule bgzip_transcriptome:
  input:
    "results/Transcriptome/{ref}/{ref}.transcriptome.fna"
  output:
    "results/Transcriptome/{ref}/{ref}.transcriptome.fna.gz"
  shell:
    "bgzip {input}"

#samtools index reference genome

rule samtools_index_genome:
  input:
    "results/Genome/{ref}/{ref}.genome.fna.gz"
  output:
    "results/Genome/{ref}/{ref}.genome.fna.gz.fai",
    "results/Genome/{ref}/{ref}.genome.fna.gz.gzi"
  shell:
    "samtools faidx {input}"

#process transcriptome to get multi-exon transcriptome, coords, junctions, introns, and intron lengths

rule process_transcriptome:
  input:
    "results/Transcriptome/{ref}/{ref}.transcriptome.fna.gz"
  output:
    temp("results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna"),
    "results/Transcriptome/{ref}/{ref}.transcriptome.coords.tsv",
    "results/Transcriptome/{ref}/{ref}.transcriptome.junctions.tsv",
    "results/Transcriptome/{ref}/{ref}.transcriptome.introns.tsv",
    "results/Transcriptome/{ref}/{ref}.transcriptome.intronlengths.tsv"
  script:
    "../scripts/process_transcriptome.py"

rule bgzip_multiexontranscriptome:
  input:
    "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna"
  output:
    "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz"
  shell:
    "bgzip {input}"
