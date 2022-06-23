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
    #expand("Reference/{ref}.genomic.fna.gz", ref=config["ref"])
    #expand("Reference/{ref}.genomic.gtf.gz", ref=config["ref"])
    #expand("Reference/{ref}.rna_from_genomic.fna.gz", ref=config["ref"])
    #get_all_short_fastq
    #get_all_hifi_fastq
    #get_all_pb_fastq
    #get_all_ont_fastq
    #expand("Reference/{ref}.rna_from_genomic.fna.gz.bwt", ref=config["ref"])
    #expand("Transcript/{ref}.transcript.coords.tsv", ref=config["ref"])
    #expand("Junctions/{ref}.shortregions.fasta", ref=config["ref"])
    #expand("Junctions/{ref}.longregions.fasta", ref=config["ref"])
    #expand("BAMS/{ref}/{sample}/short_paired_end/{ref}.{sample}.short.sorted.bam", ref=config["ref"], sample=config["fastqs"])
    #expand("BAMS/{ref}/{sample}/transcript/long/{ref}.{sample}.transcript.long.sorted.bam", ref=config["ref"], sample=config["fastqs"])
    #expand("BAMS/{ref}/{sample}/genome/long/{ref}.{sample}.genome.long.sorted.bam", ref=config["ref"], sample=config["fastqs"])
    #expand("Ambiguous/{ref}.ambiguous.long.spannedjunctions.tsv", ref=config["ref"])
    #expand("Spanned/{ref}/long/{ref}.{sample}.longreads.spannedjunctions.tsv", ref=config["ref"], sample=long_samples)
    #expand("Summary/{ref}/long/{ref}.{sample}.longreads.numspannedjunctions.tsv", ref=config["ref"], sample=long_samples)
    #expand("SFS/{ref}/long/{ref}.longreads.freqtable.tsv", ref=config["ref"])
    #expand("Spanned/{ref}/long/{ref}.{sample}.longreads.spannedjunctions.readnames.txt", ref=config["ref"], sample=long_samples)
    #expand("Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.bam", ref=config["ref"], sample=long_samples)
    #expand("AS/{ref}/SUMMARY/{ref}.SUMMARY.alldiff.txt", ref=config["ref"])
    #expand("Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sam", ref=config["ref"], sample=long_samples)
    #expand("Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sorted.sam", ref=config["ref"], sample=long_samples)
    #expand("Anchor/{ref}/{sample}/{ref}.{sample}.transcript.long.subset.geneid.sorted.sam", ref=config["ref"], sample=long_samples)
    #expand("AS/{ref}/{sample}/{ref}.{sample}.genome.AS", ref=config["ref"], sample=long_samples)
    #expand("AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.AS", ref=config["ref"], sample=long_samples)
    #expand("AS/{ref}/SUMMARY/{ref}.SUMMARY.alldiff.txt", ref=config["ref"], sample=long_samples)
    #expand("Analysis/{ref}/{sample}/{ref}.{sample}.introns", ref=config["ref"], sample=long_samples)
    #expand("Analysis/{ref}/{sample}/{ref}.{sample}.introns.missing_percent", ref=config["ref"], sample=long_samples)
    expand("Analysis/{ref}/{sample}/{ref}.{sample}.intronsummary", ref=config["ref"], sample=long_samples)

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

def get_short_fastq(wildcards):
  return config["fastqs"][wildcards.sample]["short_paired_end"][int(wildcards.i)][int(wildcards.j)]

rule link_short:
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

rule bwa_index_transcript:
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

rule minimap2_hifi_index_transcript:
  input:
    "Reference/{ref}.rna_from_genomic.fna.gz"
  output:
    "Reference/{ref}.rna_from_genomic.pacbio_hifi.mmi"
  threads: 3
  shell:
    "minimap2 -x map-hifi -d {output} {input}"

rule minimap2_pb_index_transcript:
  input:
    "Reference/{ref}.rna_from_genomic.fna.gz"
  output:
    "Reference/{ref}.rna_from_genomic.pacbio_clr.mmi"
  threads: 3
  shell:
    "minimap2 -x map-pb -d {output} {input}"

rule minimap2_ont_index_transcript:
  input:
    "Reference/{ref}.rna_from_genomic.fna.gz"
  output:
    "Reference/{ref}.rna_from_genomic.ont.mmi"
  threads: 3
  shell:
    "minimap2 -x map-ont -d {output} {intput}"

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

def get_scripts(wildcards):
  return config["scripts_directory"]

rule get_transcript_junctions:
  input:
    "Reference/{ref}.rna_from_genomic.fna.gz"
  output:
    "Junctions/{ref}.junctions.tsv",
    "Junctions/{ref}.shortregions.txt",
    "Transcript/{ref}.transcript.coords.tsv",
    "Junctions/{ref}.introns.tsv",
    "Junctions/{ref}.longregions.txt"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/get_transcript_junctions.py"

#this step requires the reference to be zipped with bgzip
rule get_junction_fasta_short:
  input:
    ref="Reference/{ref}.genomic.fna.gz",
    regions="Junctions/{ref}.shortregions.txt"
  output:
    "Junctions/{ref}.shortregions.fasta"
  shell:
    "samtools faidx {input.ref} -r {input.regions} > {output}"

#this step requires the reference to be zipped with bgzip
rule get_junction_fasta_long:
  input:
    ref="Reference/{ref}.genomic.fna.gz",
    regions="Junctions/{ref}.longregions.txt"
  output:
    "Junctions/{ref}.longregions.fasta"
  shell:
    "samtools faidx {input.ref} -r {input.regions} > {output}"

rule bwa_aln_transcript:
  input:
    "Reference/{ref}.rna_from_genomic.fna.gz.amb",
    "Reference/{ref}.rna_from_genomic.fna.gz.ann",
    "Reference/{ref}.rna_from_genomic.fna.gz.bwt",
    "Reference/{ref}.rna_from_genomic.fna.gz.pac",
    "Reference/{ref}.rna_from_genomic.fna.gz.sa",
    ref="Reference/{ref}.rna_from_genomic.fna.gz",
    reads="FASTQS/{sample}/short_paired_end/{sample}.short.{i}.R{j}.fastq.gz"
  output:
    temp("BAMS/{ref}/{sample}/short_paired_end/{ref}.{sample}.short.{i}.{j}.sai")
  threads: 16
  shell:
    "bwa aln -t {threads} {input.ref} {input.reads} > {output}"

rule bwa_sampe_transcript:
  input:
    ref="Reference/{ref}.rna_from_genomic.fna.gz",
    sai1="BAMS/{ref}/{sample}/short_paired_end/{ref}.{sample}.short.{i}.1.sai",
    sai2="BAMS/{ref}/{sample}/short_paired_end/{ref}.{sample}.short.{i}.2.sai",
    read1="FASTQS/{sample}/short_paired_end/{sample}.short.{i}.R1.fastq.gz",
    read2="FASTQS/{sample}/short_paired_end/{sample}.short.{i}.R2.fastq.gz"
  output:
    temp("BAMS/{ref}/{sample}/short_paired_end/{ref}.{sample}.short.{i}.sam")
  shell:
    "bwa sampe {input.ref} {input.sai1} {input.sai2} {input.read1} {input.read2} > {output}"

rule sam_to_bam_short:
  input:
    "BAMS/{ref}/{sample}/short_paired_end/{ref}.{sample}.short.{i}.sam"
  output:
    temp("BAMS/{ref}/{sample}/short_paired_end/{ref}.{sample}.short.{i}.bam")
  threads: 16
  shell:
    "samtools view -F 4 -@ {threads} -b -o {output} {input}"

rule sort_bam_short:
  input:
    "BAMS/{ref}/{sample}/short_paired_end/{ref}.{sample}.short.{i}.bam"
  output:
    temp("BAMS/{ref}/{sample}/short_paired_end/{ref}.{sample}.short.{i}.sorted.bam")
  threads: 16
  shell:
    "samtools sort -@ {threads} -o {output} {input}"

def get_sorted_bam_short(wildcards):
  files = []
  if "short_paired_end" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["short_paired_end"]
    for i in i_s:
      file = f"BAMS/{wildcards.ref}/{wildcards.sample}/short_paired_end/{wildcards.ref}.{wildcards.sample}.short.{i}.sorted.bam"
      files.append(file)
  return files

rule merge_bam_short:
  input:
    get_sorted_bam_short
  output:
    "BAMS/{ref}/{sample}/short_paired_end/{ref}.{sample}.short.sorted.bam"
  threads: 16
  shell:
    "samtools merge -@ {threads} -o {output} {input}"

rule minimap_transcript_long:
	input:
		mmi="Reference/{ref}.rna_from_genomic.{tech}.mmi",
		reads="FASTQS/{sample}/{tech}/{sample}.{tech}.{i}.fastq.gz"
	output:
		temp("BAMS/{ref}/{sample}/transcript/{tech}/{ref}.{sample}.transcript.{tech}.{i}.sam")
	threads: 15
	shell:
		"minimap2 -t {threads} -a {input.mmi} {input.reads} > {output}"

rule sam_to_bam_transcript_long:
	input:
		"BAMS/{ref}/{sample}/transcript/{tech}/{ref}.{sample}.transcript.{tech}.{i}.sam"
	output:
		temp("BAMS/{ref}/{sample}/transcript/{tech}/{ref}.{sample}.transcript.{tech}.{i}.bam")
	threads: 16
	shell:
		"samtools view -F 4 -@ {threads} -b -o {output} {input}"

rule sort_bam_transcript_long:
	input:
		"BAMS/{ref}/{sample}/transcript/{tech}/{ref}.{sample}.transcript.{tech}.{i}.bam"
	output:
		temp("BAMS/{ref}/{sample}/transcript/{tech}/{ref}.{sample}.transcript.{tech}.{i}.sorted.bam")
	threads: 16
	shell:
		"samtools sort -@ {threads} -o {output} {input}"

def get_sorted_bam_transcript_long(wildcards):
  files = []
  if "pacbio_hifi" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["pacbio_hifi"]
    for i in i_s:
      file = f"BAMS/{wildcards.ref}/{wildcards.sample}/transcript/pacbio_hifi/{wildcards.ref}.{wildcards.sample}.transcript.pacbio_hifi.{i}.sorted.bam"
      files.append(file)
  if "pacbio_clr" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["pacbio_clr"]
    for i in i_s:
      file = f"BAMS/{wildcards.ref}/{wildcards.sample}/transcript/pacbio_clr/{wildcards.ref}.{wildcards.sample}.transcript.pacbio_clr.{i}.sorted.bam"
      files.append(file)
  if "ont" in config["fastqs"][wildcards.sample]:
    i_s = config["fastqs"][wildcards.sample]["ont"]
    for i in i_s:
      file = f"BAMS/{wildcards.ref}/{wildcards.sample}/transcript/ont/{wildcards.ref}.{wildcards.sample}.transcript.ont.{i}.sorted.bam"
      files.append(file)
  return files

rule merge_bam_transcript_long:
	input:
		get_sorted_bam_transcript_long
	output:
		"BAMS/{ref}/{sample}/transcript/long/{ref}.{sample}.transcript.long.sorted.bam"
	threads: 16
	shell:
		"samtools merge -@ {threads} -o {output} {input}"

rule minimap_genome_long:
	input:
		mmi="Reference/{ref}.genomic.{tech}.mmi",
		reads="FASTQS/{sample}/{tech}/{sample}.{tech}.{i}.fastq.gz"
	output:
		temp("BAMS/{ref}/{sample}/genome/{tech}/{ref}.{sample}.genome.{tech}.{i}.sam")
	threads: 15
	shell:
		"minimap2 -t {threads} -a {input.mmi} {input.reads} > {output}"

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

rule minimap2_junctionreads_to_transcript:
	input:
		mmi="Reference/{ref}.rna_from_genomic.pacbio_hifi.mmi",
		reads="Junctions/{ref}.longregions.fasta"
	output:
		"BAMS/{ref}/junction/{ref}.junctionreads.long.transcript.sam"
	threads: 15
	shell:
		"minimap2 --eqx -t {threads} -a {input.mmi} {input.reads} > {output}"

rule find_junctionreads_alignments_long:
  input:
    junctions="Junctions/{ref}.junctions.tsv",
    bam="BAMS/{ref}/junction/{ref}.junctionreads.long.transcript.sam"
  output:
    spanningalignments="Ambiguous/{ref}.ambiguous.long.spanningalignments.sam",
    spannedjunctions="Ambiguous/{ref}.ambiguous.long.spannedjunctions.tsv",
    flankregions="Flanks/{ref}/junction/{ref}.junction.flankregions.tsv"
  params:
    scripts=get_scripts
  shell:
    "samtools view -h {input.bam} | python {params.scripts}/find_longreads_intronless_junction_spanning_alignments.py {input.junctions} {output.spanningalignments} {output.spannedjunctions} {output.flankregions}"

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

rule find_longreads_alignments:
  input:
    junctions="Junctions/{ref}.junctions.unambiguous.long.tsv",
    bam="BAMS/{ref}/{sample}/transcript/long/{ref}.{sample}.transcript.long.sorted.bam"
  output:
    spanningalignments="Spanned/{ref}/long/{ref}.{sample}.longreads.spanningalignments.sam",
    spannedjunctions="Spanned/{ref}/long/{ref}.{sample}.longreads.spannedjunctions.tsv",
    flankregions="Flanks/{ref}/long/{ref}.{sample}.longreads.flankregions.tsv"
  params:
    scripts=get_scripts
  shell:
    "samtools view -h {input.bam} | python {params.scripts}/find_longreads_intronless_junction_spanning_alignments.py {input.junctions} {output.spanningalignments} {output.spannedjunctions} {output.flankregions}"

rule summarize_alignments_long:
  input:
    "Transcript/{ref}.transcript.coords.tsv",
    "Junctions/{ref}.junctions.tsv",
    "Spanned/{ref}/long/{ref}.{sample}.longreads.spannedjunctions.tsv"
  output:
    "Summary/{ref}/long/{ref}.{sample}.longreads.coverage.across.junctions.tsv",
    "Summary/{ref}/long/{ref}.{sample}.longreads.numspannedjunctions.tsv"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/summarize_spanned_junctions.py"

rule make_freqtable_long:
  input:
    expand("Summary/{{ref}}/long/{{ref}}.{sample}.longreads.numspannedjunctions.tsv", sample=long_samples)
  output:
    "SFS/{ref}/long/{ref}.longreads.freqtable.tsv",
    "SFS/{ref}/long/{ref}.longreads.freqtable.singlyspanned.tsv",
    "SFS/{ref}/long/{ref}.longreads.freqtable.multiplyspanned.tsv"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/make_freqtable_longreads.py"

rule get_supporting_readnames:
  input:
    "Spanned/{ref}/long/{ref}.{sample}.longreads.spannedjunctions.tsv"
  output:
    "Spanned/{ref}/long/{ref}.{sample}.longreads.spannedjunctions.readnames.txt"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/get_supporting_readnames.py"

rule subset_genome_bam:
  input:
    reads="Spanned/{ref}/long/{ref}.{sample}.longreads.spannedjunctions.readnames.txt",
    bam="BAMS/{ref}/{sample}/genome/long/{ref}.{sample}.genome.long.sorted.bam"
  output:
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.sam"
  shell:
    "samtools view -N {input.reads} {input.bam} > {output}"

rule add_geneid_to_subset_sam:
  input:
    "SFS/{ref}/long/{ref}.longreads.freqtable.tsv",
    "Transcript/{ref}.transcript.coords.tsv",
    "Spanned/{ref}/long/{ref}.{sample}.longreads.spanningalignments.sam",
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.sam"
  output:
    "Anchor/{ref}/{sample}/{ref}.{sample}.transcript.long.subset.geneid.sam",
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sam"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/add_geneid_to_subset_sam.py"

rule sort_subset_geneid_sam:
  input:
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sam"
  output:
    "Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sorted.sam"
  shell:
    "sort -k1,1 {input} > {output}"

rule sort_transcript_geneid_sam:
  input:
    "Anchor/{ref}/{sample}/{ref}.{sample}.transcript.long.subset.geneid.sam"
  output:
    "Anchor/{ref}/{sample}/{ref}.{sample}.transcript.long.subset.geneid.sorted.sam"
  shell:
    "sort -k1,1 {input} > {output}"

rule analyze_anchors:
  input:
    "Anchor/{ref}/{sample}/{ref}.{sample}.transcript.long.subset.geneid.sorted.sam",
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
    "sort -k1,1 -k6,6n -k7,7n {input} > {output}"

rule cluster_anchors:
  input:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.sorted.AS"
  output:
    "AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.AS"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/cluster_anchors.py"

rule categorize_AS:
  input:
    expand("AS/{{ref}}/{sample}/{{ref}}.{sample}.genome.clustered.AS", sample=long_samples),
    "Transcript/{ref}.transcript.coords.tsv"
  output:
    "AS/{ref}/SUMMARY/{ref}.SUMMARY.alldiff.txt",
    "AS/{ref}/SUMMARY/{ref}.SUMMARY.allsame.txt",
    "AS/{ref}/SUMMARY/{ref}.SUMMARY.somediffsomesame.txt"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/categorize_AS.py"

rule get_GeneID_introns:
  input:
    "Summary/{ref}/long/{ref}.{sample}.longreads.coverage.across.junctions.tsv",
    "Junctions/{ref}.introns.tsv"
  output:
    "Analysis/{ref}/{sample}/{ref}.{sample}.introns"
  params:
    scripts=get_scripts
  script:
    "{params.scripts}/get_GeneID_introns.py"

rule missing_intron_percent:
  input:
    "Analysis/{ref}/{sample}/{ref}.{sample}.introns",
    "Spanned/{ref}/long/{ref}.{sample}.longreads.spannedjunctions.tsv",
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
#    "Transcript/{ref}.transcript.coords.tsv",
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
