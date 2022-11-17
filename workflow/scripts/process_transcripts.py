import gzip
import subprocess

with gzip.open("/nobackup/rogers_research/Rhyker/rogers_temp/Rhyker/RetroGenes_dmel/Reference/ncbi_dmel.rna_from_genomic.fna.gz", "rt") as input_file:
	for line in input_file:
		if line.startswith(">"):
			transcript = line.strip().split()[0][1:]
			subprocess.run(f"mkdir '/nobackup/rogers_research/Rhyker/rogers_temp/Rhyker/RetroGenes_dmel/BAMS/ncbi_dmel/ISO1/transcriptome/short_paired_end/broken/{transcript}'", shell=True)
			subprocess.run(f"bwa mem -t 63 -B 3 'Transcriptome/ncbi_dmel_transcripts/{transcript}.fa' FASTQS/ISO1/short_paired_end/ISO1.short.1.R1.fastq.gz FASTQS/ISO1/short_paired_end/ISO1.short.1.R2.fastq.gz > 'BAMS/ncbi_dmel/ISO1/transcriptome/short_paired_end/broken/{transcript}/ncbi_dmel.ISO1.transcriptome.illumina.1.sam'", shell=True)
			subprocess.run(f"samtools view -b -F 4 '/nobackup/rogers_research/Rhyker/rogers_temp/Rhyker/RetroGenes_dmel/BAMS/ncbi_dmel/ISO1/transcriptome/short_paired_end/broken/{transcript}/ncbi_dmel.ISO1.transcriptome.illumina.1.sam' -o '/nobackup/rogers_research/Rhyker/rogers_temp/Rhyker/RetroGenes_dmel/BAMS/ncbi_dmel/ISO1/transcriptome/short_paired_end/broken/{transcript}/ncbi_dmel.ISO1.transcriptome.illumina.1.bam'", shell=True)
			subprocess.run(f"rm '/nobackup/rogers_research/Rhyker/rogers_temp/Rhyker/RetroGenes_dmel/BAMS/ncbi_dmel/ISO1/transcriptome/short_paired_end/broken/{transcript}/ncbi_dmel.ISO1.transcriptome.illumina.1.sam'", shell=True)
			subprocess.run(f"samtools sort '/nobackup/rogers_research/Rhyker/rogers_temp/Rhyker/RetroGenes_dmel/BAMS/ncbi_dmel/ISO1/transcriptome/short_paired_end/broken/{transcript}/ncbi_dmel.ISO1.transcriptome.illumina.1.bam' -o '/nobackup/rogers_research/Rhyker/rogers_temp/Rhyker/RetroGenes_dmel/BAMS/ncbi_dmel/ISO1/transcriptome/short_paired_end/broken/{transcript}/ncbi_dmel.ISO1.transcriptome.illumina.1.sorted.bam'", shell=True)
			subprocess.run(f"rm '/nobackup/rogers_research/Rhyker/rogers_temp/Rhyker/RetroGenes_dmel/BAMS/ncbi_dmel/ISO1/transcriptome/short_paired_end/broken/{transcript}/ncbi_dmel.ISO1.transcriptome.illumina.1.bam'", shell=True)
