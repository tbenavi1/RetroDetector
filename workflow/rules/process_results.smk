rule get_retrogene_consensus_sequences:
  input:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.best.AS.diff",
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supportingalignments.genome.sorted.bam",
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supportingalignments.genome.sorted.bam.bai"
  output:
    "results/retrogenes/{ref}/{sample}/long/consensus/{ref}.{sample}.retrogenes.consensus.txt"
  script:
    "../scripts/get_retrogene_consensus_sequences.py"

rule get_retrogene_needles:
  input:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.best.AS.diff",
    "results/Transcriptome/{ref}/{ref}.transcriptome.fna.gz",
    "results/retrogenes/{ref}/{sample}/long/consensus/{ref}.{sample}.retrogenes.consensus.txt"
  output:
    "results/retrogenes/{ref}/{sample}/long/consensus/{ref}.{sample}.retrogenes.needles.txt"
  script:
    "../scripts/get_retrogene_needles.py"

rule final_results:
  input:
    "results/Transcriptome/{ref}/{ref}.transcriptome.coords.tsv",
    "results/Transcriptome/{ref}/{ref}.transcriptome.fna.gz",
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.best.AS.diff",
    "results/retrogenes/{ref}/{sample}/long/consensus/{ref}.{sample}.retrogenes.needles.txt"
  output:
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.final.results.tsv"
  script:
    "../scripts/final_results.py"
