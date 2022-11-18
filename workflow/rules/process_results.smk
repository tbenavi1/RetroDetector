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
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.final.results.tsv",
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.duplicate.results.tsv",
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.lowconfidence.results.tsv"
  script:
    "../scripts/final_results.py"

rule filter_mainchroms:
  input:
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.final.results.tsv"
  output:
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.final.mainchrom.results.tsv"
  script:
    "../scripts/filter_mainchroms.py"

rule chrom_movement:
  input:
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.final.mainchrom.results.tsv"
  output:
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.chrommovement.results.tsv"
  script:
    "../scripts/chrom_movement.py"
