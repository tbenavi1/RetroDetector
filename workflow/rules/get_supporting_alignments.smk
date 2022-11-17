rule final_genome_alignments:
  input:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.best.AS.diff",
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.sam",
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.AS"
  output:
    temp("results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supportingalignments.genome.sam")
  script:
    "../scripts/final_genome_alignments.py"
    
rule bam_final_genome_alignments:
  input:
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supportingalignments.genome.sam"
  output:
    temp("results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supportingalignments.genome.bam")
  threads: 32
  shell:
    "samtools view -@ {threads} -b -o {output} {input}"

rule sort_final_genome_alignments:
  input:
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supportingalignments.genome.bam"
  output:
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supportingalignments.genome.sorted.bam"
  threads: 32
  shell:
    "samtools sort -@ {threads} -o {output} {input}"

rule index_final_genome_alignments:
  input:
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supportingalignments.genome.sorted.bam"
  output:
    "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.retrogenes.supportingalignments.genome.sorted.bam.bai"
  threads: 32
  shell:
    "samtools index -@ {threads} {input}"
