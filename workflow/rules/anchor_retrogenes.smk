#Analyze and sort the genome alignments to find clusters of genome locations with potential retrogenes

rule analyze_anchors:
  input:
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.geneid.sorted.sam",
    "results/Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sorted.sam",
  output:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.AS"
  script:
    "../scripts/analyze_anchors.py"

rule sort_anchors:
  input:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.AS"
  output:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.sorted.AS"
  shell:
    "cut -f1-13 {input} | sort -k1,1 -k2,2 -k9,9 -k10,10n -k11,11n > {output}"

rule cluster_anchors:
  input:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.sorted.AS"
  output:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.AS"
  script:
    "../scripts/cluster_anchors.py"

#Find clusters with the best read support, and categorize by whether it overlaps parent gene

rule sort_clustered_anchors:
  input:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.AS"
  output:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.sorted.AS"
  shell:
    "sort -k1,1 -k4,4nr {input} > {output}"

rule best_AS:
  input:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.sorted.AS"
  output:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.best.AS"
  script:
    "../scripts/best_AS.py"

rule categorize_AS:
  input:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.best.AS",
    "results/Transcriptome/{ref}/{ref}.transcriptome.coords.tsv"
  output:
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.best.AS.diff",
    "results/AS/{ref}/{sample}/{ref}.{sample}.genome.best.AS.same"
  script:
    "../scripts/categorize_AS.py"
