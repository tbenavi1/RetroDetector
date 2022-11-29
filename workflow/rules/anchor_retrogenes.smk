# Analyze and sort the genome alignments to find clusters of genome locations with potential retrogenes


rule analyze_anchors:
    input:
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.transcriptome.long.subset.geneid.sorted.sam",
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.genome.long.subset.geneid.sorted.sam",
    output:
        "results/AS/{ref}/{sample}/{ref}.{sample}.genome.AS",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/analyze_anchors/{ref}.{sample}.log",
    script:
        "../scripts/analyze_anchors.py"


rule sort_anchors:
    input:
        "results/AS/{ref}/{sample}/{ref}.{sample}.genome.AS",
    output:
        "results/AS/{ref}/{sample}/{ref}.{sample}.genome.sorted.AS",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/sort_anchors/{ref}.{sample}.log",
    shell:
        "cut -f1-13 {input} | sort -k1,1 -k2,2 -k9,9 -k10,10n -k11,11n > {output}"


rule cluster_anchors:
    input:
        "results/AS/{ref}/{sample}/{ref}.{sample}.genome.sorted.AS",
    output:
        "results/AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.AS",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/cluster_anchors/{ref}.{sample}.log",
    script:
        "../scripts/cluster_anchors.py"


# Find clusters with the best read support, and categorize by whether it overlaps parent gene


rule sort_clustered_anchors:
    input:
        "results/AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.AS",
    output:
        "results/AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.sorted.AS",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/sort_clustered_anchors/{ref}.{sample}.log",
    shell:
        "sort -k1,1 -k4,4nr {input} > {output}"


rule best_AS:
    input:
        "results/AS/{ref}/{sample}/{ref}.{sample}.genome.clustered.sorted.AS",
    output:
        "results/AS/{ref}/{sample}/{ref}.{sample}.genome.best.AS",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/best_AS/{ref}.{sample}.log",
    script:
        "../scripts/best_AS.py"


rule categorize_AS:
    input:
        "results/AS/{ref}/{sample}/{ref}.{sample}.genome.best.AS",
        "results/Transcriptome/{ref}/{ref}.transcriptome.coords.tsv",
    output:
        "results/AS/{ref}/{sample}/{ref}.{sample}.genome.best.AS.diff",
        "results/AS/{ref}/{sample}/{ref}.{sample}.genome.best.AS.same",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/categorize_AS/{ref}.{sample}.log",
    script:
        "../scripts/categorize_AS.py"
