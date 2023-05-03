rule short_supporting_transcriptome:
    input:
        bam="results/BAMS/{ref}/{sample}/transcriptome/short/{ref}.{sample}.transcriptome.short.sorted.bam",
        bai="results/BAMS/{ref}/{sample}/transcriptome/short/{ref}.{sample}.transcriptome.short.sorted.bam.bai",
    output:
        bam="results/final_outputs/{ref}/{sample}/{sample}.shortread.supportingalignments.transcriptome.bam",
        bai="results/final_outputs/{ref}/{sample}/{sample}.shortread.supportingalignments.transcriptome.bam.bai",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/short_supporting_transcriptome/{ref}.{sample}.log",
    shell:
        "cp {input.bam} {output.bam}; "
        "cp {input.bai} {output.bai}"


rule long_supporting_genome:
    input:
        bam=f"results/retrogenes/{{ref}}/{{sample}}/long/{{ref}}.{{sample}}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.retrogenes.supportingalignments.genome.sorted.bam",
        bai=f"results/retrogenes/{{ref}}/{{sample}}/long/{{ref}}.{{sample}}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.retrogenes.supportingalignments.genome.sorted.bam.bai",
    output:
        bam="results/final_outputs/{ref}/{sample}/{sample}.longread.supportingalignments.genome.bam",
        bai="results/final_outputs/{ref}/{sample}/{sample}.longread.supportingalignments.genome.bam.bai",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/long_supporting_genome/{ref}.{sample}.log",
    shell:
        "mv {input.bam} {output.bam}; "
        "mv {input.bai} {output.bai}"


rule long_supporting_transcriptome:
    input:
        bam=f"results/Anchor/{{ref}}/{{sample}}/{{ref}}.{{sample}}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.transcriptome.long.subset.sorted.bam",
        bai=f"results/Anchor/{{ref}}/{{sample}}/{{ref}}.{{sample}}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.transcriptome.long.subset.sorted.bam.bai",
    output:
        bam="results/final_outputs/{ref}/{sample}/{sample}.longread.supportingalignments.transcriptome.bam",
        bai="results/final_outputs/{ref}/{sample}/{sample}.longread.supportingalignments.transcriptome.bam.bai",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/long_supporting_transcriptome/{ref}.{sample}.log",
    shell:
        "mv {input.bam} {output.bam}; "
        "mv {input.bai} {output.bai}"


rule retrogene_consensus_sequences:
    input:
        f"results/retrogenes/{{ref}}/{{sample}}/long/{{ref}}.{{sample}}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.retrogenes.final.mainchrom.results.tsv",
    output:
        "results/final_outputs/{ref}/{sample}/{sample}.retrogene_consensus_sequences.fa",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/retrogene_consensus_sequences/{ref}.{sample}.log",
    script:
        "../scripts/retrogene_consensus_sequences.py"


# rule short_retrogene_final_output:
#  input:
#
#  output:
#    "results/final_outputs/{ref}/{sample}/{sample}.shortread.retrogene_final_output.tsv"
#  conda:
#    "../envs/samtools.yaml"
#  log:
#    "logs/short_retrogene_final_output/{ref}.{sample}.log"
#  script:
#    "../scripts/short_retrogene_final_output"
# rule long_retrogene_final_output:
#  input:
#
#  output:
#    "results/final_outputs/{ref}/{sample}/{sample}.longread.retrogene_final_output.tsv"
