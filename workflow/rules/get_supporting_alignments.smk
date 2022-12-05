rule final_genome_alignments:
    input:
        "results/AS/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.genome.best.AS.diff",
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.genome.long.subset.sam",
        "results/AS/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.genome.AS",
    output:
        temp(
            "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.retrogenes.supportingalignments.genome.sam"
        ),
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/final_genome_alignments/{ref}.{sample}.log",
    script:
        "../scripts/final_genome_alignments.py"


rule bam_final_genome_alignments:
    input:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.retrogenes.supportingalignments.genome.sam",
    output:
        temp(
            "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.retrogenes.supportingalignments.genome.bam"
        ),
    threads: 32
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/bam_final_genome_alignments/{ref}.{sample}.log",
    shell:
        "samtools view -@ {threads} -b -o {output} {input}"


rule sort_final_genome_alignments:
    input:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.retrogenes.supportingalignments.genome.bam",
    output:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.retrogenes.supportingalignments.genome.sorted.bam",
    threads: 32
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/sort_final_genome_alignments/{ref}.{sample}.log",
    shell:
        "samtools sort -@ {threads} -o {output} {input}"


rule index_final_genome_alignments:
    input:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.retrogenes.supportingalignments.genome.sorted.bam",
    output:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.retrogenes.supportingalignments.genome.sorted.bam.bai",
    threads: 32
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/index_final_genome_alignments/{ref}.{sample}.log",
    shell:
        "samtools index -@ {threads} {input}"
