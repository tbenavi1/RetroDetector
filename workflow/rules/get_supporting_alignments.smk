rule final_genome_alignments:
    input:
        "results/AS/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.genome.best.AS.diff",
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.genome.long.subset.sam",
        "results/AS/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.genome.AS",
    output:
        temp(
            "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.retrogenes.supportingalignments.genome.sam"
        ),
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/final_genome_alignments/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.log",
    script:
        "../scripts/final_genome_alignments.py"


rule bam_final_genome_alignments:
    input:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.retrogenes.supportingalignments.genome.sam",
    output:
        temp(
            "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.retrogenes.supportingalignments.genome.bam"
        ),
    threads: 32
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/bam_final_genome_alignments/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.log",
    shell:
        "samtools view -@ {threads} -b -o {output} {input}"


rule sort_final_genome_alignments:
    input:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.retrogenes.supportingalignments.genome.bam",
    output:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.retrogenes.supportingalignments.genome.sorted.bam",
    threads: 32
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/sort_final_genome_alignments/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.log",
    shell:
        "samtools sort -@ {threads} -o {output} {input}"


rule index_final_genome_alignments:
    input:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.retrogenes.supportingalignments.genome.sorted.bam",
    output:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.retrogenes.supportingalignments.genome.sorted.bam.bai",
    threads: 32
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/index_final_genome_alignments/{ref}.{sample}.junctover{junction_overhang}.insertthresh{max_insertions}.totalsupport{read_support}.log",
    shell:
        "samtools index -@ {threads} {input}"
