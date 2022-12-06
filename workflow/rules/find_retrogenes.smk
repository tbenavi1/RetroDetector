# Find long read transcriptome alignments which span exon-exon junctions but don't have intronic sequence


rule find_retrogene_alignments_long:
    input:
        junctions="results/Transcriptome/{ref}/{ref}.transcriptome.junctions.tsv",
        bam="results/BAMS/{ref}/{sample}/transcriptome/long/{ref}.{sample}.transcriptome.long.sorted.bam",
    output:
        spanningalignments="results/Spanned/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.longreads.spanningalignments.sam",
        spanningalignmentstxt="results/Spanned/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.longreads.spanningalignments.txt",
        spannedjunctions="results/Spanned/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.longreads.spannedjunctions.tsv",
        flankregions="results/Flanks/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.longreads.flankregions.tsv",
    params:
        junction_overhang=config["junction_overhang"],
        insertions_threshold=config["insertions_threshold"],
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/find_retrogene_alignments_long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.log",
    shell:
        "samtools view -h {input.bam} | python workflow/scripts/find_retrogene_alignments_long.py {params.junction_overhang} {params.insertions_threshold} {input.junctions} {output.spanningalignments} {output.spanningalignmentstxt} {output.spannedjunctions} {output.flankregions}"


# Only include transcriptome alignments that support a spanned exon-exon junction with sufficient read coverage


rule subset_transcriptome_bam:
    input:
        "results/Spanned/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.longreads.spannedjunctions.tsv",
        "results/Spanned/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.longreads.spanningalignments.sam",
        "results/Spanned/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.longreads.spanningalignments.txt",
    output:
        "results/Spanned/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.longreads.spannedjunctions.readnames.txt",
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.transcriptome.long.subset.sam",
    params:
        junction_total_read_support_threshold=config[
            "junction_total_read_support_threshold"
        ],
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/subset_transcriptome_bam/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.log",
    script:
        "../scripts/subset_transcriptome_bam.py"


# Only include genome alignments that have read names corresponding to the transcriptome alignments above


rule subset_genome_bam:
    input:
        reads="results/Spanned/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.longreads.spannedjunctions.readnames.txt",
        bam="results/BAMS/{ref}/{sample}/genome/long/{ref}.{sample}.genome.long.sorted.bam",
    output:
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.genome.long.subset.sam",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/subset_genome_bam/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.log",
    shell:
        "samtools view -h -N {input.reads} {input.bam} > {output}"


# Add parent GeneIDs to the alignments so that we can sort and process everything one parent gene at a time


rule add_geneid_to_subset_sams:
    input:
        "results/Transcriptome/{ref}/{ref}.transcriptome.coords.tsv",
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.transcriptome.long.subset.sam",
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.genome.long.subset.sam",
    output:
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.transcriptome.long.subset.geneid.sam",
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.genome.long.subset.geneid.sam",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/add_geneid_to_subset_sams/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.log",
    script:
        "../scripts/add_geneid_to_subset_sam.py"


# sort by parent GeneID


rule sort_genome_geneid_sam:
    input:
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.genome.long.subset.geneid.sam",
    output:
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.genome.long.subset.geneid.sorted.sam",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/sort_genome_geneid_sam/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.log",
    shell:
        "sort -k1,1 {input} > {output}"


rule sort_transcriptome_geneid_sam:
    input:
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.transcriptome.long.subset.geneid.sam",
    output:
        "results/Anchor/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.transcriptome.long.subset.geneid.sorted.sam",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/sort_transcriptome_geneid_sam/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.log",
    shell:
        "sort -k1,1 {input} > {output}"


# Find alignments of real short reads that overlap genes


rule find_short_alignments_overlapping_genes:
    input:
        "results/Transcriptome/{ref}/{ref}.transcriptome.coords.tsv",
        "results/BAMS/{ref}/{sample}/genome/short/{ref}.{sample}.genome.short.sorted.bam",
        "results/BAMS/{ref}/{sample}/genome/short/{ref}.{sample}.genome.short.sorted.bam.bai",
    output:
        "results/Spanned/{ref}/{sample}/short/{ref}.{sample}.shortreads.overlapping.genes.tsv",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/find_short_alignments_overlapping_genes/{ref}.{sample}.log",
    script:
        "../scripts/find_short_alignments_overlapping_genes.py"


# Find short read transcriptome alignments which span exon-exon junctions but don't have intronic sequence
rule find_shortreads_alignments:
    input:
        bam="results/BAMS/{ref}/{sample}/transcriptome/short/{ref}.{sample}.transcriptome.short.sorted.bam",
        junctions="results/Transcriptome/{ref}/{ref}.transcriptome.junctions.tsv",
        intronlengths="results/Transcriptome/{ref}/{ref}.transcriptome.intronlengths.tsv",
        overlapping="results/Spanned/{ref}/{sample}/short/{ref}.{sample}.shortreads.overlapping.genes.tsv",
    output:
        spanningalignments="results/Spanned/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.shortreads.spanningalignments.sam",
        spannedjunctions="results/Spanned/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.shortreads.spannedjunctions.tsv",
    params:
        junction_overhang=config["junction_overhang"],
        insertions_threshold=config["insertions_threshold"],
        short_read_spanning_alignment_minimum_matching_bp=config[
            "short_read_spanning_alignment_minimum_matching_bp"
        ],
        short_read_spanning_alignment_minimum_transcriptomic_insertion_size=config[
            "short_read_spanning_alignment_minimum_transcriptomic_insertion_size"
        ],
        short_read_spanning_alignment_maximum_transcriptomic_insertion_size=config[
            "short_read_spanning_alignment_maximum_transcriptomic_insertion_size"
        ],
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/find_shortreads_alignments/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.log",
    shell:
        "samtools view -h {input.bam} | python workflow/scripts/find_retrogene_alignments_short.py {params.junction_overhang} {params.insertions_threshold} {params.short_read_spanning_alignment_minimum_matching_bp} {params.short_read_spanning_alignment_minimum_transcriptomic_insertion_size} {params.short_read_spanning_alignment_maximum_transcriptomic_insertion_size} {input.junctions} {input.intronlengths} {input.overlapping} {output.spanningalignments} {output.spannedjunctions}"
