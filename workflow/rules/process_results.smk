rule get_retrogene_consensus_sequences:
    input:
        "results/AS/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.genome.best.AS.diff",
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.supportingalignments.genome.sorted.bam",
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.supportingalignments.genome.sorted.bam.bai",
    output:
        "results/retrogenes/{ref}/{sample}/long/consensus/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.consensus.txt",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/get_retrogene_consensus_sequences/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.log",
    script:
        "../scripts/get_retrogene_consensus_sequences.py"


rule get_retrogene_needles:
    input:
        "results/AS/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.genome.best.AS.diff",
        "results/Transcriptome/{ref}/{ref}.transcriptome.fna.gz",
        "results/retrogenes/{ref}/{sample}/long/consensus/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.consensus.txt",
    output:
        "results/retrogenes/{ref}/{sample}/long/consensus/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.needles.txt",
    conda:
        "../envs/emboss.yaml"
    log:
        "logs/get_retrogene_needles/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.log",
    script:
        "../scripts/get_retrogene_needles.py"


rule final_results:
    input:
        "results/Transcriptome/{ref}/{ref}.transcriptome.coords.tsv",
        "results/Transcriptome/{ref}/{ref}.transcriptome.fna.gz",
        "results/AS/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.genome.best.AS.diff",
        "results/retrogenes/{ref}/{sample}/long/consensus/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.needles.txt",
    output:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.final.results.tsv",
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.duplicate.results.tsv",
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.lowconfidence.results.tsv",
    params:
        junction_total_read_support_threshold=config[
            "junction_total_read_support_threshold"
        ],
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/final_results/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.log",
    script:
        "../scripts/final_results.py"


rule filter_mainchroms:
    input:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.final.results.tsv",
    output:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.final.mainchrom.results.tsv",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/filter_mainchroms/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.log",
    script:
        "../scripts/filter_mainchroms.py"


rule chrom_movement:
    input:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.final.mainchrom.results.tsv",
    output:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.chrommovement.results.tsv",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/chrom_movement/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.log",
    script:
        "../scripts/chrom_movement.py"


rule summarize_alignments_short:
    input:
        "results/Transcriptome/{ref}/{ref}.transcriptome.junctions.tsv",
        "results/Spanned/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.shortreads.spannedjunctions.tsv",
    output:
        "results/Summary/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.shortreads.coverage.across.junctions.tsv",
        "results/Summary/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.shortreads.numspannedjunctions.tsv",
        "results/Summary/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.shortreads.no_non_overlapping.coverage.across.junctions.tsv",
        "results/Summary/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.shortreads.no_non_overlapping.numspannedjunctions.tsv",
        "results/Summary/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.shortreads.no_alternate.coverage.across.junctions.tsv",
        "results/Summary/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.shortreads.no_alternate.numspannedjunctions.tsv",
    params:
        short_read_strong_spanning_alignment_expected_genomic_insertion_size_threshold=config[
            "short_read_strong_spanning_alignment_expected_genomic_insertion_size_threshold"
        ],
        junction_total_read_support_threshold=config[
            "junction_total_read_support_threshold"
        ],
        junction_strong_short_read_support_threshold=config[
            "junction_strong_short_read_support_threshold"
        ],
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/summarize_alignments_short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.log",
    script:
        "../scripts/summarize_spanned_junctions_short.py"


rule long_truepositives:
    input:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.final.mainchrom.results.tsv",
    output:
        "results/Summary/{ref}/{sample}/junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.long_truepositives.txt",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/long_truepositives/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.log",
    shell:
        "cut -f1 {input} > {output}"


rule long_truepositives2:
    input:
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.final.results.tsv",
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.duplicate.results.tsv",
        "results/retrogenes/{ref}/{sample}/long/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.retrogenes.lowconfidence.results.tsv",
    output:
        "results/Summary/{ref}/{sample}/junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.long_truepositives2.txt",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/long_truepositives2/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.log",
    shell:
        "cat {input} | grep -v 'putative retrogenes' | cut -f1 > {output}"


rule paper_results:
    input:
        "results/Summary/{ref}/{sample}/junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.long_truepositives.txt",
        "results/Summary/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.shortreads.numspannedjunctions.tsv",
        "results/Summary/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.shortreads.no_non_overlapping.numspannedjunctions.tsv",
        "results/Summary/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.shortreads.no_alternate.numspannedjunctions.tsv",
    output:
        "results/Summary/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.results.txt",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/paper_results/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.log",
    script:
        "../scripts/paper_results.py"


rule paper_results2:
    input:
        "results/Summary/{ref}/{sample}/junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.long_truepositives2.txt",
        "results/Summary/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.shortreads.numspannedjunctions.tsv",
        "results/Summary/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.shortreads.no_non_overlapping.numspannedjunctions.tsv",
        "results/Summary/{ref}/{sample}/short/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.shortreads.no_alternate.numspannedjunctions.tsv",
    output:
        "results/Summary/{ref}/{sample}/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.results2.txt",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/paper_results2/{ref}.{sample}.junctover{junction_overhang}.insertthresh{insertions_threshold}.totalsupport{junction_total_read_support_threshold}.{strongthreshold}.log",
    script:
        "../scripts/paper_results.py"
