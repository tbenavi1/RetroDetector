# bgzip reference genome


rule unzip_genome:
    input:
        lambda wildcards: config["ref"][wildcards.ref]["genome"],
    output:
        "results/Genome/{ref}/{ref}.genome.fna",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/unzip_genome/{ref}.log",
    shell:
        "gunzip -c {input} > {output}"


rule bgzip_genome:
    input:
        "results/Genome/{ref}/{ref}.genome.fna",
    output:
        "results/Genome/{ref}/{ref}.genome.fna.gz",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/bgzip_genome/{ref}.log",
    shell:
        "bgzip {input}"


# bgzip reference transcriptome


rule unzip_transcriptome:
    input:
        lambda wildcards: config["ref"][wildcards.ref]["transcriptome"],
    output:
        "results/Transcriptome/{ref}/{ref}.transcriptome.fna",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/unzip_transcriptome/{ref}.log",
    shell:
        "gunzip -c {input} > {output}"


rule bgzip_transcriptome:
    input:
        "results/Transcriptome/{ref}/{ref}.transcriptome.fna",
    output:
        "results/Transcriptome/{ref}/{ref}.transcriptome.fna.gz",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/bgzip_transcriptome/{ref}.log",
    shell:
        "bgzip {input}"


# samtools index reference genome


rule samtools_index_genome:
    input:
        "results/Genome/{ref}/{ref}.genome.fna.gz",
    output:
        "results/Genome/{ref}/{ref}.genome.fna.gz.fai",
        "results/Genome/{ref}/{ref}.genome.fna.gz.gzi",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools_index_genome/{ref}.log",
    shell:
        "samtools faidx {input}"


# process transcriptome to get multi-exon transcriptome, coords, junctions, introns, and intron lengths


rule process_transcriptome:
    input:
        "results/Transcriptome/{ref}/{ref}.transcriptome.fna.gz",
    output:
        temp("results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna"),
        "results/Transcriptome/{ref}/{ref}.transcriptome.coords.tsv",
        "results/Transcriptome/{ref}/{ref}.transcriptome.junctions.tsv",
        "results/Transcriptome/{ref}/{ref}.transcriptome.introns.tsv",
        "results/Transcriptome/{ref}/{ref}.transcriptome.intronlengths.tsv",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/process_transcriptome/{ref}.log",
    script:
        "../scripts/process_transcriptome.py"


rule bgzip_multiexontranscriptome:
    input:
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna",
    output:
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/bgzip_multiexontranscriptome/{ref}.log",
    shell:
        "bgzip {input}"
