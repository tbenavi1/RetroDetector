# Index transcriptome according to short-read/long-read sequencing technology


rule bwa_index_transcriptome:
    input:
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz",
    output:
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz.amb",
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz.ann",
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz.bwt",
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz.pac",
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz.sa",
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_index_transcriptome/{ref}.log",
    shell:
        "bwa index {input}"


# we use more pemissive map-ont when mapping to transcriptome for all long-read technologies
rule minimap_index_transcriptome:
    input:
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz",
    output:
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.mmi",
    threads: 3
    conda:
        "../envs/minimap.yaml"
    log:
        "logs/minimap_index_transcriptome/{ref}.log",
    shell:
        "minimap2 -x map-ont -d {output} {input}"


# Index genome according to short-read/long-read sequencing technology


rule bwa_index_genome:
    input:
        "results/Genome/{ref}/{ref}.genome.fna.gz",
    output:
        "results/Genome/{ref}/{ref}.genome.fna.gz.amb",
        "results/Genome/{ref}/{ref}.genome.fna.gz.ann",
        "results/Genome/{ref}/{ref}.genome.fna.gz.bwt",
        "results/Genome/{ref}/{ref}.genome.fna.gz.pac",
        "results/Genome/{ref}/{ref}.genome.fna.gz.sa",
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_index_genome/{ref}.log",
    shell:
        "bwa index {input}"


rule minimap_index_genome:
    input:
        "results/Genome/{ref}/{ref}.genome.fna.gz",
    output:
        "results/Genome/{ref}/{ref}.genome.{tech}.mmi",
    threads: 3
    params:
        preset=get_minimap_preset,
    conda:
        "../envs/minimap.yaml"
    log:
        "logs/minimap_index_genome/{ref}.{tech}.log",
    shell:
        "minimap2 -x {params.preset} -d {output} {input}"


# Map real short and long reads to transcriptome


rule bwa_mem_transcriptome_short:
    input:
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz.amb",
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz.ann",
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz.bwt",
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz.pac",
        "results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz.sa",
        ref="results/Transcriptome/{ref}/{ref}.multiexontranscriptome.fna.gz",
        read1=lambda wildcards: config["fastqs"][wildcards.sample]["short"][
            int(wildcards.i)
        ][1],
        read2=lambda wildcards: config["fastqs"][wildcards.sample]["short"][
            int(wildcards.i)
        ][2],
    output:
        temp(
            "results/BAMS/{ref}/{sample}/transcriptome/short/{ref}.{sample}.transcriptome.illumina.{i,[0-9]+}.sam"
        ),
    threads: 32
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_mem_transcriptome_short/{ref}.{sample}.{i}.log",
    shell:
        "bwa mem -B 3 -D 0.0 -X 10.0 -h 1000 -t {threads} {input.ref} {input.read1} {input.read2} > {output}"


rule minimap_transcriptome_long:
    input:
        mmi="results/Transcriptome/{ref}/{ref}.multiexontranscriptome.mmi",
        reads=lambda wildcards: config["fastqs"][wildcards.sample][wildcards.tech][
            int(wildcards.i)
        ],
    output:
        temp(
            "results/BAMS/{ref}/{sample}/transcriptome/{tech}/{ref}.{sample}.transcriptome.{tech}.{i,[0-9]+}.sam"
        ),
    threads: 31
    conda:
        "../envs/minimap.yaml"
    log:
        "logs/minimap_transcriptome_long/{ref}.{sample}.{tech}.{i}.log",
    shell:
        "minimap2 -x map-ont --sam-hit-only -P -t {threads} -a {input.mmi} {input.reads} > {output}"


# Map real short and long reads to genome


rule bwa_mem_genome_short:
    input:
        "results/Genome/{ref}/{ref}.genome.fna.gz.amb",
        "results/Genome/{ref}/{ref}.genome.fna.gz.ann",
        "results/Genome/{ref}/{ref}.genome.fna.gz.bwt",
        "results/Genome/{ref}/{ref}.genome.fna.gz.pac",
        "results/Genome/{ref}/{ref}.genome.fna.gz.sa",
        ref="results/Genome/{ref}/{ref}.genome.fna.gz",
        read1=lambda wildcards: config["fastqs"][wildcards.sample]["short"][
            int(wildcards.i)
        ][1],
        read2=lambda wildcards: config["fastqs"][wildcards.sample]["short"][
            int(wildcards.i)
        ][2],
    output:
        temp(
            "results/BAMS/{ref}/{sample}/genome/short/{ref}.{sample}.genome.illumina.{i,[0-9]+}.sam"
        ),
    threads: 32
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/bwa_mem_genome_short/{ref}.{sample}.{i}.log",
    shell:
        "bwa mem -t {threads} {input.ref} {input.read1} {input.read2} > {output}"


rule minimap_genome_long:
    input:
        mmi="results/Genome/{ref}/{ref}.genome.{tech}.mmi",
        reads=lambda wildcards: config["fastqs"][wildcards.sample][wildcards.tech][
            int(wildcards.i)
        ],
    output:
        temp(
            "results/BAMS/{ref}/{sample}/genome/{tech}/{ref}.{sample}.genome.{tech}.{i,[0-9]+}.sam"
        ),
    threads: 31
    params:
        preset=get_minimap_preset,
    conda:
        "../envs/minimap.yaml"
    log:
        "logs/minimap_genome_long/{ref}.{sample}.{tech}.{i}.log",
    shell:
        "minimap2 -x {params.preset} -t {threads} -a {input.mmi} {input.reads} > {output}"
