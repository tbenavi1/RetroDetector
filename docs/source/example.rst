Example Dataset
===============

Drosophila yakuba
-----------------

In this example, we run the entire RetroDetector pipeline on nine *Drosophila yakuba* samples from the Island of Mayotte. Each sample has paired-end Illumina and PacBio HiFi sequencing data. The first step is to download the necessary reference files:

.. code-block:: console
   
   mkdir Reference
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/746/365/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_genomic.fna.gz -O Reference/ncbi_dyak.genomic.fna.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/746/365/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_genomic.gtf.gz -O Reference/ncbi_dyak.genomic.gtf.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/746/365/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1/GCF_016746365.2_Prin_Dyak_Tai18E2_2.1_rna_from_genomic.fna.gz -O Reference/ncbi_dyak.rna_from_genomic.fna.gz

The next step is to download the short-read paired-end sequencing data. Make sure you have `sra-tools <https://github.com/ncbi/sra-tools>`_ installed. If you need more assistance please follow the instructions `here <https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump>`_.

.. code-block:: console

   mkdir FASTQS
   prefetch SRR001
   prefetch SRR002
   prefetch SRR003
   prefetch SRR004
   prefetch SRR005
   prefetch SRR006
   prefetch SRR007
   prefetch SRR008
   prefetch SRR009

   fasterq-dump SRR001 --outdir FASTQS
   fasterq-dump SRR002 --outdir FASTQS
   fasterq-dump SRR003 --outdir FASTQS
   fasterq-dump SRR004 --outdir FASTQS
   fasterq-dump SRR005 --outdir FASTQS
   fasterq-dump SRR006 --outdir FASTQS
   fasterq-dump SRR007 --outdir FASTQS
   fasterq-dump SRR008 --outdir FASTQS
   fasterq-dump SRR009 --outdir FASTQS
 
The next step is to download the long-read sequencing data:

.. code-block:: console

   prefetch SRR001
   prefetch SRR002
   prefetch SRR003
   prefetch SRR004
   prefetch SRR005
   prefetch SRR006
   prefetch SRR007
   prefetch SRR008
   prefetch SRR009

   fasterq-dump SRR001 --outdir FASTQS
   fasterq-dump SRR002 --outdir FASTQS
   fasterq-dump SRR003 --outdir FASTQS
   fasterq-dump SRR004 --outdir FASTQS
   fasterq-dump SRR005 --outdir FASTQS
   fasterq-dump SRR006 --outdir FASTQS
   fasterq-dump SRR007 --outdir FASTQS
   fasterq-dump SRR008 --outdir FASTQS
   fasterq-dump SRR009 --outdir FASTQS

The next step is to set up the configuration file config.yaml which tells RetroDetector where to find the necessary input files. For this example, you can download `config.yaml <https://raw.githubusercontent.com/tbenavi1/RetroDetector/main/example/config.yaml>`_ or copy the following code to a new file called config.yaml:

.. code-block:: console

   ref:
     ncbi_dyak:
       genomic_fna: Reference/ncbi_dyak.genomic.fna.gz
       genomic_gtf: Reference/ncbi_dyak.genomic.gtf.gz
       rna_from_genomic: Reference/ncbi_dyak.rna_from_genomic.fna.gz

   fastqs:
     short_paired_end:
       BE_1:
         1: [FASTQS/BE_1_R1.fastq.gz, FASTQS/BE_1_R2.fastq.gz]
       BE_13:
         1: [FASTQS/BE_13_R1.fastq.gz, FASTQS/BE_13_R2.fastq.gz]
       BE_18:
         1: [FASTQS/BE_18_R1.fastq.gz, FASTQS/BE_18_R2.fastq.gz]
       BE_8:
         1: [FASTQS/BE_8_R1.fastq.gz, FASTQS/BE_8_R2.fastq.gz]
       MTS_22:
         1: [FASTQS/MTS_22_R1.fastq.gz, FASTQS/MTS_22_R2.fastq.gz]
       MTS_25:
         1: [FASTQS/MTS_25_R1.fastq.gz, FASTQS/MTS_25_R2.fastq.gz]
       MTS_26:
         1: [FASTQS/MTS_26_R1.fastq.gz, FASTQS/MTS_26_R2.fastq.gz]
       SOU_15:
         1: [FASTQS/SOU_15_R1.fastq.gz, FASTQS/SOU_15_R2.fastq.gz]
       SOU_9:
         1: [FASTQS/SOU_9_R1.fastq.gz, FASTQS/SOU_9_R2.fastq.gz]
     long:
       BE_1:
         1: FASTQS/BE_1_long.fastq.gz
       BE_13:
         1: FASTQS/BE_13_long.fastq.gz
       BE_18:
         1: FASTQS/BE_18_long.fastq.gz
       BE_8:
         1: FASTQS/BE_8_long.fastq.gz
       MTS_22:
         1: FASTQS/MTS_22_long.fastq.gz
       MTS_25:
         1: FASTQS/MTS_25_long.fastq.gz
       MTS_26:
         1: FASTQS/MTS_26_long.fastq.gz
       SOU_15:
         1: FASTQS/SOU_15_long.fastq.gz
       SOU_9:
         1: FASTQS/SOU_9_long.fastq.gz
