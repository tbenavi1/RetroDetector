Example Dataset
===============

Drosophila melanogaster
-----------------------

In this example, we run the entire RetroDetector pipeline on a *Drosophila melanogaster* sample with Illumina paired-end and PacBio HiFi sequencing data. To begin, make a working directory and change into that directory. Please make sure to run all of the following code while you are in your working directory. The first step is to download the necessary reference files from NCBI:

.. code-block:: console
   
   mkdir Reference
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz -O Reference/ncbi_dmel.genomic.fna.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gtf.gz -O Reference/ncbi_dmel.genomic.gtf.gz
   wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_rna_from_genomic.fna.gz -O Reference/ncbi_dmel.rna_from_genomic.fna.gz

The next step is to download and gzip the sequencing data. Make sure you have `sra-tools <https://github.com/ncbi/sra-tools>`_ installed. If you need more assistance please follow the instructions for `downloading <https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit>`_, `installing <https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit>`_, and using `fasterq-dump <https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump>`_.

.. code-block:: console

   mkdir FASTQS
   prefetch SRR10238607 --max-size 21G
   fasterq-dump SRR10238607 --outdir FASTQS
   gzip FASTQS/SRR10238607.fastq
   prefetch SRR10728584
   fasterq-dump SRR10728584 --outdir FASTQS
   gzip FASTQS/SRR10728584_1.fastq
   gzip FASTQS/SRR10728584_2.fastq
 
The next step is to set up the configuration file config.yaml which tells RetroDetector where to find the necessary input files. For this example, you can download `config.yaml <https://raw.githubusercontent.com/tbenavi1/RetroDetector/main/example/config.yaml>`_ to your working directory, or you can copy the following code to a new file in your working directory called config.yaml:

.. code-block:: console

   ref:
     ncbi_dmel:
       genomic_fna: Reference/ncbi_dmel.genomic.fna.gz
       genomic_gtf: Reference/ncbi_dmel.genomic.gtf.gz
       rna_from_genomic: Reference/ncbi_dmel.rna_from_genomic.fna.gz

   fastqs:
     short_paired_end:
       ISO1:
         1: [FASTQS/SRR10728584_1.fastq.gz, FASTQS/SRR10728584_2.fastq.gz]
     long:
       ISO1:
         1: FASTQS/SRR10238607.fastq.gz
