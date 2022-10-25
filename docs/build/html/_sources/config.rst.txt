Configuration
=============

Whenever you want to run RetroDetector to analyze the data for a new project, you need to set up a configuration file to tell RetroDetector where to find the reference genome, reference transcriptome, and whole genome sequencing files. 

Make sure you are in the working directory for your project and then copy the example config.yaml file from the RetroDetector installation location to your current directory. For example, if RetroDetector is installed at /user/software/RetroDetector, you would run:

.. code-block:: console

  cp /user/software/RetroDetector/example/config.yaml .

The first several lines of the example config.yaml file look like this:

.. code-block:: console

  #Please update this directory to the scripts directory where RetroDetector is installed
  scripts_directory:
  
  #Please update the filenames of the genome fasta file,
  #the genome gtf file, and the transcriptome file
  ref:
    "ncbi_dmel":
      genomic_fna: Downloads/ncbi_dmel.genomic.fna.gz
      genomic_gtf: Downloads/ncbi_dmel.genomic.gtf.gz
      rna_from_genomic: Downloads/ncbi_dmel.rna_from_genomic.fna.gz
  
  #Please update the filenames of any sequencing data you have.
  #Feel free to add as many samples as you want to analyze.
  fastqs:
    "ISO1":
      short_paired_end:
        1:
          1: Downloads/SRR10728584_1.fastq.gz
          2: Downloads/SRR10728584_2.fastq.gz
      pacbio_hifi:
        1: Downloads/SRR10238607.fastq.gz

The first step is to edit config.yaml to add the scripts directory location for your RetroDetector installation. For example, if RetroDetector is installed at /user/software/RetroDetector, you would edit the scripts directory line of config.yaml as follows:

.. code-block:: console

  scripts_directory: /user/software/RetroDetector/scripts

The next step is to edit the "ref" section of config.yaml with the locations for the gzipped reference genome, the gzipped reference gtf file, and the gzipped reference transcriptome. You can also rename the reference. If you change the reference name in the future because there are updated reference files or if you add a reference because you want to analyze your data against a closely related species, RetroDetector will automatically rerun the pipeline once the config file is updated.

Finally, edit the "fastqs" section of config.yaml with any samples you have and their corresponding sequencing data. Each sample may have "short_paired_end", "pacbio_hifi", "pacbio_clr", or "ont" sections for Illumina paired-end, PacBio HiFi, PacBio CLR, and Oxford Nanopore sequencing data respectively. For each sequencing technology, you may include as many files as you have by numbering them as above with 1, 2, etc. You may also include as many samples as you have.

To continue with setting up RetroDetector for your project, return to :doc:`start`.
