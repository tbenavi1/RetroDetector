FAQs
====

What do I do if I have more than one set of short-read FASTQ files or more than one long-read FASTQ file for a given sample?
----------------------------------------------------------------------------------------------------------------------------

All you need to do is to adjust the config.yaml file accordingly. For example, if you have two sets of short-read FASTQ files and three long-read FASTQ files for sample1, then the fastqs section of config.yaml would look something like this:

.. code-block:: console

   fastqs:
     short_paired_end:
       sample1:
         1: [sample1_short1_R1.fastq, sample1_short1_R2.fastq]
         2: [sample1_short2_R1.fastq, sample1_short2_R2.fastq]
     long:
       sample1:
         1: sample1_long1.fastq
         2: sample1_long2.fastq
         3: sample1_long3.fastq
