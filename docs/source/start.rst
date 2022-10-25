Quick Start
===========

The first step is to install RetroDetector and to ensure that you have all the necessary dependencies. Please see the :doc:`installation` section for further information.

RetroDetector is built on the `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ workflow management system. The benefit of this is that Snakemake keeps track of all the steps in the pipeline and automatically takes care of resource management.

Any time you want to run RetroDetector to analyze the data for a new project, all you have to do is set up a new working directory and a new configuration file with the locations of all the data files.

So, make a working directory and change into that directory. For the rest of these instructions, please make sure you are in your working directory.

Now you can set up the configuration file according to your reference and sequencing files. Please see the :doc:`config` section for further information.

Once you have set up the configuration file, you can finally run the RetroDetector pipeline. First, activate the snakemake environment and perform a dry run:

.. code-block:: console

   conda activate snakemake
   snakemake -np 

This should print to the console something like this:

.. code-block:: console

  this is a test

If there are no errors, then you can run RetroDetector. To run the pipeline, you have to tell RetroDetector how many cores you would like to use. For example, if you want to use 32 cores you would run:

.. code-block:: console

  snakemake --cores 32

RetroDetector will then map the sequencing reads to the transcriptome and genome and then identify any retrogenes. For more information about the pipeline, please see the :doc:`pipeline` section. For more information about the final results produced by RetroDetector, please see the :doc:`results` section.

If you want to see how to run RetroDetector on an example dataset, please see the :doc:`example` section. For FAQs, please see the :doc:`faqs` section. For advanced users, if you want information about adjusting advanced parameters in the pipeline, please see the :doc:`parameters` section.

If you encounter any issues while running RetroDetector, please submit a GitHub issue `here <https://github.com/tbenavi1/RetroDetector/issues>`_.
