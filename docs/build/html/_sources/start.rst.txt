Quick Start
===========

The first step is to install RetroDetector and to ensure that you have all the necessary dependencies. Please see the :doc:`installation` section for further information.

The next step is to make a working directory and to change into that directory. For the rest of these instructions, please make sure you are in your working directory.

Now you can set up the configuration file according to your reference and sequencing files. Please see the :doc:`config` section for further information. One option is to download the example config.yaml to your working directory from `here <https://raw.githubusercontent.com/tbenavi1/RetroDetector/main/example/config.yaml>`_ and make the necessary edits. Alternatively, you can make a new config.yaml in your working directory.

Then, you should copy the 

Now we can finally run the RetroDetector pipeline. First, activate the snakemake environment and perform a dry run:

.. code-block:: console

   conda activate snakemake
   snakemake -np 


