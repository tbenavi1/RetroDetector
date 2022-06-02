Installation and Dependencies
=============================

Installation
------------

To install RetroDetector, clone the GitHub repository:

.. code-block:: console

   git clone https://github.com/tbenavi1/RetroDetector.git

Dependencies
------------

RetroDetector is built on the `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ workflow management system. For full details on how to install Snakemake, please visit `here <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_. To summarize, the first step is to install mambaforge. If you have an x86_64 Linux machine then you can run the following commands:

.. code-block:: console

   wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
   bash Mambaforge-Linux-x86_64.sh

You will then need to close and restart your terminal. If you want to prevent Conda from activating the base environment by default, you can run the following:

.. code-block:: console

   conda config --set auto_activate_base false

Now, you can install Snakemake:

.. code-block:: console

   conda activate base
   mamba create -c conda-forge -c bioconda -n snakemake snakemake


Finally, you can activate the snakemake environment using:

.. code-block:: console

   conda activate snakemake

Next, you will need to install `BWA <https://github.com/lh3/bwa>`_ and `minimap2 <https://github.com/lh3/minimap2>`_ in order to align the sequencing reads to the reference files. Full instructions can be found in the corresponding GitHub pages. To summarize, run the following commands in the folder where you normally install software:

.. code-block:: console

   git clone https://github.com/lh3/bwa.git
   cd bwa
   make
   cd ..
   git clone https://github.com/lh3/minimap2
   cd minimap2
   make
   cd ..

Then, make sure to add the paths to BWA and minimap2 to your $PATH environmental variable in order to make sure that 

