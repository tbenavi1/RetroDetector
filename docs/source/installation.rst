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

If you don't have this setup, you can download and install the file that corresponds to your OS and architecture, according to the table `here <https://github.com/conda-forge/miniforge#mambaforge>`_.

You will then need to close and restart your terminal. If you want to prevent Conda from activating the base environment by default, you can run the following:

.. code-block:: console

   conda config --set auto_activate_base false

Now, you can install Snakemake:

.. code-block:: console

   conda activate base
   mamba create -c conda-forge -c bioconda -n snakemake snakemake


Now, you can activate the snakemake environment using:

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

Then, make sure to add the paths to BWA and minimap2 to your $PATH environment variable. For example, if your software directory is /user/software, you could run:

.. code-block:: console

   export PATH=/user/software/minimap2:$PATH
   export PATH=/user/software/bwa:$PATH

Next, you will need to install SAMtools. Please follow the instructions `here <http://www.htslib.org/download/>`_ to download and install the latest version of SAMtools and HTSlib. HTSlib is required in order to use bgzip. For example, if your software directory is /user/software, you could run the following (for SAMtools 1.15.1):

.. code-block:: console

   wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2
   tar -xf samtools-1.15.1.tar.bz2
   cd samtools-1.15.1
   ./configure --prefix=/user/software/samtools-1.15.1
   make
   make install
   cd ..
   wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2
   tar -xf htslib-1.15.1.tar.bz2
   cd htslib-1.15.1
   ./configure --prefix=/user/software/htslib-1.15.1
   make
   make install
   cd ..

Finally, make sure to add the paths to SAMtools and HTSlib to your $PATH environment variable. For example: 

.. code-block:: console

   export PATH=/user/software/samtools-1.15.1/bin:$PATH
   export PATH=/user/software/htslib-1.15.1/bin:$PATH
