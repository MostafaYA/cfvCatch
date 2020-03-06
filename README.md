# General

What it is. What it does. How to cite


# Install CfwCatch

To install WGSBAC, you need to install Miniconda, Snakemake and then download the source code.
All ohter dependencies are downloaded durign run-time.


Miniconda:

    # Download conda installer
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

    # Set permissions to execute
    chmod +x Miniconda3-latest-Linux-x86_64.sh 	

    # Execute. Make sure to "yes" to add the conda to your PATH
    ./Miniconda3-latest-Linux-x86_64.sh 		

    # Add channels
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

Snakemake (Version 5.2.2 or newer):

    conda install -c bioconda -c conda-forge snakemake


 Install GIT    

    https://gist.github.com/derhuerst/1b15ff4652a867391f03


Download WGSBAC:

    git clone https://gitlab.com/FLI_Bioinfo/WGSBAC.git

Add WGSBAC to PATH:

    export PATH=$PATH:<path to WGSBAC folder>






# Running cfwCatch

Running only with obligatory parameters, WGSBAC will run FASTQ, calculate coverage, assembles genomes,run Quast on genomes and identify the genus and species using Kraken and MASH.


    wgsbac.pl  --(t)able <metatdata.csv>   --(k)raken <path to Kraken2 DB>  <Optional Parameters>
    
    
# Expert options



