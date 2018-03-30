PHATCAT: PHAge Therapy Candidate assesmenT
Version 0.1.0

Abstract
=========

With a sharp increase in multi-drug resistant bacterial infections, there has been a renewed interest in phage therapy.
Large and diverse phage libraries (collections of phages whose host range and genomes have been characterized) are
essential for effective use of phage therapy. In order to develop these libraries and gain Federal Drug Administration
(FDA) approval, precise and detailed genomic characterization is necessary. The foundational step in doing so is
producing a coding complete genome. Assembling phages to coding complete status requires several special considerations:
i) Phages have small genome sizes (ranging from ~2.5 kb to ~500 kb) and modern high-throughput sequencing (HTS)
technologies are capable of producing a vast amount of data. This can lead to over-sequencing which impedes the ability
of genome assembly software to produce full and accurate contigs. ii) Many phages have redundant genomic ends which can
lead to assembly artifacts in the form of artificial overlaps at the ends of phage contigs. iii) Phages have a variety
of packaging mechanisms which lead to different types of genomic termini. Different types of packaging mechanisms
require different methods to resolve the termini. iv) Finally, an accurate coding complete assembly for therapeutic
phages is critical because errors in assemblies can result in incorrect sequences that alter outcomes of genomic
screening for bacterial virulence factors and antibiotic resistance genes, both of which present severe consequences for
therapeutic applications in a patient. Here we present an automated pipeline for bringing phages to coding complete from
raw paired-end sequencing data. The pipeline deals with the considerations listed above by: i) performing both full and
subsampled assemblies, ii) automatically removing artificial repeats at the ends of phage contigs, iii) using PhageTerm
to resolve genomic termini, iv) and performing two assembly pipelines in parallel and making sure different tools
agree on a consensus phage contig.  This pipeline is implemented in Python, can be run locally or in a Sun Grid Engine
(SGE) cluster computing environment, and will be publicly available at https://github.com/BDRD-Genomics.

Overview
=========

Raw paired-end illumina data is used as input.

This pipeline runs on the Luigi framework. Please refer to Luigi documentation
(luigi.readthedocs.io/en/) for instructions on how to set up the luigi daemon and
monitoring jobs with the web interface.


To Install:
===========

For Mac and Linux only.

If running on a SGE cluster you make nodes available to phatcat by adding
them to the queue 'all.q'

CLC Assembly Cell (commercial software) is required to run this pipeline.
https://www.qiagenbioinformatics.com/products/clc-assembly-cell/
Add directory containing Assembly Cell tools to your PATH.
If running on a cluster, you need to add nodes which are licensed to run
Assembly Cell to the SGE queue 'clc_q'

Download the latest version of Anaconda for python 3 from their website.
https://www.anaconda.com/download/
Make sure to download to a location that all nodes in all.q can access.
Make sure the anaconda bin directory is in your PATH.

Edit install.sh so the condabin variable contains the path to the bin
directory of the anaconda version you just downloaded.
Run install.sh

Install PhageTerm
https://sourceforge.net/projects/phageterm/


To Configure
=============

Add path to this directory to PYTHONPATH. This can be done by adding this to your .bashrc/.bashprofile
> export PYTHONPATH="${PYTHONPATH}:<PATH TO PHATCAT>"

Copy luigi.cfg.tmpl to luigi.cfg
> cp luigi.cfg.tmpl luigi.cfg

Create an environment variable called 'LUIGI_CONFIG_PATH' and point this to the
luigi.cfg file you just created.

Edit luigi.cfg to desired settings. A description of each required setting is found below

GLOBALS - Parameters read into the pipeline which apply to more than one task.
------------------------------------------------------------------------------
NUM_THREADS: Number of threads tasks will run with. Set to the number of cores
available in your computing environment.

OUTPUT_DIR: Location where the pipeline will output files to. Make sure users
running pipeline have write permissions here.

PRIMARY_CONDA_ENV: The conda environment which contains the dependencies. Leave
as phatcat unless you change the env in the install script.

SGE: Can be True or False. This will tell phatcat if it should attempt to submit
tasks to an SGE cluster. If false, it will run locally.

NUM_NODES: This needs to be filled out if SGE is set to True. Set to an integer
that describes the number of nodes in the cluster you would like phatcat to run on.


PhageTerm - Parameters specific to PhageTerm task
-------------------------------------------------
phageterm_loc: Location of the phageterm script you downloaded early. Should be
full path, including filename

python2env: A python environment which runs python 2 and has dependencies for phageterm


To Run
======

Make sure you are using the Python interpreter downloaded with anaconda (may require
activating phatcat conda environment).

Start the luigi daemon (luigid) - For further information see luigi documentation
(luigi.readthedocs.io/en/)

Create a json file holding the parameters you would like to pass to the pipeline.
The json should be structured as an array of objects, where each object represents
a sample and key-value pairs in the object are parameters for that sample.
Each sample should have 4 parameters:
  - r1: full path to forward reads
  - r2: full path to reverse reads
  - phage_name: The name of the sample (what you want listed on output)
  - sample_size: The number of reads to subsample to (Should sample to 80-100x coverage over genome)
For an example of how this json file should look, see the example_parameter_list.json file.

Start the pipeline by running the phatcat.py script with path to the json file you created
as a parameter. Ex:
> python phatcat.py example_parameter_list.json





