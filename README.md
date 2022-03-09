# ExempliPhi 

ExempliPhi is an automated genome assembly pipeline for phages. It takes raw paired-end illumina reads as input and attempts to give a 
confident, correctly arranged genome as output. It does this in a modular way by executing a series of tasks. 

This series of tasks includes:
1. Read quality control trimming and filtering 
   - EDGE QC (FaQCs)
   - CLC QC
2. Read subsampling
3. Assembly
   - SPAdes
   - CLC Assembly Cell
4. Assembly artifact removal
5. Taxonomic Profiling
6. Comparing contigs
7. Comparing phage contig to terminase and integrase database
8. Resolving genomic termini
9. Generating a report

---

**_Abstract_**

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

---


## Getting Started

These instructions will help you get ExempliPhi running on your environment. 


### Installing

For Mac and Linux only.

CLC Assembly Cell (commercial software) is required to run this pipeline.
https://www.qiagenbioinformatics.com/products/clc-assembly-cell/
Add directory containing Assembly Cell tools to your PATH.
If running on a cluster, you add nodes which are licensed to run Assembly Cell 
to a queue (SGE) or partition (SLURM). The name of this partition will be 
entered in the configuration file. 

Most dependencies are downloaded and installed through the conda package management system. 
Download the latest version of Anaconda/miniconda for python 3 from their website.
https://www.anaconda.com/download/
Make sure the anaconda bin directory is in your PATH and that the CONDA_EXE environment variable is set.

After installing conda, you can acquire dependencies by running **install.sh.** This will create two conda environments:
exempliphi and py2. The 'exempliphi' conda environment will be used to run exempliphi and most tasks within the pipeline.
The 'py2' conda environment will be used to run PhageTerm.

PhageTerm is a software tool used to fix and reorient the termini of phage genomes. Install Phage Term from:
https://sourceforge.net/projects/phageterm/
Make sure the location where you install it is accessible from all cluster nodes if you are using a cluster. The location
of PhageTerm must be given as a configuration later.

ExempliPhi performs taxonomic profiling on contigs and unmapped reads (See latest manuscript for details on how this 
works). This functionality depends on NCBI's nt (non-redundant nucleotide) BLAST database. The easiest way to download nt
is with the update_blastdb.pl tool (this tool is available from within the exempliphi conda environment). Find details on
how to run it *[here](https://www.ncbi.nlm.nih.gov/books/NBK537770/)*. The path to the location of this database will be
given as a configuration later. In addition, the taxonomic profiling module requires a sqlite database which is created
by update_taxonomy_database.sh. Prior to running this script, the 'taxdb_loc' configuration must be set so 
update_taxonomy_database.sh knows where to save the database to (see following section for configuration instructions).
Taxonomic profiling also uses the tool qualex, which is in the 3rd_party directory. Qualex must be compiled by going into
3rd_party/qualex-ms-1.2 and running make. In addition, you need to run make in 3rd_party/qualex-ms-1.2/converter to get the 
tool asc2bin. You may need to install blas, lapack, and fortran libraries for this to work.
If on ubuntu, you can do this with the following command.
`apt-get install -y libblas-dev liblapack-dev gfortran`
For taxonomic profiling to work on headless servers, it must be run with xvfb. If this is the case install xvfb.
`apt-get -y install xvfb`

There are various other dynamically linked libraries that tools which are used by ExempliPhi mat try to link. You can ensure you have all of them with the following command:
`apt-get update && apt-get install -y curl ca-certificates libapr1-dev libaprutil1-dev build-essential iputils-ping openssh-client ffmpeg libsm6 libxext6 libncurses5 `

### Configuring

This pipeline runs on the Luigi framework and uses Luigi's configuration system for defining certain parameters.
Please read Luigi's documentation on this topic for more detailed information 
(luigi.readthedocs.io/en/stable/configuration.html). 

A template configuration file is provided as luigi.cfg.tmpl. 
To make this file active, copy it to the same directory as luigi.cfg
```
cp luigi.cfg.tmpl luigi.cfg
```

Now edit luigi.cfg to desired settings. Each setting is described in the template. Settings which are commented out are
optional. Settings which are not commented out must be set for ExempliPhi to run properly.


### Installation and configuration checklist

1. Install CLC Assembly Cell 
2. Install Conda
3. Run install.sh
4. Install PhageTerm
5. Download NCBI's nt BLAST database
6. Compile qualex and asc2bin
7. Install xvfb if needed
8. Install additional libraries
9. Copy luigi.cfg.tmpl to luigi.cfg
10. Edit luigi.cfg to desired configurations
11. Run update_taxonomy_database.sh 


### Running

Make sure you are using the Python interpreter downloaded with anaconda (may require
activating exempliphi conda environment).

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

Start the pipeline by running the exempliphi.py script with path to the json file you created
as a parameter. Ex:
```
python exempliphi.py example_parameter_list.json
```

## Citation

Philipson, C.W.; Voegtly, L.J.; Lueder, M.R.; Long, K.A.; Rice, G.K.; Frey, K.G.; Biswas, B.; Cer, R.Z.; Hamilton, T.; Bishop-Lilly, K.A.	Characterizing Phage Genomes for Therapeutic Applications. Viruses 2018, 10, 188.

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details

## Contact

logan.j.voegtly.ctr@mail.mil

## Disclaimer

This work was funded by Naval Medical Research Center’s Advanced Medical Development Program, work unit number A1704. Casandra W. Philipson, Kimberly A. Bishop-Lilly are employees of the US government and LCDR Theron Hamilton is a military service member. This work was prepared as a part of official duties. Title 17 U.S.C. 105 provides that ‘Copyright protection under this title is not available for any work of the United States Government.’ Title 17 U.S.C. 101 defines a U.S. Government work as a work prepared by a military service member or employee of the U.S. Government as part of a person’s official duties. The views expressed in this article are those of the authors and do not necessarily reflect the official policy or position of the Department of the Navy, Defense Threat Reduction Agency, Department of Defense, nor the U.S. Government.


Version 2