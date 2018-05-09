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
5. Comparing contigs
6. Resolving genomic termini
7. Generating a report

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

If running on a SGE cluster you make nodes available to exempliphi by adding
them to the queue 'all.q'

CLC Assembly Cell (commercial software) is required to run this pipeline.
https://www.qiagenbioinformatics.com/products/clc-assembly-cell/
Add directory containing Assembly Cell tools to your PATH.
If running on a cluster, you need to add nodes which are licensed to run
Assembly Cell to the SGE queue 'clc_q'

Dependencies are downloaded and installed as conda packages. Download the latest version of Anaconda for python 3 from their website.
https://www.anaconda.com/download/
Make sure the anaconda bin directory is in your PATH.

Edit install.sh so the condabin variable contains the path to the bin
directory of the anaconda version you just downloaded.
Run install.sh

Install PhageTerm
https://sourceforge.net/projects/phageterm/

Add path to this directory to PYTHONPATH. This can be done by adding this to your .bashrc/.bashprofile
```
export PYTHONPATH="${PYTHONPATH}:<PATH TO EXEMPLIPHI>"
```

### Configuring

This pipeline runs on the Luigi framework and uses Luigi's configuration system for defining certain parameters.
Please read Luigi's documentation on this topic for more detailed information 
(luigi.readthedocs.io/en/stable/configuration.html). 

A template configuration file is provided as luigi.cfg.tmpl. This file contains a base set of required configurations. 
However, any luigi parameter defined in the code can be defined with the configuration file. To make this file active, 
copy it to the same directory as luigi.cfg
```
cp luigi.cfg.tmpl luigi.cfg
```

To ensure that the configurations will be read regardless of name or location of the configuration file, you can set an 
environment variable pointing to the file.
```
export LUIGI_CONFIG_PATH=<Full path to configuration file>
```

Edit luigi.cfg to desired settings. A description of each required setting is found below

#### GLOBALS 

##### _Parameters read into the pipeline which apply to more than one task._

---

**NUM_THREADS -** _Number of threads tasks will run with. Set to the number of cores available in your computing environment._

**OUTPUT_DIR -** _Location where the pipeline will output files to. Make sure users running pipeline have write permissions here._

**PRIMARY_CONDA_ENV -** _The conda environment which contains the dependencies. Leave as exempliphi unless you install dependencies to a different environment._

**SGE -** _Can be True or False. This will tell exempliphi if it should attempt to submit tasks to an SGE cluster. If false, it will run locally._

**NUM_NODES -** _This needs to be filled out if SGE is set to True. Set to an integer that describes the number of nodes in the cluster you would like exempliphi to run on._

#### PhageTerm

##### _Parameters specific to PhageTerm task_

---

**phageterm_loc -** _Location of the phageterm script you downloaded early. Should be full path, including filename_

**python2env -** _A conda environment which runs python 2 and has dependencies for phageterm_

#### core

##### _Core Luigi specific parameters_

---

**log_level -** _Verbosity of output. Possible values: DEBUG, INFO, WARNING, ERROR, or CRITICAL._


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

matthew.r.lueder.ctr@mail.mil 

## Disclaimer

This work was funded by Naval Medical Research Center’s Advanced Medical Development Program, work unit number A1704. Casandra W. Philipson, Kimberly A. Bishop-Lilly are employees of the US government and LCDR Theron Hamilton is a military service member. This work was prepared as a part of official duties. Title 17 U.S.C. 105 provides that ‘Copyright protection under this title is not available for any work of the United States Government.’ Title 17 U.S.C. 101 defines a U.S. Government work as a work prepared by a military service member or employee of the U.S. Government as part of a person’s official duties. The views expressed in this article are those of the authors and do not necessarily reflect the official policy or position of the Department of the Navy, Defense Threat Reduction Agency, Department of Defense, nor the U.S. Government.


_Version 0.1.1 beta_