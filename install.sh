# EDIT FOLLOWING LINE BEFORE RUNNING
condabin="<path to anaconda bin directory here>"
${condabin}/conda create -y -n phatcat
source activate phatcat
${condabin}/conda config --add channels conda-forge
${condabin}/conda config --add channels anaconda
${condabin}/conda config --add channels bioconda
${condabin}/conda config --add channels r
${condabin}/conda config --add channels auto
${condabin}/conda update -y conda
${condabin}/conda install -y -c anaconda luigi pigz lxml networkx
${condabin}/conda install -y -c conda-forge bzip2
${condabin}/conda install -y -c bioconda biopython spades seqtk blast samtools perl-string-approx perl-parallel-forkmanager
${condabin}/conda install -y -c r r-base
pip install docx-mailmerge
source deactivate
${condabin}/conda create -y --name py27 python=2.7
source activate py27
${condabin}/conda install -y matplotlib numpy pandas scikit-learn scipy statsmodels reportlab
source deactivate

