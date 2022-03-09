# Copyright 2018-2019 Leidos Inc. Naval Medical Research Center Biological Defense Research Directorate
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
echo "${CONDA_EXE}"

${CONDA_EXE} update -y conda
${CONDA_EXE} create -y -n exempliphi
${CONDA_EXE} install -y -n exempliphi -c anaconda luigi pigz lxml networkx sqlite
${CONDA_EXE} install -y -n exempliphi -c conda-forge bzip2
${CONDA_EXE} install -y -n exempliphi -c bioconda biopython spades seqtk blast samtools perl-string-approx perl-parallel-forkmanager
${CONDA_EXE} install -y -n exempliphi -c etetoolkit ete3
${CONDA_EXE} install -y -n exempliphi -c r r-base
${CONDA_EXE} run -n exempliphi python -m pip install docx-mailmerge
${CONDA_EXE} run -n exempliphi pip install igraph
${CONDA_EXE} create -y --name py2 python=2 matplotlib numpy pandas scikit-learn scipy statsmodels reportlab

cd "$(dirname "$0")/blastx_dbs"
${CONDA_EXE} run -n exempliphi makeblastdb -dbtype prot -input_type fasta -in phage_integrases.fasta -title "Phage Integrase Database" -out phage_integrase_db
${CONDA_EXE} run -n exempliphi makeblastdb -dbtype prot -input_type fasta -in phage_terminases.fasta -title "Phage Terminase Database" -out phage_terminase_db

