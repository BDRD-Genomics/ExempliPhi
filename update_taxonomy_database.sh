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

#!/bin/bash
# This script downloads/updates the ExempliPhi taxonomy database.

# Files needed
# taxdump.tar.gz (ETE can fetch)

# Configs needed
# Database file location
set -e

curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz | gzip -d > nucl_gb.accession2taxid
curl -L http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz | gzip -d > taxdump

dbfile=$(${CONDA_EXE} run -n exempliphi python -c '
import sys, os, configparser
from ete3 import NCBITaxa

def blockStdout():
    sys.stdout = open(os.devnull, "w")
def enableStdout():
    sys.stdout = sys.__stdout__

blockStdout()
luigi_config = configparser.ConfigParser()
luigi_config.read("luigi.cfg")
dbfile = luigi_config["Taxonomic_Profiling_Classify"]["taxdb_loc"]
enableStdout()
print(dbfile, end="")
blockStdout()
ncbi = NCBITaxa(dbfile=dbfile, taxdump_file="taxdump")
ncbi.update_taxonomy_database()
')

cat <<EOF > add_accession_mapping_to_database.sql
DROP TABLE IF EXISTS acc2taxid;
CREATE TABLE acc2taxid (
  accession TEXT PRIMARY KEY,
  accession_version UNIQUE,
  taxid integer NOT NULL,
  gi integer
);
.separator "\t"
.import nucl_gb.accession2taxid acc2taxid
DELETE FROM acc2taxid WHERE accession_version="accession.version";
EOF

${CONDA_EXE} run -n exempliphi sqlite3 -init add_accession_mapping_to_database.sql ${dbfile}

rm nucl_gb.accession2taxid
rm taxdump
rm add_accession_mapping_to_database.sql

echo "Update Complete";
