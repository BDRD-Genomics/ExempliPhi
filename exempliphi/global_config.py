import luigi
import os
from datetime import datetime


class Globals(luigi.Config):
    '''
    Global variables. Set using luigi configuration file.
    '''
    # Path to installation of exempliphi
    PIPELINE_ROOT = luigi.Parameter(default=os.path.split(os.path.dirname(os.path.realpath(__file__)))[0])
    # Directory to write output to
    OUTPUT_DIR = luigi.Parameter(default=os.path.join(os.path.split(os.path.dirname(os.path.realpath(__file__)))[0], 'output'))
    # Number of cores available in computing environment
    NUM_THREADS = luigi.IntParameter(default=1)
    # Coda environment holding dependencies
    PRIMARY_CONDA_ENV = luigi.Parameter(default='exempliphi')
    # Path to conda executable. Default expects it to be in PATH.
    CONDA_EXE = luigi.Parameter(default='conda')
    # Lower boundary of insert size
    INSERT_SIZE_LB = luigi.IntParameter(default=250)  # Defaults were recommendations from Dr. Ken Frey
    # Upper boundary of insert size
    INSERT_SIZE_UB = luigi.IntParameter(default=500)
    # NCBI nt blast database
    NT = luigi.Parameter()
    # Number of workers to use when running pipeline
    NUM_WORKERS = luigi.IntParameter(default=10)
    # What cluster scheduling software do you use? SGE, SLURM, or False for single node
    CLUSTER = luigi.Parameter(default=False)
    # Date of run
    RUN_DATE = luigi.DateSecondParameter(default=datetime.now())


global_config = Globals()
