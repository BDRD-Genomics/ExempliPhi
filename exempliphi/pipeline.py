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
import luigi
import subprocess
import os
import shutil
import logging
import re
import types
import sys
import pickle
import platform
from collections import OrderedDict, namedtuple
from datetime import datetime
from io import StringIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


from luigi_cluster.sge import SGEJobTask
from luigi_cluster.slurm import SLURMJobTask
from exempliphi.global_config import global_config


'''
************************************************************
                     GLOBAL FUNCTIONS
************************************************************
'''
def get_stats_from_assembly(fasta):
    '''
    Calculate:
        - # contigs over 700bp
        - Largest contig coverage
        - Largest contig size
        - N50
        - Total assembly size
    :param fasta: assembly file to read in
    '''
    # Look through contigs file and find longest sequence
    record_list = []
    num_contigs_gt700 = 0
    for record in SeqIO.parse(fasta, 'fasta'):
        record_list.append(record)
        if len(record.seq) > 700:
            num_contigs_gt700 += 1

    record_list.sort(key=lambda record: len(record.seq), reverse=True)

    # Calculate N50 and N90
    assembly_size = sum(len(record.seq) for record in record_list)
    running_sum = 0
    n50 = None
    n90 = None
    for i in range(len(record_list)):
        current_contig_len = len(record_list[i].seq)
        running_sum += current_contig_len

        # N50
        if running_sum >= assembly_size * 0.5 and n50 is None:
            n50 = current_contig_len

        # N90
        if running_sum >= assembly_size * 0.9:
            n90 = current_contig_len
            break

    longest_contig_length = len(record_list[0].seq)

    return n50, n90, longest_contig_length, num_contigs_gt700, assembly_size


'''
************************************************************
                          CLASSES
************************************************************
'''
class ClusterTaskParameters(luigi.Config):
    """
    A Generic class to set parameters for the SLURMJobTask
    The SLURMJobTask is not reading in parameters from luigi.cfg
    """
    shared_tmp_dir = luigi.Parameter(default=global_config.OUTPUT_DIR)
    no_tarball = luigi.BoolParameter(default=True)
    default_sge_queue = luigi.Parameter(default='all.q')
    default_slurm_partition = luigi.Parameter(default='normal')
    dont_remove_tmp_dir = luigi.BoolParameter(default=False)
    poll_time = luigi.IntParameter(default=5)


if global_config.CLUSTER == 'SGE':
    BaseTask = SGEJobTask
elif global_config.CLUSTER == 'SLURM':
    BaseTask = SLURMJobTask
else:
    BaseTask = luigi.Task


class ExempliTask(BaseTask):
    '''
    This is a base class for other classes in the pipeline to inherit from. It helps set the structure for a task.
    '''
    n_cpu = luigi.IntParameter(default=global_config.NUM_THREADS)
    shared_tmp_dir = ClusterTaskParameters().shared_tmp_dir
    no_tarball = ClusterTaskParameters().no_tarball
    queue = ClusterTaskParameters().default_sge_queue
    partition = ClusterTaskParameters().default_slurm_partition
    dont_remove_tmp_dir = ClusterTaskParameters().dont_remove_tmp_dir
    poll_time = ClusterTaskParameters().poll_time
    run_date = luigi.DateSecondParameter(default=global_config.RUN_DATE)
    base_dir = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super(ExempliTask, self).__init__(*args, **kwargs)
        self.job_name = "%s__%s" % (self.task_family, self.phage_name)
        self.failed = False

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task. Each task must have a unique output folder.
        '''
        raise NotImplementedError

    def out_file_path(self, temp=False):
        '''
        Returns location of output files. Must return a dict (keys identify file, values are full paths to file).
        Must have keys for 'merge_values' and 'report'!
        '''
        raise NotImplementedError

    def output(self):
        output_files = self.out_file_path()
        target_dict = {}

        for key in output_files.keys():
            target_dict[key] = luigi.LocalTarget(output_files[key])

        return target_dict

    def complete(self):
        if os.path.isdir(self.out_dir()):
            return True
        else:
            return False

    def work(self):
        self.start_task()

    def run(self):
        if global_config.CLUSTER == "SGE" or global_config.CLUSTER == "SLURM":
            return super().run()

        else:
            return self.start_task()

    def start_task(self):
        # Generate task id
        param_map = {}
        for param in self.param_kwargs.keys():
            param_map[param] = str(self.param_kwargs[param])

        self.task_uid = '%s_%s' % (
            luigi.task.task_id_str(self.get_task_family(), param_map),
            datetime.now().strftime('%Y_%m_%d_%H_%M_%S_%f')
        )

        # Create output folder
        self._create_output_folder(self.out_dir(temp=True))

        # Set up logger
        self.logger = logging.getLogger(self.task_uid)
        self.logger.setLevel(logging.DEBUG)
        self.fh = logging.FileHandler(os.path.join(self.out_dir(temp=True), '%s.log' % self.task_uid))
        self.fh.setFormatter(logging.Formatter('\n%(asctime)s level=%(levelname)s:\n%(message)s'))
        self.logger.addHandler(self.fh)

        # Show task info in log
        if global_config.CLUSTER is not False:
            init_info_message = '##### Running %s on %s #####\nParameters:' % (type(self).__name__, platform.node())
        else:
            init_info_message = '##### Running %s #####\nParameters:' % type(self).__name__
        for key in param_map.keys():
            init_info_message = '%s\n\t%s = %s' % (init_info_message, key, param_map[key])
        self.logger.info(init_info_message)

        self.logger.info('Task uid = {}'.format(self.task_uid))

        # Values to be added to report are put into this dictionary
        self.merge_values = {}

        try:
            dep_gen = self.do_task()
            if isinstance(dep_gen, types.GeneratorType):
                self.logger.info('Dynamic dependencies encountered. Attempting to resolve.')
                return dep_gen

        except Exception:
            self.logger.exception('Task %s failed!' % type(self).__name__)
            raise

        finally:
            self._update_merge_vals()

        assert self._check_output()
        self.logger.info('Task completed without errors!')

    def do_task(self):
        '''
        Implement task logic here!
        '''
        raise NotImplementedError

    def on_failure(self, exception):
        self.failed = True
        super(ExempliTask, self).on_failure(exception)

    def _run_command(self, command_params, condaenv='', **kwargs):
        '''
        Wrapper to run commands using subprocess given a string
        :param command_params: string of commands to be run
        :param stdout: file path to stdout
        :param stderr: file path to stderr
        :return: Nothing
        '''
        if condaenv:  # Prepend conda activate condaenv && to command
            command_params.insert(0, condaenv)
            command_params.insert(0, '-n')
            command_params.insert(0, 'run')
            command_params.insert(0, global_config.CONDA_EXE)
        elif global_config.PRIMARY_CONDA_ENV:
            command_params.insert(0, global_config.PRIMARY_CONDA_ENV)
            command_params.insert(0, '-n')
            command_params.insert(0, 'run')
            command_params.insert(0, global_config.CONDA_EXE)

        # Ensure all params are strings
        command_params = [str(x) for x in command_params]

        # Prepare command_string
        command_string = ' '.join(command_params)

        self.logger.info('Running command:\n%s' % command_string.strip())

        # Execute command
        cmd_result = subprocess.run(command_string, shell=True, stderr=subprocess.PIPE, **kwargs)

        if cmd_result.stderr:
            self.logger.warning('Command produced output to stderr:\n%s' % cmd_result.stderr.decode('unicode_escape'))

        if cmd_result.returncode != 0:
            self.logger.error('Command exited with error code %i' % cmd_result.returncode)

        return cmd_result

    def _check_output(self):
        '''
        Checks the temporary dir to make sure all required files are there and not empty
        :param task: A reference to the task which is calling this function
        :return: True if successful, False otherwise
        '''
        self.logger.info('Checking temp folder for task %s' % self.task_id)
        for key in self.out_file_path(temp=True).keys():
            file_path = self.out_file_path(temp=True)[key]
            self.logger.info('Checking for existence of {}'.format(file_path))
            if not os.path.isfile(file_path):
                self.logger.error('OUTPUT FILE MISSING: file %s was not created' % file_path)
                return False
            else:
                self.logger.info('{} exists. Now checking if file is empty.'.format(file_path))
                if os.path.getsize(file_path) == 0:
                    self.logger.error('OUTPUT FILE EMPTY: file %s is empty' % file_path)
                    return False

        # All files exist and are not empty, so now we will rename the folder
        self.logger.info('Task passed all checks. Renaming directory...')
        os.rename(self.out_dir(temp=True), self.out_dir(temp=False))
        # time.sleep(2)
        return True

    def _update_merge_vals(self):
        '''
        Collect all merge values from previous steps (all tasks before this in dependency tree) and then integrate new
        merge values saved in self.merge_values. This results in a new dictionary containing all merge values up to and
        including this task. This dictionary is used to generate a report and then saved in a pickle.
        '''
        merge_vals = {}
        if 'merge_values' in self.input():
            with open(self.input()['merge_values'].path, 'rb') as fh:
                merge_vals.update(pickle.load(fh))
        else:
            for key in self.requires():
                with open(self.input()[key]['merge_values'].path, 'rb') as fh:
                    merge_vals.update(pickle.load(fh))
        merge_vals.update(self.merge_values)
        with open(self.out_file_path(temp=True)['merge_values'], "wb") as fh:
            pickle.dump(merge_vals, fh)

    def _create_output_folder(self, path):
        head, tail = os.path.split(path)
        if head and not os.path.isdir(head):
            self._create_output_folder(head)

        if not os.path.isdir(path):
            try:
                os.mkdir(path)

            except FileExistsError as e:
                # Race condition encountered. Shouldn't cause further problems
                pass



class SeqRecord_Header_Handler():
    '''
    This class is for adding/extracting/and deleting key-value information in fasta headers parsed by biopython
    '''
    def __init__(self, SeqRecord):
        self.record = SeqRecord

    def get_header_features(self):
        '''
        Returns header Key-value pairs as a dict
        '''
        return {kv.split('=')[0]: kv.split('=')[1] if len(kv.split('=')) > 1 else None for
                kv in self.record.description.split()[1:]}

    def set_header_features(self, **kwargs):
        '''
        Sets header key-values to represent input kwargs
        :param feature_dict: python dict storing KV pairs to be written into header
        '''
        description = self.record.id

        for key in kwargs.keys():
            description = '%s %s=%s' % (description, str(key), str(kwargs[key]))

        self.record.description = description

    def delete_header_feature(self, key):
        '''
        Delete a key-value pair from header
        :param key: key of pair to remove
        '''
        current_features = self.get_header_features()
        if key in current_features:
            del current_features[key]
            self.set_header_features(**current_features)
        else:
            raise ValueError('Key %s not found in features' % str(key))

    def add_header_feature(self, key, value):
        '''
        Add key-value feature to a header
        '''
        current_features = self.get_header_features()
        current_features[str(key)] = str(value)
        self.set_header_features(**current_features)


'''
************************************************************
                     TASKS BEGIN HERE
************************************************************
'''
class Copy_Reads(ExempliTask):
    '''
    This task copies reads from illumina folder to a local directory
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # SLURM specific
    n_cpu = 1

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task.
        '''
        folder = 'raw-temp' if temp else 'raw'
        return os.path.join(self.base_dir, self.phage_name, 'reads', folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        file_extension = 'fastq.gz' if self.r1.endswith('.gz') else 'fastq'
        return {
            'r1': os.path.join(self.out_dir(temp), f'%s_raw.1.{file_extension}' % self.phage_name),
            'r2': os.path.join(self.out_dir(temp), f'%s_raw.2.{file_extension}' % self.phage_name),
            'merge_values': os.path.join(self.out_dir(temp), '%s_Copy_Reads_merge_values.p' % self.phage_name)
        }

    def do_task(self):
        self.logger.info('output file paths: {}'.format(self.out_file_path()))
        shutil.copyfile(self.r1, self.out_file_path(temp=True)['r1'])
        shutil.copyfile(self.r2, self.out_file_path(temp=True)['r2'])
        self.merge_values = self._get_merge_values()

    def _get_merge_values(self):
        '''
        Collect merge values for report which were made available at this step.
        :return: dict of merge values
        '''
        merge_values = {}

        # Phage name
        merge_values['Phage_Name'] = self.phage_name

        # Name of the sequencing run
        self._sequencing_run_name(merge_values)

        # Date data received
        merge_values['Date_Data_Received'] = datetime.fromtimestamp(
            os.path.getctime(self.r1)).strftime('%B %d, %Y %H:%M:%S')

        # Date pipeline ran
        merge_values['DateAnalysisRan'] = datetime.now().strftime('%B %d, %Y %H:%M:%S')

        # Total number of reads
        merge_values['TotalNumReads'] = self._count_fastq(self.out_file_path(temp=True)['r1']) * 2

        # Location of output
        merge_values['OutputLocation'] = os.path.join(self.base_dir, self.phage_name)

        return merge_values

    def _sequencing_run_name(self, merge_values):
        '''
        Determine sequencing run name from path to raw reads if possible
        '''
        head, tail = os.path.split(os.path.split(self.r1)[0])
        if tail == "BaseCalls":
            head, tail = os.path.split(head)
            if tail == "Intensities":
                head, tail = os.path.split(head)
                if tail == "Data":
                    merge_values['Sequencing_Run_Name'] = os.path.split(head)[1]

    def _count_fastq(self, file_path):
        '''
        Counts number of reads in a fastq file
        '''
        root, ext = os.path.splitext(file_path)
        if ext == '.gz':
            command_params = [
                'pigz', '-dc', file_path, '|', 'wc', '-l'
            ]
        else:
            command_params = [
                'wc', '-l', file_path
            ]

        cmd_out = self._run_command(command_params, stdout=subprocess.PIPE)
        num_reads = int(int(cmd_out.stdout.split()[0]) / 4)

        return num_reads


class EDGE_QC(ExempliTask):
    '''
    Perform Quality filtering/trimming using the same script used by EDGE
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # Trim quality level
    q = luigi.IntParameter(default=30)

    # Average quality cutoff
    avg_q = luigi.IntParameter(default=30)

    # Minimum read length
    min_L = luigi.IntParameter(default=50)

    # "N" base cutoff
    n = luigi.IntParameter(default=2)

    # Low complexity filter ratio
    lc = luigi.FloatParameter(default=0.85)

    # Adapter FASTA
    adapter = luigi.Parameter(default="")

    # Cut #bp from ends
    fiveEnd = luigi.IntParameter(default=0)
    threeEnd = luigi.IntParameter(default=0)

    def _check_output(self):
        '''
        Override check output so it doesn't fail when no trimmed reads are created
        '''
        self.logger.info('Checking temp folder for task %s' % self.task_uid)
        for key in self.out_file_path(temp=True).keys():

            file_path = self.out_file_path(temp=True)[key]
            if not os.path.isfile(file_path):
                self.logger.error('OUTPUT FILE MISSING: file %s was not created' % file_path)
                return False

            elif key != 's':
                if os.path.getsize(file_path) == 0:
                    self.logger.error('OUTPUT FILE EMPTY: file %s is empty' % file_path)
                    return False

        # All files exist and are not empty, so now we will move the folder
        os.rename(self.out_dir(temp=True), self.out_dir(temp=False))
        return True

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        folder = 'edge-temp' if temp else 'edge'
        return os.path.join(self.base_dir, self.phage_name, 'reads', folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        return {
            'r1': os.path.join(self.out_dir(temp), '%s.1.trimmed.fastq' % self.phage_name),
            'r2': os.path.join(self.out_dir(temp), '%s.2.trimmed.fastq' % self.phage_name),
            's': os.path.join(self.out_dir(temp), '%s.unpaired.trimmed.fastq' % self.phage_name),
            'stats': os.path.join(self.out_dir(temp), '%s.stats.txt' % self.phage_name),
            'qc_pdf': os.path.join(self.out_dir(temp), '%s_qc_report.pdf' % self.phage_name),
            'merge_values': os.path.join(self.out_dir(temp), '%s_EDGE_QC_merge_values.p' % self.phage_name)
        }

    def requires(self):
        return Copy_Reads(r1=self.r1, r2=self.r2, phage_name=self.phage_name, base_dir=self.base_dir)

    def do_task(self):
        command_params = [
            'perl', os.path.join(global_config.PIPELINE_ROOT, 'illumina_fastq_QC.pl'),
            '-p', self.input()['r1'].path, self.input()['r2'].path,
            '-q', self.q,
            '-avg_q', self.avg_q,
            '-min_L', self.min_L,
            '-n', self.n,
            '-lc', self.lc,
            '-5end', self.fiveEnd,
            '-3end', self.threeEnd,
            '-d', self.out_dir(temp=True),
            '-prefix', self.phage_name,
            '-t', self.n_cpu
        ]

        if self.adapter:
            command_params.extend(['-adapter', self.adapter])

        self._run_command(command_params)

        self.merge_values = self._get_merge_values()

    def _get_merge_values(self):
        '''
        Collect merge values for report which were made available at this step.
        :return: dict of merge values
        '''
        merge_values = {}

        # Average quality before QC
        qual_matrix = os.path.join(self.out_dir(temp=True), 'qa.%s.quality.matrix' % self.phage_name)
        merge_values['AvgQual'] = self._get_avg_qual(qual_matrix)
        os.remove(qual_matrix)

        # Parse number/percent of reads which passed EDGE QC from stats file
        after_trimming = False
        for line in open(self.out_file_path(temp=True)['stats'], 'r'):
            if line.startswith('Reads #: ') and not after_trimming:
                merge_values['TotalNumReads'] = int(line.strip()[9:])

            elif line.startswith('After Trimming'):
                after_trimming = True

            elif line.startswith('Reads #: ') and after_trimming:
                merge_values['NumQCReadsEDGE'] = int(line.strip()[9:])
                merge_values['PercentQCReadsEDGE'] = '{:.2%}'.format(
                    float(line.strip()[9:]) / float(merge_values['TotalNumReads'])
                )

        return merge_values


    def _get_avg_qual(self, qual_matrix):
        '''
        Calculate the average quality from quality matrix (read counts for nucleotide pos. x q-score)
        '''
        num_bases = 0
        total_score = 0
        for line in open(qual_matrix):
            read_counts = line.split('\t')
            for i in range(0, len(read_counts)):
                count = int(read_counts[i])
                num_bases += count
                total_score += count * i

        return total_score/num_bases


class CLC_QC(ExempliTask):
    '''
    Perform Quality filtering/trimming using the CLC Assembly Cell's script
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # Set the minimum quality for a good nucleotide
    cutoff = luigi.IntParameter(default=30)

    # Set the maximum fraction of bad nucleotides to define a good quality region
    badfraction = luigi.FloatParameter(default=0.1)

    # Set the fraction of the read that must be of good quality
    lengthfraction = luigi.FloatParameter(default=0)

    # Set the minimum length of output reads
    minlength = luigi.IntParameter(default=50)

    # Cluster specific
    n_cpu = 1

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        folder = 'clc-temp' if temp else 'clc'
        return os.path.join(self.base_dir, self.phage_name, 'reads', folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        return {
            'r1': os.path.join(self.out_dir(temp), '%s.1.clc-trimmed.fastq' % self.phage_name),
            'r2': os.path.join(self.out_dir(temp), '%s.2.clc-trimmed.fastq' % self.phage_name),
            'stats': os.path.join(self.out_dir(temp), '%s.clc-qc.stats.txt' % self.phage_name),
            'merge_values': os.path.join(self.out_dir(temp), '%s_CLC_QC_merge_values.p' % self.phage_name)
        }

    def requires(self):
        return Copy_Reads(r1=self.r1, r2=self.r2, phage_name=self.phage_name, base_dir=self.base_dir)

    def do_task(self):
        interleaved_path = os.path.join(self.out_dir(temp=True), '%s.clc-trimmed.interleaved.fastq' % self.phage_name)

        command_params = [
            'clc_quality_trim',
            '-c', self.cutoff,
            '-m', self.minlength,
            '-l', self.lengthfraction,
            '-b', self.badfraction,
            '-p', interleaved_path,
            '-r', '-i', self.input()['r1'].path, self.input()['r2'].path
        ]

        with open(self.out_file_path(temp=True)['stats'], 'w') as stats:
            self._run_command(command_params, stdout=stats)

        # Pull apart interleaved reads
        i = 0
        with open(self.out_file_path(temp=True)['r1'], 'w') as r1, \
                open(self.out_file_path(temp=True)['r2'], 'w') as r2, open(interleaved_path) as handle:

            for title, seq, qual in SeqIO.QualityIO.FastqGeneralIterator(handle):
                if i % 2 == 0:
                    r1.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                else:
                    r2.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
                i += 1

        os.remove(interleaved_path)

        self.merge_values = self._get_merge_values()

    def _get_merge_values(self):
        '''
        Collect merge values for report which were made available at this step.
        :return: dict of merge values
        '''
        merge_values = {}

        # Parse number/percent of reads which passed CLC QC from stats file
        for line in open(self.out_file_path(temp=True)['stats'], 'r'):
            if line.startswith('Output reads:'):
                line_list = line.split()
                merge_values['NumQCReadsCLC'] = int(line_list[-3])
                merge_values['PercentQCReadsCLC'] = '%s%%' % line_list[-2]
                break

        return merge_values


class Subsample_Reads(ExempliTask):
    '''
    Randomly sample a certain number of reads. Uses seqtk.
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # Number of reads to sample
    sample_size = luigi.IntParameter()

    # Seed for random number generation
    seed = luigi.IntParameter(default=5)

    # 'CLC' or 'EDGE'
    assembly_branch = luigi.Parameter()

    run_locally = True

    # def __init__(self, *args, **kwargs):
    #     super(Subsample_Reads, self).__init__(*args, **kwargs)
    #     self.job_name = "%s__%s__%s" % (self.task_family, self.phage_name, self.assembly_branch)

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        if self.assembly_branch == 'EDGE':
            folder = 'edge-temp' if temp else 'edge'
        elif self.assembly_branch == 'CLC':
            folder = 'clc-temp' if temp else 'clc'
        else:
            raise ValueError('Assembly branch "%s" not recognized' % self.assembly_branch)

        return os.path.join(self.base_dir, self.phage_name, 'reads', 'subsampled', folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        r1 = '%s_Sampled.%i.1.fastq' % (self.phage_name, self.sample_size)
        r2 = '%s_Sampled.%i.2.fastq' % (self.phage_name, self.sample_size)
        merge_values = '%s_%s_Subsample_Reads_merge_values.p' % (self.phage_name, self.assembly_branch)
        return {
            'r1': os.path.join(self.out_dir(temp), r1),
            'r2': os.path.join(self.out_dir(temp), r2),
            'merge_values': os.path.join(self.out_dir(temp), merge_values)
        }

    def requires(self):
        if self.assembly_branch == 'EDGE':
            return EDGE_QC(r1=self.r1, r2=self.r2, phage_name=self.phage_name, base_dir=self.base_dir)
        elif self.assembly_branch == 'CLC':
            return CLC_QC(r1=self.r1, r2=self.r2, phage_name=self.phage_name, base_dir=self.base_dir)

    def do_task(self):
        command_1_params = [
            'seqtk', 'sample',
            '-s%i' % self.seed,
            self.input()['r1'].path,
            int(self.sample_size/2)
        ]

        command_2_params = [
            'seqtk', 'sample',
            '-s%i' % self.seed,
            self.input()['r2'].path,
            int(self.sample_size/2)
        ]

        # Create files which will hold reads
        with open(self.out_file_path(temp=True)['r1'], 'w') as r1_file:
            self._run_command(command_1_params, stdout=r1_file)

        with open(self.out_file_path(temp=True)['r2'], 'w') as r2_file:
            self._run_command(command_2_params, stdout=r2_file)

        self.merge_values = self._get_merge_values()

    def _get_merge_values(self):
        '''
        Collect merge values for report which were made available at this step.
        :return: dict of merge values
        '''
        merge_values = {}
        if self.sample_size:
            merge_values['SSNumReads'] = self.sample_size
        return merge_values


class SPAdes(ExempliTask):
    '''
    Run SPAdes assembler on a set of reads.
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # Number of reads sampled
    sample_size = luigi.IntParameter(default=0)

    # Run with meta-spades?
    meta = luigi.BoolParameter(default=False)

    def __init__(self, *args, **kwargs):
        super(SPAdes, self).__init__(*args, **kwargs)
        if self.sample_size:
            self.job_name = "%s__%s__subsampled" % (self.task_family, self.phage_name)

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        base = os.path.join(self.base_dir, self.phage_name, 'assemblies')
        if self.sample_size:
            base = os.path.join(base, 'subassemblies')

        folder = 'spades-temp' if temp else 'spades'
        return os.path.join(base, folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        if self.sample_size:
            prefix = '%s_SPAdes_subsampled_%i' % (self.phage_name, self.sample_size)
        else:
            prefix = '%s_SPAdes' % self.phage_name

        return {
            'contigs': os.path.join(self.out_dir(temp), '%s_contigs.fasta' % prefix),
            'merge_values': os.path.join(self.out_dir(temp), '%s_merge_values.p' % prefix)
        }

    def requires(self):
        if self.sample_size:
            return Subsample_Reads(r1=self.r1, r2=self.r2, sample_size=self.sample_size,
                                   phage_name=self.phage_name, assembly_branch='EDGE', base_dir=self.base_dir)
        else:
            return EDGE_QC(r1=self.r1, r2=self.r2, phage_name=self.phage_name, base_dir=self.base_dir)

    def do_task(self):
        command_params = [
            'spades.py',
            '-1', self.input()['r1'].path,
            '-2', self.input()['r2'].path,
            '-t', self.n_cpu,
            '-o', self.out_dir(temp=True)
        ]

        if not self.sample_size:
            if os.path.isfile(self.requires().out_file_path()['s']) and \
                    os.stat(self.requires().out_file_path()['s']).st_size > 0:
                command_params.insert(5, '-s')
                command_params.insert(6, self.requires().out_file_path()['s'])

        if self.meta:
            command_params.insert(1, '--meta')

        self._run_command(command_params)

        self.logger.info('[SPAdes] Reformatting header') # Doesn't work!
        records = []
        header_name = "%s_SPAdes_assembly" % self.phage_name
        if self.sample_size:
            header_name = "%s_subsampled_at_%i" % (header_name, self.sample_size)

        with open(os.path.join(self.out_dir(temp=True), 'contigs.fasta')) as contigs_file:
            i = 1
            for record in SeqIO.parse(contigs_file, 'fasta'):
                record.description = "%s_contig_%i" % (header_name, i)
                record.id = record.description
                record.name = record.description
                records.append(record)
                i += 1

        SeqIO.write(records, self.out_file_path(temp=True)['contigs'], 'fasta')
        os.remove(os.path.join(self.out_dir(temp=True), 'contigs.fasta'))

        self.merge_values = self._get_merge_values()

    def get_max_kmer_size(self, temp=False):
        out_dir = self.out_dir(temp)
        kmer_sizes = []
        sub_dirs = [x for x in os.listdir(out_dir) if os.path.isdir(os.path.join(out_dir, x))]
        for sub_dir in sub_dirs:
            m = re.fullmatch('K([0-9]+)', sub_dir)
            if m:
                kmer_sizes.append(int(m.group(1)))
        return max(kmer_sizes)

    def _get_merge_values(self):
        '''
        Collect merge values for report which were made available at this step.
        :return: dict of merge values
        '''
        merge_values = {}
        n50, n90, longest_contig_length, num_contigs_gt700, assembly_size = \
            get_stats_from_assembly(self.out_file_path(temp=True)['contigs'])

        merge_values['SPAdesMaxKmer'] = self.get_max_kmer_size(temp=True)

        if self.sample_size:
            merge_values['NumSPAdesSSContigs'] = num_contigs_gt700
            merge_values['SPAdesSSLargestContigSize'] = longest_contig_length
            merge_values['SPAdesSSN50'] = n50
            merge_values['SPAdesSSN90'] = n90
            merge_values['SPAdesSSAssemblySize'] = assembly_size

        else:
            merge_values['NumSPAdesContigs'] = num_contigs_gt700
            merge_values['SPAdesLargestContigSize'] = longest_contig_length
            merge_values['SPAdesN50'] = n50
            merge_values['SPAdesN90'] = n90
            merge_values['SPAdesAssemblySize'] = assembly_size

        return merge_values


class CLC_Assembly(ExempliTask):
    """
    Run CLC assembler on a set of reads
    """
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # Number of reads sampled
    sample_size = luigi.IntParameter(default=0)

    # Try to estimate insert size?
    estimate_distances = luigi.BoolParameter(default=True)

    word_size = luigi.IntParameter(default=64)
    bubble_size = luigi.IntParameter(default=200)

    # Cluster specific
    partition = luigi.Parameter(default='clc')
    queue = luigi.Parameter(default='clc.q')

    def __init__(self, *args, **kwargs):
        super(CLC_Assembly, self).__init__(*args, **kwargs)
        if self.sample_size:
            self.job_name = "%s__%s__subsampled" % (self.task_family, self.phage_name)

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        base = os.path.join(self.base_dir, self.phage_name, 'assemblies')
        if self.sample_size:
            base = os.path.join(base, 'subassemblies')

        folder = 'clc-temp' if temp else 'clc'

        return os.path.join(base, folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        if self.sample_size:
            prefix = '%s_CLC_Assembly_subsampled_%i' % (self.phage_name, self.sample_size)
        else:
            prefix= '%s_CLC_Assembly' % self.phage_name

        return {
            'contigs': os.path.join(self.out_dir(temp), '%s_contigs.fasta' % prefix),
            'merge_values': os.path.join(self.out_dir(temp), '%s_merge_values.p' % prefix)
        }

    def requires(self):
        if self.sample_size:
            return Subsample_Reads(r1=self.r1, r2=self.r2, sample_size=self.sample_size,
                                   phage_name=self.phage_name, assembly_branch='CLC', base_dir=self.base_dir)
        else:
            return CLC_QC(r1=self.r1, r2=self.r2, phage_name=self.phage_name, base_dir=self.base_dir)

    def do_task(self):
        command_params = [
            'clc_assembler',
            '-o', self.out_file_path(temp=True)['contigs'],
            '--cpus', self.n_cpu,
            '-w', self.word_size,
            '-b', self.bubble_size,
            '-p', 'fb', 'ss', global_config.INSERT_SIZE_LB, global_config.INSERT_SIZE_UB,
            '-q', '-i', self.input()['r1'].path, self.input()['r2'].path
        ]

        if self.estimate_distances:
            command_params.insert(-4, '--estimatedistances %s_estimated_dist.txt' %
                                  os.path.join(self.out_dir(temp=True), self.phage_name))

        self._run_command(command_params)

        self.merge_values = self._get_merge_values()

    def _get_merge_values(self):
        '''
        Collect merge values for report which were made available at this step.
        :return: dict of merge values
        '''
        merge_values = {}
        n50, n90, longest_contig_length, num_contigs_gt700, assembly_size = \
            get_stats_from_assembly(self.out_file_path(temp=True)['contigs'])

        if self.sample_size:
            merge_values['NumCLCSSContigs'] = num_contigs_gt700
            merge_values['CLCSSLargestContigSize'] = longest_contig_length
            merge_values['CLCSSN50'] = n50
            merge_values['CLCSSN90'] = n90
            merge_values['CLCSSAssemblySize'] = assembly_size

        else:
            merge_values['NumCLCContigs'] = num_contigs_gt700
            merge_values['CLCLargestContigSize'] = longest_contig_length
            merge_values['CLCN50'] = n50
            merge_values['CLCN90'] = n90
            merge_values['CLCAssemblySize'] = assembly_size

        return merge_values


class Map_Reads_To_Assembly(ExempliTask):
    '''
    Use CLC's read mapper to map reads to an assembly and output:
     - coverage metrics
     - BAM alignment file
     - Corrected contigs (from calling consensus)
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # Number of reads sampled
    sample_size = luigi.IntParameter(default=0)

    # 'CLC' or 'EDGE'
    assembly_branch = luigi.Parameter()

    # Cluster specific
    partition = luigi.Parameter(default='clc')
    queue = luigi.Parameter(default='clc.q')

    def __init__(self, *args, **kwargs):
        super(Map_Reads_To_Assembly, self).__init__(*args, **kwargs)
        self.job_name = "%s__%s__%s" % (self.task_family, self.phage_name, self.assembly_branch)
        if self.sample_size:
            self.job_name = "%s__subsampled" % self.job_name

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        folder = 'read_mapping-temp' if temp else 'read_mapping'
        return os.path.join(self.requires()['assembly'].out_dir(), folder)

    def out_file_base(self):
        if self.sample_size:
            return '%s_%s_subsampled_%i_assembly_read_mapping' % \
                   (self.phage_name, self.assembly_branch, self.sample_size)

        else:
            return '%s_%s_assembly_read_mapping' % (self.phage_name, self.assembly_branch)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        out_file_base = self.out_file_base()

        return {
            'mapping': os.path.join(self.out_dir(temp), '%s.cas' % out_file_base),
            'contigs': os.path.join(self.out_dir(temp), '%s.corrected_contigs.fasta' % out_file_base),
            'unmapped_reads': os.path.join(self.out_dir(temp), '%s.unmapped_reads.fasta' % out_file_base),
            'merge_values': os.path.join(self.out_dir(temp), '%s_merge_values.p' % out_file_base)
        }

    def requires(self):
        dependencies = {}

        # SPAdes assembly
        if self.assembly_branch == 'EDGE':
            dependencies['assembly'] = SPAdes(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                              sample_size=self.sample_size, base_dir=self.base_dir)
            if self.sample_size:
                dependencies['reads'] = Subsample_Reads(r1=self.r1, r2=self.r2, sample_size=self.sample_size,
                                                        phage_name=self.phage_name, assembly_branch='EDGE', base_dir=self.base_dir)
            else:
                dependencies['reads'] = EDGE_QC(r1=self.r1, r2=self.r2, phage_name=self.phage_name, base_dir=self.base_dir)

        # CLC assembly
        elif self.assembly_branch == 'CLC':
            dependencies['assembly'] = CLC_Assembly(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                                    sample_size=self.sample_size, base_dir=self.base_dir)
            if self.sample_size:
                dependencies['reads'] = Subsample_Reads(r1=self.r1, r2=self.r2, sample_size=self.sample_size,
                                                        phage_name=self.phage_name, assembly_branch='CLC', base_dir=self.base_dir)
            else:
                dependencies['reads'] = CLC_QC(r1=self.r1, r2=self.r2, phage_name=self.phage_name, base_dir=self.base_dir)

        return dependencies

    def do_task(self):
        # Perform mapping
        self.logger.info('[Map_Reads_To_Assembly] Performing Mapping')
        self._run_command([
            'clc_mapper',
            '-a', 'global',
            '-d', '-z', self.input()['assembly']['contigs'].path,
            '-p', 'fb', 'ss', global_config.INSERT_SIZE_LB, global_config.INSERT_SIZE_UB,
            '-q', '-i', self.input()['reads']['r1'].path, self.input()['reads']['r2'].path,
            '-o', self.out_file_path(temp=True)['mapping']
        ])

        result = self._run_command([
            'clc_mapping_info',
            self.out_file_path(temp=True)['mapping']
        ], stdout=subprocess.PIPE)
        mapping_stats = StringIO(result.stdout.decode())
        percent_reads_mapped, per_contig_stats = self._parse_mapping_stats(mapping_stats)

        self.logger.info('[Map_Reads_To_Assembly] Extracting corrected contigs')
        self._run_command([
            'clc_extract_consensus',
            '-a', self.out_file_path(temp=True)['mapping'],
            '-o', self.out_file_path(temp=True)['contigs']
        ])

        self.logger.info('[Map_Reads_To_Assembly] Extracting unmapped reads')
        self._run_command([
            'clc_unmapped_reads',
            '-a', self.out_file_path(temp=True)['mapping'],
            '-o', self.out_file_path(temp=True)['unmapped_reads']
        ])

        longest_contig_cov = self._update_header_info(per_contig_stats)

        self.merge_values = self._get_merge_values(percent_reads_mapped, longest_contig_cov)

    def _parse_mapping_stats(self, mapping_stats):
        '''
        Takes output from clc_mapping_info command and parses it.
        :return: Returns a python dict mapping contig # to per-contig stats (# reads mapped and average coverage). Also
                 returns total % of reads mapped to assembly.
        '''
        ContigStats = namedtuple('ContigStats', ['reads', 'coverage'])
        per_contig_stats = {} # { contig # : ContigStats }

        in_per_contig_stats = False # Keeps track if we reached per-contig stats section of output
        for line in mapping_stats:
            line_split = line.split()
            if len(line_split) >= 4:
                if line_split[0] == 'Mapped' and line_split[1] == 'reads':
                    percent_reads_mapped = float(line_split[3])

                # Parse per-contig stats here
                if in_per_contig_stats:
                    per_contig_stats[int(line_split[0])] = ContigStats(
                        reads=int(line_split[2]),
                        coverage=float(line_split[3])
                    )

                if line_split[0] == 'Contig' and line_split[1] == 'Sites' and line_split[2] == 'Reads':
                    in_per_contig_stats = True

        return percent_reads_mapped, per_contig_stats

    def _update_header_info(self, per_contig_stats):
        '''
        Add coverage stats and lengths to header
        :return: average coverage of longest contig (finding this here so we don't have to iterate contigs again)
        '''
        longest_contig_cov = 0
        longest_contig_len = 0
        records = []
        with open(self.out_file_path(temp=True)['contigs']) as contigs_file:
            for idx, record in enumerate(SeqIO.parse(contigs_file, 'fasta')):
                h_handler = SeqRecord_Header_Handler(record)
                h_handler.set_header_features(
                    length=len(record.seq),
                    num_reads=per_contig_stats[idx+1].reads,
                    avg_cov=per_contig_stats[idx+1].coverage
                )
                records.append(record)

                if len(record.seq) > longest_contig_len:
                    longest_contig_cov = per_contig_stats[idx+1].coverage
                    longest_contig_len = len(record.seq)

        SeqIO.write(records, self.out_file_path(temp=True)['contigs'], 'fasta')
        return longest_contig_cov

    def _get_merge_values(self, percent_reads_mapped, longest_contig_cov):
        '''
        Collect merge values for report which were made available at this step.
        :return: dict of merge values
        '''
        merge_values = {}

        if self.assembly_branch == 'CLC':
            if self.sample_size:
                merge_values['CLCSSLargestContigCov'] = longest_contig_cov
                merge_values['CLCSSPercReads'] = '{:.2%}'.format(percent_reads_mapped/100)

            else:
                merge_values['CLCLargestContigCov'] = longest_contig_cov
                merge_values['CLCPercReads'] = '{:.2%}'.format(percent_reads_mapped / 100)

        elif self.assembly_branch == 'EDGE':
            if self.sample_size:
                merge_values['SPAdesSSLargestContigCov'] = longest_contig_cov
                merge_values['SPAdesSSPercReads'] = '{:.2%}'.format(percent_reads_mapped / 100)
            else:
                merge_values['SPAdesLargestContigCov'] = longest_contig_cov
                merge_values['SPAdesPercReads'] = '{:.2%}'.format(percent_reads_mapped / 100)

        return merge_values


class Taxonomic_Profiling_MegaBLAST(ExempliTask):
    '''
    Perform taxonomic profiling on sample to classify all contigs from assembly and all unmapped reads.
    Abundances are reported by number of reads. (Contigs are valued by the # of reads mapped to them)
    Step 1
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # CLC, EDGE, or AUTO: Select the assembly to use manually or automatically select quickest option
    assembly_branch = luigi.Parameter(default='AUTO')

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        folder = 'MegaBLAST-temp' if temp else 'MegaBLAST'
        return os.path.join(self.base_dir, self.phage_name, 'taxonomy', folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        return {
            'blast_results_contigs': os.path.join(self.out_dir(temp), '%s_megablast_contigs_results.xml' % self.phage_name),
            'blast_results_reads': os.path.join(self.out_dir(temp), '%s_megablast_reads_results.xml' % self.phage_name),
            'merge_values': os.path.join(self.out_dir(temp), '%s_Taxonomic_Profiling_MegaBLAST_merge_values.p' % self.phage_name)
        }

    def requires(self):
        clc = Map_Reads_To_Assembly(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                    sample_size=0, assembly_branch='CLC', base_dir=self.base_dir)
        edge = Map_Reads_To_Assembly(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                         sample_size=0, assembly_branch='EDGE', base_dir=self.base_dir)

        if self.assembly_branch == 'CLC':
            return {'CLC': clc}
        elif self.assembly_branch == 'EDGE':
            return {'EDGE': edge}
        elif self.assembly_branch == 'AUTO':
            return {
                'CLC': clc,
                'EDGE': edge,
            }

    def do_task(self):
        if self.assembly_branch == 'AUTO':
            # if auto, decide what assembly is best based on N50, SPAdes is preferred in event of tie
            with open(self.input()['CLC']['merge_values'].path, 'rb') as fh:
                clc_merge_values = pickle.load(fh)
            # clc_size = clc_merge_values['CLCAssemblySize']
            clc_n50 = clc_merge_values['CLCN50']

            with open(self.input()['EDGE']['merge_values'].path, 'rb') as fh:
                spades_merge_values = pickle.load(fh)
            # spades_size = spades_merge_values['SPAdesAssemblySize']
            spades_n50 = spades_merge_values['SPAdesN50']

            # clc_size += get_stats_from_assembly(self.input()['CLC']['unmapped_reads'].path)[4]
            # spades_size += get_stats_from_assembly(self.input()['EDGE']['unmapped_reads'].path)[4]

            # Time to choose which data to continue with
            # if clc_size < spades_size:
            if clc_n50 > spades_n50:
                selected_branch = 'CLC'
                self.logger.info('Using CLC assembly for taxonomic profiling')
                contigs = self.input()['CLC']['contigs'].path
                unmapped_reads = self.input()['CLC']['unmapped_reads'].path

            else:
                selected_branch = 'EDGE'
                self.logger.info('Using EDGE assembly for taxonomic profiling')
                contigs = self.input()['EDGE']['contigs'].path
                unmapped_reads = self.input()['EDGE']['unmapped_reads'].path

        else:
            selected_branch = self.assembly_branch
            self.logger.info('Using %s assembly for taxonomic profiling' % self.assembly_branch)
            contigs = self.input()[self.assembly_branch]['contigs'].path
            unmapped_reads = self.input()[self.assembly_branch]['unmapped_reads'].path

        self.logger.info('[Taxonomic_Profiling_MegaBLAST] Blasting Contigs')
        self._run_command([
            'blastn',
            '-task', 'megablast',
            '-db', global_config.NT,
            '-query', contigs,
            '-num_threads', self.n_cpu,
            '-outfmt', '5',  # XML output format
            '-evalue', '0.01',
            '-max_target_seqs', '50',
            '-max_hsps', 50,
            '-out', self.out_file_path(temp=True)['blast_results_contigs']
        ])

        self.logger.info('[Taxonomic_Profiling_MegaBLAST] Blasting Unmapped Reads')
        self._run_command([
            'blastn',
            '-task', 'megablast',
            '-db', global_config.NT,
            '-query', unmapped_reads,
            '-num_threads', self.n_cpu,
            '-outfmt', 5,  # XML output format
            '-evalue', 0.01,
            '-max_target_seqs', 50,
            '-max_hsps', 10,
            '-out', self.out_file_path(temp=True)['blast_results_reads']
        ])

        self.merge_values = {'SelectedBranch': selected_branch}


class Taxonomic_Profiling_Classify(ExempliTask):
    '''
    Parse BLAST results and use them to classify sequences
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # CLC, EDGE, or AUTO: Select the assembly to use manually or automatically select quickest option
    assembly_branch = luigi.Parameter(default='AUTO')

    # [0-1] Abundance of reads a taxonomic classification needs to be considered significant. (0.50 = 50%)
    significance_threshold = luigi.FloatParameter(default=0.002)

    # Use xvfb instead of X server to render images.
    use_xvfb = luigi.BoolParameter(default=False)

    # Location of 'taxa.sqlite'. Include database filename and extension.
    taxdb_loc = luigi.Parameter(default=global_config.PIPELINE_ROOT)

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        folder = 'classification-temp' if temp else 'classification'
        return os.path.join(self.base_dir, self.phage_name, 'taxonomy', folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        return {
            'full_tree': os.path.join(self.out_dir(temp), '%s_read_classification_full_tree.png' % self.phage_name),
            'simple_tree': os.path.join(self.out_dir(temp), '%s_read_classification_simplified_tree.png' % self.phage_name),
            'merge_values': os.path.join(self.out_dir(temp), '%s_Taxonomic_Profiling_Classify_merge_values.p' % self.phage_name)
        }

    def requires(self):
        clc = Map_Reads_To_Assembly(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                    sample_size=0, assembly_branch='CLC', base_dir=self.base_dir)
        edge = Map_Reads_To_Assembly(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                     sample_size=0, assembly_branch='EDGE', base_dir=self.base_dir)

        required = {
            'megablast': Taxonomic_Profiling_MegaBLAST(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                                       assembly_branch=self.assembly_branch, base_dir=self.base_dir),
        }

        if self.assembly_branch == 'CLC':
            required['mapping'] = clc
            return required

        elif self.assembly_branch == 'EDGE':
            required['mapping'] = edge
            return required

        elif self.assembly_branch == 'AUTO':
            required['CLC'] = clc
            required['EDGE'] = edge
            return required

    def do_task(self):
        if self.assembly_branch == 'AUTO':
            # Dependency on read mapping dependent on selected branch
            with open(self.input()['megablast']['merge_values'].path, 'rb') as fh:
                selected_branch = pickle.load(fh)['SelectedBranch']
            mapping = self.input()[selected_branch]

        else:
            mapping = self.input()['mapping']

        command_args = [
            'python', os.path.join(global_config.PIPELINE_ROOT, 'exempliphi', 'taxonomy.py'),
            '--taxdb_loc', self.taxdb_loc,
            '--blast_results_contigs', self.input()['megablast']['blast_results_contigs'].path,
            '--blast_results_reads', self.input()['megablast']['blast_results_reads'].path,
            '--output', self.out_dir(temp=True),
            '--full_tree_name', os.path.split(self.out_file_path(temp=True)['full_tree'])[1],
            '--simplified_tree_name', os.path.split(self.out_file_path(temp=True)['simple_tree'])[1],
            '--sig_threshold', self.significance_threshold,
            '--contigs', mapping['contigs'].path,
            '--unmapped_reads', mapping['unmapped_reads'].path,
            '--num_procs', self.n_cpu
        ]

        # At this time (Jan 2022), the -d argument is not implemented for xvfb-run on ubuntu
        if self.use_xvfb:
            result = subprocess.run(['xvfb-run -d echo'], shell=True)
            if result.returncode != 0:
                flag = '-a'
            else:
                flag = '-d'
            command_args.insert(0, f'xvfb-run {flag}')

        cmd_output = self._run_command(command_args, stdout=subprocess.PIPE)
        self.logger.info('Taxonomy command stdout:\n%s' % cmd_output.stdout.decode('unicode_escape'))


class Extract_Longest_Contig(ExempliTask):
    '''
    Pull out the longest contig from fasta and put it in a separate file
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # Number of reads sampled
    sample_size = luigi.IntParameter(default=0)

    # 'CLC' or 'EDGE'
    assembly_branch = luigi.Parameter()

    # Cluster specific
    n_cpu = 1

    def __init__(self, *args, **kwargs):
        super(Extract_Longest_Contig, self).__init__(*args, **kwargs)
        self.job_name = "%s__%s__%s" % (self.task_family, self.phage_name, self.assembly_branch)
        if self.sample_size:
            self.job_name = "%s__subsampled" % self.job_name

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        folder = 'longest_contig-temp' if temp else 'longest_contig'
        return os.path.join(os.path.dirname(self.input()['contigs'].path), folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        if self.sample_size:
            prefix = '%s_%s_longest_contig_subsampled_%i' % \
                          (self.phage_name, self.assembly_branch, self.sample_size)
        else:
            prefix = '%s_%s_longest_contig' % (self.phage_name, self.assembly_branch)

        return {
            'contig': os.path.join(self.out_dir(temp), '%s.fa' % prefix),
            'merge_values': os.path.join(self.out_dir(temp), '%s_merge_values.p' % prefix)
        }

    def requires(self):
        if self.assembly_branch == 'EDGE':
            return SPAdes(r1=self.r1, r2=self.r2, phage_name=self.phage_name, sample_size=self.sample_size, base_dir=self.base_dir)
        elif self.assembly_branch == 'CLC':
            return Map_Reads_To_Assembly(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                         sample_size=self.sample_size, assembly_branch='CLC', base_dir=self.base_dir)

    def do_task(self):
        # Look through contigs file and find longest sequence
        longest_record = None
        for record in SeqIO.parse(self.input()['contigs'].path, 'fasta'):
            if longest_record is None or len(record.seq) > len(longest_record.seq):
                longest_record = record

        if longest_record is None:
            raise FileNotFoundError('Error: Unable to parse assembly output at path: %s (branch=%s, sample_size=%i)'
                                    % (self.input()['contigs'].path, self.assembly_branch, self.sample_size))

        # Write contig to file
        SeqIO.write(longest_record, self.out_file_path(temp=True)['contig'], 'fasta')


class Remove_Overlap(ExempliTask):
    '''
    Remove the 127 bp SPAdes overlap from a sequence if it exists
        - input = fasta file
        - ouput = new fasta file with overlap removed. If no overlap detected,
        a copy of the input contig is put in output directory
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # Number of reads sampled
    sample_size = luigi.IntParameter(default=0)

    # 'CLC' or 'EDGE'
    assembly_branch = luigi.Parameter()

    # Cluster specific
    n_cpu = 1

    def __init__(self, *args, **kwargs):
        super(Remove_Overlap, self).__init__(*args, **kwargs)
        self.job_name = "%s__%s__%s" % (self.task_family, self.phage_name, self.assembly_branch)
        if self.sample_size:
            self.job_name = "%s__subsampled" % self.job_name

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        folder = 'overlap_removed-temp' if temp else 'overlap_removed'
        return os.path.join(os.path.dirname(self.input()['contig'].path), folder)

    def out_file_base(self):
        if self.sample_size:
            return '%s_%s_longest_contig_overlap_removed_subsampled_%i' % \
                   (self.phage_name, self.assembly_branch, self.sample_size)
        else:
            return '%s_%s_longest_contig_overlap_removed' % (self.phage_name, self.assembly_branch)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        out_file_base = self.out_file_base()
        return {
            'contig': os.path.join(self.out_dir(temp), '%s.fa' % out_file_base),
            'merge_values': os.path.join(self.out_dir(temp), '%s_merge_values.p' % out_file_base)
        }

    def requires(self):
        return Extract_Longest_Contig(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                      sample_size=self.sample_size, assembly_branch=self.assembly_branch, base_dir=self.base_dir)

    def do_task(self):
        # Read input fasta file to get sequence of contig
        record = SeqIO.read(self.input()['contig'].path, 'fasta')
        # save overlap length for merge values
        overlap_size = 0

        if self.assembly_branch == "EDGE":
            overlap_size = self.requires().requires().get_max_kmer_size()
            record.id = "%s_EDGE_Longest_Contig_Overlap_Removed" % self.phage_name

            # Overlaps in SPAdes are 127 bp long, so compare first 127 bases to last 127 bases

            if record.seq[:overlap_size] == record.seq[-overlap_size:]:
                record.seq = record.seq[:-overlap_size]
                record.description = "length_of_overlap=%i" % overlap_size
                self.logger.info('[Remove_Overlap] Overlap detected in SPAdes assembly. '
                                 'Outputting new file with overlap removed.\n')

            else:
                record.description = "length_of_overlap=%i" % 0
                self.logger.warning('[Remove_Overlap] Largest contig from SPAdes does not contain a 127 bp overlap.')

        elif self.assembly_branch == "CLC":
            '''
            CLC overlaps can vary in length and may not be present even if the genome is circular. To deal with this, 
            we remove any overlap regardless of size, but keep the original longest contig to compare 
            in downstream tasks (in case bases on ends were the same due to chance).
            '''
            record.id = "%s_CLC_Longest_Contig_Overlap_Removed" % self.phage_name
            repeat_found = False
            mid_base = int(len(record.seq) / 2)
            for repeat_size in reversed(range(1, mid_base + 1)):
                if record.seq[:repeat_size] == record.seq[-repeat_size:]:
                    record.seq = record.seq[:-repeat_size]
                    self.logger.info('[Remove_Overlap] Overlap detected in CLC assembly. '
                                     'Outputting new file with overlap removed.\n')
                    overlap_size = repeat_size
                    record.description = "length_of_overlap=%i" % repeat_size
                    repeat_found = True
                    break

            if not repeat_found:
                record.description = "length_of_overlap=%i" % 0
                self.logger.info('[Remove_Overlap] No overlap found in largest contig from CLC assembly')

        # Write contig to file
        record.description = "length=%i %s" % (len(record.seq), record.description)
        SeqIO.write(record, self.out_file_path(temp=True)['contig'], 'fasta')

        # Update merge values and generate report
        self.merge_values = self._get_merge_values(overlap_size)

    def _get_merge_values(self, ol_len):
        '''
        Collect merge values for report which were made available at this step.
        :return: dict of merge values
        '''
        merge_values = {}
        if self.assembly_branch == 'CLC':
            if self.sample_size:
                merge_values['CLCSSOverlap'] = ol_len
            else:
                merge_values['CLCOverlap'] = ol_len

        elif self.assembly_branch == 'EDGE':
            if self.sample_size:
                merge_values['SPAdesSSOverlap'] = ol_len
            else:
                merge_values['SPAdesOverlap'] = ol_len

        return merge_values


class Compare_Sequences(ExempliTask):
    '''
    Compare the largest contigs produced by all assemblies and show which are identical
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # Number of reads sampled
    sample_size_clc = luigi.IntParameter(default=0)
    sample_size_edge = luigi.IntParameter(default=0)

    # Cluster specific
    n_cpu = 1

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        folder = 'comparison-temp' if temp else 'comparison'
        return os.path.join(self.base_dir, self.phage_name, 'assemblies', folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        return {
            'results': os.path.join(self.out_dir(temp), '%s_assembly_sequence_comparison.txt' % self.phage_name),
            'representative_contig': os.path.join(self.out_dir(temp), "%s_representative_contig.fasta" % self.phage_name),
            'merge_values': os.path.join(self.out_dir(temp), "%s_Compare_Sequences_merge_values.p" % self.phage_name)
        }

    def requires(self):
        dependencies = {
            'EDGE full assembly': Remove_Overlap(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                                 sample_size=0, assembly_branch='EDGE', base_dir=self.base_dir),
            'CLC full assembly': Remove_Overlap(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                                sample_size=0, assembly_branch='CLC', base_dir=self.base_dir),
            'CLC full assembly (without overlap removal)': Extract_Longest_Contig(r1=self.r1, r2=self.r2,
                                                                                  phage_name=self.phage_name,
                                                                                  sample_size=0, assembly_branch='CLC', base_dir=self.base_dir)
        }
        if self.sample_size_edge:
            dependencies['EDGE subsampled assembly'] = Remove_Overlap(r1=self.r1, r2=self.r2,
                                                                      phage_name=self.phage_name,
                                                                      sample_size=self.sample_size_edge,
                                                                      assembly_branch='EDGE', base_dir=self.base_dir)
        if self.sample_size_clc:
            dependencies['CLC subsampled assembly'] = Remove_Overlap(r1=self.r1, r2=self.r2,
                                                                     phage_name=self.phage_name,
                                                                     sample_size=self.sample_size_clc,
                                                                     assembly_branch='CLC', base_dir=self.base_dir)
            dependencies['CLC subsampled assembly (without overlap removal)'] = \
                Extract_Longest_Contig(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                       sample_size=self.sample_size_clc, assembly_branch='CLC', base_dir=self.base_dir)

        return dependencies

    def do_task(self):
        # Collect sequences to compare
        sequences = [] # List of tuples
        for key in self.requires():
            record = SeqIO.read(self.input()[key]['contig'].path, 'fasta')
            sequences.append((key, record.seq))

        # Bin sequences into sets based off sequence similarity
        list_of_sets = []
        for sequence in sequences:
            sequence_classified = False

            for seq_set in list_of_sets:
                if self._compare_sequences(seq_set[0][1], sequence[1]):
                    seq_set.append(sequence)
                    sequence_classified = True
                    break

            if not sequence_classified:
                list_of_sets.append([sequence])

        # Write results to output file
        with open(self.out_file_path(temp=True)['results'], 'w') as output:
            for seq_set in range(len(list_of_sets)):
                for seq in range(len(list_of_sets[seq_set])):
                    output.write(list_of_sets[seq_set][seq][0])

                    if seq != len(list_of_sets[seq_set]) - 1:
                        output.write(" == ")

                if seq_set != len(list_of_sets) - 1:
                    output.write(" != ")

        # Update merge values and generate report
        self.merge_values = self._get_merge_values()

        # Select a representative contig if possible (Contig which CLC and SPAdes agree on)
        possible_reps = []
        for seq_set in list_of_sets:
            containsEdgeContig = False
            containsCLCContig = False

            for seq in seq_set:
                if 'EDGE' in seq[0]:
                    containsEdgeContig = True

                elif 'CLC' in seq[0]:
                    containsCLCContig = True

                if containsCLCContig and containsEdgeContig:
                    possible_reps.append(seq_set)
                    break

        if len(possible_reps) > 1:
            self.logger.error("[Compare_Sequences] More than one possible representative contigs found! Can not proceed. See %s for details"
                         % self.out_file_path(temp=True)['results'])
            sys.exit(100)

        elif len(possible_reps) == 0:
            self.logger.error("[Compare_Sequences] CLC and SPAdes could not agree upon a contig. Manual intervention required.")
            sys.exit(101)

        # Create SeqRecord object for representative contig
        rep_record = SeqRecord(
            possible_reps[0][0][1],
            id="%s_representative_sequence" % self.phage_name,
            description="length=%i" % len(possible_reps[0][0][1])
        )
        SeqIO.write(rep_record, self.out_file_path(temp=True)['representative_contig'], 'fasta')

    @staticmethod
    def _compare_sequences(seq1, seq2):
        '''
        Compare two sequences to determine if the are the same. This is more involved than just comparing two strings.
        This function will return true if sequences are the same regardless of:
            a) Starting position
            b) Orientation (it looks at reverse complement)
        :seq1/seq2: Sequences to compare. Must be Bio.Seq objects.
        :return: True if sequences are the same, false otherwise
        '''
        # Try to find start of match
        loc = (seq1 + seq1).find(seq2)
        if loc == -1:
            seq2 = seq2.reverse_complement()
            loc = (seq1 + seq1).find(seq2)

            if loc == -1:
                return False

        if (seq1 + seq1)[loc:len(seq1) + loc] == seq2:
            return True

        else:
            return False

    def _get_merge_values(self):
        '''
        Collect merge values for report which were made available at this step.
        :return: dict of merge values
        '''
        merge_values = {}
        with open(self.out_file_path(temp=True)['results'], 'r') as results:
            merge_values['AssemblyComparison'] = results.read()

        return merge_values


class PhageTerm(ExempliTask):
    '''
    Run PhageTerm on a fasta file
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # Number of reads sampled
    sample_size_clc = luigi.IntParameter(default=0)
    sample_size_edge = luigi.IntParameter(default=0)

    # Location of PhageTerm script
    phageterm_loc = luigi.Parameter()

    # Python 2 conda environment
    python2env = luigi.Parameter(default='py2')

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        folder = 'PhageTerm-temp' if temp else 'PhageTerm'
        return os.path.join(self.base_dir, self.phage_name, 'coding_complete', 'annotations', folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        return {
            'phageterm_report': os.path.join(self.out_dir(temp), '%s_PhageTerm_report.pdf' % self.phage_name),
            # 'type': os.path.join(self.out_dir(temp), 'Type'),
            'merge_values': os.path.join(self.out_dir(temp), '%s_Phageterm_merge_values.p' % self.phage_name)
        }

    def requires(self):
        return {
            'contig': Compare_Sequences(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                        sample_size_clc=self.sample_size_clc, sample_size_edge=self.sample_size_edge, base_dir=self.base_dir),
            'reads': EDGE_QC(r1=self.r1, r2=self.r2, phage_name=self.phage_name, base_dir=self.base_dir)
        }

    def do_task(self):
        command_params = [
            self.phageterm_loc,
            '-f', self.input()['reads']['r1'].path,
            '-p', self.input()['reads']['r2'].path,
            '-r', self.input()['contig']['representative_contig'].path,
            '-c', self.n_cpu,
            '-n', self.phage_name
        ]

        # Change working directory so output ends up in output folder (Script has no parameter to set this)
        old_wd = os.getcwd()
        os.chdir(self.out_dir(temp=True))
        self._run_command(command_params, condaenv=self.python2env)
        os.chdir(old_wd)


class BlastX(ExempliTask):
    '''
    Search sequence against specialized blast databases to look for the presense of various types of proteins.
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # Number of reads sampled
    sample_size_clc = luigi.IntParameter(default=0)
    sample_size_edge = luigi.IntParameter(default=0)

    # Location of blast DBs
    terminase_blast_db = luigi.Parameter(
        os.path.join(
            os.path.join(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'blastx_dbs'
            ), 'phage_terminase_db'
        )
    )
    integrase_blast_db = luigi.Parameter(
        os.path.join(
            os.path.join(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'blastx_dbs'
            ), 'phage_integrase_db'
        )
    )

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task
        '''
        folder = 'blastx-temp' if temp else 'blastx'
        return os.path.join(self.base_dir, self.phage_name, 'coding_complete', 'annotations', folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files
        '''
        return {
            'terminase_results': os.path.join(self.out_dir(temp), '%s_terminase_blastx_results.html' % self.phage_name),
            'integrase_results': os.path.join(self.out_dir(temp), '%s_integrase_blastx_results.html' % self.phage_name),
            'merge_values': os.path.join(self.out_dir(temp), '%s_blastx_merge_values.p' % self.phage_name)
        }

    def requires(self):
        return {
            'compare_sequences': Compare_Sequences(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
                                        sample_size_clc=self.sample_size_clc, sample_size_edge=self.sample_size_edge, base_dir=self.base_dir)
            # 'phageterm': PhageTerm(r1=self.r1, r2=self.r2, phage_name=self.phage_name,
            #                             sample_size_clc=self.sample_size_clc, sample_size_edge=self.sample_size_edge, base_dir=self.base_dir)
        }

    def do_task(self):
        # Determine termini type
        # with open(self.input()['phageterm']['type'].path) as type_fh:
        # line = type_fh.readline()
        #
        # # DTR and pac phages use PhageTerm output for blastx query
        # if 'DTR' in line or 'pac' in line:
        #     self.logger.info('Using PhageTerm output sequence as query for blastx')
        #     genome_path = os.path.join(self.requires()['phageterm'].out_dir(),
        #                                '%s_sequence.fasta' % self.phage_name)
        #
        # else:
        self.logger.info('Using Compare_Sequences output sequence as query for blastx')
        genome_path = os.path.join(self.requires()['compare_sequences'].out_dir(),
                                   '%s_representative_contig.fasta' % self.phage_name)

        self.logger.info('Blasting for integrase')
        self._run_command([
            'blastx',
            '-query', genome_path,
            '-db', self.integrase_blast_db,
            '-html',
            '-out', self.out_file_path(temp=True)['integrase_results'],
            '-num_threads', self.n_cpu,
            '-evalue', '0.01'
        ])

        self.logger.info('Blasting for terminase')
        self._run_command([
            'blastx',
            '-query', genome_path,
            '-db', self.terminase_blast_db,
            '-html',
            '-out', self.out_file_path(temp=True)['terminase_results'],
            '-num_threads', self.n_cpu,
            '-evalue', '0.01'
        ])


class Coding_Complete_Pipeline(luigi.WrapperTask):
    '''
    Generate final word document. Starting point for coding complete pipeline.
    '''
    # Paired reads (in fastq format)
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    # Name of the phage
    phage_name = luigi.Parameter()

    # Number of reads sampled from each set of reads
    sample_size = luigi.IntParameter(default=50000)

    perform_taxonomic_profiling = luigi.BoolParameter(default=True)

    base_dir = luigi.Parameter()

    logger = logging.getLogger('luigi-interface')

    def requires(self):
        return {
            'EDGE': EDGE_QC(r1=self.r1, r2=self.r2, phage_name=self.phage_name, base_dir=self.base_dir),
            'CLC': CLC_QC(r1=self.r1, r2=self.r2, phage_name=self.phage_name, base_dir=self.base_dir)
        }

    def additional_depends(self):
        # Need dynamic dependencies so we don't try to subsample a greater number of reads than we have left after QC
        edge_num_reads = pickle.load(open(self.input()['EDGE']['merge_values'].path, 'rb'))['NumQCReadsEDGE']
        clc_num_reads = pickle.load(open(self.input()['CLC']['merge_values'].path, 'rb'))['NumQCReadsCLC']

        sample_size_edge = 0
        if edge_num_reads > self.sample_size:
            sample_size_edge = self.sample_size

        sample_size_clc = 0
        if clc_num_reads > self.sample_size:
            sample_size_clc = self.sample_size

        dependencies = {
            'edge_full_assembly_mapping': Map_Reads_To_Assembly(
                r1=self.r1, r2=self.r2,
                phage_name=self.phage_name,
                sample_size=0,
                assembly_branch='EDGE',
                base_dir=self.base_dir
            ),
            'phage_term': PhageTerm(
                r1=self.r1, r2=self.r2,
                phage_name=self.phage_name,
                sample_size_edge=sample_size_edge,
                sample_size_clc=sample_size_clc,
                base_dir=self.base_dir
            ),
            'blastx': BlastX(
                r1=self.r1, r2=self.r2,
                phage_name=self.phage_name,
                sample_size_edge=sample_size_edge,
                sample_size_clc=sample_size_clc,
                base_dir=self.base_dir
            ),
            'compare_sequences': Compare_Sequences(
                r1=self.r1, r2=self.r2,
                phage_name=self.phage_name,
                sample_size_edge=sample_size_edge,
                sample_size_clc=sample_size_clc,
                base_dir=self.base_dir
            )
        }
        if sample_size_edge:
            dependencies['edge_subsampled_assembly_mapping'] = Map_Reads_To_Assembly(
                r1=self.r1, r2=self.r2,
                phage_name=self.phage_name,
                sample_size=sample_size_edge,
                assembly_branch='EDGE',
                base_dir=self.base_dir
            )

        if self.perform_taxonomic_profiling:
            dependencies['taxonomy'] = Taxonomic_Profiling_Classify(r1=self.r1, r2=self.r2, phage_name=self.phage_name, base_dir=self.base_dir)

        return dependencies

    def complete(self):
        return (all(r.complete() for r in luigi.task.flatten(self.requires())) and all(r.complete() for r in luigi.task.flatten(self.additional_depends()))) or any(r.failed for r in luigi.task.flatten(self.requires()))

    def run(self):
        yield self.additional_depends()


