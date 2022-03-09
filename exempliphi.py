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
from ete3 import NCBITaxa

import argparse
import json
import pickle
import os
from mailmerge import MailMerge


def parse_args():
    arg_parser = argparse.ArgumentParser(
        description='Run ExempliPhi pipeline with data specified in input json file')

    arg_parser.add_argument('json_path', type=str, help='Input json file')

    arg_parser.add_argument('-l', dest='local', action='store_true', help='Run locally (without SGE)')

    arg_parser.add_argument(
        '-o', dest='output_folder', required=False, type=str,
        help='Relative to "OUTPUT_DIR" defined in [Globals] of luigi.cfg')

    return arg_parser.parse_args()


def generate_report(output_path, PIPELINE_ROOT, merge_values):
    '''
    Generate a report by mail merging merge_values into Phage_Genomics_Report_Template.docx
    :param output_path: Full file path to output report (ex. foo/bar/reportname.docx)
    :param merge_values: Values to merge into template. Open template to see fields.
    '''
    # Create a document from template
    template = os.path.join(PIPELINE_ROOT, 'Phage_Genomics_Report_Template.docx')
    document = MailMerge(template)

    # Convert all values to formatted strings
    for key in merge_values.keys():
        if type(merge_values[key]) == int:
            merge_values[key] = '{:,}'.format(merge_values[key])

        elif type(merge_values[key]) == float:
            merge_values[key] = '{:,.2f}'.format(merge_values[key])

        else:
            merge_values[key] = str(merge_values[key])

    # Add collected values to template
    document.merge(**merge_values)

    # Save document
    document.write(output_path)


def main(args):
    print('Using ExempliPhi version 2')
    # Set luigi.cfg
    os.environ['LUIGI_CONFIG_PATH'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'luigi.cfg')
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    import luigi
    from exempliphi.pipeline import Coding_Complete_Pipeline
    from exempliphi.global_config import global_config

    # Load in parameters for each phage from json file
    with open(args.json_path) as params:
        parameter_list = json.load(params)

    if args.output_folder:
        out_dir = os.path.join(global_config.OUTPUT_DIR, args.output_folder)
    else:
        out_dir = global_config.OUTPUT_DIR

    job_array = []
    for phage_params in parameter_list:
        job_array.append(Coding_Complete_Pipeline(**phage_params, base_dir=out_dir))

    if args.local:
        luigi.build(job_array, no_lock=False, local_scheduler=True)
    else:
        # NUM_NODES = cfg.get('Globals', 'NUM_NODES', default=8)
        luigi.build(job_array, no_lock=True, workers=global_config.NUM_WORKERS)

    for phage_params in parameter_list:
        # Create folder for report - if old report exists, delete it
        report_dir = os.path.join(out_dir, phage_params['phage_name'], 'report')
        report_path = os.path.join(report_dir, '%s_Phage_Genomics_Report.docx' % phage_params['phage_name'])
        if os.path.isfile(report_path):
            os.remove(report_path)
        elif not os.path.isdir(report_dir):
            os.mkdir(report_dir)

        # Collect merge values
        merge_values = {}
        for root, dirs, files in os.walk(os.path.join(out_dir, phage_params['phage_name'])):
            for file in files:
                if file[-14:] == 'merge_values.p':
                    merge_values.update(pickle.load(open(os.path.join(root, file), 'rb')))

        pickle.dump(merge_values, open(os.path.join(out_dir, phage_params['phage_name'],
                                                    '%s_merge_values.p' % phage_params['phage_name']), "wb"))
        generate_report(report_path, global_config.PIPELINE_ROOT, merge_values)


if __name__ == '__main__':
    main(parse_args())

