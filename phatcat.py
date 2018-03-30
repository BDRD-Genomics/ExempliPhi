import luigi
from luigi.configuration import LuigiConfigParser
import json
import pickle
from mailmerge import MailMerge

from phatcat.pipeline import *


def generate_report(output_path, PIPELINE_ROOT, merge_values):
    '''
    Generate a report by mail merging merge_values into Phage_Genomics_Report_Template.docx
    :param output_path: Full file path to output report (ex. foo/bar/reportname.docx)
    :param merge_values: Values to merge into template. Open template to see fields.
    '''
    # Create a document from template
    template = os.path.join(PIPELINE_ROOT, 'Phage_Genomics_Report_Template.docx')
    document = MailMerge(template)

    # Convert all values to strings
    for key in merge_values.keys():
        merge_values[key] = str(merge_values[key])

    # Add collected values to template
    document.merge(**merge_values)

    # Save document
    document.write(output_path)


if __name__ == '__main__':
    print('Using PHATCAT version 0.1.0 (beta)')

    if not len(sys.argv) > 1:
        # TODO: print help
        raise RuntimeError('Need phage parameter json as input')

    # Load in parameters for each phage from json file
    with open(sys.argv[1]) as params:
        parameter_list = json.load(params)

    # Load configurations from luigi.config
    cfg = LuigiConfigParser().instance()
    OUTPUT_DIR = cfg.get('Globals', 'OUTPUT_DIR',
                         default=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'output'))
    PIPELINE_ROOT = cfg.get('Globals', 'PIPELINE_ROOT', default=os.path.dirname(os.path.realpath(__file__)))

    job_array = []
    for phage_params in parameter_list:
        job_array.append(Coding_Complete_Pipeline(**phage_params))

    if 'local' in sys.argv:
        luigi.build(job_array, no_lock=False, local_scheduler=True)
    else:
        NUM_NODES = cfg.get('Globals', 'NUM_NODES', default=8)
        luigi.build(job_array, no_lock=True, workers=NUM_NODES)

    for phage_params in parameter_list:
        # Create folder for report - if old report exists, delete it
        report_dir = os.path.join(OUTPUT_DIR, phage_params['phage_name'], 'report')
        report_path = os.path.join(report_dir, '%s_Phage_Genomics_Report.docx' % phage_params['phage_name'])
        if os.path.isfile(report_path):
            os.remove(report_path)
        elif not os.path.isdir(report_dir):
            os.mkdir(report_dir)

        # Collect merge values
        merge_values = {}
        for root, dirs, files in os.walk(os.path.join(OUTPUT_DIR, phage_params['phage_name'])):
            for file in files:
                if file[-14:] == 'merge_values.p':
                    merge_values.update(pickle.load(open(os.path.join(root, file), 'rb')))

        generate_report(report_path, PIPELINE_ROOT, merge_values)

    # Collect
    print('Here')
