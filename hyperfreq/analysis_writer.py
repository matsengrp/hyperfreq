"""This module provides the functionality for writing out analysis results"""

import csv
import os
import logging
import copy


BASE_COLNAMES = ['sequence', 'hm_pos', 'cutoff_cdf', 'map', 'ltmap', 'fisher_pvalue',
                'focus_pos', 'control_pos', 'focus_neg', 'control_neg']


def quant_namer(quant):
    n_part = str(quant).split('.')[1]
    if len(n_part) == 1:
        n_part += "0"
    return "q" + n_part


class AnalysisWriter(object):
    def __init__(self, prefix, pattern_name, quants, cdfs):

        # Go through and make the appropriate colnames
        colnames = copy.copy(BASE_COLNAMES)
        if pattern_name == "call":
            colnames.insert(2, 'call_pattern')
        colnames += ["q_{0}".format(q) for q in quants]
        colnames += ["cdf_{0}".format(c) for c in cdfs]

        # Create the file magicks
        filename = '{0}.{1}.csv'.format(prefix, pattern_name)
        self.pattern_name = pattern_name
        self.handle = file(filename, 'w')
        self.writer = csv.DictWriter(self.handle, colnames, extrasaction='ignore')
        self.writer.writeheader()

    def writerow(self, data):
        self.writer.writerow(data[self.pattern_name])

    def close(self):
        self.handle.close()


def write_analysis(analysis, prefix, patterns_names, quants, cdfs, call_only=True):
    """This function orchestrates the entire analysis result dump. By default, it only writes out call pattern
    analysis results and the sites file, containing information about mutated columns. Optionally, it can also
    output analysis results for each of the contexts separately."""


    log_filename = prefix + '.log'
    # Don't want these getting all hugelike for now
    try:
        os.remove(log_filename)
    except OSError:
        pass
    logging.basicConfig(filename=log_filename)
    logging.captureWarnings(True)

    def handler_builder(pattern_name):
        # little helper to tidy things
        return AnalysisWriter(prefix, pattern_name, quants, cdfs)

    # Create the analysis writers
    anal_writers = [handler_builder('call')]
    if not call_only:
        anal_writers += [handler_builder(p) for p in patterns_names]

    sites_handle = file(prefix + '.sites.csv', 'w')
    sites_writer = csv.writer(sites_handle)
    sites_writer.writerow(['sequence', 'call_pattern', 'context', 'column'])

    for result in analysis:
        for anal_writer in anal_writers:
            anal_writer.writerow(result)

        base_site_row = [result['call'][x] for x in ('sequence', 'call_pattern')]
        mut_columns = result['call']['mut_columns']
        mut_contexts = result['call']['mut_contexts']

        for column, context in zip(mut_columns, mut_contexts):
            sites_writer.writerow(base_site_row + [context, column])

    sites_handle.close()
    for handler in anal_writers:
        handler.close()


