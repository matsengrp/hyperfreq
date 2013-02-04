import csv
from hyperfreq import HyperfreqAlignment


BASE_ROWNAMES = ['sequence', 'cluster', 'br_left', 'br_median', 'br_right', 'fisher_pvalue', 'hm_pos',
                'n_focus_pos', 'n_control_pos', 'n_focus_neg', 'n_control_neg'] 

translations = dict(sequence='name')


class AnalysisWriter(object):
    def __init__(self, writeable):
        self.writeable = writeable


    def gross_header(self, gross_writer):
        sorted_contexts = list(set(self.contexts)) if self.contexts else []
        sorted_contexts.sort(cmp=HyperfreqAlignment.context_sorter)

        by_seq_writer.writerow(HyperfreqAlignment.BASE_ROWNAMES + sorted_contexts)

    def by_seq_header(self, by_seq_writer):
        pass

    def write_gross(self, write_to):
        writer = 
        writer.writerow(['cluster', 'sequence', 'column','context'])
        pass

    def write_by_seq(self, write_to):
        pass


# For sets
def write_analysis(self, gross_handle, by_seq_handle):
    """Once analyzed, write the results to files (once global, or by site, the other by sequence)."""
    sorted_contexts = list(self.contexts)
    sorted_contexts.sort(cmp=HyperfreqAlignment.context_sorter)
    gross_writer = csv.writer(gross_handle)
    by_seq_writer = csv.writer(by_seq_handle)

    gross_writer.writerow(['cluster', 'sequence', 'column', 'context'])

    by_seq_writer.writerow(HyperfreqAlignment.BASE_ROWNAMES + sorted_contexts)

    for cluster in self.cluster_alns:
        aln = self.cluster_alns[cluster]
        aln.write_analysis(gross_writer=gross_writer, by_seq_writer=by_seq_writer, cluster=cluster,
                header=False, contexts=sorted_contexts)

# For alignments
def write_analysis(analysis,
        gross_handle=None, gross_writer=None,
        by_seq_handle=None, by_seq_writer=None,
        cluster=None, header=True, contexts=None):
    """This method can be either used by itself or by the HyperfreqAlign.Set class instances to write for
    many sets. As such, it is possible either to pass in a gross_handle file for csv writer, or a
    csvwriter object which already has column names written to it."""
    if not gross_writer:
        gross_writer = csv.writer(gross_handle)

    if not by_seq_writer:
        by_seq_writer = csv.writer(by_seq_handle)

    if header:
        write_header(gross_writer, by_seq_writer)

    else:
        sorted_contexts = contexts

    for seq in self:
        # Gross writer now does by seq
        if seq.hm_pos:
            for i in seq.focus_pos_indices:
                row = [cluster, seq.name, i+1, self.context(i)]
                gross_writer.writerow(row)

        row = [seq.name, cluster, seq.br_left, seq.br_median, seq.fisher_pvalue, seq.hm_pos, seq.n_focus_pos, seq.n_control_pos,
                seq.n_focus_neg, seq.n_control_neg]
        def get_ctxt(c):
            try:
                return seq.contexts[c]
            except KeyError:
                return 0

        row += [get_ctxt(ctxt) for ctxt in sorted_contexts]
        by_seq_writer.writerow(row)

