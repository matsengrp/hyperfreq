from Bio import Align
from Bio.Align import AlignInfo
import warnings
import sekhon
import csv


class HyperfreqAlignment(Align.MultipleSeqAlignment):

    RESIDUES = ['A', 'C', 'T', 'G']

    MUT_PATTERNS = [(x,y) for x in RESIDUES for y in RESIDUES]
    MUT_PATTERNS.sort()


    @staticmethod
    def context_sorter(c1, c2):
        l1, l2 = len(c1), len(c2)
        if l1 != l2:
            return l1 - l2
        else:
            return 2 * int(c1 > c2) - 1

    @staticmethod
    def trans_col_name(trans_tuple):
        return "{}_to_{}".format(*trans_tuple)


    def __res_ratio__(self, residue, index):
        return float(self[:,index].count(residue)) / float(self.__len__())

    def __res_sites__(self, residue, consensus_threshold=0.5):
        aln_length = self.get_alignment_length()
        return [i for i in xrange(0, aln_length) if
                self.__res_ratio__(residue,i) >= consensus_threshold]

    def __mutation_sites__(self, orig_residue, mut_residue,
            consensus_threshold=0.5, trans_threshold=1):
        possible_sites = self.__res_sites__(orig_residue, consensus_threshold)
        return [i for i in possible_sites if
                self[:,i].count(mut_residue) >= trans_threshold]


    def __init__(self, *args, **kwargs):
        """This inherits from biopythons MultipleSeqAlignment, and does basically the same stuff, but also
        lets one specify a reference sequence for comparison istead of a consensus."""
        #import pdb; pdb.set_trace()
        reference_sequence = kwargs['reference_sequence']
        del kwargs['reference_sequence']
        super(HyperfreqAlignment, self).__init__(*args, **kwargs)
        self.reference_sequence = reference_sequence.seq if reference_sequence else None


    def analyze_hypermuts(self, consensus_threshold=0.5, mut_trans=('G','A'), control_trans=('C', 'T'),
            pvalue_cutoff=0.05, prob_diff=0.0):
        """This lovely bunch of coconuts does all of the grunt work, running through the alignment to find
        hypermutation on a gross and by_seq basis.
        """
        aln_info = AlignInfo.SummaryInfo(self)
        self.mut_trans, self.control_trans = mut_trans, control_trans
        ref_seq = self.reference_sequence if self.reference_sequence else aln_info.dumb_consensus(threshold=consensus_threshold)

        def findall(ref_seq, residue):
            return [i for i in range(0, len(ref_seq)) if ref_seq[i] == residue]

        def context(i):
            try:
                return str(ref_seq[i:i+2])
            except IndexError:
                return str(ref_seq[i:i+1])

        ref_seq_indices = {}
        for residue in HyperfreqAlignment.RESIDUES:
            ref_seq_indices[residue] = findall(ref_seq, residue)

        for seq in self:
            seq.mut_indices = {}
            for cons_res in ref_seq_indices.keys():
                for seq_res in HyperfreqAlignment.RESIDUES:
                    seq.mut_indices[(cons_res, seq_res)] = [i for i in ref_seq_indices[cons_res] if
                            seq[i] == seq_res]

            seq.n_muts = len(seq.mut_indices[mut_trans])
            seq.n_controls = len(seq.mut_indices[control_trans])

            seq.n_mut_ctxt = len(ref_seq_indices[mut_trans[0]])
            seq.n_control_ctxt = len(ref_seq_indices[control_trans[0]])

            seq.pvalue = 1.0 - sekhon.test(
                seq.n_muts, seq.n_mut_ctxt - seq.n_muts,
                seq.n_controls, seq.n_control_ctxt - seq.n_controls,
                prob_diff=prob_diff)
            seq.hm_pos = seq.pvalue <= pvalue_cutoff

            seq.contexts = {}
            for i in seq.mut_indices[mut_trans]:
                try:
                    seq.contexts[context(i)] += 1
                except KeyError:
                    seq.contexts[context(i)] = 1

        self.mut_seqs = [s for s in self if s.hm_pos]
        self.mut_aln = Align.MultipleSeqAlignment(self.mut_seqs)
        self.mut_indices = list(set([i for s in self.mut_seqs for i in s.mut_indices[mut_trans]]))
        self.mut_indices.sort()
        # That is, the 1-based index positions
        self.mut_columns = [i+1 for i in self.mut_indices]
        self.mut_contexts = [context(i) for i in self.mut_indices]
        self.muts_per_site = [self.mut_aln[:,i].count(mut_trans[1]) for i in self.mut_indices]

        return self


    def split_hypermuts(self, hm_columns=None):
        '''Produce the hypermut positive and hypermut negative alignments'''
        if hm_columns:
            hm_indices = list(set(map(lambda n: n - 1, hm_columns)))
            hm_indices.sort()
        else:
            # if mut_indices is not specified, use the analysis results
            hm_indices = self.mut_indices

        # soi is either a seq or index - handle appropriately
        def hyp_reducer(soi, i):
            seq1 = self[:,soi:soi+1] if type(soi) == int else soi
            seq2 = self[:,i:i+1]
            return seq1 + seq2
        self.hm_pos_aln = reduce(hyp_reducer, hm_indices)

        self.hm_neg_aln = self[:,:hm_indices[0]]
        n_hypermut = len(hm_indices)
        for i in range(0,n_hypermut-1):
            start_i = hm_indices[i] + 1
            stop_i = hm_indices[i+1]
            self.hm_neg_aln += self[:,start_i:stop_i]

        self.hm_neg_aln += self[:,hm_indices[-1]:]

        return self

    def write_analysis(self, gross_handle=None, gross_writer=None, by_seq_handle=None, by_seq_writer=None,
            cluster=None, header=True):
        """This method can be either used by itself or by the HyperfreqAlign.Set class instances to write for
        many sets. As such, it is possible either to pass in a gross_handle file for csv writer, or a
        csvwriter object which already has column names written to it."""
        if not gross_writer:
            gross_writer = csv.writer(gross_handle)

        if not by_seq_writer:
            by_seq_writer = csv.writer(by_seq_handle)

        if header:
            sorted_contexts = list(set(self.contexts))
            sorted_contexts.sort(cmp=HyperfreqAlignment.context_sorter)
            gross_writer.writerow(['cluster','column','context','count'])

            by_seq_writer.writerow(['sequence', 'cluster', 'pvalue', 'hm_pos',
                    'n_muts', 'n_controls', 'n_mut_ctxt', 'n_control_ctxt'] +
                    [HyperfreqAlignment.trans_col_name(trans) for trans in HyperfreqAlignment.MUT_PATTERNS] +
                    sorted_contexts)

        for i in xrange(0, len(self.mut_indices)):
            row = [cluster, self.mut_columns[i], self.mut_contexts[i], self.muts_per_site[i]]
            gross_writer.writerow(row)

        for seq in self:
            row = [seq.name, cluster, seq.pvalue, seq.hm_pos, seq.n_muts, seq.n_controls, seq.n_mut_ctxt,
                    seq.n_control_ctxt]
            row += [len(seq.mut_indices[trans]) for trans in HyperfreqAlignment.MUT_PATTERNS]
            def get_ctxt(c):
                try:
                    return seq.contexts[c]
                except KeyError:
                    return 0

            row += [get_ctxt(ctxt) for ctxt in sorted_contexts]
            by_seq_writer.writerow(row)


    class Set:
        """This class organizes the logic of dealing with clustered alignments in a cohesive fashion, so that
        each cluster alignmnet can be evaluated with respect to it's own consensus, and the results compiled
        for the entire alignment."""
        def __init__(self, seq_records, cluster_map=None, clusters=None, reference_sequences=None):
            self.seq_records = seq_records
            self.cluster_map = cluster_map if cluster_map else {'all': seq_records.keys()}
            self.clusters = clusters if clusters else self.cluster_map.keys()
            self.cluster_alns = {}
            for cluster in self.clusters:
                cluster_seqs = [self.seq_records[x] for x in self.cluster_map[cluster]]
                missing_ref_seqs = []
                try:
                    ref_seq = reference_sequences[cluster] if reference_sequences else None
                except KeyError:
                    missing_ref_seqs.append(cluster)
                    ref_seq = None
                self.cluster_alns[cluster] = HyperfreqAlignment(cluster_seqs, reference_sequence=ref_seq)
            if len(missing_ref_seqs) > 0:
                warnings.warn("Clusters missing representatives: " + ', '.join(missing_ref_seqs))
            self.contexts = set()


        def analyze_hypermuts(self, consensus_threshold=0.5, mut_trans=('G','A'), control_trans=('C', 'T'),
            pvalue_cutoff=0.05, prob_diff=0.0):
            """Run the analysis for each cluster's HyperfreqAlignment."""
            for aln in self.cluster_alns.values():
                aln.analyze_hypermuts(consensus_threshold,
                        mut_trans=mut_trans,
                        control_trans=control_trans,
                        prob_diff=prob_diff)
                self.contexts.update(aln.mut_contexts)


        def write_analysis(self, gross_handle, by_seq_handle):
            """Once analyzed, write the results to files (once global, or by site, the other by sequence)."""
            sorted_contexts = list(self.contexts)
            sorted_contexts.sort(cmp=HyperfreqAlignment.context_sorter)
            gross_writer = csv.writer(gross_handle)
            by_seq_writer = csv.writer(by_seq_handle)

            gross_writer.writerow(['cluster','column','context','count'])

            by_seq_writer.writerow(['sequence', 'cluster', 'pvalue', 'hm_pos',
                    'n_muts', 'n_controls', 'n_mut_ctxt', 'n_control_ctxt'] +
                    [HyperfreqAlignment.trans_col_name(trans) for trans in HyperfreqAlignment.MUT_PATTERNS] +
                    sorted_contexts)

            for cluster in self.cluster_alns:
                aln = self.cluster_alns[cluster]
                aln.write_analysis(gross_writer=gross_writer, by_seq_writer=by_seq_writer, cluster=cluster,
                        header=False)


