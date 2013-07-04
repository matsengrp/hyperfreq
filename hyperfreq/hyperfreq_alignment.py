from Bio import Align, SeqRecord, Seq
from Bio.Align import AlignInfo
from betarat import BetaRat
from time import time

import warnings
import csv
import fisher

VERBOSE = False


class HyperfreqAlignment(Align.MultipleSeqAlignment):

    RESIDUES = ['A', 'C', 'T', 'G']

    MUT_PATTERNS = [(x,y) for x in RESIDUES for y in RESIDUES]
    MUT_PATTERNS.sort()

    # For keeping output code clean. May need to move this to an abstraction layer to make more flexible.
    BASE_ROWNAMES = ['sequence', 'cluster', 'br_left', 'br_median', 'br_max', 'fisher_pvalue', 'hm_pos',
                    'n_focus_pos', 'n_control_pos', 'n_focus_neg', 'n_control_neg'] 

    @staticmethod
    def context_sorter(c1, c2):
        """ Sorts contexts first by length, then by actual string content. """
        l1, l2 = len(c1), len(c2)
        if l1 != l2:
            return l1 - l2
        else:
            return 2 * int(c1 > c2) - 1

    @staticmethod
    def mut_col_name(mut_tuple):
        """ Little helper function for making pretty mutation pattern names."""
        return "{}_to_{}".format(*mut_tuple)


    def __res_ratio__(self, residue, index):
        """ Finds the fraction of sequences in the alignment that, at the given index/column, have the
        specified residue. """
        return float(self[:,index].count(residue)) / float(self.__len__())

    def __res_sites__(self, residue, consensus_threshold=0.5):
        """ Finds the indices of sites in the alignment where the specified residue is at least as abundant as
        the specified consensus_threshold."""
        aln_length = self.get_alignment_length()
        return [i for i in xrange(0, aln_length) if
                self.__res_ratio__(residue,i) >= consensus_threshold]

    def __consensus__(self, consensus_threshold):
        """Returns a Bio.Align.AlignInfo dumb_consensus of this alignment at the given threshold. """
        aln_info = AlignInfo.SummaryInfo(self)
        return aln_info.dumb_consensus(threshold=consensus_threshold)


    def __init__(self, *args, **kwargs):
        """This inherits from biopythons MultipleSeqAlignment, and does basically the same stuff, but also
        lets one specify a reference sequence (Bio.Seq type) for comparison istead of a consensus."""

        def strip_option(option_name, default=None):
            try:
                value = kwargs[option_name]
                del kwargs[option_name]
            except KeyError:
                value = None
            return value

        self.consensus_threshold = strip_option('consensus_threshold', 0.5)
        ref_seq = strip_option('reference_sequence')
        super(HyperfreqAlignment, self).__init__(*args, **kwargs)
        self.reference_sequence = ref_seq if ref_seq else self.__consensus__(self.consensus_threshold)


    def context(self, i):
        """Context for a given index i (0-based)"""
        try:
            return str(self.reference_sequence[i:i+2])
        except IndexError:
            return str(self.reference_sequence[i:i+1])


    def analyze(self, mutation_pattern, consensus_threshold=None,
            br_left_cutoff=1.8, significance_level=0.05, prior=(0.5, 1.0), pos_quants_only=True, quants=None):
        """This is where all of the grunt work happens; running through the alignment to find
        hypermutation on a gross and by_seq basis. Consensus threshold defaults to that of the hyperfreq
        alingment initialization (by passing None) but can be overridden herek.
        """

        # Want to make it possible to override the consensus_threshold established on initialization. Making
        # sure that it's created on initialization first though, makes things safter, as you know there is
        # always a reference seq
        if consensus_threshold and consensus_threshold != self.consensus_threshold:
            self.consensus_threshold = consensus_threshold
            self.reference_sequence = self.__consensus__(consensus_threshold)

        focus_pattern, control_pattern = mutation_pattern

        focus_ref_indices = focus_pattern.ref_match_indices(self.reference_sequence)
        control_ref_indices = control_pattern.ref_match_indices(self.reference_sequence)

        for seq in self:
            hm_data = dict()

            focus_pos_indices = focus_pattern.pos_indices(seq, focus_ref_indices)

            focus_pos = len(focus_pos_indices)
            focus_neg = len(focus_pattern.neg_indices(seq, focus_ref_indices))
            control_pos = len(control_pattern.pos_indices(seq, control_ref_indices))
            control_neg = len(control_pattern.neg_indices(seq, control_ref_indices))

            counts = [focus_pos, control_pos, focus_neg, control_neg]
            beta_rat = BetaRat(*counts, prior=prior)

            if VERBOSE:
                t = time()
                print "On sequence", seq.name, beta_rat,
            
            # XXX - need to make name cdf1 reflect chosen br_left_cutoff (or whatever we end up callig it)
            cdf1 = beta_rat.cdf(br_left_cutoff)
            hm_pos = cdf1 < significance_level
            br_map = beta_rat.map()

            fisher_pvalue = fisher.pvalue(*counts).right_tail

            hm_data = dict(
                    sequence=seq.name,
                    hm_pos=hm_pos,
                    cdf1=cdf1,
                    map=br_map,
                    fisher=fisher_pvalue,
                    focus_pos=focus_pos,
                    focus_neg=focus_neg,
                    control_pos=control_pos,
                    control_neg=control_neg)

            # If this flag is set, we only want to compute the quantiles if the sequence is gonna be positive
            if hm_pos or not pos_quants_only:
                for quant in quants:
                    hm_data['q_{}'.format(quant)] = beta_rat.ppf(quant)

            #for cdf in cdfs:
                #hm_data['cdf_{}'.format(quant)] = beta_rat.cdf(quant)

            if VERBOSE:
                print "Time:", time() - t

            yield hm_data



    def split_hypermuts(self, hm_columns=None):
        '''Produce the hypermut positive and hypermut negative alignments'''
        if hm_columns or hm_columns == []:
            hm_indices = list(set(map(lambda n: n - 1, hm_columns)))
            hm_indices.sort()
        else:
            # if hm_columns is not specified, use the analysis results
            hm_indices = self.hm_pos_indices

        # soi is either a seq or index - handle appropriately
        def hyp_reducer(soi, i):
            seq1 = self[:,soi:soi+1] if type(soi) == int else soi
            seq2 = self[:,i:i+1]
            return seq1 + seq2

        init = type(self)([SeqRecord.SeqRecord(Seq.Seq(''), id=self[i].id) for i in xrange(len(self))])
        self.hm_pos_aln = reduce(hyp_reducer, hm_indices, init)

        if hm_indices:
            self.hm_neg_aln = self[:,:hm_indices[0]]
            n_hypermut = len(hm_indices)
            for i in range(0, n_hypermut - 1):
                start_i = hm_indices[i] + 1
                stop_i = hm_indices[i+1]
                self.hm_neg_aln += self[:,start_i:stop_i]

            self.hm_neg_aln += self[:,hm_indices[-1]+1:]

        else:
            self.hm_neg_aln = self

        return self


    def write_analysis(self,
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
            gross_writer.writerow(['cluster', 'sequence', 'column','context'])
            # XXX - make sure we fix column names up to be dynamic for quants, etc
            by_seq_writer.writerow(HyperfreqAlignment.BASE_ROWNAMES)

        for seq in self:
            # Gross writer now does by seq
            if seq.hm_pos:
                # XXX - going to have to make sure that focus_pos_indices can come from somewhere here
                for i in seq.focus_pos_indices:
                    row = [cluster, seq.name, i+1, self.context(i)]
                    gross_writer.writerow(row)

            row = [seq.name, cluster, seq.br_left, seq.br_median, seq.br_max, seq.fisher_pvalue, seq.hm_pos, seq.n_focus_pos, seq.n_control_pos,
                    seq.n_focus_neg, seq.n_control_neg]

            by_seq_writer.writerow(row)



    class Set:
        """This class organizes the logic of dealing with clustered alignments in a cohesive fashion, so that
        each cluster alignmnet can be evaluated with respect to it's own consensus, and the results compiled
        for the entire alignment."""
        def __init__(self, seq_records, cluster_map=None, clusters=None, reference_sequences=None,
                consensus_threshold=0.5):
            """Options one can specify here are reference sequence and consensus_threshold, both of which are
            passed along to HyperfreqAlignment instantiations. Note that one can override the consensus
            threshold (or reference_seqs) passed here by passing in a consensus_threshold to analyze
            hypermuts."""
            self.seq_records = seq_records
            self.cluster_map = cluster_map if cluster_map else {'all': seq_records.keys()}
            self.clusters = clusters if clusters else self.cluster_map.keys()
            self.cluster_alns = {}
            for cluster in self.clusters:
                cluster_seqs = [self.seq_records[x] for x in self.cluster_map[cluster]]
                missing_ref_seqs = []
                try:
                    ref_seq = reference_sequences[cluster].seq if reference_sequences else None
                except KeyError:
                    missing_ref_seqs.append(cluster)
                    ref_seq = None
                self.cluster_alns[cluster] = HyperfreqAlignment(cluster_seqs, reference_sequence=ref_seq,
                        consensus_threshold=consensus_threshold)
            if len(missing_ref_seqs) > 0:
                warnings.warn("Clusters missing representatives: " + ', '.join(missing_ref_seqs))
            self.contexts = set()


        def analyze(self, focus_pattern, control_pattern, consensus_threshold=None,
                br_left_cutoff=2.0, prior=(0.5, 1.0), pos_quants_only=True):
            # XXX - Update doc
            """Run the analysis for each cluster's HyperfreqAlignment. It is possible to specify the mutation
            transition, the control transition here, and well as the probability difference and pvalue cutoff
            which for the decision procedure.  Note that if you wish to change the consensus_threshold used to
            instantiate the Set (or override the reference_sequences), that can be done here."""
            for aln in self.cluster_alns.values():
                aln.analyze(focus_pattern, control_pattern, consensus_threshold,
                        br_left_cutoff=br_left_cutoff, prior=prior, pos_quants_only=pos_quants_only)


        def write_analysis(self, gross_handle, by_seq_handle):
            """Once analyzed, write the results to files (once global, or by site, the other by sequence)."""
            sorted_contexts = list(self.contexts)
            sorted_contexts.sort(cmp=HyperfreqAlignment.context_sorter)
            gross_writer = csv.writer(gross_handle)
            by_seq_writer = csv.writer(by_seq_handle)

            gross_writer.writerow(['cluster', 'sequence', 'column', 'context'])

            by_seq_writer.writerow(HyperfreqAlignment.BASE_ROWNAMES + sorted_contexts)
                    #[HyperfreqAlignment.mut_col_name(mut) for mut in HyperfreqAlignment.MUT_PATTERNS] +

            for cluster in self.cluster_alns:
                aln = self.cluster_alns[cluster]
                aln.write_analysis(gross_writer=gross_writer, by_seq_writer=by_seq_writer, cluster=cluster,
                        header=False, contexts=sorted_contexts)


