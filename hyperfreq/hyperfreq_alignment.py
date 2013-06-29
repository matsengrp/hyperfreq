from Bio import Align
from Bio import SeqRecord
from Bio import Seq
from Bio.Align import AlignInfo
from beta_rat import BetaRat
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


    def single_nt_mut_analysis(self):
        def findall(ref_seq, residue):
            return [i for i in range(0, len(self.reference_sequence)) if
                    self.reference_sequence[i] == residue]

        ref_seq_indices = {}
        #if type(self.reference_sequence) == 
        for residue in HyperfreqAlignment.RESIDUES:
            ref_seq_indices[residue] = findall(self.reference_sequence, residue)

        for seq in self:
            seq.mut_indices = {}
            for cons_res in ref_seq_indices.keys():
                for seq_res in HyperfreqAlignment.RESIDUES:
                    # XXX - interesting. this could be a good place to try and catch people tring to use a
                    # sequence of the wrong length for a reference sequence. Leads to IndexError in seq[i] if
                    # ref is longer than it should be
                    seq.mut_indices[(cons_res, seq_res)] = [i for i in ref_seq_indices[cons_res] if
                            seq[i] == seq_res]

            seq.contexts = {}
            for i in seq.mut_indices[self.focus_pattern.mutation]:
                try:
                    seq.contexts[self.context(i)] += 1
                except KeyError:
                    seq.contexts[self.context(i)] = 1

    def analyze_hypermuts(self, focus_pattern, control_pattern, consensus_threshold=None,
            br_left_cutoff=1.8, significance_level=0.05, prior=(0.5, 1.0), pos_quants_only=False):
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

        self.focus_pattern, self.control_pattern = focus_pattern, control_pattern

        focus_ref_indices = self.focus_pattern.ref_match_indices(self.reference_sequence)
        control_ref_indices = self.control_pattern.ref_match_indices(self.reference_sequence)
        
        # XXX - Needed for doing the context; be aware that taking this out would lead to an error later where
        # seq.context is comiled for other analysis. Could fix this at some point if we want to make this
        # optional
        self.single_nt_mut_analysis()

        for seq in self:
            seq.focus_pos_indices = self.focus_pattern.pos_indices(seq, focus_ref_indices)
            seq.focus_neg_indices = self.focus_pattern.neg_indices(seq, focus_ref_indices)
            seq.control_pos_indices = self.control_pattern.pos_indices(seq, control_ref_indices)
            seq.control_neg_indices = self.control_pattern.neg_indices(seq, control_ref_indices)

            seq.n_focus_pos = len(seq.focus_pos_indices)
            seq.n_control_pos = len(seq.control_pos_indices)
            seq.n_focus_neg = len(seq.focus_neg_indices)
            seq.n_control_neg = len(seq.control_neg_indices)

            counts = [seq.n_focus_pos, seq.n_control_pos, seq.n_focus_neg, seq.n_control_neg]
            seq.beta_rat = BetaRat(*counts, prior=prior)

            if VERBOSE:
                t = time()
                print seq.name, seq.beta_rat,

            seq.hm_pos = seq.beta_rat.cdf(br_left_cutoff) < significance_level
            # If this flag is set, we only want to compute the quantiles if the sequence is gonna be positive
            if pos_quants_only and not seq.hm_pos:
                seq.br_left = None
                seq.br_median = None
            else:
                seq.br_left = seq.beta_rat.ppf(0.05)
                seq.br_median = seq.beta_rat.ppf(0.5)

            seq.br_max = seq.beta_rat.pdf_max()

            if VERBOSE:
                print "Time:", time() - t

            seq.fisher_pvalue = fisher.pvalue(*counts).right_tail


        self.hm_pos_seqs = [s for s in self if s.hm_pos]
        self.hm_pos_aln = Align.MultipleSeqAlignment(self.hm_pos_seqs)
        self.hm_pos_indices = list(set([i for s in self.hm_pos_seqs for i in s.focus_pos_indices]))
        self.hm_pos_indices.sort()
        # That is, the 1-based index positions
        self.mut_columns = [i+1 for i in self.hm_pos_indices]
        self.mut_contexts = [self.context(i) for i in self.hm_pos_indices]
        # XXX - don't seem to have actually been using this. Probably just wanted it for stats. May throw back
        # in later
        #self.muts_per_site = [self.mut_aln[:,i].count(self.focus_pattern.mutation[1]) for i in self.hm_pos_indices]

        return self


    def split_hypermuts(self, hm_columns=None):
        '''Produce the hypermut positive and hypermut negative alignments'''
        # Come one python... Is [] really false?
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
            sorted_contexts = list(set(self.contexts)) if self.contexts else []
            sorted_contexts.sort(cmp=HyperfreqAlignment.context_sorter)
            gross_writer.writerow(['cluster', 'sequence', 'column','context'])

            by_seq_writer.writerow(HyperfreqAlignment.BASE_ROWNAMES + sorted_contexts)
                    #[HyperfreqAlignment.mut_col_name(trans) for trans in HyperfreqAlignment.MUT_PATTERNS] +

        else:
            sorted_contexts = contexts

        for seq in self:
            # Gross writer now does by seq
            if seq.hm_pos:
                for i in seq.focus_pos_indices:
                    row = [cluster, seq.name, i+1, self.context(i)]
                    gross_writer.writerow(row)

            row = [seq.name, cluster, seq.br_left, seq.br_median, seq.br_max, seq.fisher_pvalue, seq.hm_pos, seq.n_focus_pos, seq.n_control_pos,
                    seq.n_focus_neg, seq.n_control_neg]
            # XXX - again, flaggify
            #row += [len(seq.mut_indices[trans]) for trans in HyperfreqAlignment.MUT_PATTERNS]
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


        def analyze_hypermuts(self, focus_pattern, control_pattern, consensus_threshold=None,
                br_left_cutoff=2.0, prior=(0.5, 1.0), pos_quants_only=False):
            # XXX - Update doc
            """Run the analysis for each cluster's HyperfreqAlignment. It is possible to specify the mutation
            transition, the control transition here, and well as the probability difference and pvalue cutoff
            which for the decision procedure.  Note that if you wish to change the consensus_threshold used to
            instantiate the Set (or override the reference_sequences), that can be done here."""
            for aln in self.cluster_alns.values():
                aln.analyze_hypermuts(focus_pattern, control_pattern, consensus_threshold,
                        br_left_cutoff=br_left_cutoff, prior=prior, pos_quants_only=pos_quants_only)
                # XXX - Should come up with something smart here in case we don't compute the context
                self.contexts.update(aln.mut_contexts)


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


