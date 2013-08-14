from Bio import Align, SeqRecord, Seq
from Bio.Align import AlignInfo
from betarat import BetaRat
from time import time

import warnings
import fisher

VERBOSE = False


"Analysis defaults, across entire code base (including cli)"
analysis_defaults = dict(consensus_threshold=None,
        rpr_cutoff=1.0,
        significance_level=0.05,
        prior=(0.5, 1.0),
        pos_quants_only=True,
        quants=[],
        cdfs=[])

def apply_analysis_defaults(func):
    "Decorator for applying analysis defaults to various analysis functions"
    def decorated(*args, **kw_args):
        new_kw_args = analysis_defaults.copy()
        new_kw_args.update(kw_args)
        return func(*args, **new_kw_args)
    decorated.func_name = func.func_name
    return decorated


class HyperfreqAlignment(Align.MultipleSeqAlignment):

    RESIDUES = ['A', 'C', 'T', 'G']

    MUT_PATTERNS = [(x,y) for x in RESIDUES for y in RESIDUES]
    MUT_PATTERNS.sort()

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


    @apply_analysis_defaults
    def analyze(self, mutation_pattern, **kw_args):
        """This is where all of the grunt work happens; running through the alignment to find
        hypermutation on a gross and by_seq basis. Consensus threshold defaults to that of the hyperfreq
        alingment initialization (by passing None) but can be overridden herek.
        """
        # Want to make it possible to override the consensus_threshold established on initialization. Making
        # sure that it's created on initialization first though, makes things safter, as you know there is
        # always a reference seq
        consensus_threshold = kw_args['consensus_threshold']
        if consensus_threshold and consensus_threshold != self.consensus_threshold:
            self.consensus_threshold = consensus_threshold
            self.reference_sequence = self.__consensus__(consensus_threshold)

        focus_pattern, control_pattern = mutation_pattern

        focus_indices = focus_pattern.ref_match_indices(self.reference_sequence)
        control_indices = control_pattern.ref_match_indices(self.reference_sequence)

        for seq in self:
            yield self.analyze_sequence(seq, mutation_pattern, focus_indices, control_indices, **kw_args)


    @apply_analysis_defaults
    def analyze_sequence(self, seq, mutation_pattern, focus_indices=None, control_indices=None, **kw_args):
        """This function does the work of hypermutation analysis on a sequence by sequence basis. Focus and
        control indices can be computed and passed in once for efficiency. """
        focus_pattern, control_pattern = mutation_pattern

        # This makes it possible to not pass in focus indices and control indices if lazy - probably doesn't improve performance
        # that much so may just set to always do this in time.
        focus_indices = focus_indices if focus_indices else focus_pattern.ref_match_indices(self.reference_sequence)
        control_indices = control_indices if control_indices else control_pattern.ref_match_indices(self.reference_sequence)

        # Compute contingency table
        focus_pos_indices = focus_pattern.pos_indices(seq, focus_indices)
        focus_pos = len(focus_pos_indices)
        focus_neg = len(focus_pattern.neg_indices(seq, focus_indices))
        control_pos = len(control_pattern.pos_indices(seq, control_indices))
        control_neg = len(control_pattern.neg_indices(seq, control_indices))

        counts = [focus_pos, control_pos, focus_neg, control_neg]
        beta_rat = BetaRat(*counts, prior=kw_args['prior'])

        if VERBOSE:
            t = time()
            print "On sequence", seq.name, beta_rat

        # Start running stats
        cutoff_cdf = beta_rat.cdf(kw_args['rpr_cutoff'])
        hm_pos = cutoff_cdf < kw_args['significance_level']
        br_map = beta_rat.map()
        br_ltmap = beta_rat.ltmap()

        fisher_pvalue = fisher.pvalue(*counts).right_tail

        mut_columns = [i + 1 for i in focus_pos_indices] if hm_pos else []

        # Construct the dict we'll be returning for this sequence and pattern (will need to add cdfs and so on
        # to it)
        hm_data = dict(
                sequence=seq.name,
                hm_pos=hm_pos,
                cutoff_cdf=cutoff_cdf,
                map=br_map,
                ltmap=br_ltmap,
                fisher=fisher_pvalue,
                focus_pos=focus_pos,
                focus_neg=focus_neg,
                control_pos=control_pos,
                control_neg=control_neg,
                mut_columns=mut_columns)

        # If this flag is set, we only want to compute the quantiles if the sequence is gonna be positive
        for quant in kw_args['quants']:
            key = 'q_{}'.format(quant)
            if hm_pos or not kw_args['pos_quants_only']:
                hm_data[key] = beta_rat.ppf(quant)
            else:
                hm_data[key] = None

        # Always compute whatever cdfs requested, since they are fairly cheap to process
        for cdf in kw_args['cdfs']:
            hm_data['cdf_{}'.format(cdf)] = beta_rat.cdf(cdf)

        if VERBOSE:
            print "Time:", time() - t

        return hm_data


    def split_hypermuts(self, hm_columns):
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


        @apply_analysis_defaults
        def analyze(self, pattern, **kw_args):
            # XXX - Update doc
            """Run the analysis for each cluster's HyperfreqAlignment. It is possible to specify the mutation
            transition, the control transition here, and well as the probability difference and pvalue cutoff
            which for the decision procedure.  Note that if you wish to change the consensus_threshold used to
            instantiate the Set (or override the reference_sequences), that can be done here."""
            for aln in self.cluster_alns.values():
                aln.analyze(pattern, **kw_args)


