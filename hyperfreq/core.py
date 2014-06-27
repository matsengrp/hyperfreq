from Bio import Align, SeqRecord, Seq
from Bio.Align import AlignInfo
import betarat
from betarat import BetaRat
from time import time
import logging
import itertools
import copy
import re
import warnings
import fisher

VERBOSE = False


"Analysis defaults, across entire code base (including cli)"
analysis_defaults = dict(
        rpr_cutoff=1.0,
        significance_level=0.05,
        prior=(0.5, 1.0),
        pos_quants_only=True,
        quants=[],
        caller="map",
        cdfs=[])
analysis_defaults.update(betarat.defaults)


def apply_analysis_defaults(func):
    "Decorator for applying analysis defaults to various analysis functions"
    def decorated(*args, **kw_args):
        new_kw_args = analysis_defaults.copy()
        new_kw_args.update(kw_args)
        return func(*args, **new_kw_args)
    decorated.func_name = func.func_name
    return decorated


def make_pattern_call(analyses, caller):
    """Takes a list of analysis results (as dicts), presumably each for a different context, and returns the
    list index and entry identified as having the strongest hypermutation signal. Actual evaluation depends on
    1) whether the sequence is identified in a given context as actually being positive, and 2) which context
    has the more extreme value of the `caller` statistic. Larger values of the caller statistic will be seen
    as more extreme, except for `caller` values matching "cutoff_cdf" or "cdf_*", for which smaller values are
    considered more extreme."""
    # This takes care of the reverse sorting for CDF values
    sign = -1 if caller == 'cutoff_cdf' or re.match('cdf_', caller) else 1
    def sort_tuple(result):
        try:
            caller_result = result[caller]
        except KeyError:
            raise ValueError, """The call statistic {caller} does not exist for some sequence. If you are
            using a quantile (such as q_0.05) for calling, try specifying -Q instead of -q from the command
            line (or positive_quants_only=False via the API); otherwise quantiles are only computed for
            sequences evaluated as positive.""".format(caller=caller)
        return (result['hm_pos'], sign * caller_result)
    # Where all the work is - a sort that prioritizes hm_pos, then the caller statistic value
    _, index = max((sort_tuple(result), i) for i, result in enumerate(analyses))
    return index, analyses[index]



class Alignment(Align.MultipleSeqAlignment):

    def __init__(self, *args, **kwargs):
        """This inherits from biopythons MultipleSeqAlignment, and does basically the same stuff, but also
        lets one either specify a reference sequence (Bio.Seq type) or a consensus threshold for computing the
        consensus of the alignment as the reference sequence. Only one of these options should be specified.
        If nothing is specified, the consensus is taken with a threshold of None (See Bio.Align.AlignInfo's
        dumb_consensus function for details)."""

        def strip_option(option_name, default=None):
            try:
                value = kwargs[option_name]
                del kwargs[option_name]
            except KeyError:
                value = None
            return value

        consensus_threshold = strip_option('consensus_threshold')
        ref_seq = strip_option('reference_sequence')
        super(Alignment, self).__init__(*args, **kwargs)
        if consensus_threshold and ref_seq:
            raise ValueError, "Only one of reference_sequence or consensus_threshold should be specified"
        # reference sequence should be manually re-specified if this need to be changed
        self.reference_sequence = ref_seq if ref_seq else self.__consensus__(consensus_threshold)


    def __consensus__(self, consensus_threshold):
        """Returns a Bio.Align.AlignInfo dumb_consensus of this alignment at the given threshold. """
        aln_info = AlignInfo.SummaryInfo(self)
        return aln_info.dumb_consensus(threshold=consensus_threshold)


    def context(self, i):
        """Context for a given index i (0-based)"""
        try:
            return str(self.reference_sequence[i:i+2])
        except IndexError:
            return str(self.reference_sequence[i:i+1])


    @apply_analysis_defaults
    def multiple_context_analysis(self, mutation_patterns, **kw_args):
        """Run hypermutation analysis for multiple contexts and cull together all of the data, including
        calling of the pattern with the highest evidence of hypermutation. kw_arg `caller` specifies the
        metric use in making the call pattern decision."""
        # The implementation here is a little weird. We go through and create a list tuples of pattern and analysis generator
        # We do this because we need to iterative one by one through eadh seqeunce in each of the context analyses
        analyzers = (self.analyze(x, **kw_args) for x in mutation_patterns)
        caller = kw_args['caller']
        for seq_analyses in itertools.izip(*analyzers):
            # Iterate through each sequence and get a tuple of results, one for each context, as `seq_analyses`.
            # Then go through and figure out the call information and return a map from which we can fetch all
            # data for each sequence easily
            call_index, call_analysis = make_pattern_call(seq_analyses, caller)
            call_pattern = mutation_patterns[call_index]
            call_data = copy.copy(call_analysis)
            call_data['call_pattern'] = call_pattern.name
            sequence_results = dict((pattern.name, seq_analyses[i]) for i, pattern in enumerate(mutation_patterns))
            sequence_results['call'] = call_data
            yield sequence_results


    @apply_analysis_defaults
    def analyze(self, mutation_pattern, **kw_args):
        """Returns an iterator of hypermutation analysis results for the given mutation pattern."""
        focus_pattern, control_pattern = mutation_pattern

        focus_indices = focus_pattern.ref_match_indices(self.reference_sequence)
        control_indices = control_pattern.ref_match_indices(self.reference_sequence)

        for seq in self:
            yield self.analyze_sequence(seq, mutation_pattern, focus_indices, control_indices, **kw_args)


    @apply_analysis_defaults
    def analyze_sequence(self, seq, mutation_pattern, focus_indices=None, control_indices=None, **kw_args):
        """Hypermutation analysis of a specific sequence. Arguments focus_indices and control_indices make it
        possible to compute these values once, outside of this function, and reused for multiple evaluations."""
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
        logging.info(" On sequence: " + seq.name + "; with counts: " + str(counts))

        # Start running stats
        cutoff_cdf = beta_rat.cdf(kw_args['rpr_cutoff'], quadr_maxiter=kw_args['quadr_maxiter'])
        hm_pos = cutoff_cdf < kw_args['significance_level']
        br_map = beta_rat.map()
        br_ltmap = beta_rat.lt_map()

        fisher_pvalue = fisher.pvalue(*counts).right_tail

        mut_columns = [i + 1 for i in focus_pos_indices] if hm_pos else []
        mut_contexts = [self.context(i) for i in focus_pos_indices] if hm_pos else []

        # Construct the dict we'll be returning for this sequence and pattern (will need to add cdfs and so on
        # to it)
        hm_data = dict(
                sequence=seq.name,
                hm_pos=hm_pos,
                cutoff_cdf=cutoff_cdf,
                map=br_map,
                ltmap=br_ltmap,
                fisher_pvalue=fisher_pvalue,
                focus_pos=focus_pos,
                focus_neg=focus_neg,
                control_pos=control_pos,
                control_neg=control_neg,
                mut_columns=mut_columns,
                mut_contexts=mut_contexts)

        # If this flag is set, we only want to compute the quantiles if the sequence is gonna be positive
        for quant in kw_args['quants']:
            key = 'q_{}'.format(quant)
            if hm_pos or not kw_args['pos_quants_only']:
                hm_data[key] = beta_rat.ppf(quant,
                        optim_maxiter=kw_args['optim_maxiter'],
                        quadr_maxiter=kw_args['quadr_maxiter'])
            else:
                hm_data[key] = None

        # Always compute whatever cdfs requested, since they are fairly cheap to process
        for cdf in kw_args['cdfs']:
            hm_data['cdf_{}'.format(cdf)] = beta_rat.cdf(cdf, quadr_maxiter=kw_args['quadr_maxiter'])

        if VERBOSE:
            print "Time:", time() - t

        return hm_data


    def split_hypermuts(self, hm_columns):
        """Produce the hypermut positive and hypermut negative alignments"""
        hm_indices = list(set(map(lambda n: n - 1, hm_columns)))
        hm_indices.sort()

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


class AlignmentSet:
    """This class organizes the logic of dealing with clustered alignments in a cohesive fashion, so that
    each cluster alignmnet can be evaluated with respect to it's own consensus, and the results compiled
    for the entire alignment."""
    def __init__(self, seq_records, cluster_map=None, clusters=None, reference_sequences=None,
            consensus_threshold=None):
        """Only required argument is seq_records, which should be a dictionary mapping sequence names to
        SeqRecord objects. The optional cluster map specifie how these sequences partition into
        clusters; if unspecified, everything is take as one cluster named "all". Other options:

        * clusters: specify which of the clusters to keep. None defaults to all clusters present.
        * [reference_sequences | consensus_threshold]: Only one of these should be specified.
              - reference_sequences: a dict mapping cluster names to sequences.
              - consensus_threshold: a float or None; passed to AlignInfo's dumb_consensus function.
              By default consensus sequences are used with no threshold if both options are unset.
        """
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
            self.cluster_alns[cluster] = Alignment(cluster_seqs, reference_sequence=ref_seq,
                    consensus_threshold=consensus_threshold)
        if len(missing_ref_seqs) > 0:
            warnings.warn("Clusters missing representatives: " + ', '.join(missing_ref_seqs))
        self.contexts = set()


    @apply_analysis_defaults
    def multiple_context_analysis(self, patterns, **kw_args):
        """Run the analysis for each cluster's Alignment. Any keyword arguments passed to this
        function get passed along to the individual function."""
        return itertools.chain(*(aln.multiple_context_analysis(patterns, **kw_args) for aln in
                self.cluster_alns.values()))


