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
    def mut_col_name(mut_tuple):
        return "{}_to_{}".format(*mut_tuple)


    def __res_ratio__(self, residue, index):
        return float(self[:,index].count(residue)) / float(self.__len__())

    def __res_sites__(self, residue, consensus_threshold=0.5):
        aln_length = self.get_alignment_length()
        return [i for i in xrange(0, aln_length) if
                self.__res_ratio__(residue,i) >= consensus_threshold]

    # XXX
    # Hmm... not even sure this was getting used anymore...
    #def __mutation_sites__(self, orig_residue, mut_residue,
            #consensus_threshold=0.5, trans_threshold=1):
        #possible_sites = self.__res_sites__(orig_residue, consensus_threshold)
        #return [i for i in possible_sites if
                #self[:,i].count(mut_residue) >= trans_threshold]


    def __consensus__(self, consensus_threshold):
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


    # XXX - may want to leave this here for doing stats
    def single_nt_mut_analysis(self):
        def findall(ref_seq, residue):
            return [i for i in range(0, len(self.reference_sequence)) if
                    self.reference_sequence[i] == residue]

        ref_seq_indices = {}
        for residue in HyperfreqAlignment.RESIDUES:
            ref_seq_indices[residue] = findall(self.reference_sequence, residue)

        for seq in self:
            seq.mut_indices = {}
            for cons_res in ref_seq_indices.keys():
                for seq_res in HyperfreqAlignment.RESIDUES:
                    seq.mut_indices[(cons_res, seq_res)] = [i for i in ref_seq_indices[cons_res] if
                            seq[i] == seq_res]

            for i in seq.mut_indices[self.focus_pattern.mutation]:
                try:
                    seq.contexts[self.context(i)] += 1
                except KeyError:
                    seq.contexts[self.context(i)] = 1


    def analyze_hypermuts(self, focus_pattern, control_pattern, consensus_threshold=None,
            pvalue_cutoff=0.05, prob_diff=0.0):
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

        for seq in self:
            seq.focus_pos_indices = self.focus_pattern.pos_indices(seq, focus_ref_indices)
            seq.focus_neg_indices = self.focus_pattern.neg_indices(seq, focus_ref_indices)
            seq.control_pos_indices = self.control_pattern.pos_indices(seq, control_ref_indices)
            seq.control_neg_indices = self.control_pattern.neg_indices(seq, control_ref_indices)

            seq.n_focus_pos = len(seq.focus_pos_indices)
            seq.n_control_pos = len(seq.control_pos_indices)
            seq.n_focus_neg = len(seq.focus_neg_indices)
            seq.n_control_neg = len(seq.control_neg_indices)

            seq.pvalue = 1.0 - sekhon.test(seq.n_focus_pos, seq.n_control_pos, seq.n_focus_neg, seq.n_control_neg,
                    prob_diff=prob_diff)
            seq.hm_pos = seq.pvalue <= pvalue_cutoff

            seq.contexts = {}


        self.hm_pos_seqs = [s for s in self if s.hm_pos]
        self.hm_pos_aln = Align.MultipleSeqAlignment(self.hm_pos_seqs)
        self.hm_pos_indices = list(set([i for s in self.hm_pos_seqs for i in s.focus_pos_indices]))
        self.hm_pos_indices.sort()
        # That is, the 1-based index positions
        self.mut_columns = [i+1 for i in self.hm_pos_indices]
        # XXX - need to replace this by a more flexible site probing system
        self.mut_contexts = [self.context(i) for i in self.hm_pos_indices]
        # XXX - don't seem to have actually been using this. Probably just wanted it for stats. May throw back
        # in later
        #self.muts_per_site = [self.mut_aln[:,i].count(mut_trans[1]) for i in self.mut_indices]

        return self


    def split_hypermuts(self, hm_columns=None):
        '''Produce the hypermut positive and hypermut negative alignments'''
        if hm_columns:
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
        self.hm_pos_aln = reduce(hyp_reducer, hm_indices)

        self.hm_neg_aln = self[:,:hm_indices[0]]
        n_hypermut = len(hm_indices)
        for i in range(0, n_hypermut - 1):
            start_i = hm_indices[i] + 1
            stop_i = hm_indices[i+1]
            self.hm_neg_aln += self[:,start_i:stop_i]

        self.hm_neg_aln += self[:,hm_indices[-1]+1:]

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

            by_seq_writer.writerow(['sequence', 'cluster', 'pvalue', 'hm_pos',
                    'n_focus_pos', 'n_controls_pos', 'n_focus_neg', 'n_control_neg'] +
                    # XXX - throw in a flag for putting this stuff in if so desired
                    #[HyperfreqAlignment.mut_col_name(trans) for trans in HyperfreqAlignment.MUT_PATTERNS] +
                    sorted_contexts)

        else:
            sorted_contexts = contexts

        for seq in self:
            # Gross writer now does by seq
            if seq.hm_pos:
                # XXX - how is this not breaking? We'll have to see what output is actually atm
                for i in seq.mut_indices[self.mut_trans]:
                    row = [cluster, seq.name, i+1, self.context(i)]
                    gross_writer.writerow(row)

            row = [seq.name, cluster, seq.pvalue, seq.hm_pos, seq.n_focus_pos, seq.n_control_pos,
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
                pvalue_cutoff=0.05, prob_diff=0.0):
            """Run the analysis for each cluster's HyperfreqAlignment. It is possible to specify the mutation
            transition, the control transition here, and well as the probability difference and pvalue cutoff
            which for the decision procedure.  Note that if you wish to change the consensus_threshold used to
            instantiate the Set (or override the reference_sequences), that can be done here."""
            for aln in self.cluster_alns.values():
                aln.analyze_hypermuts(focus_pattern, control_pattern, consensus_threshold, prob_diff=prob_diff)
                # XXX - shit. Looks like I do have to just go through and fetch these or I'll get any array
                self.contexts.update(aln.mut_contexts)


        def write_analysis(self, gross_handle, by_seq_handle):
            """Once analyzed, write the results to files (once global, or by site, the other by sequence)."""
            sorted_contexts = list(self.contexts)
            sorted_contexts.sort(cmp=HyperfreqAlignment.context_sorter)
            gross_writer = csv.writer(gross_handle)
            by_seq_writer = csv.writer(by_seq_handle)

            gross_writer.writerow(['cluster', 'sequence', 'column', 'context'])

            by_seq_writer.writerow(['sequence', 'cluster', 'pvalue', 'hm_pos',
                    'n_focus_pos', 'n_control_pos', 'n_focus_neg', 'n_control_neg'] +
                    #[HyperfreqAlignment.mut_col_name(mut) for mut in HyperfreqAlignment.MUT_PATTERNS] +
                    sorted_contexts)

            for cluster in self.cluster_alns:
                aln = self.cluster_alns[cluster]
                aln.write_analysis(gross_writer=gross_writer, by_seq_writer=by_seq_writer, cluster=cluster,
                        header=False, contexts=sorted_contexts)


