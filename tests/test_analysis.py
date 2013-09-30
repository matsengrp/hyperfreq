import unittest
import helpers
from Bio import Seq
from hyperfreq.hyperfreq_alignment import HyperfreqAlignment
from hyperfreq.cluster import load_cluster_map
from hyperfreq import mut_pattern

old_focus_pattern = mut_pattern.MutPattern(('G','A'), '')
old_control_pattern = mut_pattern.MutPattern(('C','T'), '')

old_pattern = mut_pattern.MutPatternPair(old_focus_pattern, old_control_pattern, 'old')

class TestBasicAnalysis(unittest.TestCase):
    def setUp(self):
        aln_string = """
        >seq1
        GTCAGTCAGTCAGTCACCCC
        >seq2
        GTCAGTCAGTCAGTCACCCC
        >seq3
        ATCAATCAGTCAATCACCCC
        """
        self.seqs = helpers.parse_fasta(aln_string)
        self.aln = HyperfreqAlignment(self.seqs.values())

    def test_consensus_reference(self):
        self.assertEqual(str(self.aln.reference_sequence), 'GTCAGTCAGTCAGTCACCCC')

    def test_hm_pos(self):
        for result in self.aln.analyze((old_focus_pattern, old_control_pattern)):
            hm_pos = result['hm_pos']
            if result['sequence'] in ['seq1', 'seq2']:
                self.assertFalse(hm_pos)
            else:
                self.assertTrue(hm_pos)

    def test_hm_pos_indices(self):
        hm_pos_indices = [1, 5, 13]
        for result in self.aln.analyze((old_focus_pattern, old_control_pattern)):
            if result['hm_pos']:
                self.assertEqual(result['mut_columns'], hm_pos_indices)
            else:
                self.assertEqual(result['mut_columns'], [])


class TestMutCounts(unittest.TestCase):
    def setUp(self):
        aln_string = """
        >seq1
        GGTGACGCT
        >seq2
        AGTAACGCT
        >seq3
        GGTAACACT
        """
        ref_seq = Seq.Seq('GGTGACGCT')
        self.aln = HyperfreqAlignment(helpers.parse_fasta(aln_string).values(), reference_sequence=ref_seq)

    def __test_counts__(self, pattern, real_counts):
        for result in self.aln.analyze(pattern):
            seq_counts = [result[x] for x in ('focus_pos', 'control_pos', 'focus_neg', 'control_neg')]
            seq_real_counts = real_counts[result['sequence']]
            self.assertEqual(seq_counts, seq_real_counts)

    def test_ga_counts(self):
        self.__test_counts__(mut_pattern.GA, dict(
                seq1=[0, 0, 1, 3],
                seq2=[1, 1, 0, 2],
                seq3=[1, 1, 0, 2]))

    def test_gg_counts(self):
        self.__test_counts__(mut_pattern.GG, dict(
                seq1=[0, 0, 1, 3],
                seq2=[1, 1, 0, 2],
                seq3=[0, 2, 1, 1]))

    def test_gr_counts(self):
        self.__test_counts__(mut_pattern.GR, dict(
                seq1=[0, 0, 2, 2],
                seq2=[2, 0, 0, 2],
                seq3=[1, 1, 1, 1]))

    def test_gm_counts(self):
        self.__test_counts__(mut_pattern.GM, dict(
                seq1=[0, 0, 2, 2],
                seq2=[1, 1, 1, 1],
                seq3=[2, 0, 0, 2]))

    def test_gv_counts(self):
        self.__test_counts__(mut_pattern.GV, dict(
                seq1=[0, 0, 3, 1],
                seq2=[2, 0, 1, 1],
                seq3=[2, 0, 1, 1]))


class TestCallPatterns(unittest.TestCase):
    """The following motivates these tests:
        BetaRat(5.5, 0.5, 11.0, 16.0) cdf: 0.00392333077990058 map: 3.27321043801743 ppf: 2.47393902966
        BetaRat(4.5, 1.5, 6.0, 21.0) cdf: 0.00570220397191901 map: 3.475475215347 ppf: 1.98195901549
    As you can see, different methods can make different calls."""

    def setUp(self):
        ref_seq = Seq.Seq('GA'*9 + 'GC'*6 + 'GT'*15)
        query_seq = 'AA'*4 + 'GA'*5 + 'AC' + 'GC'*5 + 'GT'*15
        aln_string = '>seq1\n{0}'.format(query_seq)
        self.aln = HyperfreqAlignment(helpers.parse_fasta(aln_string).values(), reference_sequence=ref_seq)
        self.mutation_patterns=[mut_pattern.GA, mut_pattern.GM]

    def test_map_caller(self):
        for result in self.aln.multiple_context_analysis(self.mutation_patterns, caller='map'):
            call = result['call']
            self.assertEqual(call['call_pattern'], 'GA')

    def test_ppf_caller(self):
        for result in self.aln.multiple_context_analysis(self.mutation_patterns, caller='q_0.05',
                quants=[0.05], pos_quants_only=False):
            call = result['call']
            self.assertEqual(call['call_pattern'], 'GM')

    def test_cdf_caller(self):
        """Note that CDF should be evaluated as smaller => more extreme, whereas the other statistics are the
        inverse"""
        for result in self.aln.multiple_context_analysis(self.mutation_patterns, caller='cutoff_cdf'):
            call = result['call']
            self.assertEqual(call['call_pattern'], 'GM')

    def test_ppf_caller_without_all_quants(self):
        """Should raise if we don't have the given statistic to test for call comparison."""
        analysis = self.aln.multiple_context_analysis(self.mutation_patterns, caller='q_0.05')
        with self.assertRaises(ValueError):
            analysis.next()

    def test_not_calling_negatives(self):
        """We should only call sequences when they are actually positive, even if there are non-positive
        sequences with more extreme call statistics."""
        for result in self.aln.multiple_context_analysis(self.mutation_patterns, caller='map',
                significance_level=0.005):
            call = result['call']
            self.assertEqual(call['call_pattern'], 'GM')


class TestAlignmentSet(unittest.TestCase):
    def setUp(self):
        aln_string = """
        >seq1.1
        GTCAGTCGGTCGGTCAGCCC
        >seq1.2
        GTCAGTCGGTCGGTCAGCCC
        >seq1.3
        AGCAATCAGGCAATCAGCCC
        >seq2.1
        GGTTAACCGGTTAACCGCCC
        >seq2.2
        GGTTAACCGGTTAACCGCCC
        >seq2.3
        AATTAACCAATTAACCACCC
        """
        self.cluster_map = {'cluster1': ['seq1.1', 'seq1.2', 'seq1.3'],
                'cluster2': ['seq2.1', 'seq2.2', 'seq2.3']}
        self.seqs = helpers.parse_fasta(aln_string)

    def test_without_cluster_map(self):
        aln_set = HyperfreqAlignment.Set(self.seqs)
        self.assertEqual(len(aln_set.clusters), 1)
        self.assertEqual(len(aln_set.cluster_alns), 1)
        self.assertIn('all', aln_set.cluster_alns)
        for result in aln_set.multiple_context_analysis([old_pattern]):
            hm_pos = result['call']['hm_pos']
            if result['call']['sequence'] == 'seq2.3':
                self.assertTrue(hm_pos)
            else:
                self.assertFalse(hm_pos)

    def test_with_cluster_map(self):
        aln_set = HyperfreqAlignment.Set(self.seqs, cluster_map=self.cluster_map)
        self.assertEqual(len(aln_set.clusters), 2)
        self.assertEqual(len(aln_set.cluster_alns), 2)
        for result in aln_set.multiple_context_analysis([old_pattern]):
            hm_pos = result['call']['hm_pos']
            if result['call']['sequence'] in ['seq1.3', 'seq2.3']:
                self.assertTrue(hm_pos)
            else:
                self.assertFalse(hm_pos)

    def test_with_reference_seqs(self):
        ref_seqs = helpers.parse_fasta("""
        >cluster1
        AAAAAAAAAAAAAAAAACCC
        >cluster2
        CCTTGGCCGGTTGGCCGCCC
        """)
        aln_set = HyperfreqAlignment.Set(self.seqs, cluster_map=self.cluster_map,
            reference_sequences=ref_seqs)
        self.assertEqual(len(aln_set.clusters), 2)
        for result in aln_set.multiple_context_analysis([old_pattern]):
            hm_pos = result['call']['hm_pos']
            if result['call']['sequence'] in ['seq2.{}'.format(i) for i in [1,2,3]]:
                self.assertTrue(hm_pos)
            else:
                self.assertFalse(hm_pos)


class TestLoadClusterMap(unittest.TestCase):
    def test_loading(self):
        cluster_string = """cluster,sequence
        cluster1,seq1
        cluster1,seq2
        cluster2,seq3
        cluster2,seq4
        """
        cluster_map = {'cluster1': ['seq1', 'seq2'], 'cluster2': ['seq3', 'seq4']}
        loaded_cluster_map = load_cluster_map(helpers.fake_file(cluster_string), cluster_col='cluster')
        self.assertEqual(loaded_cluster_map, cluster_map)


class TestContextBasedEvaluation(unittest.TestCase):
    def setUp(self):
        ref_seq = helpers.parse_fasta("""
        >all
        GGGGGGGGGTGTGTGTGT""")
        self.aln = HyperfreqAlignment(helpers.parse_fasta("""
        >seq1
        GGGGGGGGGTGTGTGTGT
        >seq2
        AGAGAGAGGTGTGTGTGT
        >seq3
        GGGGGGGGGTATATATAT
        """).values(), reference_sequence = ref_seq['all'])

    def test_gg(self):
        for result in self.aln.analyze(mut_pattern.GG):
            hm_pos = result['hm_pos']
            if result['sequence'] == 'seq2':
                self.assertTrue(hm_pos)
            else:
                self.assertFalse(hm_pos)


