import unittest
import helpers
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


