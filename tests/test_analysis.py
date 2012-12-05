import unittest
import helpers
from hyperfreq.hyperfreq_alignment import HyperfreqAlignment
from hyperfreq.cluster import load_cluster_map
from hyperfreq import mut_pattern

old_focus_pattern = mut_pattern.MutPattern(('G','A'), '')
old_control_pattern = mut_pattern.MutPattern(('C','T'), '')

class TestBasicAnalysis(unittest.TestCase):
    def setUp(self):
        aln_string = """
        >seq1
        GTCAGTCAGTCAGTCA
        >seq2
        GTCAGTCAGTCAGTCA
        >seq3
        ATCAATCAGTCAATCA"""
        self.seqs = helpers.parse_fasta(aln_string)
        self.aln = HyperfreqAlignment(self.seqs.values())
        self.aln.analyze_hypermuts(focus_pattern=old_focus_pattern, control_pattern=old_control_pattern)

    def test_hm_pos(self):
        self.assertTrue(helpers.find_seq(self.aln, 'seq3').hm_pos)
        for seq in ['seq1', 'seq2']:
            self.assertFalse(helpers.find_seq(self.aln, seq).hm_pos)

    def test_hm_pos_indices(self):
        hm_pos_indices = [0, 4, 12]
        self.assertEqual(helpers.find_seq(self.aln, 'seq3').focus_pos_indices, hm_pos_indices)
        self.assertEqual(self.aln.hm_pos_indices, hm_pos_indices)
        self.assertEqual(self.aln.mut_columns, [i+1 for i in hm_pos_indices])


class TestAlignmentSet(unittest.TestCase):
    def setUp(self):
        aln_string = """
        >seq1.1
        GTCAGTCAGTCAGTCAG
        >seq1.2
        GTCAGTCAGTCAGTCAG
        >seq1.3
        ATCAATCAGTCAATCAG
        >seq2.1
        GGTTAACCGGTTAACCG
        >seq2.2
        GGTTAACCGGTTAACCG
        >seq2.3
        AGTTAACCAGTTAACCA
        """
        self.cluster_map = {'cluster1': ['seq1.1', 'seq1.2', 'seq1.3'],
                'cluster2': ['seq2.1', 'seq2.2', 'seq2.3']}
        self.seqs = helpers.parse_fasta(aln_string)

    def test_without_cluster_map(self):
        aln_set = HyperfreqAlignment.Set(self.seqs)
        aln_set.analyze_hypermuts(focus_pattern=old_focus_pattern, control_pattern=old_control_pattern)
        self.assertEqual(len(aln_set.clusters), 1)
        self.assertEqual(len(aln_set.cluster_alns), 1)
        self.assertIn('all', aln_set.cluster_alns)
        all_seqs = [seq for aln in aln_set.cluster_alns.values() for seq in aln]
        # Constructed to be false for this one, since not enough consensus with whole align
        self.assertFalse(helpers.find_seq(all_seqs, 'seq1.3').hm_pos)
        # XXX - why is this failing now?
        self.assertTrue(helpers.find_seq(all_seqs, 'seq2.3').hm_pos)
        self.assertEqual(sum([s.hm_pos for s in all_seqs]), 1)

    def test_with_cluster_map(self):
        aln_set = HyperfreqAlignment.Set(self.seqs, cluster_map=self.cluster_map)
        aln_set.analyze_hypermuts(focus_pattern=old_focus_pattern, control_pattern=old_control_pattern)
        self.assertEqual(len(aln_set.clusters), 2)
        self.assertEqual(len(aln_set.cluster_alns), 2)
        all_seqs = [seq for aln in aln_set.cluster_alns.values() for seq in aln]
        self.assertTrue(helpers.find_seq(all_seqs, 'seq1.3').hm_pos)
        self.assertTrue(helpers.find_seq(all_seqs, 'seq2.3').hm_pos)
        self.assertEqual(sum([s.hm_pos for s in all_seqs]), 2)

    def test_with_reference_seqs(self):
        ref_seqs = helpers.parse_fasta("""
        >cluster1
        AAAAAAAAAAAAAAAAA
        >cluster2
        CCTTGGCCGGTTGGCCG
        """)
        aln_set = HyperfreqAlignment.Set(self.seqs, cluster_map=self.cluster_map,
            reference_sequences=ref_seqs)
        aln_set.analyze_hypermuts(focus_pattern=old_focus_pattern, control_pattern=old_control_pattern)
        all_seqs = [seq for aln in aln_set.cluster_alns.values() for seq in aln]
        for seq in ['seq1.{}'.format(i) for i in [1,2,3]]:
            self.assertFalse(helpers.find_seq(all_seqs, seq).hm_pos)
        for seq in ['seq2.{}'.format(i) for i in [1,2,3]]:
            self.assertTrue(helpers.find_seq(all_seqs, seq).hm_pos)


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


