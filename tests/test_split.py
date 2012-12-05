import unittest
import helpers
from hyperfreq.hyperfreq_alignment import HyperfreqAlignment
from test_analysis import old_focus_pattern, old_control_pattern

class TestBasicSplit(unittest.TestCase):
    def assertSeqsEqual(self, seq_record, string):
        self.assertEqual(str(seq_record.seq), string)

    def setUp(self):
        aln_string = """
        >seq1
        GTCAGTCAGTCAGTCA
        >seq2
        GTCAGTCAGTCAGTCA
        >seq3
        ATCAATCAGTCAATCG"""
        self.seqs = helpers.parse_fasta_list(aln_string)
        self.aln = HyperfreqAlignment(self.seqs)

    def test_split_from_analysis_indices(self):
        # hm_pos should == [0, 4, 12]
        self.aln.analyze_hypermuts(focus_pattern=old_focus_pattern, control_pattern=old_control_pattern)
        self.aln.split_hypermuts()
        neg, pos = self.aln.hm_neg_aln, self.aln.hm_pos_aln
        self.assertEqual(neg.get_alignment_length(), 13)
        self.assertEqual(pos.get_alignment_length(), 3)
        self.assertEqual(neg[:,0], 'TTT')
        self.assertEqual(neg[:,12], 'AAG')
        self.assertSeqsEqual(neg[0,:], 'TCATCAGTCATCA')
        self.assertSeqsEqual(neg[2,:], 'TCATCAGTCATCG')
        self.assertEqual(pos[:,0], 'GGA')
        self.assertSeqsEqual(pos[0,:], 'GGG')
        self.assertSeqsEqual(pos[2,:], 'AAA')

    def test_manual_split(self):
        columns = [1, 2, 3, 5]
        self.aln.split_hypermuts(hm_columns=columns)
        neg, pos = self.aln.hm_neg_aln, self.aln.hm_pos_aln
        self.assertEqual(neg.get_alignment_length(), 12)
        self.assertEqual(pos.get_alignment_length(), 4)
        self.assertEqual(neg[:,0], 'AAA')
        self.assertEqual(neg[:,11], 'AAG')
        self.assertSeqsEqual(neg[0,:], 'ATCAGTCAGTCA')
        self.assertEqual(pos[:,1], 'TTT')

    def test_splitting_final_col(self):
        columns = [3, 7, 16]
        self.aln.split_hypermuts(hm_columns=columns)
        neg, pos = self.aln.hm_neg_aln, self.aln.hm_pos_aln
        self.assertEqual(neg.get_alignment_length(), 13)
        self.assertEqual(pos.get_alignment_length(), 3)
        self.assertEqual(neg[:,0], 'GGA')
        self.assertEqual(neg[:,12], 'CCC')
        self.assertSeqsEqual(neg[0,:], 'GTAGTAGTCAGTC')
        self.assertEqual(pos[:,1], 'CCC')
        self.assertEqual(pos[:,2], 'AAG')


