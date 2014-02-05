import unittest
import helpers
from hyperfreq.core import Alignment

class TestBasicSplit(unittest.TestCase):
    def assertSeqsEqual(self, seq_record, string):
        self.assertEqual(str(seq_record.seq), string)

    def setUp(self):
        aln_string = """
        >seq1
        GTCAGTCAGTCAGTCA
        GTCAGTCAGTCAGTCA
        >seq2
        GTCAGTCAGTCAGTCA
        GTCAGTCAGTCAGTCA
        >seq3
        ATCAATCAGTCAATCG
        ATCAATCAGTCAATCG"""
        self.seqs = helpers.parse_fasta_list(aln_string)
        self.aln = Alignment(self.seqs)

    def test_manual_split(self):
        columns = [1, 2, 3, 5]
        self.aln.split_hypermuts(hm_columns=columns)
        neg, pos = self.aln.hm_neg_aln, self.aln.hm_pos_aln
        self.assertEqual(neg.get_alignment_length(), 28)
        self.assertEqual(pos.get_alignment_length(), 4)
        self.assertEqual(neg[:,0], 'AAA')
        self.assertEqual(neg[:,11], 'AAG')
        self.assertSeqsEqual(neg[0,:], 'ATCAGTCAGTCAGTCAGTCAGTCAGTCA')
        self.assertEqual(pos[:,1], 'TTT')

    def test_splitting_final_col(self):
        columns = [3, 7, 32]
        self.aln.split_hypermuts(hm_columns=columns)
        neg, pos = self.aln.hm_neg_aln, self.aln.hm_pos_aln
        self.assertEqual(neg.get_alignment_length(), 29)
        self.assertEqual(pos.get_alignment_length(), 3)
        self.assertEqual(neg[:,0], 'GGA')
        self.assertEqual(neg[:,12], 'CCC')
        self.assertSeqsEqual(neg[0,:], 'GTAGTAGTCAGTCAGTCAGTCAGTCAGTC')
        self.assertEqual(pos[:,1], 'CCC')
        self.assertEqual(pos[:,2], 'AAG')

    def test_splitting_on_no_hm(self):
        columns = []
        self.aln.split_hypermuts(hm_columns=columns)
        neg, pos = self.aln.hm_neg_aln, self.aln.hm_pos_aln
        self.assertEqual(neg.get_alignment_length(), 32)
        self.assertEqual(pos.get_alignment_length(), 0)
        self.assertEqual(neg[:,0], 'GGA')
        self.assertSeqsEqual(neg[0,:], 'GTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCA')
        self.assertSeqsEqual(pos[0,:], '')


