import unittest
#import helpers
from hyperfreq import mut_pattern


class BasicTests(unittest.TestCase):
    def test_creation(self):
        pattern = mut_pattern.MutPattern(('G','A'), 'GT')
        self.assertIsInstance(pattern, mut_pattern.MutPattern)

    def test_creation_with_upstream(self):
        pattern = mut_pattern.MutPattern(('G','A'), 'GT', 'C')
        self.assertIsInstance(pattern, mut_pattern.MutPattern)

    def test_negation_without_upstream(self):
        pattern = mut_pattern.MutPattern(('G','A'), 'GT')
        neg_pattern = pattern.context_negation()
        self.assertIsInstance(neg_pattern, mut_pattern.MutPattern)
        self.assertEqual(neg_pattern.downstream_context, '(?!GT)')
        self.assertEqual(neg_pattern.upstream_context, '')

    def test_negation_with_upstream(self):
        pattern = mut_pattern.MutPattern(('G','A'), 'GT', 'C')
        neg_pattern = pattern.context_negation()
        self.assertIsInstance(neg_pattern, mut_pattern.MutPattern)
        self.assertEqual(neg_pattern.upstream_context, '(?!C)')


class FunctionalTests(unittest.TestCase):
    def setUp(self):
        self.pattern = mut_pattern.MutPattern(('G','A'), 'A')

    def test_ref_match_indices_without_context(self):
        seq = 'ACGTACGT'
        ref_indices = self.pattern.ref_match_indices(seq, enforce_context=False)
        self.assertEqual([2,6], ref_indices)

    def test_ref_match_indices_with_context(self):
        seq = 'TTGATTGC'
        ref_indices = self.pattern.ref_match_indices(seq, enforce_context=True)
        self.assertEqual([2], ref_indices)

    def test_ref_match_indices_with_overlap(self):
        seq = 'TCGGGAT'
        mp = mut_pattern.MutPattern(('G','A'), 'G')
        ref_indices = mp.ref_match_indices(seq, enforce_context=True)
        self.assertEqual([2,3], ref_indices)

    def test_ref_match_indices_with_massive_overlap(self):
        seq = 'GGGGG'
        mp = mut_pattern.MutPattern(('G','A'), 'G', 'G')
        ref_indices = mp.ref_match_indices(seq, enforce_context=True)
        self.assertEqual([1,2,3], ref_indices)

    def test_mut_pos_indices(self):
        mp = mut_pattern.MutPattern(('G','A'), '[AG]')
        query_seq = 'GAGAAGAT'
        ref_indices = mp.mut_pos_indices(query_seq)
        self.assertEqual([1,3,4], ref_indices)

    def test_pos_indices(self):
        mp = mut_pattern.MutPattern(('G','A'), '[AG]')
        query_seq = 'GAGAAGAT'
        ref_indices = mp.pos_indices(query_seq, [1,2,3])
        self.assertEqual([1,3], ref_indices)

    def test_mut_neg_indices(self):
        mp = mut_pattern.MutPattern(('G','A'), '[AG]')
        query_seq = 'GAGAAGAT'
        ref_indices = mp.mut_neg_indices(query_seq)
        self.assertEqual([0,2,5], ref_indices)


if __name__ == '__main__':
    a3g_focus = mut_pattern.A3G_FOCUS
    a3g_control = mut_pattern.A3G_CONTROL

    print a3g_focus.mut_pos_indices('AGAGGGG')
    print a3g_focus.mut_neg_indices('AGAGGGG')
    print a3g_control.mut_pos_indices('ATATGGCG')
    print a3g_control.mut_neg_indices('ATATGGCG')

