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
        self.assertEqual(neg_pattern.downstream_context, '[^(GT)]')
        self.assertEqual(neg_pattern.upstream_context, '')

    def test_negation_with_upstream(self):
        pattern = mut_pattern.MutPattern(('G','A'), 'GT', 'C')
        neg_pattern = pattern.context_negation()
        self.assertIsInstance(neg_pattern, mut_pattern.MutPattern)
        self.assertEqual(neg_pattern.upstream_context, '[^(C)]')


class FunctionalTests(unittest.TestCase):
    def setUp(self):
        self.pattern = mut_pattern.MutPattern(('G','A'), 'GT')
