import re

def regexp_negation(string):
    return "[^({})]".format(string)

class MutPattern(object):

    def build_regexp(self, i, negate_mut=False, context=True):
        """ Here `i` effectively controls whether this is the reference_seq regex (i=0) or the query_seq
        (i=1). Negate decides whether we want to capture the pos mutation or the negative (only sensible for
        i=1)."""
        mut = regexp_negation(self.mutation[i]) if negate_mut else self.mutation[i]
        (us, ds) = self.upstream_context, self.downstream_context if context else ('', '')
        string = "{}(?P<mut_index>{}){}".format(self.upstream_context, mut, self.downstream_context)
        return re.compile(string)

    def __init__(self, mutation, downstream_context, upstream_context=''):
        self.mutation = mutation
        self.downstream_context = downstream_context
        self.upstream_context = upstream_context
        self.ref_context_regexp = self.build_regexp(0)
        self.ref_wo_context_regexp = self.build_regexp(0, context=False)
        self.mut_pos_regexp = self.build_regexp(1)
        self.mut_neg_regexp = self.build_regexp(1, negate_mut=True)

    def context_negation(self):
        upstream_context = '' if self.upstream_context == '' else regexp_negation(self.upstream_context)
        downstream_context = regexp_negation(self.upstream_context)
        return MutPattern(self.mutation, downstream_context, upstream_context)

    def ref_match_indices(self, seq, enforce_context=False):
        seq = str(seq)
        regexp = self.ref_context_regexp if enforce_context else self.ref_wo_context_regexp
        #import pdb; pdb.set_trace()
        return [m.start('mut_index') for m in regexp.finditer(seq)]

    def mut_pos_indices(self, seq):
        return [m.start('mut_index') for m in self.mut_pos_regexp.finditer(seq)]

    def mut_neg_indices(self, seq):
        return [m.start('mut_index') for m in self.mut_neg_regexp.finditer(seq)]

    def pos_indices(self, seq, ref_indices):
        seq = str(seq)
        return [i for i in self.mut_pos_indices(seq) if ref_indices.count(i) > 0]

    def neg_indices(self, seq, ref_indices):
        seq = str(seq)
        return [i for i in self.mut_neg_indices(seq) if ref_indices.count(i) > 0]
    

    # XXX - need to make sure we get the matching on both vs matching on one right
    # not sure how this will work yet



A3G_FOCUS = MutPattern(('G','A'), 'G[^C]')
A3G_CONTROL = A3G_FOCUS.context_negation()

A3F_FOCUS = MutPattern(('G','A'), '[AC][^C]')
A3F_CONTROL = A3F_FOCUS.context_negation()

