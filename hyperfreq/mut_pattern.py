import re
import Bio.SeqRecord

def negate_context_regexp(string):
    return "(?!{})".format(string)

def negate_mut_regexp(string):
    return "[^{}]".format(string)

def strgfy_seqish(seqish):
    if type(seqish) == Bio.SeqRecord.SeqRecord:
        seq = seqish.seq
    else:
        seq = seqish
    return str(seq)


class MutPattern(object):

    def build_regexp(self, i, negate_mut=False, context=True):
        """ Here `i` effectively controls whether this is the reference_seq regex (i=0) or the query_seq
        (i=1). Negate decides whether we want to capture the pos mutation or the negative (only sensible for
        i=1)."""
        mut = negate_mut_regexp(self.mutation[i]) if negate_mut else self.mutation[i]
        (us, ds) = (self.upstream_context, self.downstream_context) if context else ('', '')
        # Have to be careful here to use lookaheads not to consume parts of the string so that we don't miss
        # overlapding matches
        us = '' if us == '' else "(?<={})".format(us) 
        ds = '' if ds == '' else "(?={})".format(ds)
        string = us + "(?P<mut_index>{})".format(mut) + ds
        return re.compile(string)

    def __init__(self, mutation, downstream_context, upstream_context=''):
        self.mutation = mutation
        self.downstream_context = downstream_context
        self.upstream_context = upstream_context
        self.ref_context_regexp = self.build_regexp(0)
        self.ref_wo_context_regexp = self.build_regexp(0, context=False)
        self.mut_pos_regexp = self.build_regexp(1)
        self.mut_neg_regexp = self.build_regexp(1, negate_mut=True)

    # XXX - here there be dragons! This 'sort' of works, but should really account for the length of the
    # of the context, so that's it's not biasing toward finding matches at the end of the sequence
    def context_negation(self):
        downstream_context = '' if self.downstream_context == '' else negate_context_regexp(self.downstream_context)
        upstream_context = '' if self.upstream_context == '' else negate_context_regexp(self.upstream_context)
        return MutPattern(self.mutation, downstream_context, upstream_context)

    def ref_match_indices(self, seq, enforce_context=False):
        seq = strgfy_seqish(seq)
        regexp = self.ref_context_regexp if enforce_context else self.ref_wo_context_regexp
        return [m.start('mut_index') for m in regexp.finditer(seq)]

    def mut_pos_indices(self, seq):
        seq = strgfy_seqish(seq)
        return [m.start('mut_index') for m in self.mut_pos_regexp.finditer(seq)]

    def mut_neg_indices(self, seq):
        seq = strgfy_seqish(seq)
        return [m.start('mut_index') for m in self.mut_neg_regexp.finditer(seq)]

    def pos_indices(self, seq, ref_indices):
        return [i for i in self.mut_pos_indices(seq) if ref_indices.count(i) > 0]

    def neg_indices(self, seq, ref_indices):
        return [i for i in self.mut_neg_indices(seq) if ref_indices.count(i) > 0]


class MutPatternPair(tuple):
    """This is a convenience class which gives us access to control pattern Pair as though they are tuples,
    but with a name function (for call pattern output column, etc.)"""
    def __new__(cls, focus_pattern, control_pattern, name=None):
        return super(MutPatternPair, cls).__new__(cls, tuple((focus_pattern, control_pattern)))

    def __init__(self, focus_pattern, control_pattern, name=None):
        super(MutPatternPair, self).__init__(focus_pattern, control_pattern)
        self.name = name if name else focus_pattern.mutation[0] + focus_pattern.downstream_context



GG_FOCUS = MutPattern(('G','A'), 'G')
GG_CONTROL = MutPattern(('G','A'), '[^G]')
GG = MutPatternPair(GG_FOCUS, GG_CONTROL)

GA_FOCUS = MutPattern(('G','A'), 'A')
GA_CONTROL = MutPattern(('G','A'), '[^A]')
GA = MutPatternPair(GA_FOCUS, GA_CONTROL)

GM_FOCUS = MutPattern(('G','A'), '[AC]')
GM_CONTROL = MutPattern(('G','A'), '[^AC]')
GM = MutPatternPair(GM_FOCUS, GM_CONTROL, 'GM')

GR_FOCUS = MutPattern(('G','A'), '[GA]')
GR_CONTROL = MutPattern(('G','A'), '[^GA]')
GR = MutPatternPair(GR_FOCUS, GR_CONTROL, 'GR')

GV_FOCUS = MutPattern(('G','A'), '[^T]')
GV_CONTROL = MutPattern(('G','A'), 'T')
GV = MutPatternPair(GV_FOCUS, GV_CONTROL, 'GV')

patterns = dict(GG=GG, GA=GA, GM=GM, GR=GR, GV=GV)


