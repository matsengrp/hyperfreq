import re

def regexp_negation(string):
    return "[^({})]".format(string)

class MutPattern(object):

    def __init__(self, mutation, downstream_context, upstream_context=''):
        self.mutation = mutation
        self.downstream_context = downstream_context
        self.upstream_context = upstream_context
        self.reference_regexp = re.compile(upstream_context + mutation[0] + downstream_context)
        self.mut_regexp = re.compile(upstream_context + mutation[1] + downstream_context)

    def context_negation(self):
        upstream_context = '' if self.upstream_context == '' else regexp_negation(self.upstream_context)
        downstream_context = regexp_negation(self.upstream_context)
        return MutPattern(self.mutation, downstream_context, upstream_context)
    

A3G_FOCUS = MutPattern(('G','A'), 'G[^C]')
A3G_CONTROL = A3G_FOCUS.context_negation()

A3F_FOCUS = MutPattern(('G','A'), '[AC][^C]')
A3F_CONTROL = A3F_FOCUS.context_negation()
