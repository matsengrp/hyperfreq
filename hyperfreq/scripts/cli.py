#!/usr/bin/env python

import argparse
import csv
from os import path
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from hyperfreq.cluster import load_cluster_map
from hyperfreq.hyperfreq_alignment import HyperfreqAlignment, analysis_defaults
from hyperfreq import mut_pattern, hyperfreq_alignment, __version__
from hyperfreq.analysis_writer import write_analysis


def split(args):
    hm_col_reader = csv.DictReader(args.columns)
    hm_columns = map(lambda x: int(x['column']), hm_col_reader)
    hm_columns = list(set(hm_columns))

    seq_records = SeqIO.parse(args.alignment, 'fasta')
    aln = HyperfreqAlignment(seq_records)
    aln.split_hypermuts(hm_columns = hm_columns)

    fn_base = path.join(args.out_dir, args.prefix)
    hm_pos_handle = open(fn_base + '.pos.fasta', 'w')
    hm_neg_handle = open(fn_base + '.neg.fasta', 'w')

    AlignIO.write(aln.hm_pos_aln, hm_pos_handle, 'fasta')
    AlignIO.write(aln.hm_neg_aln, hm_neg_handle, 'fasta')

    for handle in [args.alignment, args.columns, hm_pos_handle, hm_neg_handle]:
        handle.close()


def write_reference_seqs(alignments, fn_base):
    handle = open(fn_base + '.ref_seqs.fasta', 'w')
    def refseq(cluster):
        seq = alignments.cluster_alns[cluster].reference_sequence
        return SeqRecord(seq, id=cluster, name=cluster, description="")
    refseqs = (refseq(cluster) for cluster in alignments.clusters)
    SeqIO.write(refseqs, handle, 'fasta')
    handle.close()


def analyze(args):
    seq_records = SeqIO.to_dict(SeqIO.parse(args.alignment, 'fasta'))
    if args.reference_sequences:
        reference_sequences = SeqIO.to_dict(SeqIO.parse(args.reference_sequences, 'fasta'))
    else:
        reference_sequences = None

    # This lets the cluster map be aptional, so that this script can be used
    # for naive hm filtering/analysis
    cluster_map = load_cluster_map(args.cluster_map, cluster_col=args.cluster_col) if args.cluster_map else None
    alignments = HyperfreqAlignment.Set(seq_records, cluster_map,
            reference_sequences=reference_sequences)
    patterns = [mut_pattern.patterns[p] for p in args.patterns]
    
    # Create the analysis generator, and run it as we write out to files
    analysis = alignments.multiple_context_analysis(patterns, consensus_threshold=args.consensus_threshold,
            rpr_cutoff=args.rpr_cutoff, significance_level=args.significance_level, quants=args.quants,
            pos_quants_only=args.pos_quants_only, caller=args.caller, prior=args.prior, cdfs=args.cdfs)
    prefix = path.join(args.out_dir, args.prefix)
    pattern_names = [p.name for p in patterns]
    write_analysis(analysis, prefix, pattern_names, args.quants, args.cdfs, call_only=args.call_only)

    if args.write_references:
        write_reference_seqs(alignments, prefix)

    # Closing files
    args.alignment.close()
    if args.cluster_map:
        args.cluster_map.close()


def setup_common_args(subparser):
    subparser.add_argument('alignment', type=argparse.FileType('r'),
            help="""Sequence alignment on which to operate. (This argument must come before any arguments
            which take multiple inputs, such as --patterns and --cdfs)""")
    subparser.add_argument('-o', '--out-dir', default='.',
            help="Where to put files")
    subparser.add_argument('-P', '--prefix',
            help="Prefix for output files (extensions chosen automatically)")
    subparser.add_argument('-v', '--verbose', action='store_true', default=False)


class QuantAction(argparse.Action):
    """ This class is for doing some slick command line magick with specification of what quantiles to compute
    and how """
    # XXX - review: This solution to the problem leaves q around in namespace, and admittedly is a little
    # weird. Will have to poll thoughts
    default_quants = [0.05]

    def __call__(self, parser, namespace, values, option_string=None):
        if values:
            namespace.quants = values
        else:
            namespace.quants = parser.get_default('q')
            parser.set_defaults(q=None)
        namespace.pos_quants_only = True if option_string == '-q' else False

    @classmethod
    def register(cls, parser):
        """ Adds the -q/-Q argument to the parser, and sets stuff up so that the results can be accessed from
        the quants namespace, as well as pos_quants_only """
        parser.add_argument('-q', '-Q', nargs='*', type=float, action=cls,
                help="""Compute quantiles, separated by spaces. If specified with no args, default quantiles are
                %(default)s. If specified using -q, quantiles are only computed for positive sequences to save time
                (quantiles take a while). If -Q is used, specified quantiles are computed for all sequences.""")
        parser.set_defaults(q=cls.default_quants)


def setup_analyze_args(subparsers):
    def cs_arg(arg):
        return arg.split(',')
    analyze_args = subparsers.add_parser('analyze',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description='Analyze alignment for evidence of hypermuation')
    setup_common_args(analyze_args)
    # Should we get this to be smarter about whether to write out a call file when only one analysis is run?
    pattern_choices = mut_pattern.patterns.keys()
    analyze_args.add_argument('-p', '--patterns', choices=pattern_choices, default=['GG'], nargs='+',
            help="""Specify the type of apobec activity to analyze. For example, 'GG' specifies
            a focus pattern of GG to AG, characteristic of APOBEC3G activity. Multiple patterns should be
            separated by spaces. Characters M, R and V correspond to IUPAC codes.""")
    analyze_args.add_argument('--rpr-cutoff', type=float,
            help="""For hm_pos determination: if a sequence has a RPR higher than this value with confidence
            specified by --significance-level, it will be marked as hypermutation positive.""")
    analyze_args.add_argument('-s', '--significance-level', type=float,
            help="""For hm_pos determination: if, with this specified level of confidence, a sequence has a RPR
            higher than the specified rpr_cutoff, it will be marked as hypermutation positive.""")
    QuantAction.register(analyze_args)
    analyze_args.add_argument('--cdfs', nargs="+", type=float,
            help="""Specify cdfs to be computed (separated by spaces). These are computed in addition to the
            CDF of the rpr-cutoff, as described above.""")

    prior_group = analyze_args.add_mutually_exclusive_group()
    prior_group.add_argument('--prior', type=float, nargs=2,
            help="""Prior on Beta distributions. The default (%(default)s) corresponds to a belief that
            mutations are relatively rare, but unbiased with respect to context.""")
    prior_group.add_argument('--jeff', help='Use Jeffreys prior (0.5, 0.5) on Beta distributions.', action='store_const',
            dest='prior', const=(0.5, 0.5))
    prior_group.add_argument('--uniform', help='Uniform prior (1.0, 1.0) on Beta distributions.', action='store_const',
            dest='prior', const=(1.0, 1.0))

    analyze_args.add_argument('--caller',
            help="""Statistic to be used for deciding which mutation pattern has the strongest hypermutation
            signal. The choice specified should be the name of a column in the output file, such as "map",
            "cutoff_cdf" or "q_0.05". Note: the value must exist for each sequence and call pattern analyzed.
            As such, you must use the `-Q` flag (for computing all quantiles) if you wish to call based on
            quantiles.""")
    analyze_args.add_argument('-c', '--cluster-map', type=argparse.FileType('r'),
            help="CSV file mapping sequences to clusters")
    # Should make this smarter so that it guesses a few things first...
    analyze_args.add_argument('--cluster-col', default='cluster',
            help="Column in cluster_map file to be used for cluster specification")
    analyze_args.add_argument('--consensus-threshold', type=float,
            help="""Used for computing consensus sequences as reference sequences for
            HM evaluation when reference sequences are not explicity specified.""")
    # Should remove necessity for "all" and just take the first, if no matches (with warning?)
    analyze_args.add_argument('-r', '--reference-sequences', type=argparse.FileType('r'),
            help="""If specified, use the reference sequences in this file for comparison instead of consensus
            sequences. Sequence name(s) should be the names of the clusters if using a cluster map. Otherwise,
            ensure that there is a sequence named 'all' in your reference_sequences alignment.
            Clusters for which no reference sequence is specified will be compared to a computed consensus
            sequence as reference.""")
    analyze_args.add_argument('-R', '--write-references', default=False, action='store_true',
            help="""Writes reference sequences (typically consensus sequences) used for HM evalutation. If
            sequences are clustered, the reference sequence for each cluster will be given the name of the
            cluster. Consequently, the output file can be used subsequently as input to --reference-sequences""")
    analyze_args.add_argument('-F', '--full-output', default=True, action='store_false', dest='call_only',
            help="""Generate a separate file for each of the analysis patterns (default: call_only=%(default)s)""")
    for arg, default in analysis_defaults.iteritems():
        try:
            analyze_args.set_defaults(**dict([(arg, default)]))
        except KeyError:
            print "{} is not an args option".format(arg)
    analyze_args.set_defaults(prefix='hyperfreq_analysis')
    analyze_args.set_defaults(func=analyze)


def setup_split_args(subparsers):
    split_args = subparsers.add_parser('split', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    setup_common_args(split_args)
    split_args.add_argument('columns', type=argparse.FileType('r'),
            help="""File identifying hypermuted columns. Typically, you would want to use the `*.sites.csv`
            output file from running `hyperfreq analyze`, but you can contruct your own file if you like.""")
    split_args.add_argument('--column-name', default='column',
            help="""Column name in the columns file which identifies the hypermutated sites.""")
    split_args.set_defaults(prefix='hyperfreq_split')
    split_args.set_defaults(func=split)


def main():
    parser = argparse.ArgumentParser(prog='hyperfreq',
            description="""Hypermutation analysis software using BetaRat distribution for Bayesian analysis of
            the relative probability ratio (RPR) of observing mutations in two contexts.""",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))
    subparsers = parser.add_subparsers(title='subcommands', help='additional help')

    setup_analyze_args(subparsers)
    setup_split_args(subparsers)

    args = parser.parse_args()

    if args.verbose:
        hyperfreq_alignment.VERBOSE = True

    args.func(args)


if __name__ == '__main__':
    main()


