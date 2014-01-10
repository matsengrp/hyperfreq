#!/usr/bin/env python

import argparse
import csv
from hyperfreq import alnclst
from os import path
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from hyperfreq.cluster import load_cluster_map, parse_clusters
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
    import logging; logging.captureWarnings(True)
    # Fetch sequence records and analysis patterns
    seq_records = SeqIO.to_dict(SeqIO.parse(args.alignment, 'fasta'))
    patterns = [mut_pattern.patterns[p] for p in args.patterns]
    pattern_names = [p.name for p in patterns]
    prefix = path.join(args.out_dir, args.prefix)
    analysis_settings = dict(
            rpr_cutoff=args.rpr_cutoff, significance_level=args.significance_level, quants=args.quants,
            pos_quants_only=args.pos_quants_only, caller=args.caller, prior=args.prior, cdfs=args.cdfs)

    # Need to think about how best to fork things here; for instance, might make sense to let the user specify
    # the initial clusters for whatever reason... However, specifying the reference sequences shouldn't make
    # any sense there
    if args.reference_sequences:
        reference_sequences = SeqIO.to_dict(SeqIO.parse(args.reference_sequences, 'fasta'))
    else:
        reference_sequences = None

    # This lets the cluster map be aptional, so that this script can be used
    # for naive hm filtering/analysis
    cluster_map = load_cluster_map(args.cluster_map, cluster_col=args.cluster_col) if args.cluster_map else None
    alignments = HyperfreqAlignment.Set(seq_records, cluster_map, consensus_threshold=args.consensus_threshold,
            reference_sequences=reference_sequences)

    # Create the analysis generator
    analysis = alignments.multiple_context_analysis(patterns, **analysis_settings)

    if args.cluster_threshold:
        for hm_it in range(args.cluster_iterations - 1):
            print("On hm/cluster iteration", hm_it)
            # Grab the HM columns from the most recent analysis and split out the pos sites
            hm_columns = []
            for result in analysis:
                hm_columns += result['call']['mut_columns']
            hm_neg_aln = HyperfreqAlignment(seq_records.values()).split_hypermuts(hm_columns).hm_neg_aln
            # Cluster with the specified settings
            clustering = alnclst.Clustering(hm_neg_aln, args.cluster_threshold,
                    args.consensus_threshold)
            clustering = clustering.recenter(args.recentering_iterations)
            cluster_map = parse_clusters(clustering.mapping_iterator(), cluster_key=0, sequence_key=1)
            # Create the Alignment set
            clustered_alignment = HyperfreqAlignment.Set(seq_records, cluster_map,
                    consensus_threshold=args.consensus_threshold)
            analysis = clustered_alignment.multiple_context_analysis(patterns, **analysis_settings)

    # Write the final analysis to file
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
    analyze_args.add_argument('-F', '--full-output', default=True, action='store_false', dest='call_only',
            help="""Generate a separate file for each of the analysis patterns instead of just the call
            pattern (default: call_only=%(default)s)""")

    eval_group = analyze_args.add_argument_group('HM EVALUATION SETTINGS')
    # Should we get this to be smarter about whether to write out a call file when only one analysis is run?
    pattern_choices = mut_pattern.patterns.keys()
    eval_group.add_argument('-p', '--patterns', choices=pattern_choices, default=['GG'], nargs='+',
            help="""Specify the type of apobec activity to analyze. For example, 'GG' specifies
            a focus pattern of GG to AG, characteristic of APOBEC3G activity. Multiple patterns should be
            separated by spaces. Characters M, R and V correspond to IUPAC codes.""")
    eval_group.add_argument('--rpr-cutoff', type=float,
            help="""For hm_pos determination: if a sequence has a RPR higher than this value with confidence
            specified by --significance-level, it will be marked as hypermutation positive.""")
    eval_group.add_argument('-s', '--significance-level', type=float,
            help="""For hm_pos determination: if, with this specified level of confidence, a sequence has a RPR
            higher than the specified rpr_cutoff, it will be marked as hypermutation positive.""")
    QuantAction.register(eval_group)
    eval_group.add_argument('--cdfs', nargs="+", type=float,
            help="""Specify cdfs to be computed (separated by spaces). These are computed in addition to the
            CDF of the rpr-cutoff, as described above.""")
    eval_group.add_argument('--caller',
            help="""Statistic to be used for deciding which mutation pattern has the strongest hypermutation
            signal. The choice specified should be the name of a column in the output file, such as "map",
            "cutoff_cdf" or "q_0.05". Note: the value must exist for each sequence and call pattern analyzed.
            As such, you must use the `-Q` flag (for computing all quantiles) if you wish to call based on
            quantiles.""")

    prior_group = analyze_args.add_argument_group('PRIORS')
    prior_group = prior_group.add_mutually_exclusive_group()
    prior_group.add_argument('--prior', type=float, nargs=2,
            help="""Prior on Beta distributions. The default (%(default)s) corresponds to a belief that
            mutations are relatively rare, but unbiased with respect to context.""")
    prior_group.add_argument('--jeff', help='Use Jeffreys prior (0.5, 0.5) on Beta distributions.', action='store_const',
            dest='prior', const=(0.5, 0.5))
    prior_group.add_argument('--uniform', help='Uniform prior (1.0, 1.0) on Beta distributions.', action='store_const',
            dest='prior', const=(1.0, 1.0))

    refseq_group = analyze_args.add_argument_group('REFERENCE SEQS',
            """Hyperfreq requires each query sequence have a reference sequence for comparison. This
            sequence should be evolutionarily close while not exhibiting a hypermutation pattern.
            By default, hyperfreq computes this sequence as global consensus from the query
            alignment. You can also manually specify reference sequences and compute them from clusters
            of sequences (see the ITERATIVE CLUSTERING SETTINGS group for automatic identification of
            clusters).""")
    refseq_group.add_argument('-c', '--cluster-map', type=argparse.FileType('r'),
            help="""CSV file mapping sequences to clusters; cluster consensus sequences will be used as query
            sequences. Any sequences not mapped to a cluster will be implicity put in an 'all' cluster.""")
    # Should make this smarter so that it guesses a few things first...
    refseq_group.add_argument('--cluster-col', default='cluster',
            help="Column in cluster-map to be used for cluster specification")
    refseq_group.add_argument('--consensus-threshold', type=float,
            help="""For computing consensus sequences. See biopython's AlignInfo.SummmaryInfo.dumb_consensus
            method. (default: %(default)s; no threshold, most frequent base taken)""")
    # Should remove necessity for "all" and just take the first, if no matches (with warning?)
    refseq_group.add_argument('-r', '--reference-sequences', type=argparse.FileType('r'),
            help="""Manually specify reference sequences. Sequence name(s) should correspond to cluster names
            if using a cluster map. Otherwise, ensure a sequence named 'all' is included in the alignment.
            Clusters for which no reference sequence is specified will be compared to a computed consensus
            sequence.""")
    refseq_group.add_argument('-R', '--write-references', default=False, action='store_true',
            help="""Writes reference sequences used for HM evalutation. If sequences are clustered, the
            reference sequence for each cluster will be given the name of the cluster. Consequently, the
            output file can be used subsequently as input to --reference-sequences""")

    # make some mutually exclusive with refseq spec methods; other future options
    # -M --min-per-cluster-percent; -C --write-intermediate-clusters; -g --global-cons-first
    autoclst_group = analyze_args.add_argument_group('ITERATIVE CLUSTERING SETTINGS',
            """Use hyperfreq's iterative clustering strategy to find reference sequences. This involves
            iterations of hypermutation analysis, removal of hypermutated columns from alignment, and
            clustering of these 'HM free' alignments, so that clusters don't reflect hypermutation within the
            data.""")
    autoclst_group.add_argument('-t', '--cluster-threshold', type=float,
            help="""[Required for iterative clustering] If specified, triggers the iterative clustering
            algorithm with the given clustering similarity threshold.""")
    autoclst_group.add_argument('-i', '--cluster-iterations', type=int, default=5,
            help="""Number of iterations of hm analysis and clustering.""")
    autoclst_group.add_argument('-I', '--recentering-iterations', type=int, default=4,
            help="""Not to be confused with --cluster-iterations, this specifies the number of recentering
            steps to perform for each clustering step (inspired by http://goo.gl/RIoWBU).""")
    autoclst_group.add_argument('-m', '--min-per-cluster', type=int, default=5,
            help="""This value specifies the minimum number of sequences in a cluster. Clusters smaller than
            this value will be merged with the closest cluster until not small clusters are left. This avoids
            comparing a query sequence to a consensus sequence for a small enough cluster that the consensus
            reflects hypermutation patterns.""")

    # Apply analysis defaults
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


