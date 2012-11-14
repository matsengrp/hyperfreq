#!/usr/bin/env python

import argparse
import csv
from os import path
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from hyperfreq.cluster import load_cluster_map
from hyperfreq.hyperfreq_alignment import HyperfreqAlignment


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


    fn_base = path.join(args.out_dir, args.prefix)
    gross_handle = open(fn_base + '.gross.csv', 'w')
    by_seq_handle = open(fn_base + '.by_seq.csv', 'w')

    # This lets the cluster map be aptional, so that this script can be used
    # for naive hm filtering/analysis
    cluster_map = load_cluster_map(args.cluster_map, cluster_col=args.cluster_col) if args.cluster_map else None
    alignments = HyperfreqAlignment.Set(seq_records, cluster_map, args.clusters,
            reference_sequences=reference_sequences)
    alignments.analyze_hypermuts(args.consensus_threshold, control_trans=args.control_trans,
            pvalue_cutoff=0.05, prob_diff=args.prob_diff)

    alignments.write_analysis(gross_handle, by_seq_handle)

    if args.write_references:
        write_reference_seqs(alignments, fn_base)

    # Closing files
    gross_handle.close()
    by_seq_handle.close()
    args.alignment.close()
    if args.cluster_map:
        args.cluster_map.close()


def setup_common_args(subparser):
    subparser.add_argument('alignment', type=argparse.FileType('r'),
            help="Sequence alignment on which to operate")
    subparser.add_argument('--out-dir', default='.',
            help="Where to put files")
    subparser.add_argument('--prefix',
            help="Prefix for output files (extensions chosen automatically)")


def setup_analyze_args(subparsers):
    def cs_arg(arg):
        return arg.split(',')
    analyze_args = subparsers.add_parser('analyze',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description='Analyze alignment for evidence of hypermuation')
    setup_common_args(analyze_args)
    analyze_args.add_argument('--cluster-map', type=argparse.FileType('r'),
            help="CSV mapping sequences to clusters")
    analyze_args.add_argument('--cluster-col', default='cluster_id',
            help="Column in cluster map file which you would like to use for cluster specification")
    analyze_args.add_argument('--consensus-threshold', type=float,
            help="""As ratio decimal -- used for computation consensus sequences as reference sequences for
            HM evaluation when reference sequences are not explicity specified in --reference-sequences""", default=0.5)
    analyze_args.add_argument('--write-references', default=False, action='store_true',
            help="""Writes to a file the reference sequences (consensus or otherwise) used for HM evalutation
            for each cluster with cluster names as seq names (essentially in the format one would expect for
            --reference-sequences""")
    analyze_args.add_argument('--clusters', type=cs_arg,
            help='csv string - what clusters do you want to use')
    analyze_args.add_argument('--control-trans', default=('C', 'T'), type=tuple, help="Format: 'CT' for C -> T")
    analyze_args.add_argument('--prob-diff', default=0.0, type=float,
            help="Value of X in Sekhon Test of P1 - P2 > X")
    analyze_args.add_argument('--reference-sequences', type=argparse.FileType('r'),
            help="""If specified, use the reference sequences in this file for comparison instead of consensus
            sequences. Sequence name(s) should be the names of the clusters if using a cluster map. Otherwise,
            ensure that there is a reference_sequence named 'all' in your reference-sequences alignment.
            Clusters for which no reference sequence is specified will be compared to a computed consensus
            sequence as reference.""")
    analyze_args.set_defaults(prefix='hyperfreq_analysis')
    analyze_args.set_defaults(func=analyze)


def setup_split_args(subparsers):
    split_args = subparsers.add_parser('split', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    setup_common_args(split_args)
    split_args.add_argument('columns', type=argparse.FileType('r'),
            help="""File identifying hypermuted columns. Typically, you would want to use the `*.gross.csv`
            output file from running `hyperfreq analyze`, but you can contruct your own file if you like.""")
    split_args.add_argument('--column-name', default='column',
            help="""Column name in the columns file which identifies the hypermutated sites.""")
    split_args.set_defaults(prefix='hyperfreq_split')
    split_args.set_defaults(func=split)


def main():
    parser = argparse.ArgumentParser(
            description="""Hypermutation analysis software using Sekhon test""",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', help='additional help')

    setup_analyze_args(subparsers)
    setup_split_args(subparsers)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()


