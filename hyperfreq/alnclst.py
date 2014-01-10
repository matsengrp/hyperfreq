#!/usr/bin/env python

import itertools
from numpy import mean
from Bio import SeqIO, Align, SeqRecord
from Bio.Align import AlignInfo
import argparse
import csv

# TODO - get r.seq to be automatic in calls to hamming
# Make work with non-overlapping sequences....
# In particular get to work with consensus instead of reps?


def hamming_dist(seq1, seq2):
    """Normalized hamming distance that ignores deletion characters. """
    diffs = 0
    length = 0
    for x, y in zip(str(seq1), str(seq2)):
        if x == '-' or y == '-':
            continue
        elif x != y:
            diffs += 1
        length += 1
    try:
        return float(diffs) / length
    except:
        return 0.5


class Clustering(object):
    def __init__(self, seqrecords, threshold, consensus_threshold):
        """Create a new clustering, using a UCLUST-esque greedy cluster creation mechanism."""
        # iterate through each cluster and find a representative
        self.clusters = list()
        self.threshold = threshold
        self.consensus_threshold = consensus_threshold
        for record in seqrecords:
            try:
                min_dist, clst = min((c.distance(record), c) for c in self.clusters)
            except ValueError:
                min_dist = threshold + 1
            if min_dist < threshold:
                clst.add(record)
            else:
                new_cluster = Cluster(record, consensus_threshold)
                self.clusters.append(new_cluster)

    def merge_small_clusters(self, min_per_cluster):
        """Merges the smallest clusters with nearest neighbours until no cluster is smaller than min_per_cluster"""
        while any(c.size() < min_per_cluster for c in self.clusters):
            _, smallest = min((c.size(), c) for c in self.clusters)
            _, closest = min((c.distance(smallest.centroid), c) for c in self.clusters if c != smallest)
            closest.merge(smallest)
            self.clusters.remove(smallest)

    def mapping_iterator(self):
        """For spitting out to csv, essentially."""
        for i, cluster in enumerate(self.sorted_clusters()):
            for record in cluster.members:
                yield (i, record.name, cluster.distance(record))

    def sorted_clusters(self):
        """Generator for iterating over clusters in order of greatest to smallest"""
        return (c for _, c in sorted((-c.size(), c) for c in self.clusters))

    def recenter_iterator(self):
        """This iterates over cluster members of all clusters. First it yields each of the cluster centroids,
        in order of sorted_clusters. Next it yields each of the next most central members from each cluster,
        again in order of sorted clusters. It does this till there are no more cluster members."""
        cluster_iterators = (c.recenter_iterator() for c in self.sorted_clusters())
        return (r for r in itertools.chain(*itertools.izip_longest(*cluster_iterators)) if r)

    def recenter(self, n):
        """Does the requested number of recenterings."""
        clustering = self
        for i in xrange(n):
            clustering = Clustering(clustering.recenter_iterator(), threshold=clustering.threshold,
                    consensus_threshold=clustering.consensus_threshold)
        return clustering

    def write(self, handle):
        out_writer = csv.writer(handle)
        out_writer.writerow(('cluster_id', 'sequence', 'distance'))
        for row in self.mapping_iterator():
            out_writer.writerow(row)
        handle.close()


class Cluster(object):
    """This class represents a specific cluster. Initally, gets stamped out just with a representative
    sequences and empty sequence list. Sequences have to be added..."""
    def __init__(self, centroid, consensus_threshold):
        self.centroid = centroid
        self.consensus_threshold = consensus_threshold
        self.members = list()
        if self.centroid.name != 'consensus':
            self.members.append(centroid)

    def distance(self, record):
        """Simple hamming distance between rep/center and seq arg."""
        return hamming_dist(self.centroid.seq, record.seq)

    def add(self, record):
        """Add a member to this cluster"""
        if record.name != 'consensus':
            self.members.append(record)

    def merge(self, cluster):
        """Merge two clusters"""
        self.members += cluster.members

    def average_distance(self, record):
        "Average distance to other members of the cluster."
        return mean([hamming_dist(r.seq, record.seq) for r in self.members])

    def recenter(self):
        """Assign a new centroid to the cluster"""
        self.centroid = self.consensus()
        return self.centroid

    def size(self):
        "Number of members..."
        return len(self.members)

    def consensus(self):
        aln = Align.MultipleSeqAlignment(self.members)
        info = AlignInfo.SummaryInfo(aln)
        seq = info.dumb_consensus(threshold=self.consensus_threshold, ambiguous='N')
        return SeqRecord.SeqRecord(seq, name='consensus', id='consensus')

    def recenter_iterator(self):
        """This method iterates over the sequences in the cluster in order of distance from centroid,
        starting from the centroid."""
        # Note that this is where the cluster recenter is happening. Not sure if this is really the best place
        # for it. I guess it doesn't REALLY matter too much as long as it happens before this iterator is
        # called...
        self.recenter()
        yield self.centroid
        for _, record in sorted((self.distance(r), r) for r in self.members):
            yield record


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('alignment')
    parser.add_argument('output', type=argparse.FileType('w'))
    parser.add_argument('-t', '--threshold', type=float, required=True)
    parser.add_argument('-r', '--recenterings', type=int, default=0)
    parser.add_argument('-c', '--consensus-threshold', type=float, default=0.25)
    parser.add_argument('-m', '--min-per-cluster', type=int)
    return parser.parse_args()


def main():
    args = get_args()

    # Grab inputs, sort by ungapped length
    seqrecords = SeqIO.parse(args.alignment, 'fasta')
    seqrecords = (x for _, x in sorted((-len(sr.seq.ungap('-')), sr) for sr in seqrecords))

    # All of the work: cluster, recenter, and merge small clusters as necessary
    clustering = Clustering(seqrecords, args.threshold, args.consensus_threshold)
    print "Initial clustering complete. Starting recentering..."
    clustering = clustering.recenter(args.recenterings)
    print "Recenterings complete."
    if args.min_per_cluster:
        clustering.merge_small_clusters(args.min_per_cluster)
        print "Finished Merging small clusters"

    # Write output
    clustering.write(args.output)



if __name__ == '__main__':
    main()

