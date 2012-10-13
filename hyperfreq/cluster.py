import csv

def load_cluster_map(handle, cluster_col='cluster', sequence_col='sequence'):
    reader = csv.DictReader(handle)
    clusters = {}
    for mapping in reader:
        seq = mapping[sequence_col]
        try:
            clusters[mapping[cluster_col]].append(seq)
        except KeyError:
            clusters[mapping[cluster_col]] = [seq]

    return clusters

