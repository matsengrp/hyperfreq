import csv

def load_cluster_map(handle, cluster_col='cluster', sequence_col='sequence'):
    reader = csv.DictReader(handle)
    return parse_clusters(reader, cluster_col, sequence_col)


def parse_clusters(mappings, cluster_key='cluster', sequence_key='sequence'):
    clusters = {}
    for mapping in mappings:
        seq = mapping[sequence_key]
        try:
            clusters[mapping[cluster_key]].append(seq)
        except KeyError:
            clusters[mapping[cluster_key]] = [seq]

    return clusters

