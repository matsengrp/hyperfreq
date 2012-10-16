#!/usr/bin/env bash

alias anal='python -m hyperfreq.scripts.cli analyze test/data/alignment.fasta --out-dir test/out'

echo "Simple test"
anal --prefix simple
echo "Cluster test"
anal --prefix clustered --cluster-map test/data/clusters.csv
echo "Ref seqs test"
anal --prefix references --cluster-map test/data/clusters.csv --reference-sequences test/data/ref_seqs.fasta \
  --cluster-col cluster_name

alias split='python -m hyperfreq.scripts.cli split test/data/alignment.fasta --out-dir test/out'

echo "Simple split"
split test/out/simple.gross.csv --prefix simple
