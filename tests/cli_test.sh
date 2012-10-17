#!/usr/bin/env bash

alias anal='python -m hyperfreq.scripts.cli analyze tests/data/alignment.fasta --out-dir tests/out'

echo "Simple test"
anal --prefix simple
echo "Cluster test"
anal --prefix clustered --cluster-map tests/data/clusters.csv
echo "Ref seqs test"
anal --prefix references --cluster-map tests/data/clusters.csv --reference-sequences tests/data/ref_seqs.fasta \
  --cluster-col cluster_name

alias split='python -m hyperfreq.scripts.cli split tests/data/alignment.fasta --out-dir tests/out'

echo "Simple split"
split tests/out/simple.gross.csv --prefix simple
