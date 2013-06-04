#!/usr/bin/env sh

set -e

alias anal='python -m hyperfreq.scripts.cli analyze tests/data/alignment.fasta --out-dir tests/out'

echo "Simple test"
anal --prefix simple --pos-quants-only
echo "Cluster test"
anal --prefix clustered --cluster-map tests/data/clusters.csv
echo "Testing ref out"
anal --prefix ref_out --cluster-map tests/data/clusters.csv --write-references
echo "Ref seqs test"
anal --prefix references --cluster-map tests/data/clusters.csv --reference-sequences tests/data/ref_seqs.fasta \
  --cluster-col cluster_name
echo "Ref seqs test (with writing)"
anal --prefix references_ww --cluster-map tests/data/clusters.csv --reference-sequences tests/data/ref_seqs.fasta \
  --cluster-col cluster_name --write-references

alias split='python -m hyperfreq.scripts.cli split tests/data/alignment.fasta --out-dir tests/out'

echo "Simple split"
split tests/out/simple.gross.csv --prefix simple
