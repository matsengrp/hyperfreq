#!/usr/bin/env sh

set -e

anal -h
echo "Simple test"
anal --prefix simple -q 0.5
echo "Cluster test"
anal --prefix clustered --cluster-map tests/data/clusters.csv
echo "Testing ref out"
anal --prefix ref_out --cluster-map tests/data/clusters.csv --write-references
echo "Ref seqs test"
anal --prefix references --cluster-map tests/data/clusters.csv --reference-sequences tests/data/ref_seqs.fasta \
  --cluster-col cluster
echo "Ref seqs test (with writing)"
anal --prefix references_ww --cluster-map tests/data/clusters.csv --reference-sequences tests/data/ref_seqs.fasta \
  --cluster-col cluster --write-references

alias split='python -m hyperfreq.scripts.cli split tests/data/alignment.fasta --out-dir tests/out'

echo "Simple split"
split tests/out/simple.sites.csv --prefix simple
