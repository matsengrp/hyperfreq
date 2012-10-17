import tempfile
from Bio import SeqIO

def fake_file(string):
    tmp_file = tempfile.TemporaryFile()
    tmp_file.write(string.replace(' ', ''))
    tmp_file.seek(0, 0)
    return tmp_file

def parse_fasta(string):
    seqs = SeqIO.to_dict(SeqIO.parse(fake_file(string), 'fasta'))
    return seqs

def parse_fasta_list(string):
    return list(SeqIO.parse(fake_file(string), 'fasta'))

def find_seq(aln, name):
    for seq in aln:
        if seq.name == name:
            return seq


