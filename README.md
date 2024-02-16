# coursework
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 10:25:24 2024

@author: paria
"""

from Bio import SeqIO
from Bio.SeqUtils import nt_search

def remove_gaps(sequence):
    """Remove gap characters from a sequence."""
    return sequence.replace('-', '')

def parse_fasta(file_path):
    """Parse a FASTA file and return a list of SeqRecord objects."""
    return list(SeqIO.parse(file_path, "fasta"))

def simple_sequence_comparison(seq1, seq2):
    """A simplified comparison function that searches for seq1 within seq2 after removing gaps."""
    # Remove gaps from both sequences
    seq1_no_gaps = remove_gaps(str(seq1))
    seq2_no_gaps = remove_gaps(str(seq2))
    # Perform the search using nt_search from Bio.SeqUtils
    result = nt_search(seq2_no_gaps, seq1_no_gaps)
    # The result includes the pattern itself and the start positions of any matches; subtract 1 to exclude the pattern itself
    return len(result) - 1

def find_closest_match(database_sequences, test_sequence):
    """Find the database sequence that is the closest match to the test sequence."""
    best_match = None
    best_score = -1
    for db_seq in database_sequences:
        score = simple_sequence_comparison(test_sequence.seq, db_seq.seq)
        if score > best_score:
            best_match = db_seq
            best_score = score
    return best_match


database_sequences = parse_fasta('C:/Users/paria/biocomputing_prac/coursework BC/dog_breeds.fa')
test_sequences = parse_fasta('C:/Users/paria/biocomputing_prac/coursework BC/mystery.fa')


test_sequence = test_sequences[0]


closest_match = find_closest_match(database_sequences, test_sequence)

if closest_match:
    print(f"Closest match: {closest_match.description}")
else:
    print("No match found.")
