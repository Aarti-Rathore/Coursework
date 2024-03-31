# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 16:48:48 2024

@author: paria
"""

from Bio import SeqIO
from Bio.Align import PairwiseAligner

def align_sequences(seq1, seq2):
    """
    Perform alignment between two sequences and return the best alignment.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(seq1, seq2)
    # Assuming we take the first (best) alignment
    best_alignment = alignments[0]
    return best_alignment

def calculate_differences(alignment):
    """
    Calculate mismatches, insertions, and deletions from an alignment.
    """
    # Convert alignment to string format for analysis
    aligned_seq1, aligned_seq2 = alignment.format().split('\n')[:2]
    mismatches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a != b and a != '-' and b != '-')
    insertions = aligned_seq1.count('-')
    deletions = aligned_seq2.count('-')
    return mismatches, insertions, deletions

# Load the test sequence and database sequences
test_seq_record = next(SeqIO.parse('C:/Users/paria/biocomputing_prac/coursework BC/data/mystery.fa', 'fasta'))
database_seq_records = list(SeqIO.parse('C:/Users/paria/biocomputing_prac/coursework BC/data/dog_breeds.fa', 'fasta'))

# Variables to store the best match details
best_match = None
best_alignment = None
highest_score = float('-inf')

# Compare the test sequence against each sequence in the database
for db_seq_record in database_seq_records:
    alignment = align_sequences(test_seq_record.seq, db_seq_record.seq)
    if alignment.score > highest_score:
        highest_score = alignment.score
        best_match = db_seq_record
        best_alignment = alignment

# Analyze the best match
if best_match and best_alignment:
    mismatches, insertions, deletions = calculate_differences(best_alignment)
    print(f"Best match: {best_match.description}")
    print(f"Alignment score: {highest_score}")
    print(f"Mismatches: {mismatches}, Insertions: {insertions}, Deletions: {deletions}")
else:
    print("No match found.")
