# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 17:18:42 2024

@author: paria
"""

from Bio import SeqIO
from Bio.Align import PairwiseAligner
import numpy as np

def calculate_alignment_scores(database_seqs, test_seq, aligner):
    """
    Calculate alignment scores for the test sequence against each sequence in the database.
    """
    scores = []
    for db_seq in database_seqs:
        score = aligner.score(test_seq.seq, db_seq.seq)
        scores.append(score)
    return scores

# Initialize the aligner
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -1
aligner.extend_gap_score = -0.5

# Load sequences
test_seq = next(SeqIO.parse('C:/Users/paria/biocomputing_prac/coursework BC/data/mystery.fa', 'fasta'))
database_seqs = list(SeqIO.parse('C:/Users/paria/biocomputing_prac/coursework BC/data/dog_breeds.fa', 'fasta'))

# Calculate alignment scores
scores = calculate_alignment_scores(database_seqs, test_seq, aligner)

# Calculate the p-value for the highest score
max_score = max(scores)
scores_less_than_max = [score for score in scores if score < max_score]
p_value = len(scores_less_than_max) / len(scores)

print(f"Maximum alignment score: {max_score}")
print(f"P-value: {p_value}")
