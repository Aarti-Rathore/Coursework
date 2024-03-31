# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 18:28:11 2024

@author: paria
"""



from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import Bio.Phylo as Phylo

# Assuming these are the correct paths to your files
dog_breeds_path = "C:/Users/paria/biocomputing_prac/coursework BC/data/dog_breeds.fa"
mystery_path = "C:/Users/paria/biocomputing_prac/coursework BC/data/mystery.fa"

# Read sequences from the dog breeds file and the mystery sequence
dog_breeds_sequences = list(SeqIO.parse(dog_breeds_path, "fasta"))
mystery_sequence = list(SeqIO.parse(mystery_path, "fasta"))[0]  # Assuming one mystery sequence

# Combine all sequences for analysis
sequences = dog_breeds_sequences + [mystery_sequence]

# Initialize the PairwiseAligner with settings for global alignment
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -1
aligner.extend_gap_score = -0.5

# Calculate pairwise distances
names = [seq.id for seq in sequences]
distances = [[0.0 for _ in sequences] for _ in sequences]

for i, seq1 in enumerate(sequences):
    for j, seq2 in enumerate(sequences):
        if i < j:
            alignment_score = aligner.score(seq1.seq, seq2.seq)
            distance = 1 - alignment_score / max(len(seq1.seq), len(seq2.seq))
            distances[i][j] = distances[j][i] = distance

# Construct a distance matrix
distance_matrix = DistanceMatrix(names, distances)

# Construct the phylogenetic tree using the UPGMA algorithm
constructor = DistanceTreeConstructor()
tree = constructor.upgma(distance_matrix)

# Save the tree to a file in PhyloXML format for visualization
Phylo.write(tree, "phylogenetic_tree.xml", "phyloxml")

# Optionally, visualize the tree directly if you have matplotlib installed
# Uncomment the following line to visualize:
# Phylo.draw(tree)
