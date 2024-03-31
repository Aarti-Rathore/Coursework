import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

# Correctly specify the directory path only, without the filename
script_path = 'C:\\Users\\paria\\biocomputing_prac\\coursework BC\\code_src'
sys.path.append(script_path)

# Now you can import the calculate_alignment_scores function and aligner from untitled4.py
from untitled4 import calculate_alignment_scores, aligner

class TestAlignmentScores(unittest.TestCase):
    def setUp(self):
        # Setup test data
        self.test_seq = SeqRecord(Seq("ATGCTAGCTAG"), id="test")
        self.database_seqs = [SeqRecord(Seq("ATGCTAGCTAG"), id="db_seq1"),
                              SeqRecord(Seq("ATGCTAGCGTAG"), id="db_seq2")]
        # Expected scores based on the aligner's configuration in untitled4.py
        self.expected_scores = [11, 10]  # Adjust these values as necessary

    def test_calculate_alignment_scores(self):
        # Test the calculate_alignment_scores function
        scores = calculate_alignment_scores(self.database_seqs, self.test_seq, aligner)
        self.assertEqual(scores, self.expected_scores, "The alignment scores do not match the expected values.")

if __name__ == '__main__':
    unittest.main()
