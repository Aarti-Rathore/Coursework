"""
Created on Tue Mar 19 22:36:28 2024

@author: paria
"""

import unittest
import sys
from Bio.Seq import Seq

script_path = 'C:/Users/paria/biocomputing_prac/coursework BC/code_src'
if script_path not in sys.path:
    sys.path.append(script_path)

from untitled3 import align_sequences, calculate_differences

class MockAlignment:
    def format(self):
       
        return "ATGCGT\nATGCGT"

class TestDNAAlignment(unittest.TestCase):
    def setUp(self):
        
        pass

    def test_align_sequences(self):
        
        seq1 = Seq("ATGCGT")
        seq2 = Seq("ATGCGT")
        alignment = align_sequences(seq1, seq2)
       
        self.assertEqual(alignment.score, 6) 

    def test_calculate_differences(self):
        
        mismatches, insertions, deletions = calculate_differences(MockAlignment())
        self.assertEqual(mismatches, 0)
        self.assertEqual(insertions, 0)
        self.assertEqual(deletions, 0)

   

if __name__ == '__main__':
    unittest.main()
