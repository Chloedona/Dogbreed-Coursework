
import os
import unittest

import main_code

class TestDogBreedAnalysis(unittest.TestCase):
    def setUp(self):
        self.database = main_code.parse_fasta('dog_breed_database.fasta')
        self.mystery_seq = 'AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC'

    def test_percent_identity(self):
        seq1 = 'AGCTAGCT'
        seq2 = 'AGCTTGCA'
        expected_identity = 75.0
        self.assertAlmostEqual(main_code.percent_identity(seq1, seq2), expected_identity)
        print("[PASS] percent_identity returned expected value (75%)")
