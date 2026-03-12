
import os
import unittest

import main_code

class TestDogBreedAnalysis(unittest.TestCase):
    def setUp(self):
        project_root = os.path.dirname(os.path.dirname(__file__))
        database_path = os.path.join(project_root, 'Data', 'dog_breeds.fa')
        self.database = main_code.read_fasta(database_path)
        self.mystery_seq = 'AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC'

    def test_percent_identity(self):
        seq1 = 'AGCTAGCT'
        seq2 = 'AGCTTGCA'
        expected_identity = 75.0
        self.assertAlmostEqual(main_code.percent_identity(seq1, seq2), expected_identity)
        print("[PASS] percent_identity returned expected value (75%)")

    def test_percent_difference(self):
        seq1 = 'AGCTAGCT'
        seq2 = 'AGCTTGCA'
        expected_difference = 25.0
        self.assertAlmostEqual(main_code.percent_difference(seq1, seq2), expected_difference)
        print("[PASS] percent_difference returned expected value (25%)")
