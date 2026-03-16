
import os
import unittest

import main_code

class TestDogBreedAnalysis(unittest.TestCase):
    def setUp(self):
        project_root = os.path.dirname(os.path.dirname(__file__))
        database_path = os.path.join(project_root, 'Data', 'dog_breeds.fa')
        self.database = main_code.read_fasta(database_path)
        self.mystery_seq = 'AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC'

    def test_breed_name(self):
        header = " >seq1 [breed = boxer]"
        expected_breed = "boxer"
        self.assertEqual(main_code.breed_name(header), expected_breed)
        print("[PASS] breed_name extracted breed name correctly")


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
    
    def test_find_closest(self):
        best_breed, best_similarity = main_code.find_closest(self.mystery_seq, self.database)
        self.assertIsNotNone(best_breed)
        self.assertGreaterEqual(best_similarity, 0.0)
        self.assertLessEqual(best_similarity, 100.0)
        print("[PASS] find_closest returned a valid breed and similarity score")

    def test_p_values(self):
        scores = {'a': 90, 'b': 80, 'c': 90}
        pvals = main_code.p_values(scores)

        self.assertEqual(pvals['a'], 2/3)
        self.assertEqual(pvals['b'], 1.0)
        self.assertEqual(pvals['c'], 2/3)

    

