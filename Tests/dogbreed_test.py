# Import necessary modules for testing
import os # For file path manipulations
import tempfile # For creating temporary directories during tests
import unittest # For creating unit tests
# Import the main code to be tested
import main_code

class TestDogBreedAnalysis(unittest.TestCase):
    def setUp(self):
        project_root = os.path.dirname(os.path.dirname(__file__))
        database_path = os.path.join(project_root, 'Data', 'dog_breeds.fa')
        self.database = main_code.read_fasta(database_path)
        self.mystery_seq = 'AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC'

    # --- read_fasta ---
    def test_read_fasta(self):
        self.assertIsInstance(self.database, dict)
        self.assertGreater(len(self.database), 0)
        for header, seq in self.database.items():
            self.assertIsInstance(seq, str)
            self.assertGreater(len(seq), 0)
        print("[PASS] read_fasta returned a non-empty dictionary of sequences")

    # --- breed_name ---
    def test_breed_name(self):
        header = " >seq1 [breed = boxer]"
        self.assertEqual(main_code.breed_name(header), "boxer")
        print("[PASS] breed_name extracted breed name correctly")

    def test_breed_name_no_match(self):
        header = ">seq1 [organism=Canis lupus]"
        self.assertIsNone(main_code.breed_name(header))
        print("[PASS] breed_name returns None when no breed tag present")

    # --- percent_identity ---
    def test_percent_identity(self):
        seq1 = 'AGCTAGCT'
        seq2 = 'AGCTTGCA'
        self.assertAlmostEqual(main_code.percent_identity(seq1, seq2), 75.0)
        print("[PASS] percent_identity returned expected value (75%)")

    def test_percent_identity_different_lengths(self):
        # 4 matches, max length 8 → 50%
        self.assertAlmostEqual(main_code.percent_identity('AGCT', 'AGCTTTTT'), 50.0)
        print("[PASS] percent_identity handles sequences of different lengths")

    def test_percent_identity_empty(self):
        self.assertEqual(main_code.percent_identity('', ''), 0.0)
        print("[PASS] percent_identity handles empty sequences")

    # --- percent_difference ---
    def test_percent_difference(self):
        self.assertAlmostEqual(main_code.percent_difference('AGCTAGCT', 'AGCTTGCA'), 25.0)
        print("[PASS] percent_difference returned expected value (25%)")

    # --- find_closest ---
    def test_find_closest(self):
        small_db = {
            'seq_a [breed=boxer]': 'AGCTAGCT',
            'seq_b [breed=poodle]': 'TTTTTTTT',
        }
        best_breed, best_similarity = main_code.find_closest('AGCTAGCT', small_db)
        self.assertEqual(best_breed, 'seq_a [breed=boxer]')
        self.assertAlmostEqual(best_similarity, 100.0)
        print("[PASS] find_closest returned the correct closest breed")

    # --- p_values ---
    def test_p_values(self):
        scores = {'a': 90, 'b': 80, 'c': 90}
        pvals = main_code.p_values(scores)
        self.assertAlmostEqual(pvals['a'], 2/3)
        self.assertEqual(pvals['b'], 1.0)
        self.assertAlmostEqual(pvals['c'], 2/3)
        print("[PASS] p_values returned expected values")

    def test_p_values_empty(self):
        self.assertEqual(main_code.p_values({}), {})
        print("[PASS] p_values handles empty input")

    # --- build_tree ---
    def test_build_tree(self):
        sequences = {
            'seq1': 'AGCTAGCT',
            'seq2': 'AGCTTGCA',
            'seq3': 'TTCTAGCT',
        }
        tree = main_code.build_tree(sequences)
        self.assertIsNotNone(tree)
        self.assertEqual(tree.count_terminals(), 3)
        print("[PASS] build_tree returned a valid tree with the correct number of terminals")

    # --- plot_tree_image ---
    def test_plot_tree_image(self):
        sequences = {
            'seq1': 'AGCTAGCT',
            'seq2': 'AGCTTGCA',
            'seq3': 'TTCTAGCT',
        }
        tree = main_code.build_tree(sequences)
        with tempfile.TemporaryDirectory() as tmpdir:
            image_path = main_code.plot_tree_image(tree, tmpdir)
            self.assertTrue(os.path.exists(image_path))
        print("[PASS] plot_tree_image created the image file successfully")

    # --- write_report ---
    def test_write_report(self):
        top_matches = [('Boxer', 95.0, 0.05), ('Poodle', 90.0, 0.10)]
        with tempfile.TemporaryDirectory() as tmpdir:
            main_code.write_report(tmpdir, 'Boxer', 95.0, 0.05, 'tree.png', top_matches)
            report_path = os.path.join(tmpdir, 'report.txt')
            self.assertTrue(os.path.exists(report_path))
            content = open(report_path).read()
            self.assertIn('Boxer', content)
            self.assertIn('95.00', content)
            self.assertIn('0.0500', content)
        print("[PASS] write_report created the report file with expected content")

if __name__ == '__main__':
    unittest.main()
