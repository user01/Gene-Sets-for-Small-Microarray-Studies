import unittest
import pandas as pd
from .common import score_accuracy, score_recall, score_precision, score_fmeasure
from .common import pairs_to_sets


class TestScores(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame({
            'predicted': ['a','a','d','b','a','a','d'],
            'truth':     ['a','a','a','b','b','b','c'],
        })

    def test_accuracy(self):
        self.assertAlmostEqual(score_accuracy(self.df), 0.4285, 2)

    def test_recall(self):
        self.assertAlmostEqual(score_recall(self.df, 'a'), 0.666, 2)
        self.assertAlmostEqual(score_recall(self.df, 'b'), 0.333, 2)
        self.assertAlmostEqual(score_recall(self.df, 'c'), 0, 2)
        self.assertAlmostEqual(score_recall(self.df, 'd'), 0, 2)

    def test_precision(self):
        self.assertAlmostEqual(score_precision(self.df, 'a'), 0.5, 2)
        self.assertAlmostEqual(score_precision(self.df, 'b'), 1, 2)
        self.assertAlmostEqual(score_precision(self.df, 'c'), 0, 2)
        self.assertAlmostEqual(score_precision(self.df, 'd'), 0, 2)

    def test_fmeasure(self):
        self.assertAlmostEqual(score_fmeasure(self.df, 'a'), 0.571, 2)
        self.assertAlmostEqual(score_fmeasure(self.df, 'b'), 0.5, 2)
        self.assertAlmostEqual(score_fmeasure(self.df, 'c'), 0, 2)
        self.assertAlmostEqual(score_fmeasure(self.df, 'd'), 0, 2)


class TestSetCreation(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame({
            'low':  ['1','2','c','5','6','7','8','a','b'],
            'high': ['2','3','d','6','7','8','9','b','c'],
        })

    def test_basic(self):
        test_sets = pairs_to_sets(self.df, 0, 20)
        print(test_sets)
        self.assertEqual(len(test_sets), 3)
        self.assertEqual(test_sets[0], set(['1','2','3']))
        self.assertEqual(test_sets[1], set(['5','6','7','8','9']))
        self.assertEqual(test_sets[2], set(['a','b','c','d']))

    def test_limits(self):
        test_sets = pairs_to_sets(self.df, 4, 4)
        print(test_sets)
        self.assertEqual(len(test_sets), 1)
        self.assertEqual(test_sets[0], set(['a','b','c','d']))

if __name__ == '__main__':
    unittest.main()
