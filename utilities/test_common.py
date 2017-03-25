import unittest
import pandas as pd
from common import score_accuracy, score_recall, score_precision


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


if __name__ == '__main__':
    unittest.main()
