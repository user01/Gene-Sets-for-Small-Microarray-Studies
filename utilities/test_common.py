import unittest
import pandas as pd
from common import score_accuracy


class TestScores(unittest.TestCase):

    def test_accuracy(self):
        df = pd.DataFrame({
            'predicted': ['a','a','c','b','a','c'],
            'truth':     ['a','a','a','b','b','b'],
        })
        self.assertAlmostEqual(score_accuracy(df), 0.5, 2)

if __name__ == '__main__':
    unittest.main()
