import unittest
import pandas as pd
from freyja.utils import agg, prepLineageDict, prepSummaryDict


class UtilsTests(unittest.TestCase):
    def test_agg(self):
        agg_df = agg('freyja/data/outputs/')
        # simple checks, since entries are still not parsed
        self.assertTrue('mixture.tsv' in agg_df.index)
        self.assertTrue('summarized' in agg_df.columns)

    def test_prep(self):
        agg_df = pd.read_csv('freyja/data/agg_outputs.tsv',
                             skipinitialspace=True,
                             sep='\t', index_col=0)
        agg_df['linDict'] = prepLineageDict(agg_df)
        agg_df['summarized'] = prepSummaryDict(agg_df)
        self.assertTrue((agg_df['linDict'][0]['A'] > 0) &
                        (agg_df['linDict'][1]['A'] > 0))
        self.assertTrue((agg_df['summarized'][0]['A'] > 0) &
                        (agg_df['summarized'][1]['A'] > 0))
        self.assertTrue(agg_df['summarized'][1]['Delta'] > 0)


if __name__ == '__main__':
    unittest.main()
