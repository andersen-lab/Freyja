import unittest
import pandas as pd
from freyja.convert_paths2barcodes import parse_tree_paths, sortFun,\
    convert_to_barcodes, reversion_checking
import pandas.testing as pdt


class BarcodeTests(unittest.TestCase):

    def test_parse_tree_paths(self):
        df_test = pd.DataFrame({'clade': ['T0', 'T1', 'T2'],
                                'from_tree_root': ['> T1234G > A122C',
                                                   None, '> A15T > A122C']})
        df_test = parse_tree_paths(df_test)
        df_ideal = pd.DataFrame({'clade': ['T0', 'T1', 'T2'],
                                 'from_tree_root': [['T1234G', 'A122C'],
                                                    [''], ['A15T', 'A122C']]})
        df_ideal = df_ideal.set_index('clade')
        pdt.assert_frame_equal(df_test, df_ideal)

    def test_sortFun(self):
        test_string = 'A1234T'
        self.assertTrue(1234 == sortFun(test_string))

    def test_convert_to_barcodes(self):
        df_test = pd.DataFrame({'clade': ['T0', 'T1', 'T2'],
                                'from_tree_root': [['T1234G', 'A122C'],
                                                   [''], ['A15T', 'A122C']]})
        df_barcode_test = convert_to_barcodes(df_test)
        df_barcode_ideal = pd.DataFrame({'A122C': [1., 0., 1.],
                                         'A15T': [0., 0., 1.],
                                        'T1234G': [1., 0., 0.]})
        pdt.assert_frame_equal(df_barcode_test, df_barcode_ideal)

    def test_reversion_checking(self):
        df_barcode = pd.DataFrame({'A122C': [1., 0., 1.],
                                   'C122A': [0., 0., 1.],
                                   'T1234G': [1., 0., 0.],
                                   'T12346G': [0., 0., 0.]})
        df_barcode_checked = reversion_checking(df_barcode)
        df_barcode_checked_ideal = pd.DataFrame({'A122C': [1., 0., 0.],
                                                 'T1234G': [1., 0., 0.]})
        pdt.assert_frame_equal(df_barcode_checked, df_barcode_checked_ideal)

    def test_no_flip_pairs(self):
        df_barcodes = pd.read_csv('freyja/data/usher_barcodes.csv',
                                  index_col=0)
        flipPairs = [(d, d[-1] + d[1:len(d)-1]+d[0])
                     for d in df_barcodes.columns
                     if (d[-1] + d[1:len(d)-1]+d[0]) in df_barcodes.columns]
        self.assertTrue(len(flipPairs) == 0)


if __name__ == '__main__':
    unittest.main()
