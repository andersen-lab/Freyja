import unittest
from freyja.sample_deconv import *
import pandas.testing as pdt

class DeconvTests(unittest.TestCase):
    def test_buildLineageMap(self):
        mapDict = buildLineageMap()
        self.assertTrue('Alpha'==mapDict['B.1.1.7'] and 'Delta'==mapDict['AY.4'])

    # def test_build_mix_and_depth(self):
        


if __name__ == '__main__':
    unittest.main()