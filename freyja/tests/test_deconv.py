import unittest
from freyja.sample_deconv import *
import pandas.testing as pdt
import pandas.api.types as ptypes
from numpy.random import negative_binomial

class DeconvTests(unittest.TestCase):
    def test_buildLineageMap(self):
        mapDict = buildLineageMap()
        self.assertTrue('Alpha'==mapDict['B.1.1.7'] and 'Delta'==mapDict['AY.4'])

    def test_build_mix_and_depth_plus_reindex(self):
        df_barcodes = pd.read_csv('freyja/data/usher_barcodes.csv',index_col=0)
        muts = list(df_barcodes.columns)
        varFn='freyja/data/mixture.tsv'
        depthFn='freyja/data/mixture.depth'
        mix,depths = build_mix_and_depth_arrays(varFn,depthFn,muts)
        #just making sure the files we've read in are ready for use
        self.assertTrue(ptypes.is_float_dtype(mix))
        self.assertTrue(ptypes.is_float_dtype(depths)) 

        df_barcodes,mix,depths = reindex_dfs(df_barcodes,mix,depths)
        pdt.assert_index_equal(df_barcodes.columns,mix.index)
        pdt.assert_index_equal(df_barcodes.columns,depths.index)

        self.assertFalse('20C' in df_barcodes.columns)

    def test_constellation_mapping(self):
        mapDict = buildLineageMap()
        vals = [0.1,0.5,0.31,0.01,0.02,0.01]
        strains = ['B.1.617.2','Q.3','B.1.427','A.2.5','B.1.1','B.1.1.7']
        locDict = map_to_constellation(strains,vals,mapDict)
        locDict_ideal = {'Alpha':0.51,'Epsilon':0.31,'Delta':0.1,'Other':0.02,'Aaron':0.01}
        self.assertTrue(locDict==list(locDict_ideal.items()))

    def test_demixing(self):
        df_barcodes = pd.read_csv('freyja/data/usher_barcodes.csv',index_col=0)
        #build crude in silico mixture from barcodes, perform recovery
        strain1='B.1.1.7'
        strain2='B.1.427'
        mixFracs = [0.4,0.6]
        eps = 0.001 #to check accuracy
        mix =mixFracs[0]*df_barcodes.loc[strain1,]+mixFracs[1]*df_barcodes.loc[strain2,]
        depths = negative_binomial(50,0.25,size=len(mix))
        sample_strains,abundances,error = solve_demixing_problem(df_barcodes,mix,depths)
        self.assertTrue(np.abs(abundances[sample_strains.index(strain1)]-mixFracs[0]<eps))
        self.assertTrue(np.abs(abundances[sample_strains.index(strain2)]-mixFracs[1]<eps))


if __name__ == '__main__':
    unittest.main()
