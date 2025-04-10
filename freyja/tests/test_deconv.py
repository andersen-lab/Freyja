import unittest
import pandas as pd
from freyja.sample_deconv import buildLineageMap, build_mix_and_depth_arrays, \
     reindex_dfs, map_to_constellation, solve_demixing_problem, \
     perform_bootstrap
import pandas.testing as pdt
import pandas.api.types as ptypes
from numpy.random import negative_binomial


class DeconvTests(unittest.TestCase):
    def test_buildLineageMap(self):
        mapDict = buildLineageMap('')
        self.assertTrue('Alpha' == mapDict['B.1.1.7'])
        self.assertTrue('Delta' == mapDict['AY.4'])

    def test_build_mix_and_depth_plus_reindex(self):
        df_barcodes = pd.read_csv('freyja/data/usher_barcodes.csv',
                                  index_col=0)
        muts = list(df_barcodes.columns)
        varFns = ['freyja/data/mixture.tsv', 'freyja/data/mixture.vcf']
        depthFn = 'freyja/data/mixture.depth'
        covcut = 10
        autoadapt = False

        for varFn in varFns:
            mix, depths, cov, adapt = build_mix_and_depth_arrays(
                varFn, depthFn, muts, covcut, autoadapt)
            self.assertTrue(ptypes.is_float_dtype(mix))
            self.assertTrue(ptypes.is_float_dtype(depths))
            df_barcodes_reindexed, mix, depths = reindex_dfs(
                df_barcodes, mix, depths)
            pdt.assert_index_equal(df_barcodes_reindexed.columns, mix.index)
            pdt.assert_index_equal(df_barcodes_reindexed.columns, depths.index)
            self.assertFalse('20C' in df_barcodes.columns)

    def test_constellation_mapping(self):
        mapDict = buildLineageMap('')
        vals = [0.1, 0.5, 0.31, 0.01, 0.02, 0.01]
        strains = ['B.1.617.2', 'Q.3', 'B.1.427', 'A.2.5', 'B.1.1', 'B.1.1.7']
        locDict = map_to_constellation(strains, vals, mapDict)
        locDict_ideal = {'Alpha': 0.51, 'Epsilon': 0.31, 'Delta': 0.1,
                         'Other': 0.02, 'A': 0.01}
        self.assertTrue(locDict == list(locDict_ideal.items()))

    def test_demixing(self):
        df_barcodes = pd.read_csv('freyja/data/usher_barcodes.csv',
                                  index_col=0)
        # build crude in silico mixture from barcodes, perform recovery
        strain1 = 'B.1.1.7'
        strain2 = 'B.1.427'
        mixFracs = [0.4, 0.6]
        mix = mixFracs[0]*df_barcodes.loc[strain1, ]\
            + mixFracs[1]*df_barcodes.loc[strain2, ]
        # generate random sequencing depth at each position
        depths = negative_binomial(50, 0.25, size=len(mix))
        eps = 0.001
        solver = 'CLARABEL'
        sample_strains, abundances, error = solve_demixing_problem(df_barcodes,
                                                                   mix, depths,
                                                                   eps, 0.,
                                                                   1.E-8,
                                                                   solver)
        self.assertAlmostEqual(
            abundances[sample_strains.tolist().index(strain1)], mixFracs[0])
        self.assertAlmostEqual(
            abundances[sample_strains.tolist().index(strain2)], mixFracs[1])

    def test_boot(self):
        mapDict = buildLineageMap('')
        df_barcodes = pd.read_csv('freyja/data/usher_barcodes.csv',
                                  index_col=0)
        muts = list(df_barcodes.columns)
        # build crude in silico mixture from barcodes, perform recovery
        strain1 = 'B.1.1.7'
        strain2 = 'B.1.427'
        mixFracs = [0.4, 0.6]
        mix = mixFracs[0]*df_barcodes.loc[strain1, ]\
            + mixFracs[1]*df_barcodes.loc[strain2, ]
        # generate random sequencing depth at each position
        depths = pd.Series(negative_binomial(50, 0.25, size=len(mix)),
                           index=mix.index)
        eps = 0.001
        numBootstraps = 3
        n_jobs = 1
        bootseed = 10
        solver = 'CLARABEL'
        lin_df, constell_df = perform_bootstrap(df_barcodes, mix, depths,
                                                numBootstraps, eps,
                                                n_jobs, mapDict, muts,
                                                '', 'test', bootseed,
                                                solver)
        lin_out = lin_df.quantile([0.5])
        constell_out = constell_df.quantile([0.5])
        self.assertAlmostEqual(lin_out.loc[0.5, 'B.1.1.7'], 0.4, delta=0.1)
        self.assertAlmostEqual(constell_out.loc[0.5, 'Alpha'], 0.4, delta=0.1)


if __name__ == '__main__':
    unittest.main()
