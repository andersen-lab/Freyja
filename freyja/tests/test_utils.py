import unittest
import pandas as pd
from freyja.utils import agg, checkConfig, prepLineageDict, prepSummaryDict,\
    get_color_scheme, get_abundance
import os
import plotly.express as px
import yaml

df_ab_lin = pd.DataFrame.from_dict({
        "Q.3": {
            "2021-03-01": 99.86929930000001,
            "2021-03-03": 93.9596733,
            "2021-03-08": 89.77044079999999,
            "2021-03-10": 84.1962053,
            "2021-03-12": 64.68600500000001,
            "2021-03-14": 58.985296,
            "2021-03-17": 60.4915979,
            "2021-03-20": 54.0080433,
            "2021-03-25": 36.943905300000004,
            "2021-03-30": 20.987667199999997,
            "2021-03-31": 11.768960100000001,
            "2021-04-04": 0
        },
        "AY.48": {
            "2021-03-01": 0,
            "2021-03-03": 4.99648,
            "2021-03-08": 8.010876,
            "2021-03-10": 13.844148,
            "2021-03-12": 33.4055237,
            "2021-03-14": 36.6772838,
            "2021-03-17": 35.1807163,
            "2021-03-20": 41.6008049,
            "2021-03-25": 60.344586,
            "2021-03-30": 73.46381410000001,
            "2021-03-31": 87.07396659999999,
            "2021-04-04": 99.807238
        },
        "B.1.617.2": {
            "2021-03-01": 0,
            "2021-03-03": 0,
            "2021-03-08": 0,
            "2021-03-10": 0,
            "2021-03-12": 0,
            "2021-03-14": 3.37115719,
            "2021-03-17": 3.2468312,
            "2021-03-20": 2.9356701,
            "2021-03-25": 1.43804146,
            "2021-03-30": 2.22985247,
            "2021-03-31": 0,
            "2021-04-04": 0
        },
        "Other": {
            "2021-03-01": 0.13070067750000003,
            "2021-03-03": 1.0438466700499998,
            "2021-03-08": 2.2186831967230005,
            "2021-03-10": 1.9596465105119993,
            "2021-03-12": 1.9084712976769997,
            "2021-03-14": 0.966262983661,
            "2021-03-17": 1.080854620274,
            "2021-03-20": 1.4554817049518003,
            "2021-03-25": 1.27346720046,
            "2021-03-30": 3.3186662099999995,
            "2021-03-31": 1.1570733055280005,
            "2021-04-04": 0.19276375983000002
        }
    })
df_ab_sum = pd.DataFrame.from_dict({
        "Alpha": {
            "2021-03-01": 99.87789999982671,
            "2021-03-03": 93.98129492955486,
            "2021-03-08": 89.9147441699543,
            "2021-03-10": 84.22558680720658,
            "2021-03-12": 64.70828169148619,
            "2021-03-14": 59.031199999535666,
            "2021-03-17": 60.52881597979205,
            "2021-03-20": 54.05737427761939,
            "2021-03-25": 36.95459482299681,
            "2021-03-30": 21.04474169994245,
            "2021-03-31": 11.78529999719067,
            "2021-04-04": 0
        },
        "Delta": {
            "2021-03-01": 0.014058802055376244,
            "2021-03-03": 5.177708342651904,
            "2021-03-08": 9.155290001286307,
            "2021-03-10": 14.964234802582208,
            "2021-03-12": 33.61162689800541,
            "2021-03-14": 40.27064950666095,
            "2021-03-17": 38.551629204685014,
            "2021-03-20": 44.72319999997791,
            "2021-03-25": 62.0014336019168,
            "2021-03-30": 75.85518917199204,
            "2021-03-31": 87.5540522317538,
            "2021-04-04": 99.85259954102992
        },
        "Other": {
            "2021-03-01": 0.10804120463362901,
            "2021-03-03": 0.8409967276198742,
            "2021-03-08": 0.9299658256744001,
            "2021-03-10": 0.8101782067119108,
            "2021-03-12": 1.680091400472183,
            "2021-03-14": 0.6981504842771986,
            "2021-03-17": 0.9195548174455725,
            "2021-03-20": 1.2194257038786276,
            "2021-03-25": 1.0439715589713492,
            "2021-03-30": 3.1000691040175408,
            "2021-03-31": 0.6606477742925074,
            "2021-04-04": 0.14740221471037052
        }
    })
dates_to_keep = ['2021-03-01 00:00:00',
                 '2021-03-03 00:00:00',
                 '2021-03-08 00:00:00',
                 '2021-03-10 00:00:00',
                 '2021-03-12 00:00:00',
                 '2021-03-14 00:00:00',
                 '2021-03-17 00:00:00',
                 '2021-03-20 00:00:00',
                 '2021-03-25 00:00:00',
                 '2021-03-30 00:00:00',
                 '2021-03-31 00:00:00',
                 '2021-04-04 00:00:00']


class UtilsTests(unittest.TestCase):
    def setUp(self):
        agg_results = 'freyja/data/test_sweep.tsv'
        metadata = 'freyja/data/sweep_metadata.csv'
        config_file = 'freyja/data/plot_config.yml'
        self.agg_df = pd.read_csv(agg_results, skipinitialspace=True, sep='\t',
                                  index_col=0)
        self.agg_df = self.agg_df[self.agg_df['summarized'] != '[]']
        self.meta_df = pd.read_csv(metadata, index_col=0)
        self.meta_df['sample_collection_datetime'] = \
            pd.to_datetime(self.meta_df['sample_collection_datetime'])
        with open('freyja/data/lineages.yml', 'r') as f:
            try:
                lineages_yml = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                raise ValueError('lineages.yml error, run update: ' + str(exc))
        with open(config_file, "r") as f:
            try:
                config = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                raise ValueError('Error in config file: ' + str(exc))
        if config is not None:
            self.config = checkConfig(config)
        # converts lineages_yml to a dictionary where the lineage names are the
        # keys.
        self.lineage_info = {}
        for lineage in lineages_yml:
            self.lineage_info[lineage['name']] = {
                    'name': lineage['name'],
                    'children': lineage['children']
                }

    def test_agg(self):
        results = 'freyja/data/outputs/'
        agg_df = agg([results + fn
                      for fn in os.listdir('freyja/data/outputs/')])
        # simple checks, since entries are still not parsed
        self.assertTrue('mixture.tsv' in agg_df.index)
        self.assertTrue('summarized' in agg_df.columns)

    def test_prep(self):
        agg_df = pd.read_csv('freyja/data/agg_outputs.tsv',
                             skipinitialspace=True,
                             sep='\t', index_col=0)
        agg_df = prepLineageDict(agg_df, thresh=0.001)
        agg_df = prepSummaryDict(agg_df)
        # flexibility to python versions
        if isinstance(agg_df['linDict'][0], list):
            self.assertTrue((agg_df['linDict'][0][0]['A'] > 0) &
                            (agg_df['linDict'][1][0]['A'] > 0))
        else:
            self.assertTrue((agg_df['linDict'][0]['A'] > 0) &
                            (agg_df['linDict'][1]['A'] > 0))
        self.assertTrue((agg_df['summarized'][0]['A'] > 0) &
                        (agg_df['summarized'][1]['A'] > 0))
        self.assertTrue(agg_df['summarized'][1]['Delta'] > 0)

    def test_checkConfig(self):
        # config dictionaries whose format is valid.
        valid_configs = [
            {
                'Lineages': {
                    'grp_1': {
                        'name': 'grp_1',
                        'members': ['B.1.617*', 'B.1.252'],
                        'color': 'orange'
                        },
                    'grp_2': {
                        'name': 'grp_2',
                        'members': ['Q.3'],
                        'color': 'green'
                        }
                    },
                'VOC': {
                    'Delta': {
                        'name': 'Delta',
                        'color': 'default'
                        }
                    }
            },
            {
                'Lineages': {
                    'grp_1': {
                        'name': 'grp_1',
                        'members': ['B.1.617*', 'B.1.252'],
                        'color': 'orange'
                        },
                    'grp_2': {
                        'name': 'grp_2',
                        'members': ['Q.3'],
                        'color': 'green'
                        }
                    }
            },
            {
                'VOC': {
                    'Delta': {
                        'name': 'Delta',
                        'color': 'default'
                        }
                    }
            },
        ]

        output_valid_configs = [
            {
                'Lineages': {
                    'grp_1': {
                        'name': 'grp_1',
                        'members': ['B.1.617*', 'B.1.252'],
                        'color': 'orange'
                        },
                    'grp_2': {
                        'name': 'grp_2',
                        'members': ['Q.3'],
                        'color': 'green'
                        }
                    },
                'VOC': {
                    'Delta': {
                        'name': 'Delta',
                        'color': 'default'
                        }
                    }
            },
            {
                'Lineages': {
                    'grp_1': {
                        'name': 'grp_1',
                        'members': ['B.1.617*', 'B.1.252'],
                        'color': 'orange'
                        },
                    'grp_2': {
                        'name': 'grp_2',
                        'members': ['Q.3'],
                        'color': 'green'
                        }
                    },
                'VOC': {}
            },
            {
                'Lineages': {},
                'VOC': {
                    'Delta': {
                        'name': 'Delta',
                        'color': 'default'
                        }
                    }
            },
        ]

        # config dictionaries whose format is invalid.
        invalid_configs = [
            {
                'Lineages': {
                    'grp_1': {
                        'name': 'grp_1',
                        'members': ['B.1.617*', 'B.1.252'],
                        },
                    'grp_2': {
                        'name': 'grp_2',
                        'members': ['Q.3'],
                        'color': 'green'
                        }
                    },
                'VOC': {
                    'Delta': {
                        'name': 'Delta',
                        'color': 'default'
                        }
                    }
            },
            {
                'Lineages': {
                    'grp_1': {
                        'name': 'grp_1',
                        'members': ['B.1.617*', 'B.1.252'],
                        'color': 'orange'
                        },
                    'grp_2': {
                        'name': 'grp_2',
                        'color': 'green'
                        }
                    }
            },
            {
                'VOC': {
                    'Delta': {
                        'color': 'default'
                        }
                    }
            },
        ]

        #  Test valid configs
        for i, config in enumerate(valid_configs):
            self.assertEqual(checkConfig(config), output_valid_configs[i])

        # Test invalid configs
        for i, config in enumerate(invalid_configs):
            self.assertRaises(ValueError, checkConfig, config)

    def test_get_abundance(self):
        df_ab_lin.index = pd.to_datetime(df_ab_lin.index)
        df_ab_sum.index = pd.to_datetime(df_ab_sum.index)
        lin, sum, dates = get_abundance(
                self.agg_df,
                self.meta_df,
                0.01,
                False,
                {},
                self.lineage_info)
        pd.testing.assert_frame_equal(lin, df_ab_lin)
        pd.testing.assert_frame_equal(sum, df_ab_sum)
        assert type(dates) is pd.DatetimeIndex
        dates = [str(i) for i in list(dates)]
        self.assertListEqual(dates, dates_to_keep)

    def test_get_color_scheme(self):
        df_ab_lin.index = pd.to_datetime(df_ab_lin.index)
        default_color_scheme = {
            11: px.colors.qualitative.Vivid,
            24: px.colors.qualitative.Dark24
        }
        color_scheme_without_config = {
            'Q.3': 'rgb(229, 134, 6)',
            'AY.48': 'rgb(93, 105, 177)',
            'B.1.617.2': 'rgb(82, 188, 163)',
            'Other': 'rgb(153, 201, 69)'
        }

        self.assertDictEqual(get_color_scheme(
                df_ab_lin,
                default_color_scheme,
                {}
            ), color_scheme_without_config)

        df_ab_lin.rename(columns={'Q.3': 'grp_1', 'AY.48': 'grp_2'},
                         inplace=True)
        df_ab_lin.drop(columns=['B.1.617.2'], inplace=True)
        color_scheme_with_config = {
            'grp_1': 'orange',
            'grp_2': 'green',
            'Other': 'rgb(82, 188, 163)'
        }

        self.assertDictEqual(get_color_scheme(
                df_ab_lin,
                default_color_scheme,
                self.config.get('Lineages')
            ), color_scheme_with_config)


if __name__ == '__main__':
    unittest.main()
