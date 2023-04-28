import unittest
import subprocess
import os
import pandas as pd
from freyja.read_analysis_tools import extract as _extract, filter as _filter


class ReadAnalysisTests(unittest.TestCase):

    def setUp(self) -> None:
        self.input_bam = 'freyja/data/test.bam'
        self.query_file = 'freyja/data/test_query.csv'
        self.output = 'freyja/data/output.bam'
        self.refname = 'NC_045512.2'

        return super().setUp()

    def test_extract_snps(self):

        snps = 'C75T,G230A,A543C'

        reads = [
            'A01535:8:HJ3YYDSX2:4:1377:5177:13182',
            'A01535:8:HJ3YYDSX2:4:1266:9254:16329',
            'A01535:8:HJ3YYDSX2:4:1264:1642:21183',
            'A01535:8:HJ3YYDSX2:4:1318:15483:28025'
        ]

        with open(self.query_file, 'w') as f:
            f.write(snps)

        reads_found = _extract(self.query_file, self.input_bam, self.output,
                               self.refname, same_read=False)
        reads_found = [x.query_name for x in reads_found]

        self.assertFalse('A01535:8:HJ3YYDSX2:4:1235:12427:11052'
                         in reads_found)

        self.assertTrue(len(reads_found) == 7)
        self.assertTrue(len(set(reads_found)) == 4)  # 3 paired reads

        for read in reads:
            self.assertTrue(read in reads_found)
        os.remove(self.output)

    def test_extract_insertions(self):

        insertions = "(732:'TT'),(1349:'A'),(12333:'A')"

        reads = [
            'A01535:8:HJ3YYDSX2:4:1317:28763:25942',
            'A01535:8:HJ3YYDSX2:4:1141:8341:17566',
            'A01535:8:HJ3YYDSX2:4:1360:23746:26396',
            'A01535:8:HJ3YYDSX2:4:1248:27073:19820'
        ]
        with open(self.query_file, 'w') as f:
            f.write(insertions)

        reads_found = _extract(self.query_file, self.input_bam, self.output,
                               self.refname, same_read=False)
        reads_found = [x.query_name for x in reads_found]
        self.assertTrue('A01535:8:HJ3YYDSX2:4:1158:32669:6809'
                        not in reads_found)
        self.assertTrue(len(reads_found) == 8)
        self.assertTrue(len(set(reads_found)) == 4)  # 4 paired reads

        for read in reads:
            self.assertTrue(read in reads_found)
        os.remove(self.output)

    def test_extract_dels(self):

        deletions = "(1443:32),(1599:2),(2036:3)"

        reads = [
            'A01535:8:HJ3YYDSX2:4:1373:25663:31313',
            'A01535:8:HJ3YYDSX2:4:1112:19090:5713',
            'A01535:8:HJ3YYDSX2:4:1429:11550:34178',
        ]

        with open(self.query_file, 'w') as f:
            f.write(deletions)

        reads_found = _extract(self.query_file, self.input_bam, self.output,
                               self.refname, same_read=False)
        reads_found = [x.query_name for x in reads_found]
        self.assertFalse('A01535:8:HJ3YYDSX2:4:1123:4707:5165'
                         in reads_found)
        self.assertTrue(len(reads_found) == 6)
        self.assertTrue(len(set(reads_found)) == 3)  # 3 paired reads

        for read in reads:
            self.assertTrue(read in reads_found)
        os.remove(self.output)

    def test_filter_snps(self):

        snps = 'C75T,G230A,A543C'

        reads = [
            'A01535:8:HJ3YYDSX2:4:1377:5177:13182',
            'A01535:8:HJ3YYDSX2:4:1266:9254:16329',
            'A01535:8:HJ3YYDSX2:4:1264:1642:21183',
            'A01535:8:HJ3YYDSX2:4:1318:15483:28025'
        ]

        with open(self.query_file, 'w') as f:
            f.write(snps)

        reads_found = _filter(self.query_file, self.input_bam,
                              75, 600, self.output, self.refname)
        self.assertTrue('A01535:8:HJ3YYDSX2:4:1235:12427:11052'
                        in reads_found)
        for read in reads:
            self.assertFalse(read in reads_found)
        os.remove(self.output)

    def test_filter_insertions(self):

        insertions = "(732:'TT'),(1349:'A'),(12333:'A')"

        reads = [
            'A01535:8:HJ3YYDSX2:4:1317:28763:25942',
            'A01535:8:HJ3YYDSX2:4:1141:8341:17566',
            'A01535:8:HJ3YYDSX2:4:1360:23746:26396',
            'A01535:8:HJ3YYDSX2:4:1248:27073:19820'
        ]

        with open(self.query_file, 'w') as f:
            f.write(insertions)

        reads_found = _filter(self.query_file, self.input_bam, 730, 12340,
                              self.output, self.refname)

        self.assertTrue('A01535:8:HJ3YYDSX2:4:1158:32669:6809'
                        in reads_found)
        for read in reads:
            self.assertFalse(read in reads_found)
        os.remove(self.output)

    def test_filter_deletions(self):

        deletions = "(1443:32),(1599:2),(2036:3)"

        reads = [
            'A01535:8:HJ3YYDSX2:4:1373:25663:31313',
            'A01535:8:HJ3YYDSX2:4:1112:19090:5713',
            'A01535:8:HJ3YYDSX2:4:1429:11550:34178',
        ]

        with open(self.query_file, 'w') as f:
            f.write(deletions)

        reads_found = _filter(self.query_file, self.input_bam, 1400, 2100,
                              self.output, self.refname)

        self.assertTrue('A01535:8:HJ3YYDSX2:4:1123:4707:5165'
                        in reads_found)

        for read in reads:
            self.assertFalse(read in reads_found)
        os.remove(self.output)

    def test_covariants(self):
        cmd = ['freyja', 'covariants', self.input_bam,
               '21563', '25384', '--gff-file',
               'freyja/data/NC_045512_Hu-1.gff', '--output',
               'freyja/data/test_covar.tsv']
        subprocess.run(cmd)
        df = pd.read_csv('freyja/data/test_covar.tsv', sep='\t')

        patterns = []
        for c in df.iloc[:, 0]:
            patterns.append(c.split(' '))

        # Test for covariants spanning > 150bp
        self.assertTrue(['A23403G(S:D614G)', 'C23604G(S:P681R)'] in patterns)

        # Test coverage ranges
        cov_start = df.iloc[:, 2]
        cov_end = df.iloc[:, 3]

        for i in range(len(patterns)):
            for mut in patterns[i]:
                self.assertTrue(int(mut[1:6]) in range(
                    cov_start[i], cov_end[i]))

        # Test spans_region flag
        cmd = ['freyja', 'covariants', self.input_bam,
               '23550', '23700', '--output', 'freyja/data/test_covar.tsv',
               '--spans_region']
        subprocess.run(cmd)
        df = pd.read_csv('freyja/data/test_covar.tsv', sep='\t')
        patterns = []
        for c in df.iloc[:, 0]:
            patterns.append(c.split(' '))
        cov_start = df.iloc[:, 2]
        cov_end = df.iloc[:, 3]

        for i in range(len(patterns)):
            for mut in patterns[i]:
                self.assertTrue(int(mut[1:6]) >= cov_start[i] and
                                int(mut[1:6]) <= cov_end[i])
        self.assertTrue(df.shape[0] == 2)
        os.remove('freyja/data/test_covar.tsv')

    def test_plot_covariants(self):
        cmd = ['freyja', 'plot-covariants',
               'freyja/data/example_covariants0.tsv',
               '--output', 'freyja/data/test_covar_plot.png']
        subprocess.run(cmd)

        self.assertTrue('test_covar_plot.png' in os.listdir('freyja/data'))
        os.remove('freyja/data/test_covar_plot.png')


if __name__ == '__main__':
    unittest.main()
