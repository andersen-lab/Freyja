import unittest
from freyja.read_analysis_tools import filter as _filter

class FilterTests(unittest.TestCase):

    def test_snps(self):
        query_file = 'freyja/data/test_filter.csv'
        input_bam = 'freyja/data/test.bam'
        output = 'freyja/data/outputs'

        snps = 'C75T,G230A,A543C'
        reads = [
            'A01535:8:HJ3YYDSX2:4:1377:5177:13182',
            'A01535:8:HJ3YYDSX2:4:1266:9254:16329',
            'A01535:8:HJ3YYDSX2:4:1264:1642:21183',
            'A01535:8:HJ3YYDSX2:4:1318:15483:28025'
        ]

        with open(query_file, 'w') as f:
            f.write(snps)

        reads_found = _filter(query_file, input_bam, output)

        self.assertTrue(len(reads_found) == 7)
        self.assertTrue(len(set(reads_found)) == 4) # 3 paired reads

        for read in reads:
            self.assertTrue(read in reads_found)

    def test_insertions(self):
        query_file = 'freyja/data/test_filter.csv'
        input_bam = 'freyja/data/test.bam'
        output = 'freyja/data/outputs'

        insertions = "(732:'TT'),(1349:'A'),(12333:'A')"

        reads = [
            'A01535:8:HJ3YYDSX2:4:1317:28763:25942',
            'A01535:8:HJ3YYDSX2:4:1141:8341:17566',
            'A01535:8:HJ3YYDSX2:4:1360:23746:26396',
            'A01535:8:HJ3YYDSX2:4:1248:27073:19820'
        ]

        with open(query_file, 'w') as f:
            f.write('\n' + insertions)

        reads_found = _filter(query_file, input_bam, output)
        self.assertTrue(len(reads_found) == 8)
        self.assertTrue(len(set(reads_found)) == 4) # 4 paired reads

        for read in reads:
            self.assertTrue(read in reads_found)

    def test_dels(self):
        query_file = 'freyja/data/test_filter.csv'
        input_bam = 'freyja/data/test.bam'
        output = 'freyja/data/outputs'

        deletions = "(1443:32),(1599:2),(2036:3)"

        reads = [
            'A01535:8:HJ3YYDSX2:4:1373:25663:31313',
            'A01535:8:HJ3YYDSX2:4:1112:19090:5713',
            'A01535:8:HJ3YYDSX2:4:1429:11550:34178',
        ]

        with open(query_file, 'w') as f:
            f.write('\n\n' + deletions)

        reads_found = _filter(query_file, input_bam, output)
        self.assertTrue(len(reads_found) == 6)
        self.assertTrue(len(set(reads_found)) == 3) # 3 paired reads

        for read in reads:
            self.assertTrue(read in reads_found)
            
if __name__ == '__main__':
    unittest.main()