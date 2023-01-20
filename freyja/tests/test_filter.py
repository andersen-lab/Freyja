import unittest
from freyja.read_analysis_tools import filter as _filter
import os
import subprocess

class FilterTests(unittest.TestCase):
            
    def test_snps(self):
        query_file = 'freyja/data/test_filter.csv'
        bam_dir = 'freyja/data/test.bam'
        output = 'freyja/data/outputs/test_filtered.bam'


        snps = 'C75T,G230A,A543C'
        reads = [
            'A01535:8:HJ3YYDSX2:4:1377:5177:13182',
            'A01535:8:HJ3YYDSX2:4:1377:5177:13182',
            'A01535:8:HJ3YYDSX2:4:1266:9254:16329',
            'A01535:8:HJ3YYDSX2:4:1266:9254:16329'
            'A01535:8:HJ3YYDSX2:4:1264:1642:21183',
            'A01535:8:HJ3YYDSX2:4:1264:1642:21183',
            'A01535:8:HJ3YYDSX2:4:1318:15483:28025'
        ]

        with open(query_file, 'w') as f:
            f.write(snps)

        reads_found = _filter(query_file, bam_dir, output)

        self.assertTrue(len(reads_found) == 7)
        self.assertTrue(len(set(reads_found)) == 4) # 3 reverse reads

        for read in reads_found:
            self.assertTrue(read in reads)

    def test_insertions(self):
        query_file = 'freyja/data/test_filter.csv'
        bam_dir = 'freyja/data/test.bam'
        output = 'freyja/data/outputs/test_filtered.bam'


    

if __name__ == '__main__':
    unittest.main()