import unittest
import os
import sys


class VariantsTests(unittest.TestCase):
    def test_variants(self):
        filePathVariants = 'freyja/data/test.tsv'
        filePathDepth = 'freyja/data/test.depth'

        if os.path.exists(filePathDepth):
            os.remove(filePathDepth)
            os.remove(filePathVariants)
        varCmd = "freyja variants freyja/data/test.bam "\
                 "--variants freyja/data/test --depths freyja/data/test.depth"
        sys.stdout.flush()
        os.system(varCmd)
        self.assertTrue(os.path.exists(filePathVariants))
        self.assertTrue(os.path.exists(filePathDepth))


if __name__ == '__main__':
    unittest.main()
