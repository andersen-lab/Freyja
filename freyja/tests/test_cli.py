import unittest
import os


def file_exists(directory, filename):
    file_path = os.path.join(directory, filename)
    return os.path.exists(file_path)


class CommandLineTests(unittest.TestCase):

    def test_version(self):
        os.system('freyja --version')

    def test_demix(self):
        os.system('freyja demix freyja/data/test.tsv freyja/data/test.depth \
                   --output test.demixed.tsv')
        self.assertTrue(file_exists('.', "test.demixed.tsv"))

    def test_demix_with_cutoff(self):
        os.system('freyja demix freyja/data/test.tsv freyja/data/test.depth \
                   --output test.demixed.tsv --depthcutoff 100 --lineageyml \
                   freyja/data/lineages.yml')
        self.assertTrue(file_exists('.', "test_collapsed_lineages.yml"))

    def test_plot(self):
        os.system('freyja plot freyja/data/aggregated_result.tsv \
                   --output test_plot.pdf')
        self.assertTrue(file_exists('.', "test_plot.pdf"))

    def test_plot_time(self):
        os.system('freyja plot freyja/data/test_sweep.tsv \
                   --times freyja/data/sweep_metadata.csv \
                   --output test_plot_time.pdf \
                   --config freyja/data/plot_config.yml --lineageyml \
                   freyja/data/lineages.yml --interval D')
        self.assertTrue(file_exists('.', "test_plot_time.pdf"))

    def test_growth_rate(self):
        os.system('freyja relgrowthrate freyja/data/test_sweep.tsv \
                   freyja/data/sweep_metadata.csv \
                   --output test_growth_rates.csv \
                   --config freyja/data/plot_config.yml \
                   --lineageyml freyja/data/lineages.yml')
        self.assertTrue(file_exists('.', "test_growth_rates.csv"))

    def test_dash(self):
        os.system('freyja dash freyja/data/test_sweep.tsv \
                   freyja/data/sweep_metadata.csv \
                   freyja/data/title.txt \
                   freyja/data/introContent.txt \
                   --lineageyml freyja/data/lineages.yml \
                   --output test_dash.html')
        self.assertTrue(file_exists('.', "test_dash.html"))

    def test_get_lineage_def(self):
        os.system('freyja get-lineage-def B.1.1.7 '
                  '--annot freyja/data/NC_045512_Hu-1.gff '
                  '--ref freyja/data/NC_045512_Hu-1.fasta '
                  '--output lineage_def.txt')
        self.assertTrue(file_exists('.', "lineage_def.txt"))
        with open('lineage_def.txt') as f:
            self.assertEqual(len(f.readlines()), 28)

    def test_boot(self):
        os.system('freyja boot '
                  'freyja/data/test.tsv freyja/data/test.depth '
                  '--nt 10 --nb 10 --output_base boot_output '
                  '--bootseed 10')
        self.assertTrue(file_exists('.', "boot_output_lineages.csv"))


if __name__ == '__main__':
    unittest.main()
