import csv
import re
import subprocess
import unittest
import os


def file_exists(directory, filename):
    file_path = os.path.join(directory, filename)
    return os.path.exists(file_path)


def _parse_version(version_str):
    return tuple(int(p) for p in re.findall(r'\d+', version_str))


class CommandLineTests(unittest.TestCase):

    def test_version(self):
        os.system('freyja --version')

    @unittest.skipUnless(
        os.environ.get('RUN_RELEASE_CHECKS'),
        "release-only check; set RUN_RELEASE_CHECKS=1 to run "
        "(requires full git tag history)"
    )
    def test_version_bumped_for_release(self):
        # Guards against forgetting to bump the version in freyja/_cli.py
        # before cutting a new release. Only meant to be run manually via
        # the "Pre-release version check" GitHub Actions workflow, not as
        # part of normal PR/CI runs.
        with open('freyja/_cli.py') as f:
            cli_source = f.read()

        match = re.search(r"@click\.version_option\('([^']+)'\)", cli_source)
        self.assertIsNotNone(
            match, "Could not find click.version_option(...) in _cli.py")
        cli_version = _parse_version(match.group(1))

        tags = subprocess.run(
            ['git', 'tag', '--list', 'v*'],
            capture_output=True, text=True, check=True
        ).stdout.split()
        self.assertTrue(
            tags, "No git tags found - fetch full tag history "
            "(git fetch --tags) before running this check")

        latest_tag_version = max(_parse_version(t) for t in tags)

        self.assertGreater(
            cli_version, latest_tag_version,
            f"freyja/_cli.py version {match.group(1)} is not newer than "
            f"the latest released tag v{'.'.join(map(str, latest_tag_version))}"
            f". Bump click.version_option in freyja/_cli.py before "
            f"releasing."
        )

    def test_demix(self):
        os.system('freyja demix freyja/data/test.tsv freyja/data/test.depth \
                   --output test.demixed.tsv')
        self.assertTrue(file_exists('.', "test.demixed.tsv"))

        with open('test.demixed.tsv') as f:
            reader = csv.reader(f, delimiter='\t')
            result = {row[0]: row[1] for row in reader if len(row) >= 2}

        self.assertEqual(result['pathogen'], 'SARS-CoV-2')
        self.assertGreater(float(result['coverage']), 90.0)

        lineages = result['lineages'].split()
        self.assertIn('A', lineages)
        self.assertIn('AY.48', lineages)

        abundances = [float(a) for a in result['abundances'].split()]
        self.assertAlmostEqual(sum(abundances), 1.0, places=2)

    def test_demix_with_vcf(self):
        os.system('freyja demix freyja/data/test.vcf freyja/data/test.depth \
                   --output test.demixed.tsv')
        self.assertTrue(file_exists('.', "test.demixed.tsv"))

    def test_demix_with_cutoff(self):
        os.system('freyja demix freyja/data/test.tsv freyja/data/test.depth \
                   --output test_demixed.tsv --depthcutoff 100 --lineageyml \
                   freyja/data/lineages.yml')
        self.assertTrue(file_exists('.',
                                    "test_demixed_collapsed_lineages.yml"))

    def test_cov_res(self):
        os.system('freyja cov-res --region_start 22000 --region_end 25000 \
                   --output test_cov_res')
        self.assertTrue(
            file_exists('.', "test_cov_res_collapsed_lineages.yml"))

    def test_covariants_assign(self):
        os.system('freyja covariants-assign \
                   freyja/data/covariants_example.tsv \
                   --region_start 22000 --region_end 25000 \
                   --output test_covariants_assign')
        self.assertTrue(
            file_exists('.', "test_covariants_assign_collapsed_lineages.yml"))
        self.assertTrue(
            file_exists('.', "test_covariants_assign_clusters.tsv"))
        self.assertTrue(
            file_exists('.', "test_covariants_assign_abundances.tsv"))

    def test_covariants_assign_region(self):
        os.system('freyja covariants-assign \
                   freyja/data/covariants_example.tsv --collapse region \
                   --region_start 22000 --region_end 25000 \
                   --output test_covariants_assign_region')
        self.assertTrue(file_exists(
            '.', "test_covariants_assign_region_abundances.tsv"))

    def test_covariants_assign_ties(self):
        for ties in ('mrca', 'majority'):
            os.system(
                'freyja covariants-assign freyja/data/covariants_example.tsv'
                ' --region_start 22000 --region_end 25000 --ties '
                f'{ties} --collapse region '
                f'--output test_covariants_assign_ties_{ties}')
            self.assertTrue(file_exists(
                '.', f"test_covariants_assign_ties_{ties}_clusters.tsv"))

    def test_covariants_assign_groups(self):
        for assign_by in ('mrca', 'majority'):
            os.system(
                'freyja covariants-assign freyja/data/covariants_example.tsv'
                ' --region_start 22000 --region_end 25000 --groups '
                'freyja/data/example_lineage_groups.yml --assign_by '
                f'{assign_by} --max_missing 1 --max_extra 1 --unassigned drop '
                f'--output test_covariants_assign_{assign_by}')
            self.assertTrue(file_exists(
                '.', f"test_covariants_assign_{assign_by}"
                     "_group_abundances.tsv"))

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

    def test_ampstat(self):
        os.system('freyja ampliconstat '
                  '--primer freyja/data/ARTIC_V4-1.bed '
                  '--input_depth freyja/data/test.depth '
                  '--output_csv test.amplicon.csv '
                  '--output_plot test.amplicon.png')
        self.assertTrue(file_exists('.', "test.amplicon.csv"))
        self.assertTrue(file_exists('.', "test.amplicon.png"))


if __name__ == '__main__':
    unittest.main()
