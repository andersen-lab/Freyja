import click
import pandas as pd
from freyja.convert_paths2barcodes import parse_tree_paths,\
     convert_to_barcodes, reversion_checking
from freyja.sample_deconv import buildLineageMap, build_mix_and_depth_arrays,\
    reindex_dfs, map_to_constellation, solve_demixing_problem
import os
import urllib.request


@click.group()
def cli():
    pass


@cli.command()
@click.argument('filename', type=click.Path(exists=True))
def barcode(filename):
    print('Building barcodes from global phylogenetic tree')
    df = pd.read_csv(filename, sep='\t')
    df = parse_tree_paths(df)
    df_barcodes = convert_to_barcodes(df)
    df_barcodes = reversion_checking(df_barcodes)
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
    df_barcodes.to_csv(os.path.join(locDir, 'data/usher_barcodes.csv'))


@cli.command()
@click.argument('variants', type=click.Path(exists=True))
@click.argument('depths', type=click.Path(exists=True))
@click.option('--output', default='demixing_result.csv', help='Output file',
              type=click.Path(exists=False))
def demix(variants, depths, output):
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
    df_barcodes = pd.read_csv(os.path.join(locDir, 'data/usher_barcodes.csv'), index_col=0)
    muts = list(df_barcodes.columns)
    mapDict = buildLineageMap()
    print('building mix/depth matrices')
    # assemble data from of (possibly) mixed samples
    mix, depths_ = build_mix_and_depth_arrays(variants, depths, muts)
    print('demixing')
    df_barcodes, mix, depths_ = reindex_dfs(df_barcodes, mix, depths_)
    sample_strains, abundances, error = solve_demixing_problem(df_barcodes, mix, depths_)
    localDict = map_to_constellation(sample_strains, abundances, mapDict)
    # assemble into series and write.
    sols_df = pd.Series(data=(localDict, sample_strains, abundances, error),
                        index=['summarized', 'lineages', 'abundances', 'resid'],
                        name=mix.name)
    sols_df.to_csv(output, sep='\t')


@cli.command()
def update():
    print('Downloading a new global tree')
    url = 'http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz'
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
    urllib.request.urlretrieve(url, os.path.join(locDir, "data/public-latest.all.masked.pb.gz"))
    print('Downloading an updated curated lineage set from outbreak.info')
    url2 = 'https://raw.githubusercontent.com/outbreak-info/outbreak.info/master/web/src/assets/genomics/curated_lineages.json'
    urllib.request.urlretrieve(url2, os.path.join(locDir, "data/curated_lineages.json"))


if __name__ == '__main__':
    cli()
