import click
import pandas as pd
from freyja.convert_paths2barcodes import parse_tree_paths,\
     convert_to_barcodes, reversion_checking
from freyja.sample_deconv import buildLineageMap, build_mix_and_depth_arrays,\
    reindex_dfs, map_to_constellation, solve_demixing_problem,\
    perform_bootstrap
from freyja.updates import download_tree, convert_tree,\
                           get_curated_lineage_data
from freyja.utils import agg, makePlot_simple, makePlot_time
import os
import subprocess
import sys

locDir = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))


@click.group()
def cli():
    pass


@cli.command()
@click.argument('variants', type=click.Path(exists=True))
@click.argument('depths', type=click.Path(exists=True))
@click.option('--eps', default=1e-3, help='minimum abundance to include')
@click.option('--barcodes', default='-1', help='custom barcode file')
@click.option('--output', default='demixing_result.csv', help='Output file',
              type=click.Path(exists=False))
def demix(variants, depths, output, eps, barcodes):
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    # option for custom barcodes
    if barcodes != '-1':
        df_barcodes = pd.read_csv(barcodes, index_col=0)
    else:
        df_barcodes = pd.read_csv(os.path.join(locDir,
                                  'data/usher_barcodes.csv'), index_col=0)
    muts = list(df_barcodes.columns)
    mapDict = buildLineageMap()
    print('building mix/depth matrices')
    # assemble data from (possibly) mixed samples
    mix, depths_ = build_mix_and_depth_arrays(variants, depths, muts)
    print('demixing')
    df_barcodes, mix, depths_ = reindex_dfs(df_barcodes, mix, depths_)
    sample_strains, abundances, error = solve_demixing_problem(df_barcodes,
                                                               mix,
                                                               depths_,
                                                               eps)
    localDict = map_to_constellation(sample_strains, abundances, mapDict)
    # assemble into series and write.
    sols_df = pd.Series(data=(localDict, sample_strains, abundances, error),
                        index=['summarized', 'lineages',
                        'abundances', 'resid'],
                        name=mix.name)
    sols_df.to_csv(output, sep='\t')


@cli.command()
def update():
    # get data from UShER
    print('Downloading a new global tree')
    download_tree()
    print('Getting outbreak data')
    get_curated_lineage_data()
    print("Converting tree info to barcodes")
    convert_tree()  # returns paths for each lineage
    # Now parse into barcode form
    lineagePath = os.path.join(os.curdir, "lineagePaths.txt")
    print('Building barcodes from global phylogenetic tree')
    df = pd.read_csv(lineagePath, sep='\t')
    df = parse_tree_paths(df)
    df_barcodes = convert_to_barcodes(df)
    df_barcodes = reversion_checking(df_barcodes)
    df_barcodes.to_csv(os.path.join(locDir, 'data/usher_barcodes.csv'))
    # delete files generated along the way that aren't needed anymore
    print('Cleaning up')
    os.remove(lineagePath)
    os.remove(os.path.join(locDir, "data/public-latest.all.masked.pb.gz"))


@cli.command()
@click.argument('bamfile', type=click.Path(exists=True))
@click.option('--ref', help='Reference',
              default=os.path.join(locDir,
                                   'data/NC_045512_Hu-1.fasta'),
              type=click.Path())
@click.option('--variants', help='Variant call output file', type=click.Path())
@click.option('--depths', help='Sequencing depth output file',
              type=click.Path())
@click.option('--refname', help='Ref name (for bams with multiple sequences)',
              default='')
def variants(bamfile, ref, variants, depths, refname):
    if len(refname) == 0:
        bashCmd = f"samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f "\
                  f"{ref} {bamfile} | tee >(cut -f1-4 > {depths}) |"\
                  f" ivar variants -p {variants} -q 20 -t 0.0 -r {ref}"
    else:
        bashCmd = f"samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f "\
                  f"{ref} {bamfile} -r {refname} | tee >(cut -f1-4 > {depths}"\
                  f") | ivar variants -p {variants} -q 20 -t 0.0 -r {ref}"
    sys.stdout.flush()  # force python to flush
    completed = subprocess.run(bashCmd, shell=True, executable="/bin/bash",
                               stdout=subprocess.PIPE)
    sys.exit(completed.returncode)


@cli.command()
@click.argument('variants', type=click.Path(exists=True))
@click.argument('depths', type=click.Path(exists=True))
@click.option('--nb', default=100, help='number of bootstraps')
@click.option('--nt', default=1, help='max number of cpus to use')
@click.option('--eps', default=1e-3, help='minimum abundance to include')
@click.option('--barcodes', default='-1', help='custom barcode file')
@click.option('--output_base', default='test', help='Output file basename',
              type=click.Path(exists=False))
def boot(variants, depths, output_base, eps, barcodes, nb, nt):
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    # option for custom barcodes
    if barcodes != '-1':
        df_barcodes = pd.read_csv(barcodes, index_col=0)
    else:
        df_barcodes = pd.read_csv(os.path.join(locDir,
                                  'data/usher_barcodes.csv'), index_col=0)
    muts = list(df_barcodes.columns)
    mapDict = buildLineageMap()
    print('building mix/depth matrices')
    # assemble data from (possibly) mixed samples
    mix, depths_ = build_mix_and_depth_arrays(variants, depths, muts)
    print('demixing')
    df_barcodes, mix, depths_ = reindex_dfs(df_barcodes, mix, depths_)
    lin_out, constell_out = perform_bootstrap(df_barcodes, mix, depths_,
                                              nb, eps, nt, mapDict, muts)
    lin_out.to_csv(output_base + '_lineages.csv')
    constell_out.to_csv(output_base + '_summarized.csv')


@cli.command()
@click.argument('results', type=click.Path(exists=True))
@click.option('--output', default='aggregated_result.tsv', help='Output file',
              type=click.Path(exists=False))
def aggregate(results, output):
    df_demix = agg(results)
    df_demix.to_csv(output, sep='\t')


@cli.command()
@click.argument('agg_results', type=click.Path(exists=True))
@click.option('--lineages', is_flag=True)
@click.option('--times', default='-1')
@click.option('--interval', default='MS')
@click.option('--colors', default='', help='path to csv of hex codes')
@click.option('--output', default='mix_plot.pdf', help='Output file')
@click.option('--windowsize', default=14)
def plot(agg_results, lineages, times, interval, output, windowsize, colors):
    agg_df = pd.read_csv(agg_results, skipinitialspace=True, sep='\t',
                         index_col=0)
    if len(colors) > 0:
        colors0 = pd.read_csv(colors, header=None).values[0]
    else:
        colors0 = colors
    if times == '-1':
        # make basic plot, without time info
        makePlot_simple(agg_df, lineages, output, colors0)
    else:
        # make time aware plot
        times_df = pd.read_csv(times, skipinitialspace=True,
                               index_col=0)
        times_df['sample_collection_datetime'] = \
            pd.to_datetime(times_df['sample_collection_datetime'])
        makePlot_time(agg_df, lineages, times_df, interval, output,
                      windowsize, colors0)


if __name__ == '__main__':
    cli()
