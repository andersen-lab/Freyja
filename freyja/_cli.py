import glob
import os
import subprocess
import sys

import click
import pandas as pd
import yaml

from freyja.convert_paths2barcodes import (check_mutation_chain,
                                           convert_to_barcodes,
                                           parse_tree_paths,
                                           reversion_checking)
from freyja.read_analysis_tools import covariants as _covariants
from freyja.read_analysis_tools import extract as _extract
from freyja.read_analysis_tools import filter as _filter
from freyja.read_analysis_tools import plot_covariants as _plot_covariants
from freyja.sample_deconv import (build_mix_and_depth_arrays, buildLineageMap,
                                  map_to_constellation, perform_bootstrap,
                                  reindex_dfs, solve_demixing_problem)
from freyja.updates import (convert_tree, convert_tree_custom,
                            download_barcodes, download_tree, get_cl_lineages,
                            get_curated_lineage_data)
from freyja.utils import (agg, calc_rel_growth_rates, checkConfig,
                          collapse_barcodes, get_abundance, make_dashboard,
                          makePlot_simple, makePlot_time,
                          read_lineage_file, handle_region_of_interest)

locDir = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))


@click.group()
@click.version_option('1.4.9')
def cli():
    pass


def print_barcode_version(ctx, param, value):
    """
    Gets the barcode version used in the program
    """
    if not value or ctx.resilient_parsing:
        return
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    f = open(os.path.join(locDir, 'data/last_barcode_update.txt'), 'r')
    click.echo('Barcode version:')
    click.echo(f.readline())
    ctx.exit()


@cli.command()
@click.argument('variants', type=click.Path(exists=True),
                help="Path to variant.tsv file generated"
                     " using freyja variants function")
@click.argument('depths', type=click.Path(exists=True),
                help="Path to depths file generated using"
                     " freyja variants function")
@click.option('--eps', default=1e-3,
              help='minimum abundance to include for each'
                   ' lineage', show_default=True)
@click.option('--barcodes', default='-1',
              help='Path to custom barcode file',
              show_default=True)
@click.option('--meta', default='-1',
              help='custom lineage to variant metadata file',
              show_default=True)
@click.option('--output', default='output name',
              help='Output file',
              type=click.Path(exists=False), show_default=True)
@click.option('--covcut', default=10,
              help='calculate percent of sites with n or greater reads',
              show_default=True)
@click.option('--confirmedonly', is_flag=True,
              help="exclude unconfirmed lineages",
              default=False, show_default=True)
@click.option('--version', is_flag=True, callback=print_barcode_version,
              expose_value=False, is_eager=True, show_default=True)
@click.option('--depthcutoff', default=0,
              help='exclude sites with coverage depth below this value and'
              'group identical barcodes', show_default=True)
@click.option('--lineageyml', default='-1',
              help='lineage hierarchy file in a yaml format',
              show_default=True)
@click.option('--adapt', default=0.,
              help='adaptive lasso penalty parameter',
              show_default=True)
@click.option('--a_eps', default=1E-8,
              help='adaptive lasso parameter, hard threshold',
              show_default=True)
@click.option('--region_of_interest', default='-1',
              help='JSON file containing region(s) of interest'
                   ' for which to compute additional coverage estimates',
              show_default=True)
def demix(variants, depths, output, eps, barcodes, meta,
          covcut, confirmedonly, depthcutoff, lineageyml,
          adapt, a_eps, region_of_interest):
    """
    Generate prevalence of lineages per sample
    """
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    # option for custom barcodes
    if barcodes != '-1':
        df_barcodes = pd.read_csv(barcodes, index_col=0)
    else:
        df_barcodes = pd.read_csv(os.path.join(locDir,
                                  'data/usher_barcodes.csv'), index_col=0)
    if confirmedonly:
        confirmed = [dfi for dfi in df_barcodes.index
                     if 'proposed' not in dfi and 'misc' not in dfi]
        df_barcodes = df_barcodes.loc[confirmed, :]

    # drop intra-lineage diversity naming (keeps separate barcodes)
    indexSimplified = [dfi.split('_')[0] for dfi in df_barcodes.index]
    df_barcodes = df_barcodes.loc[indexSimplified, :]
    df_depth = pd.read_csv(depths, sep='\t', header=None, index_col=1)
    if depthcutoff != 0:
        df_barcodes = collapse_barcodes(df_barcodes, df_depth, depthcutoff,
                                        lineageyml, locDir, output)
    muts = list(df_barcodes.columns)
    mapDict = buildLineageMap(meta)
    print('building mix/depth matrices')
    # assemble data from (possibly) mixed samples
    mix, depths_, cov = build_mix_and_depth_arrays(variants, depths, muts,
                                                   covcut)
    print('demixing')
    df_barcodes, mix, depths_ = reindex_dfs(df_barcodes, mix, depths_)
    try:
        sample_strains, abundances, error = solve_demixing_problem(df_barcodes,
                                                                   mix,
                                                                   depths_,
                                                                   eps, adapt,
                                                                   a_eps)
    except Exception as e:
        print(e)
        print('Error: Demixing step failed. Returning empty data output')
        sample_strains, abundances = [], []
        error = -1
    # merge intra-lineage diversity if multiple hits.
    if len(set(sample_strains)) < len(sample_strains):
        localDict = {}
        for jj, lin in enumerate(sample_strains):
            if lin not in localDict.keys():
                localDict[lin] = abundances[jj]
            else:
                localDict[lin] += abundances[jj]
        # ensure descending order
        localDict = dict(sorted(localDict.items(),
                                key=lambda x: x[1],
                                reverse=True))
        sample_strains = list(localDict.keys())
        abundances = list(localDict.values())

    localDict = map_to_constellation(sample_strains, abundances, mapDict)

    # assemble into series and write.
    sols_df = pd.Series(data=(localDict, sample_strains, abundances,
                              error, cov),
                        index=['summarized', 'lineages',
                        'abundances', 'resid', 'coverage'],
                        name=mix.name)

    # Determine coverage in region(s) of interest (if specified)
    if region_of_interest != '-1':

        sols_df = handle_region_of_interest(region_of_interest, sols_df,
                                            df_depth, covcut, mix.name)

    # convert lineage/abundance readouts to single line strings
    sols_df['lineages'] = ' '.join(sols_df['lineages'])
    sols_df['abundances'] = ['%.8f' % ab for ab in sols_df['abundances']]
    sols_df['abundances'] = ' '.join(sols_df['abundances'])
    sols_df.to_csv(output, sep='\t')


@cli.command()
@click.option('--outdir', default='-1',
              help='Output directory save updated files',
              show_default=True)
@click.option('--noncl', is_flag=True, default=True,
              help='only include include lineages that are'
                   'confirmed by cov-lineages',
              show_default=True)
@click.option('--buildlocal', is_flag=True, default=False,
              help='Perform barcode building locally',
              show_default=True)
def update(outdir, noncl, buildlocal):
    """
    Updates the lineage information using the latest
    outbreak data
    """
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                                          os.pardir))
    if outdir != '-1':
        # User specified directory
        locDir = outdir
    else:
        locDir = os.path.join(locDir, 'data')

    print('Getting outbreak data')
    get_curated_lineage_data(locDir)
    get_cl_lineages(locDir)
    # # get data from UShER
    if buildlocal:
        print('Downloading a new global tree')
        download_tree(locDir)
        print("Converting tree info to barcodes")
        convert_tree(locDir)  # returns paths for each lineage
        # Now parse into barcode form
        lineagePath = os.path.join(os.curdir, "lineagePaths.txt")
        print('Building barcodes from global phylogenetic tree')
        df = pd.read_csv(lineagePath, sep='\t')
        df = parse_tree_paths(df)
        df_barcodes = convert_to_barcodes(df)
        df_barcodes = reversion_checking(df_barcodes)
        df_barcodes = check_mutation_chain(df_barcodes)

        # as default:
        # since usher tree can be ahead of cov-lineages,
        # we drop lineages not in cov-lineages
        if noncl:
            # read linages.yml file
            with open(os.path.join(locDir, 'lineages.yml'), 'r') as f:
                try:
                    lineages_yml = yaml.safe_load(f)
                except yaml.YAMLError as exc:
                    raise ValueError('Error in lineages.yml file: ' + str(exc))
            lineageNames = [lineage['name'] for lineage in lineages_yml]
            df_barcodes = df_barcodes.loc[df_barcodes.index.isin(lineageNames)]
        else:
            print("Including lineages not yet in cov-lineages.")
        df_barcodes.to_csv(os.path.join(locDir, 'usher_barcodes.csv'))
        # delete files generated along the way that aren't needed anymore
        print('Cleaning up')
        os.remove(lineagePath)
        os.remove(os.path.join(locDir, "public-latest.all.masked.pb.gz"))
    else:
        print('Downloading barcodes')
        download_barcodes(locDir)


@cli.command()
@click.option('--pb', type=click.Path(exists=True),
              help='protobuf tree', show_default=True)
@click.option('--outdir', type=click.Path(exists=True),
              help='Output directory save updated files',
              show_default=True)
@click.option('--noncl', is_flag=True, default=True,
              help='only include include lineages that are'
                   'confirmed by cov-lineages', show_default=True)
def barcode_build(pb, outdir, noncl):
    """
    Building barcodes from a global tree
    """
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                                          os.pardir))
    locDir = outdir
    print('Getting outbreak data')
    get_curated_lineage_data(locDir)
    get_cl_lineages(locDir)
    # # get data from UShER
    print('Downloading a new global tree')
    print("Converting tree info to barcodes")
    convert_tree_custom(pb)  # returns paths for each lineage
    # Now parse into barcode form
    lineagePath = os.path.join(os.curdir, "lineagePaths.txt")
    print('Building barcodes from global phylogenetic tree')
    df = pd.read_csv(lineagePath, sep='\t')
    df = parse_tree_paths(df)
    df_barcodes = convert_to_barcodes(df)
    df_barcodes = reversion_checking(df_barcodes)
    df_barcodes = check_mutation_chain(df_barcodes)

    # as default:
    # since usher tree can be ahead of cov-lineages,
    # we drop lineages not in cov-lineages
    if noncl:
        # read linages.yml file
        with open(os.path.join(locDir, 'lineages.yml'), 'r') as f:
            try:
                lineages_yml = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                raise ValueError('Error in lineages.yml file: ' + str(exc))
        lineageNames = [lineage['name'] for lineage in lineages_yml]
        df_barcodes = df_barcodes.loc[df_barcodes.index.isin(lineageNames)]
    else:
        print("Including lineages not yet in cov-lineages.")
    df_barcodes.to_csv(os.path.join(locDir, 'usher_barcodes.csv'))
    # delete files generated along the way that aren't needed anymore
    print('Cleaning up')
    os.remove(lineagePath)


@cli.command()
@click.argument('bamfile', type=click.Path(exists=True))
@click.option('--ref', help='Reference file in fasta format',
              default=os.path.join(locDir,
                                   'data/NC_045512_Hu-1.fasta'),
              type=click.Path(), show_default=True)
@click.option('--variants', help='Variant calling output file',
              type=click.Path(), show_default=True)
@click.option('--depths', help='Sequencing depth output file',
              type=click.Path(), show_default=True)
@click.option('--refname', help='Ref name (for bams with multiple sequences)',
              default='', show_default=True)
@click.option('--minq', help='Minimum base quality score',
              default=20, show_default=True)
@click.option('--annot', help='AA annotation output',
              default='', show_default=True)
def variants(bamfile, ref, variants, depths, refname, minq, annot):
    """
    Perform variant calling using samtools and ivar
    """
    if len(refname) == 0:
        bashCmd = f"samtools mpileup -aa -A -d 600000 -Q {minq} -q 0 -B -f "\
                  f"{ref} {bamfile} | tee >(cut -f1-4 > {depths}) |"\
                  f" ivar variants -p {variants} -q {minq} -t 0.0 -r {ref}"
    else:
        bashCmd = f"samtools mpileup -aa -A -d 600000 -Q {minq} -q 0 -B -f "\
                  f"{ref} {bamfile} -r {refname} | tee >(cut -f1-4 > {depths}"\
                  f") | ivar variants -p {variants} -q {minq} -t 0.0 -r {ref}"
    if len(annot) > 0:
        print('Including annotation')
        bashCmd = bashCmd + f" -g {annot}"
    print(bashCmd)
    sys.stdout.flush()  # force python to flush
    completed = subprocess.run(bashCmd, shell=True, executable="/bin/bash",
                               stdout=subprocess.PIPE)
    sys.exit(completed.returncode)


@cli.command()
@click.argument('variants',
                help="path to variants.tsv file"
                     "generated using freyja variants function",
                type=click.Path(exists=True))
@click.argument('depths',
                help="path to depths file generated"
                     " using freyja variants function",
                type=click.Path(exists=True))
@click.option('--nb', default=100,
              help='number of times bootstrapping is performed',
              show_default=True)
@click.option('--nt', default=1,
              help='max number of cpus to use',
              show_default=True)
@click.option('--eps', default=1e-3,
              help='minimum abundance to include',
              show_default=True)
@click.option('--barcodes', default='-1',
              help='custom barcode file',
              show_default=True)
@click.option('--meta', default='-1',
              help='custom lineage to variant metadata file',
              show_default=True)
@click.option('--output_base', default='test',
              help='Output file basename',
              type=click.Path(exists=False), show_default=True)
@click.option('--boxplot', default='',
              help='file format of boxplot output (e.g. pdf or png)',
              show_default=True)
@click.option('--confirmedonly', is_flag=True,
              help="exclude unconfirmed lineages",
              default=False, show_default=True)
@click.option('--rawboots', is_flag=True, default=False,
              help='return raw bootstraps', show_default=True)
@click.option('--lineageyml', default='-1',
              help='lineage hierarchy file in a yaml format',
              show_default=True)
@click.option('--depthcutoff', default=0,
              help='exclude sites with coverage depth below this value and'
              'group identical barcodes', show_default=True)
def boot(variants, depths, output_base, eps, barcodes, meta,
         nb, nt, boxplot, confirmedonly, lineageyml, depthcutoff, rawboots):
    """
    Perform bootstrapping method for freyja
    """
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    # option for custom barcodes
    if barcodes != '-1':
        df_barcodes = pd.read_csv(barcodes, index_col=0)
    else:
        df_barcodes = pd.read_csv(os.path.join(locDir,
                                  'data/usher_barcodes.csv'), index_col=0)

    if confirmedonly:
        confirmed = [dfi for dfi in df_barcodes.index
                     if 'proposed' not in dfi and 'misc' not in dfi]
        df_barcodes = df_barcodes.loc[confirmed, :]

    # drop intra-lineage diversity naming (keeps separate barcodes)
    indexSimplified = [dfi.split('_')[0] for dfi in df_barcodes.index]
    df_barcodes = df_barcodes.loc[indexSimplified, :]

    df_depths = pd.read_csv(depths, sep='\t', header=None, index_col=1)
    if depthcutoff != 0:
        df_barcodes = collapse_barcodes(
            df_barcodes, df_depths, depthcutoff,
            lineageyml, locDir, output_base)

    muts = list(df_barcodes.columns)
    mapDict = buildLineageMap(meta)
    print('building mix/depth matrices')
    # assemble data from (possibly) mixed samples
    covcut = 10  # set value, coverage estimate not returned to user from boot
    mix, depths_, cov = build_mix_and_depth_arrays(variants, depths, muts,
                                                   covcut)
    print('demixing')
    df_barcodes, mix, depths_ = reindex_dfs(df_barcodes, mix, depths_)
    lin_df, constell_df = perform_bootstrap(df_barcodes, mix, depths_,
                                            nb, eps, nt, mapDict, muts,
                                            boxplot, output_base)
    if rawboots:
        lin_df.to_csv(output_base + '_lineages_boot.csv')
        constell_df.to_csv(output_base + '_summarized_boot.csv')

    q = [0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975]
    lin_out = lin_df.quantile(q)
    constell_out = constell_df.quantile(q)

    lin_out.to_csv(output_base + '_lineages.csv')
    constell_out.to_csv(output_base + '_summarized.csv')


@cli.command()
@click.argument('results',
                help="path to results",
                type=click.Path(exists=True))
@click.option('--ext', default='-1',
              help='file extension option',
              show_default=True)
@click.option('--output',
              help= "name for aggregated results",
              default='aggregated_result.tsv',
              type=click.Path(exists=False), show_default=True)
def aggregate(results, ext, output):
    """
    Aggregates all the outputs
    """
    if ext != '-1':
        results_ = [fn for fn in glob.glob(results + '*' + ext)]
    else:
        results_ = [results + fn for fn in os.listdir(results)]
    df_demix = agg(results_)
    df_demix.to_csv(output, sep='\t')


@cli.command()
@click.argument('agg_results', type=click.Path(exists=True))
@click.option('--lineages', is_flag=True, show_default=True)
@click.option('--times', default='-1', show_default=True)
@click.option('--interval', default='MS', show_default=True)
@click.option('--config', default=None, help='path to yaml file', show_default=True)
@click.option('--mincov', default=60., help='min genome coverage included', show_default=True)
@click.option('--output', default='mix_plot.pdf', help='Output file', show_default=True)
@click.option('--windowsize', default=14, show_default=True)
@click.option('--lineageyml', default='-1', help='lineage hierarchy file', show_default=True)
@click.option('--thresh', default=0.01, help='min lineage abundance included', show_default=True)
@click.option('--writegrouped', default='-1',
              help='path to write grouped lineage data', show_default=True)
def plot(agg_results, lineages, times, interval, output, windowsize,
         config, mincov, lineageyml, thresh, writegrouped):
    """
        create plots using the outputs
    """
    agg_df = pd.read_csv(agg_results, skipinitialspace=True, sep='\t',
                         index_col=0)
    # drop poor quality samples
    if 'coverage' in agg_df.columns:
        agg_df = agg_df[agg_df['coverage'] > mincov]
    else:
        print('WARNING: Freyja should be updated \
to include coverage estimates.')

    if agg_df.shape[0] == 0:
        print('ERROR: No samples matching coverage requirements, \
so no plot will be generated. Try changing --mincov threshold.')
        exit()
    if config is not None:
        with open(config, "r") as f:
            try:
                config = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                raise ValueError('Error in config file: ' + str(exc))

    # convert lineages_yml to a dictionary where the lineage names are the
    # keys.
    lineage_info = read_lineage_file(lineageyml, locDir)

    if config is not None:
        config = checkConfig(config)
    else:
        config = {}
    agg_df['abundances'] = agg_df['abundances'].astype(str)
    agg_df['summarized'] = agg_df['summarized'].astype(str)
    agg_df = agg_df[agg_df['summarized'] != '[]']
    if times == '-1':
        # make basic plot, without time info
        makePlot_simple(agg_df, lineages, output, config, lineage_info,
                        thresh, writegrouped)
    else:
        # make time aware plot
        times_df = pd.read_csv(times, skipinitialspace=True,
                               index_col=0)
        times_df['sample_collection_datetime'] = \
            pd.to_datetime(times_df['sample_collection_datetime'])
        makePlot_time(agg_df, lineages, times_df, interval, output,
                      windowsize, config, lineage_info, thresh,
                      writegrouped)


@cli.command()
@click.argument('agg_results', type=click.Path(exists=True))
@click.argument('metadata', type=click.Path(exists=True))
@click.argument('title', type=click.Path(exists=True))
@click.argument('intro', type=click.Path(exists=True))
@click.option('--thresh', default=0.01, help='min lineage abundance included', show_default=True)
@click.option('--headerColor', default='#2fdcf5', help='color of header', show_default=True)
@click.option('--bodyColor', default='#ffffff', help='color of body', show_default=True)
@click.option('--scale_by_viral_load', is_flag=True,
              help='scale by viral load', show_default=True)
@click.option('--nboots', default=1000, help='Number of Bootstrap iterations', show_default=True)
@click.option('--serial_interval', default=5.5, help='Serial Interval', show_default=True)
@click.option('--config', default=None, help='path to yaml file', show_default=True)
@click.option('--mincov', default=60., help='min genome coverage included', show_default=True)
@click.option('--output', default='mydashboard.html', help='Output html file', show_default=True)
@click.option('--days', default=56, help='N Days used for growth calc', show_default=True)
@click.option('--grthresh', default=0.001, help='min avg prev. for growth', show_default=True)
@click.option('--lineageyml', default='-1', help='lineage hierarchy file', show_default=True)
@click.option('--keep_plot_files', is_flag=True, help='Keep separate plot', show_default=True)
def dash(agg_results, metadata, title, intro, thresh, headercolor, bodycolor,
         scale_by_viral_load, nboots, serial_interval, config, mincov, output,
         days, lineageyml, grthresh, keep_plot_files):
    agg_df = pd.read_csv(agg_results, skipinitialspace=True, sep='\t',
                         index_col=0)
    """
    create plots using the outputs
    """
    # drop poor quality samples
    if 'coverage' in agg_df.columns:
        agg_df = agg_df[agg_df['coverage'] > mincov]
    else:
        print('WARNING: Freyja should be updated ' +
              'to include coverage estimates.')
    agg_df = agg_df[agg_df['summarized'] != '[]']

    meta_df = pd.read_csv(metadata, index_col=0)
    meta_df['sample_collection_datetime'] = \
        pd.to_datetime(meta_df['sample_collection_datetime'])
    # read in inputs
    with open(title, "r") as f:
        titleText = ''.join(f.readlines())
    with open(intro, "r") as f:
        introText = ''.join(f.readlines())
    if config is not None:
        with open(config, "r") as f:
            try:
                config = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                raise ValueError('Error in config file: ' + str(exc))

    lineage_info = read_lineage_file(lineageyml, locDir)
    if config is not None:
        config = checkConfig(config)
    else:
        config = {}
    make_dashboard(agg_df, meta_df, thresh, titleText, introText,
                   output, headercolor, bodycolor, scale_by_viral_load, config,
                   lineage_info, nboots, serial_interval, days, grthresh,
                   keep_plot_files)


@cli.command()
@click.argument('agg_results', type=click.Path(exists=True))
@click.argument('metadata', type=click.Path(exists=True))
@click.option('--thresh', default=0.01, help='min lineage abundance in plot', show_default=True)
@click.option('--scale_by_viral_load', is_flag=True,
              help='scale by viral load', show_default=True)
@click.option('--nboots', default=1000, help='Number of Bootstrap iterations', show_default=True)
@click.option('--serial_interval', default=5.5, help='Serial Interval', show_default=True)
@click.option('--config', default=None, help='path to yaml file', show_default=True)
@click.option('--mincov', default=60., help='min genome coverage included', show_default=True)
@click.option('--output', default='rel_growth_rates.csv',
              help='Output html file', show_default=True)
@click.option('--days', default=56, help='N Days used for growth calc', show_default=True)
@click.option('--grthresh', default=0.001, help='min avg prev. for growth', show_default=True)
@click.option('--lineageyml', default='-1', help='lineage hierarchy file', show_default=True)
def relgrowthrate(agg_results, metadata, thresh, scale_by_viral_load, nboots,
                  serial_interval, config, mincov, output, days, grthresh,
                  lineageyml):
    """
    generates relative growth rate for each lineage
    """
    agg_df = pd.read_csv(agg_results, skipinitialspace=True, sep='\t',
                         index_col=0)

    # drop poor quality samples
    if 'coverage' in agg_df.columns:
        agg_df = agg_df[agg_df['coverage'] > mincov]
    else:
        print('WARNING: Freyja should be updated ' +
              'to include coverage estimates.')
    agg_df = agg_df[agg_df['summarized'] != '[]']

    meta_df = pd.read_csv(metadata, index_col=0)
    meta_df['sample_collection_datetime'] = \
        pd.to_datetime(meta_df['sample_collection_datetime'])
    if config is not None:
        with open(config, "r") as f:
            try:
                config = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                raise ValueError('Error in config file: ' + str(exc))

    lineage_info = read_lineage_file(lineageyml, locDir)
    if config is not None:
        config = checkConfig(config)
    else:
        config = {}
    df_ab_lin, df_ab_sum, dates_to_keep = get_abundance(agg_df, meta_df,
                                                        thresh,
                                                        scale_by_viral_load,
                                                        config, lineage_info)
    calc_rel_growth_rates(df_ab_lin.copy(deep=True), nboots,
                          serial_interval, output, daysIncluded=days,
                          thresh=grthresh)


@cli.command()
@click.argument('query_mutations', type=click.Path(exists=True))
@click.argument('input_bam', type=click.Path(exists=True))
@click.option('--output', default='extracted.bam',
              help='path to save extracted reads', show_default=True)
@click.option('--same_read', is_flag=True,
              help='include to specify that query reads must all occur on the\
                    same read', show_default=True)
def extract(query_mutations, input_bam, output, same_read):
    """
    extracts reads containing mutation of interest
    """
    _extract(query_mutations, input_bam, output, same_read)


@cli.command()
@click.argument('query_mutations', type=click.Path(exists=True))
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('min_site', default=0)
@click.argument('max_site', default=29903)
@click.option('--output', default='filtered.bam',
              help='path to save filtered reads', show_default=True)
def filter(query_mutations, input_bam, min_site, max_site, output):
    """
    excludes reads containing one or more mutations
    """
    _filter(query_mutations, input_bam, min_site, max_site, output)


@cli.command()
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('min_site', default=0)
@click.argument('max_site', default=29903)
@click.option('--output', default='covariants.tsv',
              help='path to save co-occurring mutations', show_default=True)
@click.option('--ref-genome', type=click.Path(exists=True),
              default=os.path.join(locDir, 'data/NC_045512_Hu-1.fasta'), show_default=True)
@click.option('--gff-file', type=click.Path(exists=True),
              default=None,
              help=('path to gff file corresponding to reference genome. If '
                    'included, outputs amino acid mutations in addition to '
                    'nucleotide mutations.'), show_default=True)
@click.option('--min_quality', default=20,
              help='minimum quality for a base to be considered', show_default=True)
@click.option('--min_count', default=10,
              help='minimum count for a set of mutations to be saved', show_default=True)
@click.option('--spans_region', is_flag=True,
              help=('if included, consider only reads that span the region '
                    'defined by (min_site, max_site)'), show_default=True)
@click.option('--sort_by', default='count',
              help=('method by which to sort covariants patterns. Set to '
                    '"count" or "freq" to sort patterns by count or frequency '
                    '(in descending order). Set to "site" to sort patterns by '
                    'start site (n ascending order).'), show_default=True)
def covariants(input_bam, min_site, max_site, output,
               ref_genome, gff_file, min_quality, min_count, spans_region,
               sort_by):
    """
    finds co-variants (mutations co-occurring on the same read pair)
    """
    _covariants(input_bam, min_site, max_site, output,
                ref_genome, gff_file, min_quality, min_count, spans_region,
                sort_by)


@cli.command()
@click.argument('covariants', type=click.Path(exists=True))
@click.option('--output', default='covariants_plot.pdf', show_default=True)
@click.option('--num_clusters', default=10,
              help='number of clusters to plot (sorted by frequency)', show_default=True)
@click.option('--min_mutations', default=1,
              help='minimum number of mutations in a cluster to be plotted', show_default=True)
@click.option('--nt_muts', is_flag=True,
              help='if included, include nucleotide mutations in x-labels', show_default=True)
@click.option('--vmin', default=-5, help=('minimum value for colorbar'
                                          ' (log scale)'), show_default=True)
@click.option('--vmax', default=0, help=('maximum value for colorbar'
                                         ' (log scale)'), show_default=True)
def plot_covariants(covariants, output, num_clusters,
                    min_mutations, nt_muts, vmin, vmax):
    """Plot covariants output"""
    _plot_covariants(covariants, output, num_clusters,
                     min_mutations, nt_muts, vmin, vmax)


if __name__ == '__main__':
    cli()
