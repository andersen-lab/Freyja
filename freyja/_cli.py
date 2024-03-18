import os
import sys
import yaml

import click
import pandas as pd

locDir = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))


@click.group()
@click.version_option('1.5.0')
def cli():
    pass


def print_barcode_version(ctx, param, value):
    """
    Gets the barcode version used in the program

    Arguments:
     :param ctx: used for click command interface
     :param param: used for click command interface
     :param value: used for click command interface
     :return: date stamp of the barcode file used in the analysis
    """
    if not value or ctx.resilient_parsing:
        return
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    f = open(os.path.join(locDir, 'data/last_barcode_update.txt'), 'r')
    click.echo('Barcode version:')
    click.echo(f.readline())
    ctx.exit()


@cli.command(context_settings={'show_default': True})
@click.argument('variants', type=click.Path(exists=True))
@click.argument('depths', type=click.Path(exists=True))
@click.option('--eps', default=1e-3, help='minimum abundance to include')
@click.option('--barcodes', default='', help='custom barcode file')
@click.option('--meta', default='', help='custom lineage metadata file')
@click.option('--output', default='demixing_result.csv', help='Output file',
              type=click.Path(exists=False))
@click.option('--covcut', default=10, help='depth cutoff for\
                                            coverage estimate')
@click.option('--confirmedonly', is_flag=True, default=False)
@click.option('--version', is_flag=True, callback=print_barcode_version,
              expose_value=False, is_eager=True)
@click.option('--depthcutoff', default=0,
              help='exclude sites with coverage depth below this value and'
              'group identical barcodes')
@click.option('--lineageyml', default='', help='lineage hierarchy file')
@click.option('--adapt', default=0.,
              help='adaptive lasso penalty parameter')
@click.option('--a_eps', default=1E-8,
              help='adaptive lasso parameter, hard threshold')
@click.option('--region_of_interest', default='', help='JSON file containing'
              'region(s) of interest for which to compute additional coverage'
              'estimates')
@click.option('--relaxedmrca', is_flag=True, default=False,
              help='for use with depth cutoff,'
              'clusters are assigned robust mrca to handle outliers')
@click.option('--relaxedthresh', default=0.9,
              help='associated threshold for robust mrca function')
def demix(variants, depths, output, eps, barcodes, meta,
          covcut, confirmedonly, depthcutoff, lineageyml,
          adapt, a_eps, region_of_interest,
          relaxedmrca, relaxedthresh):
    """
    Generate prevalence of lineages per sample

    Arguments:
     :param variants: used to pass variant calling file generated
      using freyja variant command (tsv)
     :param depths: used to pass depth file generated using
      freyja variant command (tsv)
     :param output: used to pass output name, the default
      output will be named demixing_result.csv
     :param eps: float, if true, it is used to define
     minimum abundance of each lineage
     :param barcodes: used to pass a custom
      barcode file (csv)
     :param meta: used to pass a custom
      lineage to variant file
     :param covcut: int,if true,it is used to
      calculate percent of sites with n or greater reads
     :param confirmedonly: used to exclude unconfirmed lineages
     :param depthcutoff: used to exclude sites with
      coverage less than the specified value
     :param lineageyml: used to pass a custom lineage hierarchy file
     :param adapt: used to set adaptive lasso penalty parameter
     :param a_eps: used to set adaptive lasso
     penalty parameter hard threshold'
     :param relaxedmrca: for use with depth cutoff,
     clusters are assigned robust mrca to handle outliers
     :param relaxedthresh: associated threshold for robust mrca function
     :return : a tsv file that includes the
     lineages present,their corresponding abundances,
      and summarization by constellation.
    """
    from freyja.sample_deconv import (build_mix_and_depth_arrays,
                                      buildLineageMap,
                                      map_to_constellation,
                                      reindex_dfs, solve_demixing_problem)
    from freyja.utils import (collapse_barcodes, handle_region_of_interest)
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    # option for custom barcodes
    if barcodes != '':
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
                                        lineageyml, locDir, output,
                                        relaxedmrca, relaxedthresh)
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
    if region_of_interest != '':

        sols_df = handle_region_of_interest(region_of_interest, sols_df,
                                            df_depth, covcut, mix.name)

    # convert lineage/abundance readouts to single line strings
    sols_df['lineages'] = ' '.join(sols_df['lineages'])
    sols_df['abundances'] = ['%.8f' % ab for ab in sols_df['abundances']]
    sols_df['abundances'] = ' '.join(sols_df['abundances'])
    sols_df.to_csv(output, sep='\t')


@cli.command(context_settings={'show_default': True})
@click.option('--outdir', default='',
              help='Output directory save updated files')
@click.option('--noncl', is_flag=True, default=True,
              help='only include lineages in cov-lineages')
@click.option('--buildlocal', is_flag=True, default=False,
              help='Perform barcode building locally')
def update(outdir, noncl, buildlocal):
    """
    Updates the lineage information using the latest
    outbreak data

    Arguments:
     :param outdir: used to define the location for the
      data to be stored
     :param noncl: include proposed lineages or unconfirmed
     lineages (not in cov-lineages.org)
     :param buildlocal: calculate barcodes locally using USHER
     global phylogenetic tree
     :return : the most recent barcodes in json format
    """
    from freyja.updates import (convert_tree, download_barcodes,
                                download_tree, get_cl_lineages,
                                get_curated_lineage_data)
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                                          os.pardir))
    if outdir != '':
        # User specified directory
        locDir = outdir
    else:
        locDir = os.path.join(locDir, 'data')

    print('Getting outbreak data')
    get_curated_lineage_data(locDir)
    get_cl_lineages(locDir)
    # # get data from UShER
    if buildlocal:
        from freyja.convert_paths2barcodes import (check_mutation_chain,
                                                   convert_to_barcodes,
                                                   parse_tree_paths,
                                                   reversion_checking)
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


@cli.command(context_settings={'show_default': True})
@click.option('--pb', type=click.Path(exists=True),
              help='protobuf tree')
@click.option('--outdir', type=click.Path(exists=True),
              help='Output directory save updated files')
@click.option('--noncl', is_flag=True, default=True,
              help='only include lineages in cov-lineages')
def barcode_build(pb, outdir, noncl):
    """
    Building barcodes from a global tree

    Arguments:
     :param pb: used to provide a protobuf tree file
     :param outdir: used to define the location for the
      barcoe data to be stored
     :param noncl: include proposed lineages or unconfirmed
     lineages (not in cov-lineages.org)
     :return : a csv file containing barcodes for lineages
    """
    from freyja.convert_paths2barcodes import (check_mutation_chain,
                                               convert_to_barcodes,
                                               parse_tree_paths,
                                               reversion_checking)
    from freyja.updates import (convert_tree_custom,
                                get_cl_lineages,
                                get_curated_lineage_data)
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


@cli.command(context_settings={'show_default': True})
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
@click.option('--minq', help='Minimum base quality score',
              default=20)
@click.option('--annot', help='provide an annotation file in gff3 format',
              default='')
@click.option('--varthresh', help='Variant frequency threshold', default=0.0)
def variants(bamfile, ref, variants, depths, refname, minq, annot, varthresh):
    """
    Perform variant calling using samtools and ivar

    Arguments:
     :param bamfile: used to pass sorted bam file
     :param ref: used to pass a reference genome file
     :param variants: used to pass additional variant calling file
     :param depths: used to pass depth files
     :param refname: used to provide reference name for bam files with
     multiple refernces
     :param minq: used to provide minimum base quality
     :param annot: used to provide amino acid annotation file

     :return : a variant calling tsv file.
    """
    import subprocess
    if len(refname) == 0:
        bashCmd = f"samtools mpileup -aa -A -d 600000 -Q {minq} -q 0 -B -f "\
                  f"{ref} {bamfile} | tee >(cut -f1-4 > {depths}) |"\
                  f" ivar variants -p {variants} -q {minq} -t {varthresh}"\
                  f" -r {ref}"
    else:
        bashCmd = f"samtools mpileup -aa -A -d 600000 -Q {minq} -q 0 -B -f "\
                  f"{ref} {bamfile} -r {refname} | tee >(cut -f1-4 > {depths}"\
                  f") | ivar variants -p {variants} -q {minq} -t {varthresh}"\
                  f" -r {ref}"
    if len(annot) > 0:
        print('Including annotation')
        bashCmd = bashCmd + f" -g {annot}"
    print(bashCmd)
    sys.stdout.flush()  # force python to flush
    completed = subprocess.run(bashCmd, shell=True, executable="/bin/bash",
                               stdout=subprocess.PIPE)
    sys.exit(completed.returncode)


@cli.command(context_settings={'show_default': True})
@click.argument('variants', type=click.Path(exists=True))
@click.argument('depths', type=click.Path(exists=True))
@click.option('--nb', default=100, help='number of bootstraps')
@click.option('--nt', default=1, help='max number of cpus to use')
@click.option('--eps', default=1e-3, help='minimum abundance to include')
@click.option('--barcodes', default='', help='custom barcode file')
@click.option('--meta', default='', help='custom lineage metadata file')
@click.option('--output_base', default='test', help='Output file basename',
              type=click.Path(exists=False))
@click.option('--boxplot', default='',
              help='file format of boxplot output (e.g. pdf or png)')
@click.option('--confirmedonly', is_flag=True, default=False)
@click.option('--rawboots', is_flag=True, default=False,
              help='return raw bootstraps')
@click.option('--lineageyml', default='', help='lineage hierarchy file')
@click.option('--depthcutoff', default=0,
              help='exclude sites with coverage depth below this value and'
              'group identical barcodes')
@click.option('--relaxedmrca', is_flag=True, default=False,
              help='for use with depth cutoff,'
              'clusters are assigned robust mrca to handle outliers')
@click.option('--relaxedthresh', default=0.9,
              help='associated threshold for robust mrca function')
def boot(variants, depths, output_base, eps, barcodes, meta,
         nb, nt, boxplot, confirmedonly, lineageyml, depthcutoff,
         rawboots, relaxedmrca, relaxedthresh):
    """
    Perform bootstrapping method for freyja

    Arguments:
     :param variants: used to pass tsv format variant calling file
     :param depths: used to pass depth files
     :param variants: used to pass additional variant calling file
     :param nb: used to pass number of times bootstrapping is performed
     :param nts: used to provide number of cpus
     multiple refernces
     :param eps: used to provide minimum abundance
     :param barcodes: used to provide custom barcode file
     :param meta: used to provide custom lineage metadata file
     :param output_base: used to provide output file name
     :param boxplot: used to provide the format of the boxplot (pdf or png)
     :param confirmedonly: used to exclude unconfirmed lineages
     :param rawboots: used to return raw bootstraps values
     :param lineageyml: used to pass a custom lineage hierarchy file
     :param depthcutoff: used to exclude sites with coverage depth
     below this value andgroup identical barcodes
     :param relaxedmrca: for use with depth cutoff,
     clusters are assigned robust mrca to handle outliers
     :param relaxedthresh: associated threshold for robust mrca function

     :return : base-name_lineages.csv and base-name_summarized.csv
    """
    from freyja.sample_deconv import (build_mix_and_depth_arrays,
                                      buildLineageMap,
                                      perform_bootstrap,
                                      reindex_dfs)
    from freyja.utils import collapse_barcodes

    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    # option for custom barcodes
    if barcodes != '':
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
            lineageyml, locDir, output_base,
            relaxedmrca, relaxedthresh)

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
@click.argument('results', type=click.Path(exists=True))
@click.option('--ext', default='', help='file extension option')
@click.option('--output', default='aggregated_result.tsv', help='Output file',
              type=click.Path(exists=False))
def aggregate(results, ext, output):
    """
        Aggregates all the outputs

        Arguments:
         :param results: used to pass result files
         :param ext: used to pass extension of the files
         :param output: used to define output file name

         :return : an aggregated tsv file
        """
    import glob
    from freyja.utils import agg

    if ext != '':
        results_ = [fn for fn in glob.glob(results + '*' + ext)]
    else:
        results_ = [results + fn for fn in os.listdir(results)]
    df_demix = agg(results_)
    df_demix.to_csv(output, sep='\t')


@cli.command(context_settings={'show_default': True})
@click.argument('agg_results', type=click.Path(exists=True))
@click.option('--lineages', is_flag=True)
@click.option('--times', default='')
@click.option('--interval', default='MS')
@click.option('--config', default=None, help='path to yaml file')
@click.option('--mincov', default=60., help='min genome coverage included')
@click.option('--output', default='mix_plot.pdf', help='Output file')
@click.option('--windowsize', default=14)
@click.option('--lineageyml', default='', help='lineage hierarchy file')
@click.option('--thresh', default=0.01, help='min lineage abundance included')
@click.option('--writegrouped', default='',
              help='path to write grouped lineage data')
def plot(agg_results, lineages, times, interval, output, windowsize,
         config, mincov, lineageyml, thresh, writegrouped):
    """
        create plots using the outputs

        Arguments:
         :param agg_results: used to pass result files
         :param lineages: used to provide a lineage specific breakdown
         :param times: used to pass ample collection time information
         :param interval: used to pass ample collection time information
         :param windowsize: used to pass width of the rolling average window
         :param config: allows users to control the colors and grouping of
          lineages in the plot
         :param mincov: used to pass minimum genome coverage
         :param output: used to specify the output name
         :param lineageyml: used to pass a custom lineage hierarchy file
         :param thresh: used to pass a minum lineage abundance
         :param writegrouped: used for path to grouped lineage data
         :return : an aggregated tsv file
        """
    from freyja.utils import (checkConfig,
                              makePlot_simple, makePlot_time,
                              read_lineage_file)
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
    if times == '':
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
@click.option('--thresh', default=0.01, help='min lineage abundance included')
@click.option('--headerColor', default='#2fdcf5', help='color of header')
@click.option('--bodyColor', default='#ffffff', help='color of body')
@click.option('--scale_by_viral_load', is_flag=True,
              help='scale by viral load')
@click.option('--nboots', default=1000, help='Number of Bootstrap iterations')
@click.option('--serial_interval', default=5.5, help='Serial Interval')
@click.option('--config', default=None, help='path to yaml file')
@click.option('--mincov', default=60., help='min genome coverage included')
@click.option('--output', default='mydashboard.html', help='Output html file')
@click.option('--days', default=56, help='N Days used for growth calc')
@click.option('--grthresh', default=0.001, help='min avg prev. for growth')
@click.option('--lineageyml', default='', help='lineage hierarchy file')
@click.option('--keep_plot_files', is_flag=True, help='Keep separate plot')
def dash(agg_results, metadata, title, intro, thresh, headercolor, bodycolor,
         scale_by_viral_load, nboots, serial_interval, config, mincov, output,
         days, lineageyml, grthresh, keep_plot_files):
    agg_df = pd.read_csv(agg_results, skipinitialspace=True, sep='\t',
                         index_col=0)
    """
        create plots using the outputs

        Arguments:
         :param agg_results: used to pass result files
         :param metadata: used to pass lineage grouping config file
         :param title: used to pass dashboard title
         :param intro: used to pass dashboard introduction line
         :param thresh: used to pass a minimum lineage abundance
         :param headercolor: used to pass header color
         :param bodycolor: used to pass dashboard body color
         :param scale_by_viral_load: used to scale by viral load
         :param nboots: used to pass number of bootstrapping
         :param serial_interval: used to pass serial intervals
         :param config: allows users to control the colors and grouping of
          lineages in the plot
         :param mincov: used to pass minimum genome coverage
         :param output: used to specify the output name
         :param days: used to specify number of days for
         growth calculation
         :param lineageyml: used to pass a custom lineage hierarchy file
         :param grthresh: used to specify minimum prevalence
         to calculate relative growth rate for
         :param keep_plot_files: used to keep the intermediate html
         for the core plot

         :return : an interactive dashboard
        """
    from freyja.utils import (checkConfig, make_dashboard, read_lineage_file)
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


@cli.command(context_settings={'show_default': True})
@click.argument('agg_results', type=click.Path(exists=True))
@click.argument('metadata', type=click.Path(exists=True))
@click.option('--thresh', default=0.01, help='min lineage abundance in plot')
@click.option('--scale_by_viral_load', is_flag=True,
              help='scale by viral load')
@click.option('--nboots', default=1000, help='Number of Bootstrap iterations')
@click.option('--serial_interval', default=5.5, help='Serial Interval')
@click.option('--config', default=None, help='path to yaml file')
@click.option('--mincov', default=60., help='min genome coverage included')
@click.option('--output', default='rel_growth_rates.csv',
              help='Output html file')
@click.option('--days', default=56, help='N Days used for growth calc')
@click.option('--grthresh', default=0.001, help='min avg prev. for growth')
@click.option('--lineageyml', default='', help='lineage hierarchy file')
def relgrowthrate(agg_results, metadata, thresh, scale_by_viral_load, nboots,
                  serial_interval, config, mincov, output, days, grthresh,
                  lineageyml):
    from freyja.utils import (calc_rel_growth_rates, checkConfig,
                              get_abundance, read_lineage_file)
    agg_df = pd.read_csv(agg_results, skipinitialspace=True, sep='\t',
                         index_col=0)
    """
    generates relative growth rate for each lineage

    Arguments:
     :param agg_results: used to pass aggregated result files
     :param metadata: used to pass sample metadata such as
     collection date and viral load
     :param thresh: used to define min lineage abundance in plot
     :param scale_by_viral_load: used to scale by viral load
     :param nboots: Number of Bootstrap iterations
     :param serial_interval: used to specify serial interval
     :param config: used to pass lineage grouping config file
     in yml format
     :param mincov: used to specify minimum coverage to include
     :param days: used to specify number of days for
     growth calculation
     :param output: used to specify the output name
     :param grthresh: used to specify minimum prevalence
      to calculate relative growth rate for

     :return : a CSV file with lineage, estimated advantage
     and bootstrap 95% interval columns.
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
              help='path to save extracted reads')
@click.option('--same_read', is_flag=True,
              help='include to specify that query reads must all occur on the\
                    same read')
def extract(query_mutations, input_bam, output, same_read):
    """
    extracts reads containing mutation of interest

    Arguments:
     :param query_mutations: used to pass a set of mutations of
     interest
     :param input_bam: used to pass aligned reads in bam format
      to the reference genome
     :param same_read: used to only include reads with mutations
      of interest on the same read
     :param output: used to specify the output name

     :return : bam formatted file including reads
     with mutation of interest mutations
    """
    from freyja.read_analysis_tools import extract as _extract
    _extract(query_mutations, input_bam, output, same_read)


@cli.command(context_settings={'show_default': True})
@click.argument('query_mutations', type=click.Path(exists=True))
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('min_site', default=0)
@click.argument('max_site', default=29903)
@click.option('--output', default='filtered.bam',
              help='path to save filtered reads')
def filter(query_mutations, input_bam, min_site, max_site, output):
    """
    excludes reads containing one or more mutations

    Arguments:
     :param query_mutations: used to pass a set of mutations of
     interest
     :param input_bam: used to pass aligned reads in bam format
      to the reference genome
     :param min_site: consider read pairs that span the entire
     genomic region defined by (min_site, max_site)
     :param max_site: consider read pairs that span the entire
     genomic region defined by (min_site, max_site)
     :param output: used to specify the output name

     :return : bam formatted file not including reads with
      mutations of interest
    """
    from freyja.read_analysis_tools import filter as _filter

    _filter(query_mutations, input_bam, min_site, max_site, output)


@cli.command(context_settings={'show_default': True})
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('min_site', default=0)
@click.argument('max_site', default=29903)
@click.option('--output', default='covariants.tsv',
              help='path to save co-occurring mutations')
@click.option('--ref-genome', type=click.Path(exists=True),
              default=os.path.join(locDir, 'data/NC_045512_Hu-1.fasta'))
@click.option('--gff-file', type=click.Path(exists=True),
              default=None,
              help=('path to gff file corresponding to reference genome. If '
                    'included, outputs amino acid mutations in addition to '
                    'nucleotide mutations.'))
@click.option('--min_quality', default=20,
              help='minimum quality for a base to be considered')
@click.option('--min_count', default=10,
              help='minimum count for a set of mutations to be saved')
@click.option('--spans_region', is_flag=True,
              help=('if included, consider only reads that span the region '
                    'defined by (min_site, max_site)'))
@click.option('--sort_by', default='count',
              help=('method by which to sort covariants patterns. Set to '
                    '"count" or "freq" to sort patterns by count or frequency '
                    '(in descending order). Set to "site" to sort patterns by '
                    'start site (n ascending order).'))
def covariants(input_bam, min_site, max_site, output,
               ref_genome, gff_file, min_quality, min_count, spans_region,
               sort_by):
    """
    finds co-variants (mutations co-occurring on the same read pair)

    Arguments:
     :param input_bam: used to pass aligned reads in bam format
      to the reference genome
     :param min_site: consider read pairs that span the entire
     genomic region defined by (min_site, max_site)
     :param max_site: consider read pairs that span the entire
     genomic region defined by (min_site, max_site)
     :param output: used to define the name and path for the
     output
     :param ref_genome: used to specify reference genome file
     default reference genome is NC_045512_Hu-1
     :param gff_file: path to corresponding genome annotation file
     :param min_quality: used to define minium quality for a base
     :param min_count: used to define minimum count for
      a set of mutations to be saved
     :param spans_region: if included, consider only reads
      that span the region defined by (min_site, max_site)
     :param sort_by: method by which to sort covariants patterns

     :return : a tsv formatted file of covariants
    """
    from freyja.read_analysis_tools import covariants as _covariants
    _covariants(input_bam, min_site, max_site, output,
                ref_genome, gff_file, min_quality, min_count, spans_region,
                sort_by)


@cli.command(context_settings={'show_default': True})
@click.argument('covariants', type=click.Path(exists=True))
@click.option('--output', default='covariants_plot.pdf')
@click.option('--num_clusters', default=10,
              help='number of clusters to plot (sorted by frequency)')
@click.option('--min_mutations', default=1,
              help='minimum number of mutations in a cluster to be plotted')
@click.option('--nt_muts', is_flag=True,
              help='if included, include nucleotide mutations in x-labels')
@click.option('--vmin', default=-5, help=('minimum value for colorbar'
                                          ' (log scale)'))
@click.option('--vmax', default=0, help=('maximum value for colorbar'
                                         ' (log scale)'))
def plot_covariants(covariants, output, num_clusters,
                    min_mutations, nt_muts, vmin, vmax):
    from freyja.read_analysis_tools import plot_covariants as _plot_covariants
    _plot_covariants(covariants, output, num_clusters,
                     min_mutations, nt_muts, vmin, vmax)


if __name__ == '__main__':
    cli()
