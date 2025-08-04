import pandas as pd
import numpy as np
import json
import sys
import re
import cvxpy as cp
import os
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from tqdm import tqdm
import matplotlib
import glob
from collections import defaultdict


def buildLineageMap(locDir):
    # Parsing curated lineage data from outbreak.info
    mapDict = {}
    if os.path.isdir(locDir):
        global_dat = defaultdict(set)
        for filename in glob.glob(os.path.join(locDir, "*.json")):
            with open(filename, "r") as f0:
                dat = json.load(f0)
                dat = sorted(dat, key=lambda x: len(
                    x.get('pango_descendants', [])), reverse=True)
                for record in dat:
                    if 'who_name' in record:
                        if record['who_name'] is not None:
                            global_dat[record['who_name']].update(
                                record.get('pango_descendants', []))
        # Sort global_dat by number of descendants (largest first)
        # rather than membership checking, more specific
        # designations will overwrite broader, ancestral ones
        sorted_dat = sorted(global_dat.items(),
                            key=lambda x: len(x[1]),
                            reverse=True)
        for clade, descendants in sorted_dat:
            for descendant in descendants:
                mapDict[descendant] = clade
    else:
        if locDir == '':
            locDir = os.path.abspath(
                os.path.join(os.path.realpath(__file__), os.pardir))
            with open(os.path.join(
                    locDir, 'data/curated_lineages.json')) as f0:
                dat = json.load(f0)
        else:
            with open(locDir) as f0:
                dat = json.load(f0)
        # Sort records based on the length of 'pango_descendants'
        dat = sorted(dat, key=lambda x: len(
            x.get('pango_descendants', [])), reverse=True)
        for ind in range(len(dat)):
            if 'who_name' in dat[ind].keys():
                for d0 in dat[ind]['pango_descendants']:
                    if dat[ind]['who_name'] is not None:
                        mapDict[d0] = dat[ind]['who_name']
    return mapDict


def get_error_rate(muts, df0, df_depths0):
    from freyja.convert_paths2barcodes import sortFun
    sites = set([sortFun(m) for m in muts])
    df_depths0 = df_depths0[~df_depths0.index.isin(sites)]
    df0 = df0[~df0['POS'].isin(sites)]
    df0 = df0[df0['ALT_FREQ'] > 0]
    error_rate = df0['ALT_FREQ'].median()
    return error_rate


def build_mix_and_depth_arrays(fn, depthFn, muts, covcut, autoadapt):
    input_is_vcf = fn.lower().endswith('vcf')
    if input_is_vcf:
        df = read_snv_frequencies_vcf(fn, depthFn, muts)
    else:
        df = read_snv_frequencies_ivar(fn, depthFn, muts)
        df = df[['REF', 'POS', 'ALT', 'ALT_FREQ', 'ALT_DP']]

    # only works for substitutions, but that's what we get from usher tree
    df_depth = pd.read_csv(depthFn, sep='\t', header=None, index_col=1)[[3]]

    if autoadapt:
        adapt = get_error_rate(muts, df, df_depth)
    else:
        adapt = 0

    df['mutName'] = df['REF'] + df['POS'].astype(str) + df['ALT']
    df = df.drop_duplicates(subset='mutName')
    df.set_index('mutName', inplace=True)
    keptInds = set(muts) & set(df.index)
    mix = df.loc[list(keptInds), 'ALT_FREQ'].astype(float)
    mix.name = fn
    depths = pd.Series({kI: df_depth.loc[int(re.findall(r'\d+', kI)[0]), 3]
                        .astype(float) for kI in muts}, name=fn)
    coverage = 100.*np.sum(df_depth.loc[:, 3] >= covcut)/df_depth.shape[0]
    return mix, depths, coverage, adapt


def read_snv_frequencies_ivar(fn, depthFn, muts):
    df = pd.read_csv(fn, sep='\t')
    return df


def read_snv_frequencies_vcf(fn, depthFn, muts):
    vcfnames = []
    with open(fn, "r") as file:
        for line in file:
            if line.startswith("#CHROM"):
                vcfnames = [x for x in line.strip("\n").split('\t')]
                break
    file.close()

    df = pd.read_csv(fn, comment='#', sep=r'\s+',
                     header=None,
                     names=vcfnames)
    df["ALT_FREQ"] = df["INFO"].str.extract(r"AF=([0-9]*\.*[0-9]+)")
    if df["ALT_FREQ"].isnull().all():
        raise ValueError(
            f"No AF (allele frequency) data found in the INFO column of "
            f"VCF file: {fn} or the VCF file is empty")
    df["ALT_FREQ"] = pd.to_numeric(df["ALT_FREQ"], downcast="float")
    return df


def reindex_dfs(df_barcodes, mix, depths):
    # first, drop Nextstrain clade names.
    nxNames = df_barcodes.index[df_barcodes.index.str[0].str.isdigit()]
    df_barcodes = df_barcodes.drop(index=nxNames)
    # reindex everything to match across the dfs
    df_barcodes = df_barcodes.reindex(sorted(df_barcodes.columns), axis=1)
    mix = mix.reindex(df_barcodes.columns).fillna(0.)

    mix_as_set = set(mix.index)
    # dropping extra sequencing depth info we don't need
    depths = depths.drop(labels=[m0 for m0 in df_barcodes.columns
                                 if m0 not in mix_as_set])
    depths = depths.reindex(df_barcodes.columns).fillna(0.)
    return df_barcodes, mix, depths


def map_to_constellation(sample_strains, vals, mapDict):
    # maps lineage names to constellations
    localDict = {}
    for jj, lin0 in enumerate(sample_strains):
        if '-like' in lin0:
            lin = lin0.split('-like')[0]
        else:
            lin = lin0
        if lin in mapDict.keys():
            if mapDict[lin] not in localDict.keys():
                localDict[mapDict[lin]] = vals[jj]
            else:
                localDict[mapDict[lin]] += vals[jj]
        elif lin.startswith('A.') or lin == 'A':
            if 'A' not in localDict.keys():
                localDict['A'] = vals[jj]
            else:
                localDict['A'] += vals[jj]
        else:
            if 'Other' not in localDict.keys():
                localDict['Other'] = vals[jj]
            else:
                localDict['Other'] += vals[jj]
    # convert to descending order
    localDict = sorted(localDict.items(), key=lambda x: x[1], reverse=True)

    # cast np.float64 to float
    localDict = [(x[0], float(x[1])) for x in localDict]

    return localDict


def solve_demixing_problem(df_barcodes, mix, depths, eps, adapt,
                           a_eps, solver, max_threads=0, verbose=False):
    # single file problem setup, solving

    dep = np.log2(depths+1)
    dep = dep/np.max(dep)  # normalize depth scaling pre-optimization

    # set up and solve demixing problem
    A = np.array((df_barcodes*dep).T)
    b = np.array(pd.to_numeric(mix)*dep)
    x = cp.Variable(A.shape[1])
    cost = cp.norm(A @ x - b, 1)
    constraints = [sum(x) == 1, x >= 0]
    prob = cp.Problem(cp.Minimize(cost), constraints)

    solver_kwargs = {}
    if solver == 'ECOS':
        solver_ = cp.ECOS
    elif solver == 'OSQP':
        solver_ = cp.OSQP
    else:
        solver_ = cp.CLARABEL
        # Setup CLARABEL thread limit if necessary
        if int(max_threads) != 0:
            solver_kwargs = {"max_threads": max_threads}

    try:
        prob.solve(verbose=verbose, solver=solver_, **solver_kwargs)
    except cp.error.SolverError:
        raise ValueError('Solver error encountered, most '
                         'likely due to insufficient sequencing depth.'
                         'Try running with --depthcutoff.')
    sol = x.value
    rnorm0 = cp.norm(A @ x - b, 1).value
    if adapt > 0.:
        # if specified, run the adaptive lasso based method
        # thresholding step
        solNz = [j for j, s in enumerate(sol) if s > a_eps]
        solFlip = [1./sol[j] for j in solNz]
        A = A[:, solNz]
        x = cp.Variable(A.shape[1])
        cost = cp.norm(A @ x - b, 1) +\
            (rnorm0*adapt)*cp.norm(cp.multiply(solFlip, x), 1)
        constraints = [sum(x) == 1, x >= 0]
        prob = cp.Problem(cp.Minimize(cost), constraints)
        try:
            prob.solve(verbose=verbose, solver=solver_, **solver_kwargs)
        except cp.error.SolverError:
            raise ValueError('Solver error encountered, most '
                             'likely due to insufficient sequencing depth.'
                             'Try running with --depthcutoff.')
        sol = np.zeros(len(sol))
        sol[solNz] = x.value

    rnorm = cp.norm(A @ x - b, 1).value
    # extract lineages with non-negligible abundance
    sol[sol < eps] = 0
    nzInds = np.nonzero(sol)[0]
    sample_strains = df_barcodes.index[nzInds].to_numpy()
    abundances = sol[nzInds]
    # sort strain/abundance lists in order of decreasing abundance
    indSort = np.argsort(abundances)[::-1]
    return sample_strains[indSort], abundances[indSort], rnorm


def bootstrap_parallel(jj, samplesDefining, fracDepths_adj, mix_grp,
                       mix, df_barcodes, eps0, muts, mapDict,
                       adapt, a_eps, bootseed, solver):
    # helper function for fast bootstrap and solve
    # get sequencing depth at the position of all defining mutations
    mix_boot = mix.copy()
    seed = bootseed+jj
    rng = np.random.default_rng()
    # get the BitGenerator used by default_rng
    BitGen = type(rng.bit_generator)
    # use the state from a fresh bit generator
    rng.bit_generator.state = BitGen(seed).state

    dps = pd.Series(rng.multinomial(samplesDefining[jj],
                    fracDepths_adj, size=1)[0, :], index=fracDepths_adj.index)
    # get number of reads of each possible nucleotide
    # only for lineage defining muts observed in the dataset
    for mp in mix_grp.index:
        if len(mix_grp[mp]) == 1:
            mut0 = mix_grp[mp][0]
            if dps[mp] > 0:
                mix_boot.loc[mut0] = rng.binomial(
                    dps[mp],
                    mix.loc[mut0])/dps[mp]
            else:
                mix_boot.loc[mut0] = 0.
        elif len(mix_grp[mp]) > 1:  # if multiple muts at a single site
            # TODO: streamline-- right now looks at all positions
            probs = [mix.loc[mut0] for mut0 in mix_grp[mp]]
            probs.append(np.max([1.-np.sum(probs), 0]))
            if np.sum(probs) > 1.0:
                # correct rounding errors
                probs = np.array(probs)/np.sum(probs)
            altCounts = rng.multinomial(dps[mp], probs,
                                        size=1)[0, 0:(len(probs)-1)]
            for j, mut0 in enumerate(mix_grp[mp]):
                if dps[mp] > 0:
                    mix_boot.loc[mut0] = \
                        float(altCounts[j])/dps[mp]
                else:
                    mix_boot.loc[mut0] = 0.
    dps_ = pd.Series({kI:
                      dps[int(kI[1:(len(kI)-1)])].astype(float)
                      for kI in muts}, name='depths')
    df_barcodes, mix_boot_, dps_ = reindex_dfs(df_barcodes,
                                               mix_boot, dps_)
    sample_strains, abundances, error = solve_demixing_problem(df_barcodes,
                                                               mix_boot_,
                                                               dps_, eps0,
                                                               adapt, a_eps,
                                                               solver)
    localDict = map_to_constellation(sample_strains, abundances, mapDict)
    return sample_strains, abundances, localDict


def perform_bootstrap(df_barcodes, mix, depths_,
                      numBootstraps, eps0, n_jobs,
                      mapDict, muts, boxplot, basename, bootseed,
                      solver, adapt=0., a_eps=1E-8):
    depths_.index = depths_.index.to_series().apply(lambda x:
                                                    int(x[1:len(x)-1]))
    depths_ = depths_[~depths_.index.duplicated(keep='first')]
    totalDepth = depths_.sum()
    fracDepths = depths_/totalDepth

    fracDefining = fracDepths.sum()
    fracDepths_adj = fracDepths/fracDefining

    seed = bootseed
    rng = np.random.default_rng()
    # get the BitGenerator used by default_rng
    BitGen = type(rng.bit_generator)
    # use the state from a fresh bit generator
    rng.bit_generator.state = BitGen(seed).state
    # get the total depth at lineage defining positions
    samplesDefining = rng.binomial(totalDepth, fracDefining, numBootstraps)
    mixPos = pd.Series(mix.index,
                       index=mix.index.to_series()
                                      .apply(lambda x:
                                             int(x[1: len(x)-1])))
    mix_grp = mixPos.groupby(level=0).apply(list)
    lin_df = pd.DataFrame()
    constell_df = pd.DataFrame()
    out = Parallel(n_jobs=n_jobs)(delayed(bootstrap_parallel)(jj0,
                                                              samplesDefining,
                                                              fracDepths_adj,
                                                              mix_grp,
                                                              mix,
                                                              df_barcodes,
                                                              eps0,
                                                              muts,
                                                              mapDict,
                                                              adapt,
                                                              a_eps,
                                                              bootseed,
                                                              solver)
                                  for jj0 in tqdm(range(numBootstraps)))
    for i in range(len(out)):
        sample_lins, abundances, localDict = out[i]

        # merge intra-lineage diversity if multiple hits.
        if len(set(sample_lins)) < len(sample_lins):
            localDict = {}
            for jj, lin in enumerate(sample_lins):
                if lin not in localDict.keys():
                    localDict[lin] = abundances[jj]
                else:
                    localDict[lin] += abundances[jj]
            # ensure descending order
            localDict = dict(sorted(localDict.items(),
                                    key=lambda x: x[1],
                                    reverse=True))
            sample_lins = list(localDict.keys())
            abundances = list(localDict.values())

        lin_df = pd.concat(
            [lin_df,
             pd.DataFrame({sample_lins[j]: abundances[j]
                           for j in range(len(sample_lins))},
                          index=[0])],
            axis=0, join='outer', ignore_index=True)
        constell_df = pd.concat(
            [constell_df,
             pd.DataFrame({localDict[j][0]: localDict[j][1]
                           for j in range(len(localDict))},
                          index=[0])],
            axis=0, join='outer', ignore_index=True)
    lin_df = lin_df.fillna(0)
    constell_df = constell_df.fillna(0)

    if len(boxplot) > 0:
        if boxplot == 'pdf':
            matplotlib.rcParams['pdf.fonttype'] = 42
            matplotlib.rcParams['ps.fonttype'] = 42
        fig, ax = plt.subplots()
        lin_df.boxplot(ax=ax, rot=90)
        ax.set_ylim([0, 1])
        ax.set_aspect(2)
        ax.set_ylabel('Variant Prevalence')
        fig.tight_layout()
        fig.savefig(basename+'_lineages.'+boxplot)

        fig, ax = plt.subplots()
        constell_df.boxplot(ax=ax, rot=90)
        ax.set_ylim([0, 1])
        ax.set_aspect(2)
        ax.set_ylabel('Variant Prevalence')
        fig.tight_layout()
        fig.savefig(basename+'_summarized.'+boxplot)
    return lin_df, constell_df


if __name__ == '__main__':
    print('loading lineage models')
    # read in  barcodes.
    df_barcodes = pd.read_csv('freyja/data/usher_barcodes.csv', index_col=0)
    muts = list(df_barcodes.columns)
    mapDict = buildLineageMap()

    # grab wastewater sample data

    variants = sys.argv[1]  # variant file
    depths = sys.argv[2]  # depth file
    # assemble data from of (possibly) mixed samples
    muts = list(df_barcodes.columns)
    eps = 0.001
    mapDict = buildLineageMap()
    print('building mix/depth matrices')
    # assemble data from of (possibly) mixed samples
    mix, depths_, cov = build_mix_and_depth_arrays(variants, depths, muts)
    print('demixing')
    df_barcodes, mix, depths_ = reindex_dfs(df_barcodes, mix, depths_)
    sample_strains, abundances, error = solve_demixing_problem(df_barcodes,
                                                               mix,
                                                               depths_, eps)
    localDict = map_to_constellation(sample_strains, abundances, mapDict)
    # assemble into series and write.
    sols_df = pd.Series(data=(localDict, sample_strains, abundances, error),
                        index=['summarized', 'lineages',
                               'abundances', 'resid'],
                        name=mix.name)
    numBootstraps = 100
    n_jobs = 10
    lin_out, constell_out = perform_bootstrap(df_barcodes, mix, depths_,
                                              numBootstraps, eps,
                                              n_jobs, mapDict, muts)
