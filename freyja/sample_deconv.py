import pandas as pd
import numpy as np
import json
import sys
import re
import cvxpy as cp
import os


def buildLineageMap():
    # Parsing curated lineage data from outbreak.info
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    f0 = open(os.path.join(locDir, 'data/curated_lineages.json'))
    dat = json.load(f0)
    f0.close()

    mapDict = {}
    for ind in range(len(dat)):
        if 'who_name' in dat[ind].keys():
            for d0 in dat[ind]['pango_descendants']:
                if dat[ind]['who_name'] is not None:
                    mapDict[d0] = dat[ind]['who_name']
    return mapDict


def build_mix_and_depth_arrays(fn, depthFn, muts):
    df = pd.read_csv(fn, sep='\t')
    df_depth = pd.read_csv(depthFn, sep='\t', header=None, index_col=1)
    # only works for substitutions, but that's what we get from usher tree
    df['mutName'] = df['REF'] + df['POS'].astype(str) + df['ALT']
    df = df.drop_duplicates(subset='mutName')
    df.set_index('mutName', inplace=True)
    keptInds = set(muts) & set(df.index)
    mix = df.loc[keptInds, 'ALT_FREQ'].astype(float)
    mix.name = fn
    depths = pd.Series({kI: df_depth.loc[int(re.findall(r'\d+', kI)[0]), 3]
                        .astype(float) for kI in muts}, name=fn)
    return mix, depths


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
    for jj, lin in enumerate(sample_strains):
        if lin in mapDict.keys():
            if mapDict[lin] not in localDict.keys():
                localDict[mapDict[lin]] = vals[jj]
            else:
                localDict[mapDict[lin]] += vals[jj]
        elif 'A.' in lin or lin == 'A':
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
    return localDict


def solve_demixing_problem(df_barcodes, mix, depths, eps):
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
    prob.solve(verbose=False)
    sol = x.value
    rnorm = cp.norm(A @ x - b, 1).value
    # extract lineages with non-negligible abundance
    sol[sol < eps] = 0
    nzInds = np.nonzero(sol)[0]
    sample_strains = df_barcodes.index[nzInds].to_numpy()
    abundances = sol[nzInds]
    # sort strain/abundance lists in order of decreasing abundance
    indSort = np.argsort(abundances)[::-1]
    return sample_strains[indSort], abundances[indSort], rnorm


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
    mapDict = buildLineageMap()
    print('building mix/depth matrices')
    # assemble data from of (possibly) mixed samples
    mix, depths_ = build_mix_and_depth_arrays(variants, depths, muts)
    print('demixing')
    df_barcodes, mix, depths_ = reindex_dfs(df_barcodes, mix, depths_)
    sample_strains, abundances, error = solve_demixing_problem(df_barcodes,
                                                               mix, depths_)
    localDict = map_to_constellation(sample_strains, abundances, mapDict)
    # assemble into series and write.
    sols_df = pd.Series(data=(localDict, sample_strains, abundances, error),
                        index=['summarized', 'lineages',
                               'abundances', 'resid'],
                        name=mix.name)
