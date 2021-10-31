import pandas as pd 
import numpy as np
import os
from tqdm import tqdm 
import datetime
import json
import sys
import re
import cvxpy as cp

### get variant files. 
def buildLineageMap():
	#### make sure to pull a new one periodically!
	f0 = open('freyja/data/curated_lineages.json')
	dat = json.load(f0)
	f0.close()

	sampleFreq = 'D'#sampleFreq
	mapDict = {}
	for l in range(len(dat)):
		if 'who_name' in dat[l].keys():
			for d0 in dat[l]['pango_descendants']:
				if dat[l]['who_name'] != None:
					mapDict[d0] = dat[l]['who_name']
	return mapDict


def build_mix_and_depth_arrays(fn, depthFn,muts):
	df = pd.read_csv(fn,sep='\t')
	df_depth  = pd.read_csv(depthFn,sep='\t',header=None,index_col=1)
	df['mutName'] = df['REF'] + df['POS'].astype(str) + df['ALT']### this only works for substitutions, but that's what we get from the usher tree
	df = df.drop_duplicates(subset='mutName')
	df.set_index('mutName',inplace=True)
        keptInds = set(muts) & set(df.index)
	mix = df.loc[keptInds, 'ALT_FREQ'].astype(float)
	mix.name = fn
	depths = pd.Series({kI:df_depth.loc[int(re.findall(r'\d+',kI)[0]),3].astype(float) for kI in muts},name=fn)
	return mix,depths


def reindex_dfs(df_barcodes,mix,depths):
	#### reindex everything to match across the dfs
	df_barcodes = df_barcodes.drop(index=['20A', '20C', '20D', '20H(Beta,V2)', '21C(Epsilon)'])### drop Nextstrain clade names. 
	df_barcodes= df_barcodes.reindex(sorted(df_barcodes.columns), axis=1)
	mix = mix.reindex(df_barcodes.columns).fillna(0.)

        mix_as_set = set(mix.index)
	depths = depths.drop(labels=[m0 for m0 in df_barcodes.columns if m0 not in mix_as_set])# dropping extra sequencing depth info we don't need
	depths = depths.reindex(df_barcodes.columns).fillna(0.)
	return df_barcodes,mix,depths




def map_to_constellation(sample_strains,vals,mapDict):
	#maps lineage names to constellations
	localDict = {}
	for jj,lin in enumerate(sample_strains):
		if lin in mapDict.keys():
			if mapDict[lin] not in localDict.keys():
				localDict[mapDict[lin]] = vals[jj]
			else: 
				localDict[mapDict[lin]] += vals[jj]
		elif 'A.' in lin or lin=='A':
			if 'Aaron' not in localDict.keys():
				localDict['Aaron'] = vals[jj]
			else: 
				localDict['Aaron'] += vals[jj]
		else:
			if 'Other' not in localDict.keys():
				localDict['Other']=vals[jj]
			else:
				localDict['Other']+=vals[jj]
	#convert to descending order
	localDict = sorted(localDict.items(), key=lambda x: x[1], reverse=True)
	return localDict


def solve_demixing_problem(df_barcodes,mix,depths):
	### single file problem setup, solving
	eps = 1e-5#minimum abundance to include---TO DO:allow user to modify
	dep= np.log2(depths+1)
	dep=dep/np.max(dep)#normalize depth scaling pre-optimization

	# set up and solve demixing problem
	A = np.array((df_barcodes*dep).T)
	b = np.array(pd.to_numeric(mix)*dep)
	x = cp.Variable(A.shape[1])
	cost = cp.norm(A @ x - b,1)
	constraints = [sum(x)==1,x>=0]
	prob = cp.Problem(cp.Minimize(cost), constraints)
	prob.solve(verbose=False)
	sol = x.value
	rnorm = cp.norm(A @ x - b, 1).value
	#extract lineages with non-negligible abundance
	sol[sol<eps] = 0
	nzInds = np.nonzero(sol)[0]
	newOrd = np.argsort(sol[nzInds])
	sample_strains = df_barcodes.index[nzInds].to_list()
	# vals = sol[nzInds]/np.sum(sol[nzInds])### that 
	abundances = sol[nzInds]
	return sample_strains,abundances,rnorm

if __name__ == '__main__':
	print('loading lineage models')
	##################################################################33
	#######read in  barcodes. 
	df_barcodes = pd.read_csv('usher_barcodes.csv',index_col=0)
	muts = list(df_barcodes.columns)
	mapDict = buildLineageMap()
	print('building mix/depth matrices')

	####################################################33
	### grab wastewater sample data
	####################

	varFn = sys.argv[1]#variant file
	depthFn = sys.argv[2]#depth file
	### assemble data from of (possibly) mixed samples
	df_mix,df_depths = build_mix_and_depth_arrays(varFn,depthFn,muts)
	print('demixing')
	sample_strains,abundances,error = solve_demixing_problem(df_barcodes,df_mix,df_depths)
	localDict = mapToConstellation(sample_strains,abundances,mapDict)

	### assemble into df and write. 
	sols_df = pd.DataFrame(data=np.stack((localDict,sample_strains,abundances,error),axis=1),\
	index=[dmi.split('.')[0] for dmi in df_mix.index],\
	columns=['summarized','lineages','abundances','resid'])
	sols_df = sols_df.sort_index()
	sols_df.to_csv(sys.argv[3],sep='\t')
