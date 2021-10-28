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
	f0 = open('curated_lineages.json')
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


def build_mix_and_depth_arrays(dirName,depthDir,muts):
	df_mix = pd.DataFrame()
	df_depths = pd.DataFrame()
	fnames = os.listdir(dirName)

	filetype =fnames[0].split('.')[-1]
	for fn in tqdm(os.listdir(dirName)):
		df = pd.read_csv(dirName+fn,sep='\t')
		try:
			df_depth  = pd.read_csv(depthDir+fn.split(filetype)[0]+'depth',sep='\t',header=None,index_col=1)
		except:
			continue
		df['mutName'] = df['REF'] + df['POS'].astype(str) + df['ALT']### this only works for substitutions, but that's what we get from the usher tree
		df = df.drop_duplicates(subset='mutName')
		df.set_index('mutName',inplace=True)
		keptInds = [dfi for jdi, dfi in  enumerate(df.index) if dfi  in muts]#meh... works but should be improved
		df_mix = df_mix.append(pd.Series({kI:df.loc[kI,'ALT_FREQ'].astype(float) for kI in keptInds},name=fn))
		df_depths = df_depths.append(pd.Series({kI:df_depth.loc[int(re.findall(r'\d+',kI)[0]),3].astype(float) for kI in muts},name=fn))
	return df_mix,df_depths


def reindex_dfs(df_barcodes,df_mix,df_depths):
	#### reindex everything to match across the dfs
	df_barcodes = df_barcodes.drop(index=['20A', '20C', '20D', '20H(Beta,V2)', '21C(Epsilon)'])### drop Nextstrain clade names. 
	df_barcodes= df_barcodes.reindex(sorted(df_barcodes.columns), axis=1)
	df_mix = df_mix.reindex(df_barcodes.columns, axis=1).fillna(0.)

	df_depths = df_depths.drop(labels=[m0 for m0 in muts if m0 not in df_mix.columns.to_list()],axis=1)# dropping extra sequencing depth info we don't need
	df_depths = df_depths.reindex(df_barcodes.columns, axis=1).fillna(0.)
	return df_barcodes,df_mix,df_depths




def map_to_constellation(vals,sample_strains,mapDict):
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

	localDict = sorted(localDict.items(), key=lambda x: x[1], reverse=True)#convert to descending order
	return localDict


def solve_demixing_problem(df_barcodes,df_mix,df_depths,mapDict):
	sols = []
	errors = []
	strains = []
	summarized = []
	abundances = []
	eps = 1e-5#minimum abundance to include---TO DO:allow user to modify
	for i in tqdm(range(df_mix.shape[0])):
		dep= np.log2(df_depths.iloc[i]+1)
		dep=dep/np.max(dep)#normalize depth scaling pre-optimization

		# set up and solve demixing problem
		A = np.array((df_barcodes*dep).T)
		b = np.array(pd.to_numeric(df_mix.iloc[i])*dep)
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
		vals = sol[nzInds]/np.sum(sol[nzInds])
		strains.append(sample_strains)

		#### summarize by constellation
		localDict = mapToConstellation(vals,sample_strains,mapDict)

		abundances.append(sol[nzInds])
		summarized.append(localDict)
		sols.append(sol)
		errors.append(rnorm)
	sols_df = pd.DataFrame(data=np.stack((summarized,strains,abundances,errors),axis=1),index=[dmi.split('.')[0] for dmi in df_mix.index],columns=['summarized','lineages','abundances','resid'])
	sols_df = sols_df.sort_index()
	return sols_df




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

	dirName = sys.argv[1]#variant files dir
	depthDir = sys.argv[2]#depth files dir
	### assemble data from of (possibly) mixed samples
	df_mix,df_depths = build_mix_and_depth_arrays(dirName,depthDir,muts)
	print('demixing')
	sols_df = solve_demixing_problem(df_barcodes,df_mix,df_depths)
	sols_df.to_csv(sys.argv[3],sep='\t')
