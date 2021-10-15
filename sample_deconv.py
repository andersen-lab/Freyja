import pandas as pd 
import numpy as np
import os

from tqdm import tqdm 
import datetime
import json
import sys

### get variant files. 

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

print('loading lineage models')
##################################################################33
#######build barcodes. 

df_barcodes = pd.read_csv('usher_barcodes.csv',index_col=0)

muts = list(df_barcodes.columns)

print('building mix matrix')
#########################################################################3

####################################################33
### grab wastewater samples
####################

dirName = sys.argv[1]
from collections import Counter
import re


### assemble matrix of (possibly) mixed samples

df_mix = pd.DataFrame()
df_depths = pd.DataFrame()


fnames = os.listdir(dirName)

filetype ='pileup.variants.tsv'
depthDir = sys.argv[2]#'test_depthfiles/'

for fn in tqdm(os.listdir(dirName)):
	df = pd.read_csv(dirName+fn,sep='\t')
	try:
		df_depth  = pd.read_csv(depthDir+fn.split(filetype)[0]+'sorted.depth',sep='\t',header=None,index_col=1)
	except:
		continue

	df['mutName'] = df['REF'] + df['POS'].astype(str) + df['ALT']### this only works for substitutions, but that's what we get from the usher tree
	df = df.drop_duplicates(subset='mutName')
	df.set_index('mutName',inplace=True)
	keptInds = [dfi for jdi, dfi in  enumerate(df.index) if dfi  in muts]#meh... works but gross
	df_mix = df_mix.append(pd.Series({kI:df.loc[kI,'ALT_FREQ'].astype(float) for kI in keptInds},name=fn))
	df_depths = df_depths.append(pd.Series({kI:df_depth.loc[int(re.findall(r'\d+',kI)[0]),3].astype(float) for kI in muts},name=fn))


#### reindex everything to match across the dfs
df_barcodes.drop(index=['20A', '20C', '20D', '20H(Beta,V2)', '21C(Epsilon)'],inplace=True)### drop Nextstrain clade names. 
df_barcodes= df_barcodes.reindex(sorted(df_barcodes.columns), axis=1)
df_mix = df_mix.groupby(axis=1,level=0).max()#drop this? 
df_mix = df_mix.reindex(df_barcodes.columns, axis=1).fillna(0.)


df_depths = df_depths.drop(labels=[m0 for m0 in muts if m0 not in df_mix.columns.to_list()],axis=1)
df_depths= df_depths.groupby(axis=1,level=0).max()#drop this? 
df_depths = df_depths.reindex(df_barcodes.columns, axis=1).fillna(0.)

print('demixing')
from scipy.optimize import nnls,lsq_linear
from sklearn import linear_model
import cvxpy as cp

sols = []
errors = []
strains = []
abundances = []
eps = 1e-5#minimum abundance to include
for i in tqdm(range(df_mix.shape[0])):
	dep= np.log2(df_depths.iloc[i]+1)
	dep=dep/np.max(dep)

	A = np.array((df_barcodes*dep).T)
	b = np.array(pd.to_numeric(df_mix.iloc[i])*dep)
	x = cp.Variable(A.shape[1])
	cost = cp.norm(A @ x - b,1)
	constraints = [sum(x)==1,x>=0]
	prob = cp.Problem(cp.Minimize(cost), constraints)

	prob.solve(verbose=False)
	sol = x.value
	rnorm = cp.norm(A @ x - b, 1).value
	sol[sol<eps] = 0
	nzInds = np.nonzero(sol)[0]
	newOrd = np.argsort(sol[nzInds])
	strains.append(df_barcodes.index[nzInds].to_list())
	#### summarize by constellation
	localDict = {}
	vals = sol[nzInds]/np.sum(sol[nzInds])


	localDict = sorted(localDict.items(), key=lambda x: x[1], reverse=True)
	abundances.append(sol[nzInds])
	summarized.append(localDict)
	sols.append(sol)
	errors.append(rnorm)



sols_df = pd.DataFrame(data=np.stack((strains,abundances),axis=1),index=[dmi.split('.')[0] for dmi in df_mix.index],columns=['lineages','abundances'])
sols_df['resid'] = errors
sols_df = sols_df.sort_index()
sols_df.to_csv(sys.argv[3])

