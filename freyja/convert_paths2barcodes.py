from copy import deepcopy
import numpy as np
import pandas as pd 
import sys

def parse_tree_paths(df):
    df = df.set_index('clade')
    df = df.drop_duplicates(keep='last')### check which lineages got trimmed.
    df['from_tree_root'] = df['from_tree_root'].fillna('')
    df['from_tree_root'] = df['from_tree_root'].apply(lambda x:x.replace(' ','').strip('>').split('>'))
    return df

def sortFun(x):
	#sort based on nuc position, ignoring nuc identities
	return int(x[1:(len(x)-1)])

def convert_to_barcodes(df):
	### builds simple barcodes, not accounting for reversions
	df_barcodes = pd.DataFrame()
	for clade in df.index:
		cladeSeries = pd.Series({c:1 for c in df.loc[clade,'from_tree_root']},name=clade) #sparse,binary encoding
		df_barcodes = df_barcodes.append(cladeSeries)

	print('separating combined splits')
	df_barcodes = df_barcodes.drop(columns='')
	temp = pd.DataFrame()
	dropList = []
	for c in df_barcodes.columns:
		### if column includes multiple mutations, split into separate columns and concatenate
		if "," in c:
			for mt in c.split(","):
				temp[mt] = df_barcodes[c]
			dropList.append(c)
	df_barcodes = df_barcodes.drop(columns=dropList)
	df_barcodes = pd.concat((df_barcodes,temp),axis=1)
	df_barcodes = df_barcodes.fillna(0)
	df_barcodes = df_barcodes.groupby(axis=1,level=0).sum()
	return df_barcodes

def reversion_checking(df_barcodes):
	print('checking for mutation pairs')

	###check if a reversion is present. 
	flips = [d[-1] +d[1:len(d)-1]+d[0] for d in df_barcodes.columns[1:]]
	flipPairs = [d for d in df_barcodes.columns[1:] if d in flips]
	flipPairs.sort(key=sortFun)
	flipPairs = [[flipPairs[j],flipPairs[j+1]] for j in np.arange(0,len(flipPairs),2)]

	### subtract lower of two pair counts from both to get overall lineage defining mutations
	for fp in flipPairs:
		df_barcodes[fp] = df_barcodes[fp].subtract(df_barcodes[fp].min(axis=1),axis=0)

	#### drop all unused mutations (i.e. paired mutations with reversions)
	df_barcodes = df_barcodes.drop(columns=df_barcodes.columns[df_barcodes.sum(axis=0)==0])
	return df_barcodes


if __name__ == '__main__':

	fn = sys.argv[1]
	df = pd.read_csv(fn,sep='\t')
	df = parse_tree_paths(df)
	df_barcodes = convert_to_barcodes(df)
	df_barcodes = reversion_checking(df_barcodes)
	df_barcodes.to_csv('usher_barcodes.csv')
