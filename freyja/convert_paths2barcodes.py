from copy import deepcopy
import numpy as np
import pandas as pd 
import sys

fn = sys.argv[1]
df = pd.read_csv(fn,sep='\t')
df.set_index('clade',inplace=True)
df = df.drop_duplicates(keep='last')### check which lineages got trimmed. 
df['from_tree_root'] = df['from_tree_root'].fillna('')
df['from_tree_root'] = df['from_tree_root'].apply(lambda x:x.replace(' ','').strip('>').split('>'))


df_barcodes = pd.DataFrame()

for clade in df.index:
	cladeSeries = pd.Series({c:1 for c in df.loc[clade,'from_tree_root']},name=clade)
	df_barcodes = df_barcodes.append(cladeSeries)

print('separating combined splits')
df_barcodes.drop(columns='',inplace=True)
temp = pd.DataFrame()
dropList = []
for c in df_barcodes.columns:
	### if column includes multiple mutations, split into separate columns and concatenate
	if "," in c:
		for mt in c.split(","):
			temp[mt] = df_barcodes[c]
		dropList.append(c)
df_barcodes.drop(columns=dropList,inplace=True)
df_barcodes = pd.concat((df_barcodes,temp),axis=1)
df_barcodes = df_barcodes.fillna(0)

df_barcodes = df_barcodes.groupby(axis=1,level=0).sum()

print('checking for mutation pairs')
def sortFun(x):
	return int(x[1:(len(x)-1)])
###check if a reversion is present. 
flips = [d[-1] +d[1:len(d)-1]+d[0] for d in df_barcodes.columns[1:]]
flipPairs = [d for d in df_barcodes.columns[1:] if d in flips]
flipPairs.sort(key=sortFun)
flipPairs = [[flipPairs[j],flipPairs[j+1]] for j in np.arange(0,len(flipPairs),2)]

### subtract lower of two pair counts from both to get overall haplotype
for fp in flipPairs:
	df_barcodes[fp] = df_barcodes[fp].subtract(df_barcodes[fp].min(axis=1),axis=0)

#### drop all unused mutations (i.e. paired mutations with reversions)
df_barcodes.drop(columns=df_barcodes.columns[df_barcodes.sum(axis=0)==0])
df_barcodes.to_csv('usher_barcodes.csv')



# #############################
# #render simple visualization
# import seaborn as sns
# import matplotlib.pyplot as plt
# fig,ax = plt.subplots()
# ### minimize for figure 

# df_barcodesMini = df_barcodes.loc[df_barcodes.sum(axis=1).sort_values(ascending=False).index[np.arange(0,600,4)],:]
# df_barcodesMini = df_barcodesMini.loc[:,df_barcodes.sum(axis=0).sort_values(ascending=False).index[0:150]].drop_duplicates()
# cg = sns.clustermap(df_barcodesMini.T.drop_duplicates(),yticklabels=False,xticklabels=False,cmap='Greens',linewidths=0.25, linecolor='grey')
# cg.ax_row_dendrogram.set_visible(False)
# cg.ax_col_dendrogram.set_visible(False)
# cg.cax.set_visible(False)
# ax =cg.ax_heatmap
# ax.set_xlabel("Haplotype")
# ax.set_ylabel("Mutation")
# fig.tight_layout()
# plt.savefig('barcode_mat.pdf')
# plt.close()