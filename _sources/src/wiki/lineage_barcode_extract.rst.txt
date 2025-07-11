Lineage Barcode Extract
-------------------------------------------------------------------------------

If you’re looking to identify lineage defining mutations for your
lineage of choice, look no further! You can extract this information
using Freyja’s existing barcoding and a bit of python.

1. Run a freyja update and extract barcodes to the directory of your
   choosing (will take a few minutes):

::

   freyja update --outdir /my/local/directory/

2. Load the barcodes using pandas and extract barcodes for your lineages
   of interest

.. code:: python

   import pandas as pd

   def sortFun(x):
       # sort based on nuc position, ignoring nuc identities
       return int(x[1:(len(x)-1)])
   # replace "/my/local/directory/" with your chosen directory
   df = pd.read_csv('/my/local/directory/usher_barcodes.csv', index_col=0)

   # specify your lineages of interest here
   lins = ['B.1.617.2','BA.1','BA.2','B.1.1.7','P.1']
   df = df.loc[lins]
   # keep only columns with at least one mutation across the lins
   keepcols = df.columns[df.sum(axis=0)>0]
   df = df.loc[:, keepcols]

Optional: 3. Load in gene coordinate information to translate nucleotide
mutations into possible AA mutations. An example file for doing this is
available
`here <https://github.com/andersen-lab/Freyja/wiki/SARS-CoV-2.json>`__.
Then append AA mutations to the corresponding nucleotide mutations.

.. code:: python

   # build dataframe with gene data
   import json
   f1 = open('/my/local/directory/SARS-CoV-2.json',)
   dat = json.load(f1)
   f1.close()
   df_genes = pd.DataFrame(columns=["gene","start","end"], dtype=object)

   for d in dat["genes"]:
       new_row = {'gene':d,'start':dat["genes"][d]['coordinates']['from'],
                  'end':dat["genes"][d]['coordinates']['to']}
       df_genes = df_genes.append(new_row,ignore_index=True)
   df_genes = df_genes.sort_values('start').set_index('gene',drop=False)

   # function for converting to gene-wise coordinate numbering 
   def getGeneNum(pos,df0):
       j=0
       while j<df0.shape[0]:
           if pos<=df0.iloc[j].loc['end'] and pos>=df0.iloc[j].loc['start']:
               return df0.iloc[j].loc['gene'], (pos - df0.iloc[j].loc['start'])//3+1
           else:
               j+=1
       return '',''

   # add AA mutation info into our dataframe
   cols = list(df.columns)
   cols.sort(key=lambda x: sortFun(x))
   df = df[cols]
   for i, mut in enumerate(cols):
       pos = sortFun(mut)
       g,pa = getGeneNum(pos,df_genes)
       if len(g)==0:
           refAA=''
           mutAA=''
       else:
           posInTriple = (pos - df_genes.loc[g,'start']) % 3
           triple = ref0[pos-posInTriple-1:pos-posInTriple+2]
           refAA = triple.translate()
           tripleMut = triple
           tripleMut[posInTriple] = mut[-1]
           mutAA = tripleMut.translate()
           cols[i] = cols[i] +'('+ g+':'+str(refAA)+str(pa)+str(mutAA)+')'

   df.columns = cols

4. Export to csv format for later viewing.

.. code:: python

   df.to_csv('mutations_lineage_subset.csv')
