import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as clr

# Define columns (mutations in samples)
colnames = []
for infile in os.listdir('data/plotting'):
    coocs = pd.read_csv(os.path.join('data/plotting', infile), sep='\t',
                        header=0)
    grouped_muts = []
    for c in coocs.iloc[:, 0]:
        for mut in c.split('),'):
            if not mut[-1] == ')':
                mut += ')'
            if mut not in colnames:
                colnames.append(mut)
colnames = sorted(colnames, key=lambda x: int(x[1:6]))

data = {}

for infile in os.listdir('data/plotting'):
    date = infile.split('.')[0]
    coocs = pd.read_csv(os.path.join('data/plotting', infile), sep='\t',
                        header=0)
    grouped_muts = []
    for c in coocs.iloc[:, 0]:
        if '),' in c:
            for mut in c.split('),'):
                if not mut[-1] == ')':
                    mut += ')'
                if mut not in grouped_muts:
                    grouped_muts.append(mut)

    # coocs = coocs.loc[~coocs['COOCCURRENCES'].isin(
    #     grouped_muts)].to_dict('split')['data']

    coocs = coocs.to_dict('split')['data']
    counts = {c: 0 for c in colnames}
    for mut in colnames:
        for pair in coocs:
            if mut in pair[0]:
                counts[mut] = pair[1]
                break
    data[date] = [counts[c] for c in colnames]


df = pd.DataFrame.from_dict(data, orient='index')
df.index = pd.to_datetime(df.index).date
df.columns = colnames

df = df.sort_index()
fig, ax = plt.subplots(figsize=(15, 10))
cmap = clr.LinearSegmentedColormap.from_list(
    'rdgray', ['#D3D3D3', '#FF6347'], N=256)

print(df)
plot = sns.heatmap(df.T, ax=ax, cbar=False, square=True,
                   fmt='', linewidths=0.5, cmap=cmap, vmin=0, vmax=500)

plt.savefig('cooccurrence_plot.pdf')