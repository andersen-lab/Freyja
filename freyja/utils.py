import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import copy
import matplotlib.dates as mdates


def agg(results):
    allResults = [pd.read_csv(results + fn, skipinitialspace=True, sep='\t',
                              index_col=0) for fn in os.listdir(results)]
    df_demix = pd.concat(allResults, axis=1).T
    df_demix.index = [x.split('/')[-1] for x in df_demix.index]
    return df_demix


def prepLineageDict(agg_d0, thresh=0.001):
    agg_d0.loc[:, 'lineages'] = agg_d0['lineages']\
          .apply(lambda x:
                 x.replace("'", "")
                  .replace("]", "")
                  .replace("[", "")
                  .replace(")", "")
                  .replace("(", "")
                  .replace("\n", "")).copy()
    agg_d0 = agg_d0[agg_d0['lineages'].apply(lambda x: len(x) > 0)].copy()
    agg_d0.loc[:, 'lineages'] = agg_d0['lineages'].apply(lambda x:
                                                         re.sub(' +', ' ', x)
                                                           .split(' ')).copy()
    agg_d0.loc[:, 'abundances'] = agg_d0['abundances']\
          .apply(lambda x:
                 x.replace("'", "")
                  .replace("]", "")
                  .replace("[", "")
                  .replace(")", "")
                  .replace("(", "")
                  .replace("\n", "")).copy()
    agg_d0.loc[:, 'abundances'] = agg_d0['abundances']\
          .apply(lambda x:
                 re.sub(' +', ' ', x)
                   .split(' ')).copy()
    # print([float(abund) for abund in agg_d0.iloc[0].loc['abundances']])
    agg_d0.loc[:, 'linDict'] = [{lin: float(abund) for lin, abund in
                                zip(agg_d0.loc[samp, 'lineages'],
                                    agg_d0.loc[samp, 'abundances'])}
                                for samp in agg_d0.index]

    # get most common lineages in our dataset
    from collections import Counter
    counter_summary = Counter()
    for line in agg_d0.index:
        counter_summary += Counter(agg_d0.loc[line, 'linDict'])

    # possibly switch to a max number of lineages
    keepers = [x[0] for x in counter_summary.most_common()
               if x[1] / agg_d0.shape[0] > thresh]
    if len(keepers) < len(counter_summary):
        for j, sampLabel in enumerate(agg_d0.index):
            linDictMod = copy.deepcopy(agg_d0.loc[sampLabel, 'linDict'])
            linDictModCopy = copy.deepcopy(
                agg_d0.loc[sampLabel, 'linDict'])
            for rInd in linDictModCopy.keys():
                if linDictModCopy[rInd] < (thresh * agg_d0.shape[0]):
                    if 'Other' in linDictMod.keys():
                        linDictMod['Other'] += linDictMod[rInd]
                        del linDictMod[rInd]
                    else:
                        linDictMod['Other'] = linDictMod[rInd]
                        del linDictMod[rInd]
            agg_d0.loc[sampLabel, 'linDict'] = [linDictMod]
    return agg_d0


def prepSummaryDict(agg_d0):
    agg_d0.loc[:, 'summarized'] = agg_d0['summarized']\
          .apply(lambda x:
                 x.replace("'", "")
                  .replace("]", "")
                  .replace("[", "")
                  .replace(")", "")
                  .replace("(", "")
                  .replace("\n", "")
                  .split(', ')).copy()
    # drop any samples with NO lineages identified from analysis
    agg_d0 = agg_d0[agg_d0['summarized'].apply(lambda x: len(x) > 1)].copy()
    agg_d0.loc[:, 'summarized'] = agg_d0['summarized']\
          .apply(lambda x:
                 dict(zip(x[0::2],
                          x[1::2]))).copy()
    agg_d0.loc[:, 'summarized'] = agg_d0['summarized']\
          .apply(lambda x:
                 {k: float(v)
                  for k, v
                  in x.items()}).copy()
    return agg_d0


def makePlot_simple(agg_df, lineages, outputFn, colors0):
    cmap = plt.cm.tab20
    cmap_dict = {}
    if lineages:
        queryType = 'linDict'
        agg_df = prepLineageDict(agg_df)

    else:
        queryType = 'summarized'
        agg_df = prepSummaryDict(agg_df)
    ct = 0
    fig, ax = plt.subplots()
    # Make abundance fraction for all samples in aggregated dataset
    for k in range(0, agg_df.shape[0]):
        loc = pd.Series(agg_df.iloc[k][queryType])
        if len(colors0) == 0:
            for i, label in enumerate(loc.index):
                if label not in cmap_dict.keys():
                    cmap_dict[label] = cmap(ct / 20.)
                    ct += 1
                    ax.bar(k, agg_df.iloc[k][queryType][label],
                           width=0.75,
                           bottom=loc.iloc[0:i].sum(), label=label,
                           color=cmap_dict[label])
                else:
                    ax.bar(k, agg_df.iloc[k][queryType][label],
                           width=0.75,
                           bottom=loc.iloc[0:i].sum(), color=cmap_dict[label])
        else:
            for i, label in enumerate(loc.index):
                if label not in cmap_dict.keys():
                    cmap_dict[label] = cmap(ct / 20.)
                    ct += 1
                    ax.bar(k, agg_df.iloc[k][queryType][label],
                           width=0.75,
                           bottom=loc.iloc[0:i].sum(), label=label,
                           color=colors0[i])
                else:
                    ax.bar(k, agg_df.iloc[k][queryType][label],
                           width=0.75,
                           bottom=loc.iloc[0:i].sum(), color=colors0[i])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 4})
    ax.set_ylabel('Variant Prevalence')
    ax.set_xticks(range(0, agg_df.shape[0]))
    ax.set_xticklabels([sd.split('_')[0] for sd in agg_df.index],
                       rotation=90, fontsize=7)
    ax.set_ylim([0, 1])
    ax.set_xlim([-0.5, agg_df.shape[0] - 0.5])
    ax.set_aspect(6)
    fig.tight_layout()
    plt.savefig(outputFn)
    plt.close()


def makePlot_time(agg_df, lineages, times_df, interval, outputFn,
                  windowSize, colors0):
    cmap = plt.cm.tab20
    if lineages:
        queryType = 'linDict'
        agg_df = prepLineageDict(agg_df)

    else:
        queryType = 'summarized'
        agg_df = prepSummaryDict(agg_df)

    df_abundances = pd.DataFrame()
    for i, sampLabel in enumerate(agg_df.index):
        df_abundances = df_abundances.append(
            pd.Series(agg_df.loc[sampLabel, queryType],
                      name=times_df.loc[sampLabel,
                                        'sample_collection_datetime']))
    df_abundances = df_abundances.fillna(0)
    df_abundances = df_abundances.groupby(pd.Grouper(freq=interval)).mean()
    fig, ax = plt.subplots()
    if interval == 'D':
        df_abundances = df_abundances.rolling(windowSize, center=True,
                                              min_periods=0).mean()
        if len(colors0) == 0:
            ax.stackplot(df_abundances.index, df_abundances.to_numpy().T,
                         labels=df_abundances.columns, cmap=cmap)
        else:
            ax.stackplot(df_abundances.index, df_abundances.to_numpy().T,
                         labels=df_abundances.columns, colors=colors0)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 4})
        ax.set_ylabel('Variant Prevalence')
        ax.set_ylim([0, 1])
        plt.setp(ax.get_xticklabels(), rotation=90)
        fig.tight_layout()
        plt.savefig(outputFn)
        plt.close()
    elif interval == 'MS':
        for i in range(0, df_abundances.shape[1]):
            label = df_abundances.columns[i]
            # color = colorDict[label]
            if len(colors0) == 0:
                ax.bar(df_abundances.index, df_abundances.iloc[:, i],
                       width=14, bottom=df_abundances.iloc[:, 0:i].sum(axis=1),
                       label=label, color=cmap(i / 20.))
            else:
                ax.bar(df_abundances.index, df_abundances.iloc[:, i],
                       width=14, bottom=df_abundances.iloc[:, 0:i].sum(axis=1),
                       label=label, color=colors0[i])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 4})
        ax.set_ylabel('Variant Prevalence')
        locator = mdates.MonthLocator(bymonthday=1)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        ax.set_ylim([0, 1])
        ax.set_aspect(150)
        fig.tight_layout()
        plt.savefig(outputFn)
        plt.close()
    else:
        print('Error. Currently we only support Daily (D) \
               and Monthly (MS) time plots.')


if __name__ == '__main__':
    agg_results = 'data/agg_outputs.tsv'
    agg_df = pd.read_csv(agg_results, skipinitialspace=True, sep='\t',
                         index_col=0)
    lineages = False
    output = 'data/test.pdf'
    # make basic plot, without time info
    # makePlot_simple(agg_df, lineages, output, [])
    times_df = pd.read_csv('data/times_metadata.csv', index_col=0)
    times_df['sample_collection_datetime'] = \
        pd.to_datetime(times_df['sample_collection_datetime'])
    interval = 'D'
    windowSize = 14
    makePlot_time(agg_df, lineages, times_df,
                  interval, output, windowSize, [])
