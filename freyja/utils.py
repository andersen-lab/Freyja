import matplotlib.pyplot as plt
import pandas as pd
import re
import copy
import matplotlib.dates as mdates
import plotly.graph_objects as go
import plotly.express as px
import os


def agg(results):
    allResults = [pd.read_csv(fn, skipinitialspace=True, sep='\t',
                              index_col=0) for fn in results]
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
        dat = agg_df.iloc[k][queryType]
        # handle parsing differences across python versions
        if isinstance(dat, list):
            loc = pd.Series(dat[0])
        else:
            loc = pd.Series(dat)

        if len(colors0) == 0:
            for i, label in enumerate(loc.index):
                if label not in cmap_dict.keys():
                    cmap_dict[label] = cmap(ct / 20.)
                    ct += 1
                    # print('hih', loc.index, label, agg_df.iloc[k][queryType])
                    ax.bar(k, loc[label],
                           width=0.75,
                           bottom=loc.iloc[0:i].sum(), label=label,
                           color=cmap_dict[label])
                else:
                    ax.bar(k, loc[label],
                           width=0.75,
                           bottom=loc.iloc[0:i].sum(), color=cmap_dict[label])
        else:
            for i, label in enumerate(loc.index):
                if label not in cmap_dict.keys():
                    cmap_dict[label] = cmap(ct / 20.)
                    ct += 1
                    ax.bar(k, loc[label],
                           width=0.75,
                           bottom=loc.iloc[0:i].sum(), label=label,
                           color=colors0[i])
                else:
                    ax.bar(k, loc[label],
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
        dat = agg_df.loc[sampLabel, queryType]
        if isinstance(dat, list):
            df_abundances = df_abundances.append(
                pd.Series(agg_df.loc[sampLabel, queryType][0],
                          name=times_df.loc[sampLabel,
                                            'sample_collection_datetime']))
        else:
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


def make_dashboard(agg_df, meta_df, thresh, title, introText,
                   outputFn, headerColor):
    # beta version of dash output
    agg_df = prepLineageDict(agg_df)
    agg_df = prepSummaryDict(agg_df)

    # collect lineage data
    df_ab_lin = pd.DataFrame()
    for i, sampLabel in enumerate(agg_df.index):
        dat = agg_df.loc[sampLabel, 'linDict']
        if isinstance(dat, list):
            df_ab_lin = df_ab_lin.append(
                pd.Series(agg_df.loc[sampLabel, 'linDict'][0],
                          name=meta_df.loc[sampLabel,
                                           'sample_collection_datetime']))
        else:
            df_ab_lin = df_ab_lin.append(
                pd.Series(agg_df.loc[sampLabel, 'linDict'],
                          name=meta_df.loc[sampLabel,
                                           'sample_collection_datetime']))

    df_ab_lin = df_ab_lin.fillna(0)
    for col in df_ab_lin.columns:
        if col != 'Other':
            if df_ab_lin[col].sum() <= thresh:
                df_ab_lin['Other'] += df_ab_lin[col]
                df_ab_lin = df_ab_lin.drop(labels=[col], axis=1)

    cols0 = [c0 for c0 in df_ab_lin.columns
             if c0 != 'Other'] + ['Other']
    df_ab_lin = df_ab_lin[cols0]
    df_ab_lin = 100. * df_ab_lin

    # collect VOC summarized data
    df_ab_sum = pd.DataFrame()
    for i, sampLabel in enumerate(agg_df.index):
        dat = agg_df.loc[sampLabel, 'summarized']
        if isinstance(dat, list):
            df_ab_sum = df_ab_sum.append(
                pd.Series(agg_df.loc[sampLabel, 'summarized'][0],
                          name=meta_df.loc[sampLabel,
                                           'sample_collection_datetime']))
        else:
            df_ab_sum = df_ab_sum.append(
                pd.Series(agg_df.loc[sampLabel, 'summarized'],
                          name=meta_df.loc[sampLabel,
                                           'sample_collection_datetime']))

    df_ab_sum = df_ab_sum.fillna(0)
    for col in df_ab_sum.columns:
        # TODO: expand to sum values not in top N lineages, etc.
        if col != 'Other':
            if df_ab_sum[col].sum() <= thresh:
                df_ab_sum['Other'] += 1.0 * df_ab_sum[col]
                df_ab_sum = df_ab_sum.drop(labels=[col], axis=1)

    cols0 = [c0 for c0 in df_ab_sum.columns
             if c0 != 'Other'] + ['Other']
    df_ab_sum = df_ab_sum[cols0]
    df_ab_sum = 100. * df_ab_sum.fillna(0)
    # df_abundances = df_abundances.groupby(pd.Grouper(freq=interval)).mean()
    # fig, ax = plt.subplots()
    # if interval == 'D':
    #     df_abundances = df_abundances.rolling(windowSize, center=True,
    #                                           min_periods=0).mean()
    meta_df = meta_df.set_index('sample_collection_datetime')
    fig = go.Figure()
    for j, col in enumerate(df_ab_lin.columns):

        fig.add_trace(go.Scatter(
            x=df_ab_lin.index, y=df_ab_lin[col],
            hoverinfo='x+y',
            name=col,
            mode='markers+lines',
            hovertemplate="%{y:.1f}%",
            line=dict(width=0.5, color=px.colors.qualitative.Vivid[j]),
            visible=False,
            stackgroup='one'
        ))

    for j, col in enumerate(df_ab_sum.columns):

        fig.add_trace(go.Scatter(
            x=df_ab_sum.index, y=df_ab_sum[col],
            hoverinfo='x+y',
            name=col,
            mode='markers+lines',
            hovertemplate="%{y:.1f}%",
            line=dict(width=0.5, color=px.colors.qualitative.Pastel[j]),
            visible=True,
            stackgroup='one',
        ))
    fig.add_trace(go.Scatter(
        x=meta_df.index, y=meta_df['viral_load'],
        hoverinfo='x+y',
        mode='markers+lines',
        hovertemplate="%{y:.1f} copies/L",
        line=dict(width=1.25, color='blue'),
        visible=False,
    ))
    # load scaled abundances
    df_ab_lin_s = df_ab_lin.multiply(meta_df.viral_load,
                                     axis=0) / 100.
    for j, col in enumerate(df_ab_lin_s.columns):

        fig.add_trace(go.Scatter(
            x=df_ab_lin_s.index, y=df_ab_lin_s[col],
            hoverinfo='x+y',
            name=col,
            mode='markers+lines',
            hovertemplate="%{y:.1f} copies/L",
            line=dict(width=0.5, color=px.colors.qualitative.Vivid[j]),
            visible=False,
            stackgroup='one'
        ))
    df_ab_sum_s = df_ab_sum.multiply(meta_df.viral_load,
                                     axis=0) / 100.
    for j, col in enumerate(df_ab_sum_s.columns):

        fig.add_trace(go.Scatter(
            x=df_ab_sum_s.index, y=df_ab_sum_s[col],
            hoverinfo='x+y',
            name=col,
            mode='markers+lines',
            hovertemplate="%{y:.1f} copies/L",
            line=dict(width=0.5, color=px.colors.qualitative.Pastel[j]),
            visible=False,
            stackgroup='one'
        ))
    fig.update_layout(updatemenus=[dict(type="buttons",
                      direction='right',
                      active=0,
                      bgcolor='lightskyblue',
                      x=0.85,
                      y=1.07,
                      buttons=list([
                                   dict(label="VOC Summary",
                                        method="update",
                                        args=[{
                                             "visible":
                                             [False] * df_ab_lin.shape[1] +
                                             [True] * df_ab_sum.shape[1] +
                                             [False] +
                                             [False] * df_ab_lin_s.shape[1] +
                                             [False] * df_ab_sum_s.shape[1]},
                                            {"yaxis": {"title":
                                             'VOC Prevalence',
                                                       "ticksuffix": '%',
                                                       "range": [0, 100]}}]),
                                   dict(label="Lineages",
                                        method="update",
                                        args=[{
                                              "visible":
                                              [True] * df_ab_lin.shape[1] +
                                              [False] * df_ab_sum.shape[1] +
                                              [False] +
                                              [False] * df_ab_lin_s.shape[1] +
                                              [False] * df_ab_sum_s.shape[1]},
                                              {"yaxis": {"title":
                                               'Lineage Prevalence',
                                                         "ticksuffix": '%',
                                                         "range": [0, 100]}}]),
                                   dict(label="Viral Load",
                                        method="update",
                                        args=[{
                                              "visible":
                                              [False] * df_ab_lin.shape[1] +
                                              [False] * df_ab_sum.shape[1] +
                                              [True] +
                                              [False] * df_ab_lin_s.shape[1] +
                                              [False] * df_ab_sum_s.shape[1]},
                                              {"yaxis": {"title":
                                                         'Virus copies/L',
                                                         "range":
                                                         [0,
                                                          meta_df['viral_load']
                                                          .max() * 1.1]}}]),
                                   dict(
                                       label="Load Scaled Summary",
                                       method="update",
                                       args=[{
                                              "visible":
                                              [False] * df_ab_lin.shape[1] +
                                              [False] * df_ab_sum.shape[1] +
                                              [False] +
                                              [False] * df_ab_lin_s.shape[1] +
                                              [True] * df_ab_sum_s.shape[1]},
                                             {"yaxis": {"title":
                                              'Variant copies/L',
                                                        "range":
                                                        [0,
                                                         meta_df['viral_load']
                                                         .max() * 1.1]}}]),
                                   dict(
                                       label="Load Scaled Lineages",
                                       method="update",
                                       args=[{
                                              "visible":
                                              [False] * df_ab_lin.shape[1] +
                                              [False] * df_ab_sum.shape[1] +
                                              [False] +
                                              [True] * df_ab_lin_s.shape[1] +
                                              [False] * df_ab_sum_s.shape[1]},
                                             {"yaxis": {"title":
                                                        'Lineage copies/L',
                                                        "range":
                                                        [0,
                                                         meta_df['viral_load']
                                                         .max() * 1.1]}}])]))],
                      template="plotly_white",
                      hovermode="x unified",
                      xaxis=dict(hoverformat="%B %d, %Y"),
                      legend=dict(yanchor="top",
                                  y=0.99,
                                  xanchor="right",
                                  x=1.1,
                                  itemsizing='constant'))

    fig.update_layout(yaxis_title='VOC Prevalence',
                      yaxis_ticksuffix="%",
                      yaxis_range=[0, 100.],
                      width=1000,
                      height=600)
    fig.update_layout(margin=dict(b=0, t=10))

    fig.update_xaxes(dtick="6.048e+8", tickformat="%b %d", mirror=True,
                     showline=False, ticks="", showgrid=False)
    # generate plot as a div, to combine with additional html/css
    fig.write_html("div-plot.html", full_html=False, default_width='50%',
                   config={'displaylogo': False, 'displayModeBar': False})
    # now assemble into a web page
    header = "<html>\n \
              <head><meta charset='utf-8' /></head>\n\
              <body>\n\
              <style>\n\
              h1 {background: " + headerColor + "; color: white; height: 62px;\n\
                  font-family:  font-family: 'Helvetica Neue', sans-serif;\n\
                   font-size: 50px; font-weight: bold;}\n\
              h2 {background: mediumpurple; font-size: 24px; color: white;\n\
                  font-family:  font-family: 'Open Sans', sans-serif;}\n\
              p {font-family: sans-serif}\n\
              </style>\n\
              <div class='header' align='center'> \n\
              <h1>" + title + "</h1> \n\
              <p>" + introText + "</p> \n\
              </div>\n"
    bottom = "</body> \n\
                <p align='right'>\
                <img src='https://tinyurl.com/freyjaLogo' \
                alt='Freyja logo' width=100 height=100> </p>\n\
                <p align='right'> Made with\
                <a href= https://github.com/andersen-lab/Freyja> Freyja </a>\
                </p>\n\
                <hr>\n\
                </html>"

    centereddiv = "<div align='center'><script type=\
                  'text/javascript'>window.PlotlyConfig =\
                  {MathJaxConfig: 'local'};</script>"

    filenames = ['div-plot.html']
    with open(outputFn, 'w') as outfile:
        outfile.write(header)
        for fname in filenames:
            with open(fname) as infile:
                for j, line in enumerate(infile):
                    if j == 0:
                        outfile.write(centereddiv)
                    else:
                        outfile.write(line)
        outfile.write(bottom)
    os.remove('div-plot.html')


if __name__ == '__main__':
    agg_results = 'data/test_sweep.tsv'
    agg_df = pd.read_csv(agg_results, skipinitialspace=True, sep='\t',
                         index_col=0)
    lineages = True
    # output = 'data/test.pdf'
    # make basic plot, without time info
    # makePlot_simple(agg_df, lineages, output, [])
    # times_df = pd.read_csv('data/times_metadata.csv', index_col=0)
    # times_df['sample_collection_datetime'] = \
    #     pd.to_datetime(times_df['sample_collection_datetime'])
    interval = 'D'
    windowSize = 14
    # # makePlot_time(agg_df, lineages, times_df,
    # #               interval, output, windowSize, [])
    meta_df = pd.read_csv('data/sweep_metadata.csv', index_col=0)
    meta_df['sample_collection_datetime'] = \
        pd.to_datetime(meta_df['sample_collection_datetime'])
    thresh = 0.01
    # read in inputs
    with open("data/title.txt", "r") as f:
        title = ''.join(f.readlines())
    with open("data/introContent.txt", "r") as f:
        introText = ''.join(f.readlines())
    outputFn = 'tester0.html'
    headerColor = 'mediumpurple'
    make_dashboard(agg_df, meta_df, thresh, title,
                   introText, outputFn, headerColor)
