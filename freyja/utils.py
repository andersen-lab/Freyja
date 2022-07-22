import matplotlib.pyplot as plt
import pandas as pd
import re
import copy
import matplotlib.dates as mdates
import plotly.graph_objects as go
import plotly.express as px
import os
import yaml


def agg(results):
    allResults = [pd.read_csv(fn, skipinitialspace=True, sep='\t',
                              index_col=0) for fn in results]
    df_demix = pd.concat(allResults, axis=1).T
    df_demix.index = [x.split('/')[-1] for x in df_demix.index]
    return df_demix


# Check the format of the config file and correct if necessary.
def checkConfig(config):
    """
    Sanitize the config file
    """
    CONFIG_KEYS = {
        'VOC': ['name', 'color'],
        'Lineages': ['name', 'members', 'color']
    }
    for key, val in CONFIG_KEYS.items():
        if key not in config.keys():
            config[key] = {}
        if config[key] != {}:
            for k in val:
                for nest_key in config[key].keys():
                    if k not in config[key][nest_key].keys():
                        raise ValueError(f'{key} key missing {k} key' +
                                         ' in the config file')
    return config


# Get value from the config dictionary.
def get_value(val, dict, get_val, match_key):
    values = []
    for key, value in dict.items():
        if val == value[match_key]:
            values.append(dict[key][get_val])
    return values[0]


# Get name key from the config dictionary for which the
# queried lineage is a member.
def get_name(val, dict):
    values = []
    for key, value in dict.items():
        if val in value['members']:
            values.append(dict[key]['name'])
    return values[0]


# Used to define the color scheme for make_dashboard function.
def get_color_scheme(df, default_color_scheme, config=None):
    if config is None:
        config = {}
    default_colors = None
    color_scheme = {}
    for i in default_color_scheme.keys():
        if df.shape[1] <= int(i):
            default_colors = default_color_scheme[i]
            break

    if default_colors is None and len(config.keys()) < df.shape[1]:
        raise Exception('No default color scheme found for this number of' +
                        ' variants.')

    for i, col in enumerate(df.columns):
        if any(col == i['name'] for i in config.values()) and config != {}:
            if get_value(col, config, 'color', 'name').lower() != "default":
                color_scheme[col] = get_value(col, config, 'color', 'name')
            else:
                color_scheme[col] = default_colors[i]
        else:
            color_scheme[col] = default_colors[i]
    return color_scheme


def prepLineageDict(agg_d0, thresh=0.001, config=None, lineage_info=None):
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

    # Aggregate the lineages specified in config file into a single
    # key/value pair.
    if config is not None:
        # Returns the first name where the queried lineage is present in the
        # config memebers.
        for i, i_config in config.items():
            for lin in i_config['members']:
                if '*' in lin:
                    config[i]['members'].extend(lineage_info
                                                [lin.replace('*', '')]
                                                ['children']
                                                )
        # aggregate all the members into 1ist
        config_members = [lin for val in config.values() for lin in val[
                          'members']]
        # Aggregate lineages into a single dict key/value pair
        for j, sampLabel in enumerate(agg_d0.index):
            # This list will contain the keys that have been created
            # in the dict. This is used to ensure that the same key is
            # not created twice or that value of the key is doubled.
            created_keys = []
            linDictMod = copy.deepcopy(agg_d0.loc[sampLabel, 'linDict'])
            linDictModCopy = copy.deepcopy(
                agg_d0.loc[sampLabel, 'linDict'])
            for rInd in linDictModCopy.keys():
                if rInd in config_members:
                    if get_name(rInd, config) in created_keys:
                        linDictMod[
                            get_name(rInd, config)
                            ] += linDictMod.pop(rInd)
                    else:
                        linDictMod[
                            get_name(rInd, config)
                            ] = linDictMod.pop(rInd)
                        created_keys.append(get_name(rInd, config))
            agg_d0.loc[sampLabel, 'linDict'] = [linDictMod]

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
                   outputFn, headerColor, scale_by_viral_load,
                   config, lineage_info):
    # beta version of dash output
    agg_df = prepLineageDict(agg_df, config=config.get('Lineages'),
                             lineage_info=lineage_info)
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

    if scale_by_viral_load:
        df_ab_lin_rescale = df_ab_lin.multiply(meta_df.viral_load,
                                               axis=0)
        df_ab_lin_rescale = df_ab_lin_rescale.groupby(level=0).sum()
        df_ab_lin_rescale = df_ab_lin_rescale.sort_index()
        meta_df = meta_df.groupby('sample_collection_datetime').sum()\
            .sort_index()
        df_ab_lin = df_ab_lin_rescale.divide(meta_df.viral_load,
                                             axis=0)
        df_ab_lin = df_ab_lin.fillna(0).sort_index()
    else:
        meta_df = meta_df.groupby('sample_collection_datetime').mean()\
            .sort_index()
        df_ab_lin = df_ab_lin.groupby(level=0).mean().sort_index()

    dates_to_keep = meta_df.index[~meta_df['viral_load'].isna()]
    df_ab_sum = df_ab_sum.groupby(level=0).mean()
    df_ab_sum = df_ab_sum.sort_index()
    fig = go.Figure()

    default_color_lin = {
        11: px.colors.qualitative.Vivid,
        24: px.colors.qualitative.Dark24
    }
    color_lin = get_color_scheme(df_ab_lin,
                                 default_color_lin,
                                 config.get('Lineages'))
    for j, col in enumerate(df_ab_lin.columns):

        fig.add_trace(go.Scatter(
            x=df_ab_lin.index, y=df_ab_lin[col],
            hoverinfo='x+y',
            name=col,
            mode='markers+lines',
            hovertemplate="%{y:.1f}%",
            line=dict(width=0.5, color=color_lin[col]),
            visible=False,
            stackgroup='one'
        ))

    default_color_sum = {
        11: px.colors.qualitative.Pastel,
        24: px.colors.qualitative.Light24
    }
    color_sum = get_color_scheme(df_ab_sum,
                                 default_color_sum,
                                 config.get('VOC'))
    for j, col in enumerate(df_ab_sum.columns):
        fig.add_trace(go.Scatter(
            x=df_ab_sum.index, y=df_ab_sum[col],
            hoverinfo='x+y',
            name=col,
            mode='markers+lines',
            hovertemplate="%{y:.1f}%",
            line=dict(width=0.5, color=color_sum[col]),
            visible=True,
            stackgroup='one',
        ))
    # if needed, drop dates with missing viral load metadata
    if len(dates_to_keep) < meta_df.shape[0]:
        meta_df = meta_df.loc[dates_to_keep]
        df_ab_sum = df_ab_sum.loc[dates_to_keep]
        df_ab_lin = df_ab_lin.loc[dates_to_keep]
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

    default_color_lin_s = {
        11: px.colors.qualitative.Vivid,
        24: px.colors.qualitative.Dark24
    }
    color_lin_s = get_color_scheme(df_ab_lin_s,
                                   default_color_lin_s,
                                   config.get('Lineages'))
    for j, col in enumerate(df_ab_lin_s.columns):

        fig.add_trace(go.Scatter(
            x=df_ab_lin_s.index, y=df_ab_lin_s[col],
            hoverinfo='x+y',
            name=col,
            mode='markers+lines',
            hovertemplate="%{y:.1f} copies/L",
            line=dict(width=0.5, color=color_lin_s[col]),
            visible=False,
            stackgroup='one'
        ))

    df_ab_sum_s = df_ab_sum.multiply(meta_df.viral_load,
                                     axis=0) / 100.

    default_color_sum_s = {
        11: px.colors.qualitative.Pastel,
        24: px.colors.qualitative.Light24
    }
    color_sum_s = get_color_scheme(df_ab_sum,
                                   default_color_sum_s,
                                   config.get('VOC'))
    for j, col in enumerate(df_ab_sum_s.columns):
        fig.add_trace(go.Scatter(
            x=df_ab_sum_s.index, y=df_ab_sum_s[col],
            hoverinfo='x+y',
            name=col,
            mode='markers+lines',
            hovertemplate="%{y:.1f} copies/L",
            line=dict(width=0.5, color=color_sum_s[col]),
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
    print("Plot saved to " + outputFn)


if __name__ == '__main__':
    agg_results = 'freyja/data/test_sweep.tsv'
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
    meta_df = pd.read_csv('freyja/data/sweep_metadata.csv', index_col=0)
    meta_df['sample_collection_datetime'] = \
        pd.to_datetime(meta_df['sample_collection_datetime'])
    thresh = 0.01
    # read in inputs
    with open("freyja/data/title.txt", "r") as f:
        title = ''.join(f.readlines())
    with open("freyja/data/introContent.txt", "r") as f:
        introText = ''.join(f.readlines())
    outputFn = 'tester0.html'
    headerColor = 'mediumpurple'
    with open('freyja/data/plot_config.yml', "r") as f:
        try:
            config = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            raise ValueError('Error in config file: ' + str(exc))
    import requests
    r = requests.get('https://raw.githubusercontent.com/cov-lineages' +
                     '/lineages-website/master/data/lineages.yml')
    if r.status_code == 200:
        with open('freyja/data/lineages.yml', 'w+') as f:
            f.write(r.text)
    elif os.path.exists('freyja/data/lineages.yml'):
        print('Using existing lineages.yml')
    else:
        raise FileNotFoundError('Could not download or find lineages.yml')

    # read linages.yml file
    with open('freyja/data/lineages.yml', 'r') as f:
        try:
            lineages_yml = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            raise ValueError('Error in lineages.yml file: ' + str(exc))
    scale_by_viral_load = False
    # converts lineages_yml to a dictionary where the lineage names are the
    # keys.
    lineage_info = {}
    for lineage in lineages_yml:
        lineage_info[lineage['name']] = {'children': lineage['children']}
    config = checkConfig(config)
    make_dashboard(agg_df, meta_df, thresh, title,
                   introText, outputFn, headerColor, scale_by_viral_load,
                   config, lineage_info)
