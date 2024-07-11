import copy
import os
import re
import tqdm

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, timedelta
import yaml
from scipy.optimize import curve_fit

# parameters to make plots illustrator friendly
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def agg(results):
    allResults = [pd.read_csv(fn, skipinitialspace=True, sep='\t',
                              index_col=0) for fn in results]
    df_demix = pd.concat(allResults, axis=1).T
    df_demix.index = [x.split('/')[-1] for x in df_demix.index]
    return df_demix


def load_barcodes(barcodes):
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    if barcodes != '':
        if barcodes.endswith('csv'):
            df_barcodes = pd.read_csv(barcodes, index_col=0)
        elif barcodes.endswith('feather'):
            df_barcodes = pd.read_feather(barcodes).set_index('index')
        else:
            raise ValueError('Only csv and feather barcode ' +
                             'formats supported')
    else:
        df_barcodes = pd.read_feather(os.path.join(locDir,
                                      'data/usher_barcodes.feather')
                                      ).set_index('index')
    return df_barcodes


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


def logistic_growth(ndays, b, r):
    return 1 / (1 + (b * np.exp(-1 * r * ndays)))


# Calcualate the relative growth rates of the lineages and return a dataFrame.
def calc_rel_growth_rates(df, nboots=1000, serial_interval=5.5,
                          outputFn='rel_growth_rates.csv', daysIncluded=56,
                          thresh=0.001):
    df.index.name = 'Date'
    df.reset_index(inplace=True)
    df['Date'] = pd.to_datetime(df['Date'])
    df.columns = [dfc.split(' (')[0] for dfc in df.columns]
    df = df.set_index('Date')

    df = df.dropna(axis=0, how='all')
    df = df.fillna(0)
    df = df / 100.

    # go as far back as we can, within daysIncluded limit
    nBack = next((x[0] + 1 for x in enumerate(df.index[::-1])
                 if (df.index[-1] - x[1]).days > daysIncluded), 0)
    rel_growth_rate = {
        'Lineage': [],
        'Estimated Advantage': [],
        'Bootstrap 95% interval': [],
    }
    # get all lineages present at >0.1% average over last 8 weeks
    lineages = df.columns[df.iloc[-nBack:].mean(axis=0) > thresh]
    rate_cal = tqdm.tqdm(enumerate(lineages), total=len(lineages),
                         desc='Rate calculations for lineages/groups')
    for k, lineage in rate_cal:
        # print(f"\nCalculating relative rate for {lineage}")
        rate_cal.set_postfix({"Calculating relative rate for": lineage})
        days = np.array([(dfi - df.index[-nBack]).days
                         for j, dfi in enumerate(df.index[-nBack:])])
        data = df[lineage][-len(days):]
        fit, covar = curve_fit(
            f=logistic_growth,
            xdata=days,
            ydata=data,
            p0=[100, 0.1],
            bounds=([0, -10], [1000, 10])
        )
        # build 95% CI by bootstrapping.
        bootSims = []
        coef_ests = []
        for j in range(nboots):
            bootInds = np.random.randint(0, len(days), len(days))
            daysBoot = days[bootInds]
            dataBoot = data.iloc[bootInds]
            try:
                fitBoot, covarBoot = curve_fit(
                    f=logistic_growth,
                    xdata=daysBoot,
                    ydata=dataBoot,
                    p0=[100, 0.1],
                    bounds=([0, -10], [1000, 10])
                )
            except RuntimeError:
                print('WARNING: Optimal parameters not found: The maximum' +
                      ' number of function evaluations is exceeded.')

            bootSims.append([logistic_growth(i, *fitBoot) for i in days])
            coef_ests.append(fitBoot[1])

        bootSims = np.array(bootSims)

        coef_lower = np.percentile(coef_ests, 2.5)
        coef_upper = np.percentile(coef_ests, 97.5)

        rate0 = fit[1]

        trans_increase = serial_interval * rate0
        rel_growth_rate['Lineage'].append(lineage)
        rel_growth_rate['Estimated Advantage'].append(f'{trans_increase:.1%}')
        rel_growth_rate['Bootstrap 95% interval'].append(
            f'[{serial_interval*coef_lower:0.2%} ' +
            f', {serial_interval*coef_upper:0.2%}]')
    if outputFn.endswith('.html'):
        outputFn = outputFn.replace('.html', '_rel_growth_rates.csv')
    pd.DataFrame.from_dict(
        rel_growth_rate,
        orient='columns').sort_values(
            by='Estimated Advantage',
            ascending=False).to_csv(outputFn, index=False)
    print("CSV file saved to " + outputFn)


# Get value from the config dictionary.
def get_value(val, dict, get_val, match_key):
    values = []
    for key, value in dict.items():
        if val == value[match_key]:
            values.append(dict[key][get_val])
    return values[0]


def read_lineage_file(lineageyml, locDir, fileOnly=False):
    if lineageyml == "":
        with open(os.path.join(locDir, 'data/lineages.yml'), 'r') as f:
            try:
                lineages_yml = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                raise ValueError('Error in lineages.yml file: ' + str(exc))
    else:
        with open(lineageyml, 'r') as f:
            try:
                lineages_yml = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                raise ValueError('Error in lineages.yml file: ' + str(exc))
    if fileOnly:
        return lineages_yml
    else:
        lineage_info = {}
        for lineage in lineages_yml:
            lineage_info[lineage['name']] = {'name': lineage['name'],
                                             'children': lineage['children']}
        return lineage_info


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
        raise Exception('Too many lineages to show. Use --config to group.')

    for i, col in enumerate(df.columns):
        if any(col == i['name'] for i in config.values()) and config != {}:
            if get_value(col, config, 'color', 'name').lower() != "default":
                color_scheme[col] = get_value(col, config, 'color', 'name')
            else:
                color_scheme[col] = default_colors[i]
        else:
            color_scheme[col] = default_colors[i]
    # print(color_scheme)
    return color_scheme


def prepLineageDict(agg_d0, thresh=0.001, config=None, lineage_info=None,
                    mergeLikes=False):

    if len(agg_d0.index[agg_d0.index.duplicated(keep=False)]) > 0:
        print('WARNING: multiple samples have the same ID/filename.')
        print(agg_d0.index[agg_d0.index.duplicated(keep=False)])
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

    if mergeLikes:
        agg_d0.loc[:, 'lineages'] = agg_d0['lineages'].apply(
            lambda x: [x0.split('-like')[0] for x0 in x])
        newLins, newAbunds = [], []
        for lins, abunds in zip(agg_d0['lineages'], agg_d0['abundances']):
            linsUnique, indicesUnique = np.unique(lins, return_inverse=True)
            newLins.append(linsUnique)
            newAbund = [0]*len(linsUnique)
            for ind0, j in enumerate(indicesUnique):
                newAbund[j] += float(abunds[ind0])
            newAbunds.append(newAbund)
        agg_d0['lineages'] = newLins
        agg_d0['abundances'] = newAbunds
        # agg_d0.loc[:,'abundances'] = agg_d0[['lineages']
        # merge abundances
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
        processed_linDictMod = []
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
            processed_linDictMod.append(linDictMod)
        agg_d0.loc[:, 'linDict'] = processed_linDictMod

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


def makePlot_simple(agg_df, lineages, outputFn, config, lineage_info,
                    thresh, writeGrouped):
    if lineages:
        queryType = 'linDict'
        config = config.get('Lineages')
        agg_df = prepLineageDict(agg_df, config=config,
                                 lineage_info=lineage_info, thresh=thresh)
        if writeGrouped != '':
            agg_df.to_csv(writeGrouped, sep='\t')

    else:
        queryType = 'summarized'
        config = config.get('VOC')
        agg_df = prepSummaryDict(agg_df)
    fig, ax = plt.subplots()
    # Make abundance fraction for all samples in aggregated dataset
    labelList = []

    df_abundances = pd.DataFrame()
    for i, sampLabel in enumerate(agg_df.index):
        dat = agg_df.loc[sampLabel, queryType]
        if isinstance(dat, list):
            df_abundances = pd.concat([
                df_abundances,
                pd.Series(
                    agg_df.loc[sampLabel, queryType][0],
                    name=sampLabel)
            ], axis=1)
        else:
            df_abundances = pd.concat([
                df_abundances,
                pd.Series(
                    agg_df.loc[sampLabel, queryType],
                    name=sampLabel)
            ], axis=1)

    df_abundances = df_abundances.T
    if 'Other' in df_abundances.columns:
        cols0 = [c0 for c0 in df_abundances.columns
                 if c0 != 'Other'] + ['Other']
        df_abundances = df_abundances[cols0]
    default_cmap_dict = {
        24: px.colors.qualitative.Dark24
    }
    cmap_dict = get_color_scheme(df_abundances,
                                 default_cmap_dict,
                                 config)

    for k in range(0, agg_df.shape[0]):
        loc = df_abundances.iloc[k]

        for i, label in enumerate(loc.index):
            if label in labelList:
                ax.bar(k, loc[label],
                       width=0.75,
                       bottom=loc.iloc[0:i].sum(),
                       color=cmap_dict[label])
            else:
                ax.bar(k, loc[label],
                       width=0.75,
                       bottom=loc.iloc[0:i].sum(), label=label,
                       color=cmap_dict[label])
                labelList.append(label)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='center left',
              bbox_to_anchor=(1, 0.5))
    ax.set_ylabel('Variant Prevalence')
    ax.set_xticks(range(0, agg_df.shape[0]))
    ax.set_xticklabels(agg_df.index,
                       rotation=90, fontsize=7)
    ax.set_ylim([0, 1])
    ax.set_xlim([-0.5, agg_df.shape[0] - 0.5])
    ax.set_aspect(6)
    fig.tight_layout()
    plt.savefig(outputFn)
    plt.close()


def makePlot_time(agg_df, lineages, times_df, interval, outputFn,
                  windowSize, config, lineage_info, thresh,
                  writeGrouped):
    if lineages:
        queryType = 'linDict'
        config = config.get('Lineages')
        agg_df = prepLineageDict(agg_df, config=config,
                                 lineage_info=lineage_info, thresh=thresh)
        if writeGrouped != '':
            agg_df.to_csv(writeGrouped, sep='\t')

    else:
        queryType = 'summarized'
        config = config.get('VOC')
        agg_df = prepSummaryDict(agg_df)

    df_abundances = pd.DataFrame()
    for i, sampLabel in enumerate(agg_df.index):
        dat = agg_df.loc[sampLabel, queryType]
        if isinstance(dat, list):
            df_abundances = pd.concat([
                df_abundances,
                pd.Series(
                    agg_df.loc[sampLabel, queryType][0],
                    name=times_df.loc[sampLabel,
                                      'sample_collection_datetime'])
            ], axis=1)
        else:
            df_abundances = pd.concat([
                df_abundances,
                pd.Series(
                    agg_df.loc[sampLabel, queryType],
                    name=times_df.loc[sampLabel,
                                      'sample_collection_datetime'])
            ], axis=1)
    df_abundances = df_abundances.T
    df_abundances = df_abundances.fillna(0)
    if interval == 'W':
        # epiweek ends on sat, starts on sun
        interval = 'W-SAT'
    df_abundances = df_abundances.groupby(pd.Grouper(freq=interval)).mean()
    if 'Other' in df_abundances.columns:
        cols0 = [c0 for c0 in df_abundances.columns
                 if c0 != 'Other'] + ['Other']
        df_abundances = df_abundances[cols0]
    default_cmap_dict = {
        24: px.colors.qualitative.Dark24
    }

    cmap_dict = get_color_scheme(df_abundances,
                                 default_cmap_dict,
                                 config)
    fig, ax = plt.subplots()
    if interval == 'D':
        df_abundances = df_abundances.rolling(windowSize, center=True,
                                              min_periods=0).mean()
        ax.stackplot(df_abundances.index, df_abundances.to_numpy().T,
                     labels=df_abundances.columns, colors=cmap_dict.values())
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='center left',
                  bbox_to_anchor=(1, 0.5))
        ax.set_ylabel('Variant Prevalence')
        ax.set_ylim([0, 1])
        ax.set_xlim([df_abundances.index.min(), df_abundances.index.max()])
        plt.setp(ax.get_xticklabels(), rotation=90)
        fig.tight_layout()
        plt.savefig(outputFn)
        plt.close()
    elif interval == 'MS':
        for i in range(0, df_abundances.shape[1]):
            label = df_abundances.columns[i]
            ax.bar(df_abundances.index, df_abundances.iloc[:, i],
                   width=14, bottom=df_abundances.iloc[:, 0:i].sum(axis=1),
                   label=label, color=cmap_dict[label])
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='center left',
                  bbox_to_anchor=(1, 0.5))
        ax.set_ylabel('Variant Prevalence')
        locator = mdates.MonthLocator(bymonthday=1)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        ax.set_ylim([0, 1])
        ax.set_aspect(150)
        ax.set_xlim([df_abundances.index.min()-timedelta(days=15),
                     df_abundances.index.max()+timedelta(days=15)])
        fig.tight_layout()
        plt.savefig(outputFn)
        plt.close()
    elif interval == 'W-SAT':
        from epiweeks import Week
        weekInfo = [Week.fromdate(dfi).weektuple()
                    for dfi in df_abundances.index]
        df_abundances.index = [str(wi[0])+'-'+str(wi[1]) for wi in weekInfo]
        for i in range(0, df_abundances.shape[1]):
            label = df_abundances.columns[i]
            ax.bar(df_abundances.index, df_abundances.iloc[:, i],
                   width=0.5,
                   bottom=df_abundances.iloc[:, 0:i].sum(axis=1),
                   label=label, color=cmap_dict[label])
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='center left',
                  bbox_to_anchor=(1, 0.5))
        labelsAx = [item.split('-')[1] for item in df_abundances.index]
        ax.set_xticks(range(0, len(labelsAx)))
        ax.set_xticklabels(labelsAx)
        ax.set_ylabel('Variant Prevalence')
        ax.set_xlabel('Epiweek')
        ax.set_ylim([0, 1])
        fig.tight_layout()
        plt.savefig(outputFn)
        plt.close()
    else:
        print('Error. Currently we only support Daily (D), Weekly (W), \
               and Monthly (MS) time plots.')


def get_abundance(agg_df, meta_df, thresh, scale_by_viral_load, config,
                  lineage_info):
    agg_df = prepLineageDict(agg_df, config=config.get('Lineages'),
                             lineage_info=lineage_info, thresh=thresh)
    agg_df = prepSummaryDict(agg_df)

    if len(meta_df.index[meta_df.index.duplicated(keep=False)]) > 0:
        raise ValueError(
            'ERROR: multiple entries for same sample in metadata:'
            f'\n{meta_df.index[meta_df.index.duplicated(keep=False)]}')

    agg_df.to_csv('agg_df.csv')
    # collect lineage data
    df_ab_lin = pd.DataFrame()
    for i, sampLabel in enumerate(agg_df.index):
        dat = agg_df.loc[sampLabel, 'linDict']
        if isinstance(dat, list):
            if i == 0:
                df_ab_lin = pd.Series(
                    agg_df.loc[sampLabel, 'linDict'][0],
                    name=meta_df.loc[sampLabel,
                                     'sample_collection_datetime'])
            else:
                df_ab_lin = pd.concat([
                    df_ab_lin,
                    pd.Series(agg_df.loc[sampLabel, 'linDict'][0],
                              name=meta_df.loc[sampLabel,
                                               'sample_collection_datetime'])
                ], axis=1)
        else:
            if i == 0:
                df_ab_lin = pd.Series(
                    agg_df.loc[sampLabel, 'linDict'],
                    name=meta_df.loc[sampLabel,
                                     'sample_collection_datetime'])
            else:
                df_ab_lin = pd.concat([
                    df_ab_lin,
                    pd.Series(agg_df.loc[sampLabel, 'linDict'],
                              name=meta_df.loc[sampLabel,
                                               'sample_collection_datetime'])
                ], axis=1)
    df_ab_lin = df_ab_lin.T
    df_ab_lin.to_csv('df_ab_lin_raw.csv')
    df_ab_lin = df_ab_lin.fillna(0)
    if 'Other' not in df_ab_lin.columns:
        df_ab_lin['Other'] = 0.
    # for col in df_ab_lin.columns:
    #     if col != 'Other':
    #         if df_ab_lin[col].mean() <= thresh:
    #             df_ab_lin['Other'] += df_ab_lin[col]
    #             df_ab_lin = df_ab_lin.drop(labels=[col], axis=1)

    cols0 = [c0 for c0 in df_ab_lin.columns
             if c0 != 'Other'] + ['Other']
    df_ab_lin = df_ab_lin[cols0]
    df_ab_lin = 100. * df_ab_lin

    # collect VOC summarized data
    # df_ab_sum = pd.DataFrame()
    for i, sampLabel in enumerate(agg_df.index):
        dat = agg_df.loc[sampLabel, 'summarized']
        if isinstance(dat, list):
            if i == 0:
                df_ab_sum = pd.Series(
                    agg_df.loc[sampLabel, 'summarized'][0],
                    name=meta_df.loc[sampLabel,
                                     'sample_collection_datetime'])
            else:
                df_ab_sum = pd.concat([
                    df_ab_sum,
                    pd.Series(agg_df.loc[sampLabel, 'summarized'][0],
                              name=meta_df.loc[sampLabel,
                                               'sample_collection_datetime'])
                ], axis=1)
        else:
            if i == 0:
                df_ab_sum = pd.Series(
                    agg_df.loc[sampLabel, 'summarized'],
                    name=meta_df.loc[sampLabel,
                                     'sample_collection_datetime'])
            else:
                df_ab_sum = pd.concat([
                    df_ab_sum,
                    pd.Series(agg_df.loc[sampLabel, 'summarized'],
                              name=meta_df.loc[sampLabel,
                                               'sample_collection_datetime'])
                ], axis=1)
    df_ab_sum = df_ab_sum.T
    df_ab_sum = df_ab_sum.fillna(0)
    if 'Other' not in df_ab_sum.columns:
        df_ab_sum['Other'] = 0.
    for col in df_ab_sum.columns:
        # TODO: expand to sum values not in top N lineages, etc.
        if col != 'Other':
            if df_ab_sum[col].mean() <= thresh:
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
    dates_to_keep = dates_to_keep.intersection(df_ab_sum.index)

    df_ab_sum = df_ab_sum.groupby(level=0).mean()
    df_ab_sum = df_ab_sum.sort_index()

    return df_ab_lin, df_ab_sum, dates_to_keep


def make_dashboard(agg_df, meta_df, thresh, title, introText,
                   outputFn, headerColor, bodyColor, scale_by_viral_load,
                   config, lineage_info, nboots, serial_interval, days,
                   grthresh, keepPlotFile):
    df_ab_lin, df_ab_sum, dates_to_keep = get_abundance(agg_df, meta_df,
                                                        thresh,
                                                        scale_by_viral_load,
                                                        config, lineage_info)

    calc_rel_growth_rates(df_ab_lin.copy(deep=True), nboots,
                          serial_interval, outputFn,
                          daysIncluded=days, thresh=grthresh)

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
    meta_df = meta_df.set_index('sample_collection_datetime')
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
                                   dict(label="Variants",
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
                                       label="Load Scaled Variants",
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
                                  xanchor="left",
                                  x=1.1,
                                  itemsizing='constant'))

    fig.update_layout(yaxis_title='Variant Prevalence',
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
    # Generate a web page with the plot by placing it in the template.
    locDir = os.path.abspath(os.path.join(os.path.realpath(__file__),
                             os.pardir))
    webpage = open(os.path.join(locDir,
                                'data/dashboard_template.html'),
                   "r").read()
    webpage = webpage.replace("{title}", title)
    webpage = webpage.replace("{introText}", introText)
    webpage = webpage.replace("{plot}", open("div-plot.html", "r").read())
    webpage = webpage.replace("{lastUpdated}",
                              str(datetime.now().strftime("%b-%d-%Y %H:%M")))
    webpage = webpage.replace("{headerColor}",
                              headerColor)
    webpage = webpage.replace("{bodyColor}",
                              bodyColor)
    webpage = webpage.replace("{table}",
                              pd.read_csv(
                                  outputFn.replace('.html',
                                                   '_rel_growth_rates.csv'))
                                .to_html(index=False)
                                .replace('dataframe',
                                         'table table-bordered table-hover' +
                                         ' table-striped table-light' +
                                         ' table-bordered'))

    with open(outputFn, 'w') as outfile:
        outfile.write(webpage)
    if not keepPlotFile:
        os.remove('div-plot.html')
    else:
        print('Keeping intermediate plot')
    print("Dashboard html file saved to " + outputFn)


def collapse_barcodes(df_barcodes, df_depth, depthcutoff,
                      lineageyml, locDir, output,
                      relaxed, relaxedthresh, input_basename):
    # drop low coverage sites
    low_cov_sites = df_depth[df_depth[3].astype(int) < depthcutoff] \
        .index.astype(str)
    low_cov_muts = [mut for mut in df_barcodes.columns
                    if mut[1:-1] in low_cov_sites]
    df_barcodes = df_barcodes.drop(low_cov_muts, axis=1)

    max_depth = df_depth[3].astype(int).max()
    min_depth = df_depth[3].astype(int).min()
    sd_depth = df_depth[3].astype(int).sd()
    

    # find lineages with identical barcodes
    try:
        duplicates = df_barcodes.groupby(df_barcodes.columns.tolist()).apply(
            lambda x: tuple(x.index) if len(x.index) > 1 else None
        ).dropna()
    except ValueError:
        raise ValueError(f'Error: --depthcutoff {depthcutoff} too high '
                         f'for data with max depth {max_depth},'
                         f' min depth {min_depth},'
                         f' standard deviation {sd_depth}')

    if len(duplicates) == 0:
        return df_barcodes

    # load lineage data
    lineage_yml = read_lineage_file(lineageyml, locDir, fileOnly=True)
    lineage_data = {lineage['name']: lineage for lineage in lineage_yml}

    alias_count = {}
    collapsed_lineages = {}

    # collapse lineages into MRCA, where possible
    for tup in duplicates:
        try:
            pango_aliases = [lineage_data[lin]['alias']
                             for lin in tup]
        except KeyError:
            print('Lineage hierarchy file is likely behind'
                  ' the selected barcode file. Try updating'
                  ' the hierarchy file.')
        # handle cases where multiple lineage classes are being merged
        # e.g. (A.5, B.12) or (XBB, XBN)
        multiple_lin_classes = len(
            set([alias.split('.')[0] for alias in pango_aliases])) > 1

        if multiple_lin_classes:
            recombs = [alias for alias in pango_aliases
                       if 'recombinant_parents' in
                       lineage_data[alias.split('.')[0]]]

            # for recombinant lineages, find the parent lineages
            startTypes = set([alias.split('.')[0] for alias in pango_aliases])
            # figure out which are the candidates for recomb merging
            # if they exist
            while len(recombs) > 0:
                parent_aliases = []
                for alias in recombs:
                    if 'recombinant_parents' in lineage_data[alias.split('.')
                                                             [0]]:
                        # trace up tree until a recombination event.
                        # grab parents of recombinant
                        parents = lineage_data[alias.split('.')
                                               [0]]['recombinant_parents'
                                                    ].replace('*', ''
                                                              ).split(',')
                        parent_aliases.append([lineage_data[lin]['alias']
                                               for lin in parents])

                distinct = []
                newRecombs = []
                mergedIn = False
                for alias, pa in zip(recombs, parent_aliases):
                    for aliasP in pa:
                        if aliasP.split('.')[0] in startTypes:
                            mergedIn = True
                            # if now using same start as others,
                            # add to list of aliases
                            pango_aliases.append(aliasP)
                            if alias in pango_aliases:
                                pango_aliases.remove(alias)
                        elif 'recombinant_parents' in lineage_data[aliasP.
                                                                   split('.')
                                                                   [0]]:
                            # check if it's a different recombinant
                            newRecombs.append(aliasP)
                        else:
                            # non-recombinant, but not in current start types.
                            distinct.append(aliasP)
                if not mergedIn:
                    # if no merges, remove the recombinants
                    # and add in the parents
                    for r in recombs:
                        pango_aliases.remove(r)
                    pango_aliases.extend(distinct+newRecombs)

                startTypes = set([alias.split('.')[0]
                                  for alias in pango_aliases])
                if len(startTypes) == 1:
                    break
                recombs = [alias for alias in pango_aliases if
                           'recombinant_parents' in
                           lineage_data[alias.split('.')[0]]]

        if not relaxed:
            mrca = os.path.commonpath(
                [lin.replace('.', '/') for lin in pango_aliases]
            ).replace('/', '.')
        else:
            j0 = 1
            groupCt = float(len(pango_aliases))
            ext_counts = np.unique([lin.split('.')[0:j0]
                                    for lin in pango_aliases],
                                   return_counts=True)
            coherentFrac = np.max(ext_counts[1]) / groupCt
            if coherentFrac < relaxedthresh:
                mrca = ''
            else:
                maxLength = np.max([len(lin.split('.'))
                                    for lin in pango_aliases])
                while coherentFrac >= relaxedthresh and j0 <= maxLength:
                    ext_counts = np.unique([lin.split('.')[0:j0]
                                            if j0 <= len(lin.split('.'))
                                            else lin.split('.') +
                                            ['']*(j0-len(lin.split('.')))
                                            for lin in pango_aliases],
                                           return_counts=True,
                                           axis=0)
                    max_ind = np.argmax(ext_counts[1])
                    coherentFrac = ext_counts[1][max_ind] / groupCt
                    mrca = '.'.join(ext_counts[0][max_ind])
                    j0 += 1

        # assign placeholder if no MRCA found
        if len(mrca) == 0:
            mrca = 'Misc'
        else:
            if mrca[-1] == '.':
                mrca = mrca[0:(len(mrca)-1)]
            # otherwise, get the shortened alias, if available
            for lineage in lineage_data:
                if lineage_data[lineage]['alias'] == mrca:
                    # add flag to indicate that this is a merged lineage
                    mrca = lineage + '-like'
                    break

        # include index for multiple barcode classes with same MRCA
        if mrca not in alias_count:
            alias_count[mrca] = 0
        else:
            alias_count[mrca] += 1
            mrca += f'({alias_count[mrca]})'

        collapsed_lineages[mrca] = list(tup)
        df_barcodes = df_barcodes.rename({lin: mrca for lin in tup})
    df_barcodes = df_barcodes.drop_duplicates()

    # defaults from demix and boot
    if output == 'demixing_result.csv' or output == 'test':
        output = (f'{os.path.dirname(output)}/'
                  f'test_collapsed_lineages.yml')
    else:
        output = (f'{os.path.dirname(output)}/'
                  f'{input_basename}_collapsed_lineages.yml')

    with open(output, 'w') as f:
        yaml.dump(collapsed_lineages, f, default_flow_style=False)

    print(f'collapsed lineages saved to {output}')

    return df_barcodes


def handle_region_of_interest(region_of_interest, output_df,
                              df_depth, covcut, title):

    roi_df = pd.read_json(region_of_interest, orient='index')

    roi_df['start'] = roi_df['start'].astype(int)
    roi_df['end'] = roi_df['end'].astype(int)

    # Ensure start < end
    roi_df['start'], roi_df['end'] = zip(*roi_df.apply(
        lambda x: (x['start'], x['end']) if x['start'] < x['end']
        else (x['end'], x['start']), axis=1))

    # Ensure start > 0 and end < 29903
    roi_df['start'] = roi_df['start'].apply(lambda x: x if x > 0 else 1)
    roi_df['end'] = roi_df['end'].apply(
        lambda x: x if x < 29903 else 29903)

    roi_df.index.name = 'name'

    # Get percent coverage in each region
    roi_cov = pd.Series()
    for _, row in roi_df.iterrows():
        roi_depths = df_depth.loc[(df_depth.index >= row['start']) &
                                  (df_depth.index <= row['end'])]
        roi_cov = pd.concat([roi_cov,
                             pd.Series(
                                 (sum(roi_depths.loc[:, 3] > covcut) /
                                     len(roi_depths)) * 100,
                                 index=[row.name])
                             ])

    # Write to output
    output_df = pd.concat([output_df, roi_cov], axis=0)
    output_df.name = title

    return output_df


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
    bodyColor = 'white'
    with open('freyja/data/plot_config.yml', "r") as f:
        try:
            config = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            raise ValueError('Error in config file: ' + str(exc))
    if os.path.exists('freyja/data/lineages.yml'):
        print('Using existing lineages.yml')
    else:
        raise FileNotFoundError('Update. Could not find lineages.yml')

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
