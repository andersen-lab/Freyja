Custom Plotting Tutorial
-------------------------------------------------------------------------------

Once you’ve been sampling wastewater for long enough, you might find
that the standard ``freyja plot`` function isn’t sufficiently flexible
to build more complicated and/or curated plots of variant/lineage
dynamics. In this tutorial, we’ll talk a little bit about how to
customize plots as you see fit, using the Freyja Python library.

1. Load in necessary packages and prepare an aggregated output file
   using the ``python freyja aggregate`` command. Note: this will be in
   tsv (tab separated value) format.

.. code:: python

   import pandas as pd
   import os
   import numpy as np
   import matplotlib.pyplot as plt
   import matplotlib
   import matplotlib.dates as mdates

   # make text illustrator 
   matplotlib.rcParams['pdf.fonttype'] = 42
   matplotlib.rcParams['ps.fonttype'] = 42

   agg_df = pd.read_csv('your-aggregated-data.tsv', skipinitialspace=True, sep='\t',index_col=0) 
   ## filter samples by coverage (generally 60-70 is decent cutoff, but can vary across samples)
   agg_df = agg_df[agg_df['coverage']>70.0] 

2. Read the aggregated data in, and convert to a more friendly format
   using the ``python prepLineageDict()`` and
   ``python prepSummarizedDict()`` functions:

.. code:: python

   from freyja.utils import prepLineageDict, prepSummaryDict

   agg_df = prepSummaryDict(agg_df)
   agg_df = prepLineageDict(agg_df)


   # specify which sort of data you're looking to get out
   lineages=True # False = variant summarized, True = lineage specific 
   if lineages:
       ## lineage focused analyses
       queryType = 'linDict'
   else:
       ## analyses summarized by variant
       queryType = 'summarized'

3. Add in dates. Then convert to a single dataframe ordered by date.

.. code:: python

   # load time metadata
   times_df = pd.read_csv(times, skipinitialspace=True,
                                  index_col=0)
           times_df['sample_collection_datetime'] = \
               pd.to_datetime(times_df['sample_collection_datetime'])

   # build a dataframe 
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

   # fill nans, group data by the appropriate interval

   interval = 'MS' # 'MS' is monthly, 'D' is daily
   df_abundances = df_abundances.fillna(0)
   df_abundances = df_abundances.groupby(pd.Grouper(freq=interval)).mean()

Optional: Add in custom groupings as desired. For example, if we want to
group all AY lineages, we can adjust

.. code:: python

   df_abundances['AY.X'] = df_abundances.loc[:,df_abundances.columns.str.startswith('AY.')].sum(axis=1)
   df_abundances = df_abundances.drop(columns= df_abundances.columns[df_abundances.columns.str.startswith('AY.') & ~(df_abundances.columns=='AY.X')])

For more complex groupings, relationships between lineages are made
available via cov-lineages
`here <https://github.com/cov-lineages/lineages-website/blob/master/data/lineages.yml>`__
and can be grouped accordingly.

4. Build your custom plot, using either an existing colormap (like tab20
   in matplotlib) or a custom colormap (see colors0 in the example
   below).

.. code:: python


   windowSize = 14 # set window size for rolling average if using daily interval 

   # replace "example" with whatever you'd like to call the output file
   outputFn = 'example_' + queryType + '.pdf'

   colors0 = []# if empty, use cmap
   #example for user specified colors
   # colors0 = [#DAF7A6,#FFC300,#FF5733,#1B4F72,#808B96,#45B39D]
   cmap = plt.cm.tab20
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
               #this should work for up to 20 colors
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
