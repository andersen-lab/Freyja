.. click:: freyja._cli:dash
    :prog: freyja dash
    :nested: full
    :commands: dash
------------

**Example Usage:**

We are now providing functionality to rapidly prepare a dashboard web
page, directly from aggregated freyja output. This can be done with the
command

::

   freyja dash [aggregated-filename-tsv] [sample-metadata.csv] [dashboard-title.txt] [introContent.txt] --output [outputname.html] --lineage.yml [path-to-lineage.yml-file]

where the metadata file should have this
`form <freyja/data/sweep_metadata.csv>`__. See example
`title <freyja/data/title.txt>`__ and
`intro-text <freyja/data/introContent.txt>`__ files as well. For samples
taken the same day, we average the freyja outputs by default. However,
averaging can be performed that takes the viral loads into account using
the ``--scale_by_viral_load`` flag. The header and body color can be
changed with the ``--headerColor [mycolorname/hexcolor]`` and
``--bodyColor [mycolorname/hexcolor]`` option respectively. The
``--mincov`` option is also available, as in ``plot``. The resulting
dashboard will look like
`this <https://htmlpreview.github.io/?https://github.com/andersen-lab/Freyja/blob/main/freyja/data/test0.html>`__.

The plot can now be configured using the
``--config [path-to-plot-config-file]`` option. The `plot config
file <freyja/data/plot_config.yml>`__ is a yaml file. More information
about the plot config file can be found in the `sample config
file <freyja/data/plot_config.yml>`__. By default, this will use the
lineage hierarchy information present in ``freyja/dash/lineages.yml``,
but a custom hierarchy can be supplied using the
``--lineageyml [path-to-hierarchy-file]`` option. The
``--keep_plot_files`` option can be used keep the intermediate html for
the core plot (will be deleted following incorporation into the main
html output by default).

A CSV file will also be created along with the html dashboard which will
contain the relative growth rates for each lineage. The lineages will be
grouped together based on the ``Lineages`` key specified in the config
file if provided.