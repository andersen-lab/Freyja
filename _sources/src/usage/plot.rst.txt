.. click:: freyja._cli:plot
    :prog: freyja plot
    :nested: full
    :commands: plot
------------

**Example Usage:**

This resulting aggregated data can analyzed directly as a tsv file, or
can be visualized using

::

   freyja plot [aggregated-filename-tsv] --output [plot-filename(.pdf,.png,etc.)]

which provides a fractional abundance estimate for all aggregated
samples. To provide a specific lineage breakdown, the ``--lineages``
flag can be used. We now provide a
``--config [path-to-plot-config-file]`` option that allows users to
control the colors and grouping of lineages in the plot. The `plot
config file <freyja/data/plot_config.yml>`__ is a yaml file. More
information about the plot config file can be found in the `sample
config file <freyja/data/plot_config.yml>`__. Example outputs:

.. |testSummary| image:: ../../../freyja/data/testSummary.png
.. |test0| image:: ../../../freyja/data/test0.png

+-----------------+------------------------------+
|  **Summarized** |     **Lineage-Specific**     |
+=================+==============================+
| |testSummary|   | |test0|                      |
+-----------------+------------------------------+

If users wish to include sample collection time information, this can be
done using

::

   freyja plot [aggregated-filename-tsv] --output [plot-filename(.pdf,.png,etc.)] --times [times_metadata.csv(note csv!)] --interval [MS or D (month/day bins)] --lineageyml [path-to-lineage.yml-file]

A custom lineage hierarchy file can be provided using ``--lineageyml``
option for visualization purposes. When using the ``--interval D``
option, the ``--windowsize NN`` should also be specified, where ``NN``
is the width of the rolling average window. See
``freyja/data/times_metadata.csv`` for an example collection time
metadata file. Example outputs:

.. |test2| image:: ../../../freyja/data/test2.png
.. |test| image:: ../../../freyja/data/test.png

+--------------------+--------------------------------------------+
|  **Month binning** |     **Daily binning (with smoothing)**     |
+====================+============================================+
| |test2|            | |test|                                     |
+--------------------+--------------------------------------------+