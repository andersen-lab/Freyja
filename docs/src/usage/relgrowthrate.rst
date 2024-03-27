.. click:: freyja._cli:relgrowthrate
    :prog: freyja relgrowthrate
    :nested: full
    :commands: relgrowthrate
------------

**Example Usage:**

This command will generate the relative growth rates for each
Lineage and output it as a CSV file. The lineages can also be grouped
together by passing the ``--config [path-to-plot-config-file]`` option.
This uses the same configuration file as the ``dash`` command. The
lineages are grouped together based on the ``Lineage:`` field in the
config file. The number of bootstraps can be specified with the
``--nboots [number-of-bootstraps]`` option and the serial interval can
be specified with the ``--serial_interval`` option. The output will be a
CSV file with the following columns:

+----------------------+---------------------+------------------------+
|        Lineage       | Estimated Advantage | Bootstrap 95% interval |
+----------------------+---------------------+------------------------+
