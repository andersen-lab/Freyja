.. click:: freyja._cli:demix
    :prog: freyja demix
    :nested: full
    :commands: demix
------------
  
**Example Usage:**

After running ``freyja variants`` we can run:
``freyja demix [variants-file] [depth-file] --output [output-file]``

This outputs to a tsv file that includes the lineages present, their
corresponding abundances, and summarization by constellation. This
method also includes a ``--eps`` option, which enables the user to
define the minimum lineage abundance returned to the user
(e.g. ``--eps 0.0001``). A custom barcode file can be provided using the
``--barcodes [path-to-barcode-file]`` option. By default, freyja uses
the lineage hierarchy file located in\ ``freyja/data`` directory which
is updated everytime the ``freyja update`` command is run. The user,
however, can define a custom lineage hierarchy file
using\ ``--lineageyml [path-to-lineage-file]``. Users can get the
historic ``lineage.yml`` file at freyja-data GitHub repository
`here <https://github.com/andersen-lab/Freyja-data/tree/main/history_lineage_hierarchy>`_.
As the UShER tree now included proposed lineages, we now offer the
``--confirmedonly`` flag which removes unconfirmed lineages from the
analysis. For additional flexibility and reproducibility of analyses, a
custom lineage-to-constellation mapping metadata file can be provided
using the ``--meta`` option. A coverage depth minimum can be specified
using the ``--depthcutoff`` option, which excludes sites with coverage
less than the specified value. An example output should have the format

+-------------+------------------------------------------------------+
|             | filename                                             |
+=============+======================================================+
| summarized  | [('Delta', 0.65),  ('Other', 0.25),  ('Alpha', 0.1)] |
+-------------+------------------------------------------------------+
| lineages    | ['B.1.617.2' 'B.1.2' 'AY.6' 'Q.3']                   |
+-------------+------------------------------------------------------+
| abundances  | "[0.5 0.25 0.15 0.1]"                                |
+-------------+------------------------------------------------------+
| resid       | 3.14159                                              |
+-------------+------------------------------------------------------+
| coverage    | 95.8                                                 |
+-------------+------------------------------------------------------+

Where ``summarized`` denotes a sum of all lineage abundances in a particular WHO designation (i.e. B.1.617.2 and AY.6 abundances are summed in the above example), otherwise they are grouped into "Other". The ``lineage`` array lists the identified lineages in descending order, and  ``abundances`` contains the corresponding abundances estimates. Using the ``--depthcutoff`` option may result in some distinct lineages now having identical barcodes, which are grouped into the format ``[lineage]-like(num)`` (based on their shared phylogeny) in the output. A summary of this lineage grouping is outputted to ``[output-file]_collapsed_lineages.yml``. The value of ``resid`` corresponds to the residual of the weighted least absolute deviation problem used to estimate lineage abundances. The ``coverage`` value provides the 10x coverage estimate (percent of sites with 10 or greater reads- 10 is the default but can be modfied using the ``--covcut`` option in ``demix``). If there is an solver error during the `demix` step (generally associated with poor data quality), an error message will be returned, along with an output empty summarized, lineages, and abundances, and with resid = -1. 

**NOTE**: The ``freyja variants`` output is stable in time, and does not need to be re-run to incorporate updated lineage designations/corresponding mutational barcodes, whereas the outputs of ``freyja demix`` will change as barcodes are updated (and thus ``demix`` should be re-run as new information is made available).
