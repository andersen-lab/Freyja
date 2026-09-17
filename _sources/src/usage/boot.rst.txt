.. click:: freyja._cli:boot
    :prog: freyja boot
    :nested: full
    :commands: boot
------------

**Example Usage:**

We provide a fast bootstrapping method for freyja, which can be run
using the command

::

   freyja boot [variants-file] [depth-file] --nt [number-of-cpus] --nb [number-of-bootstraps] --output_basename [base-name]

which results in two output files: ``base-name_lineages.csv`` and
``base-name_summarized.csv``, which contain the 0.025, 0.05, 0.25, 0.5
(median),0.75, 0.95, and 0.975 percentiles for each lineage and WHO
designated VOI/VOC, respectively, as obtained via the bootstrap. A
custom lineage hierarchy file can be provided using ``--lineageyml``
option. If the ``--rawboots`` option is used, it will return two
additional output files ``base-name_lineages_boot.csv`` and
``base-name_summarized_boot.csv``, which contain the bootstrap estimates
(rather than summary statistics). We also provide the ``--eps``,
``--barcodes``, and ``--meta`` options as in ``freyja demix``. We now
also provide a ``--boxplot`` option, which should be specified in the
form ``--boxplot pdf`` if you want the boxplot in pdf format.