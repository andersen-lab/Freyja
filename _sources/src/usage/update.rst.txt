.. click:: freyja._cli:update
    :prog: freyja update
    :nested: full
    :commands: update
------------

**Example Usage:**

By default, this method ships with an existing "data/usher_barcodes.csv" file for the barcodes, and the `outbreak.info <https://outbreak.info/>`_ curated lineage metadata file for summarizing lineages by WHO designation. To update both of these we recommend running the command:
``freyja update``
which downloads new versions of the curated lineage file and barcodes (which are now stored on the github repo to save users time). If the --buildlocal flag is used, the barcodes will calculated locally using the UShER global phylogenetic `tree <http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/>`_ and saved in "data/usher_barcodes.csv". The `--outdir` option can be used to specify a local directory to store the lineage mapping and barcode files. By default, Freyja now only includes lineages that are present on `cov-lineages.org <https://cov-lineages.org/>`_. To include proposed lineages and lineages that haven't been released via cov-lineages (usually this lag is no more than a few days), the `--noncl` flag can be used.

**NOTE:** Due to the large size of the global tree, this step can be somewhat memory intensive. Providing somewhere in the range of 10GB should be sufficient to ensure the update runs to completion.
