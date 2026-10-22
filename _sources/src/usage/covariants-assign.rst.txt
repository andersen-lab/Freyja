.. click:: freyja._cli:covariants_assign
    :prog: freyja covariants-assign
    :nested: full
    :commands: covariants-assign
------------

**Example Usage:**

``freyja covariants`` reports clusters of mutations that co-occur on the
same read pair. ``freyja covariants-assign`` takes that output and assigns
each cluster to the lineages it cannot be told apart from. Each cluster is
compared against the lineage barcodes over the span its reads cover, and is
labelled with the most recent common ancestor (MRCA) of all matching
lineages, using the same MRCA logic as ``freyja demix --depthcutoff``.

.. code-block:: bash

   freyja covariants sample.bam 21563 25384 --output sample_covariants.tsv
   freyja covariants-assign sample_covariants.tsv --output sample

The matching rules can be relaxed with ``--max_missing`` (the number of a
lineage's mutations within the cluster's span that may be absent from the
cluster) and ``--max_extra`` (the number of a cluster's mutations that may
be absent from the lineage barcode). Clusters that match no lineage are
counted as ``Other`` by default; use ``--unassigned drop`` to drop them and
renormalise the abundances over the assigned clusters.

With ``--output sample``, the following files are written:

* ``sample_clusters.tsv``: the input clusters with their assigned class and
  matching lineages.
* ``sample_abundances.tsv``: summed cluster frequency for each assigned class.
* ``sample_collapsed_lineages.yml``: the lineages behind each assigned class.
* ``sample_demix.tsv``: a ``freyja demix`` style output, written when
  ``--demix_output`` is given.
* ``sample_group_abundances.tsv``: abundances per lineage group, written when
  a lineage groups yaml is given via ``--groups`` (same ``name``/``members``
  format as ``plot_config.yml``). ``--assign_by`` selects whether a cluster is
  placed in the group of its MRCA (``mrca``) or in the group holding most of
  its lineages (``majority``).
