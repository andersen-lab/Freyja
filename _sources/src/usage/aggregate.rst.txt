.. click:: freyja._cli:aggregate
    :prog: freyja aggregate
    :nested: full
    :commands: aggregate
------------

**Example Usage:**

For rapid visualization of results, we also offer two utility methods
for manipulating the “demixed” output files. The first is an aggregation
method

::

   freyja aggregate [directory-of-output-files] --output [aggregated-filename.tsv]

By default, the minimum genome coverage is set at 60 percent. To adjust
this, the ``--mincov`` option can be used (e.g. ``--mincov 75``.We also
now allow the user to specify a file extension of their choosing, using
the ``--ext`` option (for example, for ``demix`` outputs called
``X.output``)

::

   freyja aggregate [directory-of-output-files] --output [aggregated-filename.tsv] --ext output


