.. click:: freyja._cli:filter
    :prog: freyja filter
    :nested: full
    :commands: filter
------------

**Example Usage:**

Excludes reads containing one or more mutations, akin to ``freyja extract`` but with the opposite effect.

``[min-site]`` and ``[max-site]`` specify the range of reads to
include. Formatting requirements for ``[query-mutations.csv]`` can be
found below, where SNPs, insertions, and deletions are listed in
separate lines.

::

   C75T,G230A,A543C
   (732:'TT'),(1349:'A'),(12333:'A')
   (1443:32),(1599:2),(2036:3)