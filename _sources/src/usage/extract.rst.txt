.. click:: freyja._cli:extract
    :prog: freyja extract
    :nested: full
    :commands: extract
------------

**Example Usage:**

The above command will extract reads containing one or more mutations of
interest and save them to ``[input-bam]_extracted.bam``. Additionally,
the ``--same_read`` flag can be included to specify that all query
mutations must occur on the same set of paired reads (same UMI).

Formatting requirements for ``[query-mutations.csv]`` can be
found below, where SNPs, insertions, and deletions are listed in
separate lines.

::

   C75T,G230A,A543C
   (732:'TT'),(1349:'A'),(12333:'A')
   (1443:32),(1599:2),(2036:3)
