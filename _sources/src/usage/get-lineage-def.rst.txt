.. click:: freyja._cli:get_lineage_def
    :prog: freyja get-lineage-def
    :nested: full
    :commands: get_lineage_def
    
------------

**Example Usage:**

This command fetches the lineage-defining mutations for the specified lineage. Defaults to stdout if no output file is specified.

.. code-block:: bash

    freyja get-lineage-def B.1.1.7
    C241T
    C913T
    C3037T
    C3267T
    C5388A
    (cont...)


If a gene annotation file (gff3) is provided along with a reference genome (fasta), the mutations will include the respective amino acid changes for non-synonymous mutations.

.. code-block:: bash

    freyja get-lineage-def BA.2.12 --annot freyja/data/NC_045512_Hu-1.gff --ref freyja/data/NC_045512_Hu-1.fasta
    ...
    C21618T(S:T19I)
    T22200G(S:V213G)
    G22578A(S:G339D)
    C22674T(S:S371F)
    T22679C(S:S373P)
    C22686T(S:S375F)
    A22688G(S:T376A)
    G22775A(S:D405N)
    A22786C(S:R408S)
    G22813T(S:K417N)
    (cont...)
