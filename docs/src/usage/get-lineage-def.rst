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


If a gene annotation file (gff3) is provided along with a reference genome (fasta), the mutations will include the respective amino acid changes.

.. code-block:: bash

    freyja get-lineage-def B.1.1.7 --annot freyja/data/NC_045512_Hu-1.gff --ref freyja/data/NC_045512_Hu-1.fasta
    C241T(None)
    C913T(ORF1ab:S216S)
    C3037T(ORF1ab:F924F)
    C3267T(ORF1ab:T1001I)
    C5388A(ORF1ab:A1708D)
    (cont...)
