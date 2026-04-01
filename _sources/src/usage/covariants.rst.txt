.. click:: freyja._cli:covariants
    :prog: freyja covariants
    :nested: full
    :commands: covariants
------------

**Example Usage:**

In many cases, it can be useful to study covariant mutations
(i.e. mutations co-occurring on the same read pair). This outputs to a tsv file that includes the mutations present in each
set of covariants, their absolute counts (the number of read pairs with
the mutations), their coverage ranges (the minimum and maximum position
for read-pairs with the mutations), their “maximum” counts (the number
of read pairs that span the positions in the mutations), and their
frequencies (the absolute count divided by the maximum count). The ``--reference`` argument defaults
to ``freyja/data/NC_045512_Hu-1.fasta``. If you are using a different
build to perfrom alignment, it is important to pass that file in to
``--reference`` instead. Additionally, a gff file
(e.g. ``freyja/data/NC_045512_Hu-1.gff``) must be included via the
``--annot`` option to output amino acid mutations alongside
nucleotide mutations. Inclusion thresholds for read-mapping quality and
the number of observed instances of a set of covariants can be set using
``--min_quality`` and ``--min_count`` respectively.