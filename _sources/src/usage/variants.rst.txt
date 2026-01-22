.. click:: freyja._cli:variants
    :prog: freyja variants
    :nested: full
    :commands: variants
------------

**Example Usage:**


After primer trimming in iVar, we get both variant call and sequencing depth information with the command:
``freyja variants [bamfile] --variants [variant outfile name] --depths [depths outfile name] --ref [reference.fa]``
which uses both samtools and iVar. Note that the reference should match the fasta file used for alignment. In cases where multiple reference genomes are present in the reference fasta, the user can specify the name of the desired reference genome with ``--refname [name-of-reference]``. To enable alternative variant calling methods ( such as `LoFreq <https://csb5.github.io/lofreq/>`_),  we also allow users to provide a VCF file using the ``--variants`` option (in addition to the usual depth file, which can be obtained using a command like ``samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f ref.fasta sample.bam | cut -f1-4 > sample.depth``).