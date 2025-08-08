Running Freyja on other pathogens
-------------------------------------------------------------------------------

This guide provides instructions for analyzing non-SARS-CoV-2 pathogens such as
influenza or MPox using Freyja. The process is similar to SARS-CoV-2 analysis,
but with some key differences.

Data Availability
^^^^^^^^^^^^^^^^^

Data for various pathogens can be found in the following repository:
`Freyja Barcodes <https://github.com/gp201/Freyja-barcodes>`_

Folders are organized by pathogen, with each subfolder named after the date the
barcode was generated, using the format ``YYYY-MM-DD``. Barcode files are named
``barcode.csv``, and reference genome files are named ``reference.fasta``.

.. note::
        Influenza barcodes are available upon request.

Required Files
^^^^^^^^^^^^^^

To perform these analyses, you will need the following files for the MPox pathogen:

*       `test.sorted.bam <https://github.com/andersen-lab/Freyja/blob/main/docs/data/test.sorted.bam>`_: Aligned, trimmed, and sorted BAM file
*       `reference.fasta <https://github.com/gp201/Freyja-barcodes/blob/main/MPX/2024-07-24/reference.fasta>`_: Reference genome file
*       `barcode.csv <https://github.com/gp201/Freyja-barcodes/blob/main/MPX/2024-07-24/barcode.csv>`_: Barcode file


Setting Up Output Directories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since you will likely be working with multiple wastewater samples, it is
advisable to create directories for storing output files:

.. code-block:: sh

        mkdir variants_files depth_files demix_files


Analysis Steps
^^^^^^^^^^^^^^

The first step is to generate a variant file. Use the following command to
perform this step:

.. code-block:: sh

        freyja variants test.sorted.bam --ref reference.fasta --variants variants_files/test.tsv --depths depth_files/test.depth

Please note that you will be passing the reference genome file provided in the
pathogen folder as the ``--ref`` argument. In cases where multiple reference
genomes are present in the reference fasta, you can specify the name of the
desired reference genome with ``--refname [name-of-reference]``.

Once the variant file is generated, proceed to the de-mixing step with the
following command:

.. code-block:: sh

        freyja demix variants_files/test.tsv depth_files/test.depth --barcodes barcode.csv --output demix_files/test.output

Please note that you will be passing the barcode file provided in the pathogen
folder as the ``--barcodes`` argument.

Once you’ve run demix on a bunch of samples, you can aggregate all of
the output files using the command

.. code-block:: sh

        freyja aggregate demix_files/ --output bunch_of_files.tsv

From there, it’s easy to view the output files in any standard TSV viewer
(Excel, Numbers, LibreOffice Calc, etc.). You should see something like this:

.. code-block::

                summarized      lineages        abundances      resid   coverage
        test.tsv        [('Other', 0.999999999530878)]  MPX-A.3 MPX-A.2.2       0.79798000 0.20202000   7.5952064496123075      99.94117915510955
