Running Freyja on other pathogens - Measles
----

If you want to detect other pathogens using freyja,
this is a tutorial  on how to run the pipeline for Measles samples.


**STEPS**
****

1. Download the reference genome for your pathogen.
It can be downloaded
from NCBI database or use the references provided by `Nextstrain <https://nextstrain.org>`_ .
For the purpose of this tutorial, we will be using `NC_001498.1. <https://www.ncbi.nlm.nih.gov/nuccore/NC_001498.1>`_
The reference genome file can also be found `here. <https://github.com/andersen-lab/Freyja/blob/main/docs/data/measles-reference.fasta>`_


2. Download the latest lineage barcodes from `Freyja barcodes repository. <https://github.com/gp201/Freyja-barcodes>`_
For this analysis we will be using `Measles barcodes <https://github.com/gp201/Freyja-barcodes/tree/main/MEASLESgenome>`_
which can also be found `here. <https://github.com/andersen-lab/Freyja/blob/main/docs/data/measles-wg-barcode.csv>`_

Please note that there are two sets of barcodes available for Measles. Here, we use the whole genome barcodes but
you may also use the N450 region barcodes. It is important, however, to
use the right reference genome according to the barcode used. If you are using region N450 barcodes, please make sure to
use the reference that includes N450 region only.

3. Prepare your primer file in a bed format. For this tutorial, we will be using `artic-measles V1.0.0 panel <https://labs.primalscheme.com/detail/artic-measles/400/v1.0.0/?q=measles>`_
also provided `here. <https://github.com/andersen-lab/Freyja/blob/main/docs/data/artic-measles-v1.0.0.bed>`_
This file will be used to simulate measles amplicon sequencing reads and trim the primers in the downstream analysis.

4. Prepare your reads, by assessing quality of the reads and removing the sequencing adapters.
We will be using simulated reads produced by `MixAmp, <https://github.com/andersen-lab/MixAmp>`_ an amplicon read simulator developed by our team.

For instance, we generated reads from two different measles samples using this code as following:


* Pure sample: `Reads <https://github.com/andersen-lab/Freyja/blob/main/docs/data/GCA_031128185.1-simulated.fastq>`_ simulated from a H1 lineage sample individually.

* Mixed samples: `Reads <https://github.com/andersen-lab/Freyja/blob/main/docs/data/measles-mixed-simulated.fastq>`_ with mixing proportions H1:20% and D9:80%.


*Pure sample read simulation:*

.. code::

    mixamp simulate-proportions GCA_031128185.1.fna primer.bed --outdir measles-H1-100/

*Mixed sample read simulation:*

.. code::
    
    mixamp simulate-proportions GCA_031128185.1.fna,GCA_031129565.1.fna primer.bed --outdir measles-H1-20-D9-80/ --proportions 0.2,0.8


5. Align your reads to your reference genome using an aligner of your choice. 
Here, we use ``minimap2`` with parameters set for aligning short reads to a reference genome.

.. code::

    minimap2 -ax sr reference.fasta reads.fastq > aligned.sam

^^^^

**From here on, the instruction will be the same for mixed and pure samples.
Please change the file names accordingly.**

6. Convert sam format file to a bam file using samtools.

.. code:: 

   samtools view -bS aligned.sam > aligned.bam

7. Sort your bam file.

.. code:: 

    samtools sort -o aligned_sorted.bam aligned.bam

8. Index your bam file.

.. code::

    samtools index aligned_sorted.bam

9. Remove primers from the reads. The following command will remove reads with mean
quality score `-q` less than 15 and minimum read length of 50 (after trimming).
Flag `-e` keeps reads without primers in the output and 
the flag `-x` will make sure that reads that occur at the 
specified offset positions relative to primer positions will also be trimmed.


.. code::

    ivar trim -b primer.bed -p trimmed -i aligned_sorted.bam -q 15 -m 50 -e -x 3

10. Sort and index the trimmed bam file.

.. code::

    samtools sort -o trimmed_sorted.bam trimmed.bam && samtools index trimmed_sorted.bam

11. Generate coverage depth for each single genomic location in the reference.
This will be used in freyja pipeline downstream analysis.

.. code::

    samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f reference.fasta trimmed_sorted.bam | cut -f1-4 > depths.tsv

12. Call variants using variant caller of your choice. We recommend using Lofreq or ivar and included both commands for your reference.

.. code::

    # call variants using lofreq
    lofreq call -f reference.fasta -a 1 -b 1 -o variants.vcf trimmed_sorted.bam
    # Call variants using ivar
    freyja variants trimmed_sorted.bam --variants variants.tsv --depths depths.tsv --ref reference.fasta

13. Run freyja demix to estimate lineage prevalence.

.. code::

    freyja demix variants.tsv depths.tsv --output freyja_demix.txt --barcodes barcodes.csv


The final demix outputs for the pure and mixed sample are as following:

*Mixed sample output:*

.. code::

    summarized      [('Other', 0.9999999968413253)]
    lineages        MEASLES-D9 MEASLES-H1
    abundances      0.79692605 0.20307394
    resid   214.51679168207156
    coverage        91.39927016484208

*Pure sample output*

.. code::

    summarized      [('Other', 0.999999999926792)]
    lineages        MEASLES-H1
    abundances      1.00000000
    resid   53.868769540487826
    coverage        89.52434881087203
