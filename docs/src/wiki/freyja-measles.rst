Freyja for other pathogens- Measles
-------------------------------------------------------------------------------

If you want to detect other pathogens using freyja,
this is a toturial on how to run the pipeline for Measles samples.

1. Download the reference genome for your pathogen. It can be downloaded
from NCBI database or use the references provided by `Nextstrain <https://nextstrain.orgL>`_
For the purpose of this toturial, we will be using `NC_001498.1 <https://www.ncbi.nlm.nih.gov/nuccore/NC_001498.1`_


2. Download the latest lineage barcodes from `Freyja barcodes repository <https://github.com/gp201/Freyja-barcodes>`_

3. Prepare your primer bed format file. For this toturial we will be using `artic-measles V.1.0.0 panel <https://labs.primalscheme.com/detail/artic-measles/400/v1.0.0/?q=measles>`_

4. Align your reads to your reference genome using an aligner of your choice. 
Here, we use minimap2 with parameters set for aligning short reads to a reference genome.

.. code::

    minimap2 -ax sr reference.fasta reads.fastq > aligned.sam

5. Convert sam format file to a bam file using samtools

.. code:: 

   samtools view -bS aligned.sam > aligned.bam

6. Sort your bam file using samtools.

.. code:: 

    samtools sort -o aligned_sorted.bam aligned.bam

7. Index your bam file.

.. code::

    samtools index aligned_sorted.bam

8. Remove primers from the reads

.. code::

    ivar trim -b primer_file.bed -p file_prefix -i aligned_sorted.bam -q 15 -m 50 -s 4 -e