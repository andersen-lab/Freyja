Freyja for other pathogens- Measles
-------------------------------------------------------------------------------

If you want to detect other pathogens using freyja,
this is a toturial on how to run the pipeline for Measles samples.

1. Download the reference genome for your pathogen. It can be downloaded
from NCBI database or use the references provided by `Nextstrain <https://nextstrain.orgL>`_
For the purpose of this toturial, we will be using `NC_001498.1 <https://www.ncbi.nlm.nih.gov/nuccore/NC_001498.1>`_


2. Download the latest lineage barcodes from `Freyja barcodes repository <https://github.com/gp201/Freyja-barcodes>`_
For this analysis we will be using `Measles barcodes <https://github.com/gp201/Freyja-barcodes/tree/main/MEASLESN450>`_
There are two sets of barcodes available for Measles, one may use N450 region barcodes but it is important to
use the right reference genome according to the barcode.

3. Prepare your primer bed format file. For this toturial we will be using `artic-measles V.1.0.0 panel <https://labs.primalscheme.com/detail/artic-measles/400/v1.0.0/?q=measles>`_

4. Align your reads to your reference genome using an aligner of your choice. 
Here, we use minimap2 with parameters set for aligning short reads to a reference genome.
We will be using simulated reads produced using `this amplicon read simulation pipeline <https://github.com/mariaelf97/amplicon_sequencing_simulator>_

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

8. Remove primers from the reads. The following command will remove reads with mean 
quality score `-q` less than 15 with sliding window `-s` size of 4 and minimum read 
length of 50 (after trimming). `e` will make sure reads without primers are kept in the output. 
Please note that ivar trim does not trim the sequencing adapters and users need to make sure to do 
it as a pre-processing step on the reads fastq files.

.. code::

    ivar trim -b primer.bed -p trimmed -i aligned_sorted.bam -q 15 -m 50 -s 4 -e

9. Sort and index the trimmed bam file.

.. code::

    samtools sort -o trimmed_sorted.bam trimmed.bam && samtools index trimmed_sorted.bam

10. Generate coverage depth for each single genomic location in the reference.
This will be used in freyja pipeline downstream analysis(only if you are using ivar variant caller embedded in freyja)

.. code::

    samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f reference.fasta trimmed_sorted.bam | cut -f1-4 > depths.tsv

11. Call variants using either Lofreq or ivar. We have included both commands for your reference.
`-a` sets the p-value cut-off with Bonferroni correction.

.. code::

    lofreq call -f reference.fasta -a 1 -b 1 -o variants.vcf trimmed_sorted.bam 
    freyja variants trimmed_sorted.bam  --variants variants.tsv  --depths depths.tsv --ref reference.fasta

12. Run freyja demix to estimate lineage prevalence. Please note that the barcodes.tsv file is binary mutation
matrix downloaded from the `freyja barcodes repository <https://github.com/gp201/Freyja-barcodes>`_

.. code::

    freyja demix variants.tsv depths.tsv --output freyja_demix.txt --barcodes barcodes.csv