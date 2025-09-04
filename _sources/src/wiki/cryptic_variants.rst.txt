Cryptic Variant Detection
-------------------------------------------------------------------------------

Oftentimes, it is possible to detect certain variants in wastewater
weeks before they show up in clinical sequencing (if at all). In this
tutorial, we will utilize freyja’s covariants utility to detect cryptic
mutations in wastewater, and use the `Python Outbreak
API <https://github.com/outbreak-info/python-outbreak-info>`__ in order
to compare these results against clinical sequencing data published to
GISAID.

Let’s begin by creating a fresh conda environment and installing the
requisite packages:

.. code:: bash

   conda create -n cryptic-variants
   conda activate cryptic-variants

   conda install -c bioconda freyja
   pip install python-outbreak-info==1.0.1 

For this analysis, we’re going to be looking at a sample collected from
Point Loma, San Diego on Jan. 17th, 2022, available
`here <https://www.ncbi.nlm.nih.gov/sra/?term=SRR18541029>`__ via SRA.
We’ll also need the SARS-CoV-2 `reference
genome <https://www.ncbi.nlm.nih.gov/nuccore/1798174254>`__ and `gff
file <https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz>`__.
Next, run alignment, trimming, sorting, and indexing as outlined in the
Command Line Workflow tutorial. Once finished, you should have both the
processed bam file, along with the corresponding index file
(filename.bam.bai). We are now ready to run ``freyja covariants``, which
looks for SNVs that occur together on the same read. For our purposes,
we’re mainly interested in the spike gene, so min_site is set to 21563,
and max_site is set to 25384. We’ll also pass in the same reference
genome used for alignment along with the gff file for gene-level
mutation annotations.

.. code:: bash

   freyja covariants SRR18541029_trimmed.bam 21563 25384 --ref-genome NC_045512_Hu-1.fasta --gff-file NC_045512_Hu-1.gff --output SRR18541029_covariants.tsv

Let’s inspect the output using any spreadsheet software (e.g. Excel,
Numbers). Upon first glance, we see a lot of mutations characteristic of
Omicron sequences occuring together in the sample (S:A67V, S:DEL69/70,
S:E484A, S:N501Y etc.). Of note, we see a few variants that are
relatively uncommon in most other samples, such as S:G593D. A quick
ctrl+F for this SNV tells us that it has been detected alone, and also
alongside S:D614G and S:T572N in separate instances, giving us more
confidence that the SNV we’re observing is indeed present in the sample.

For our next step, we’ll be using the Python outbreak API to look at the
clinical prevalence of this SNV and deduce which lineage it most likely
belongs to. First, we’ll need to authenticate via GISAID by running the
following:

.. code:: bash

   python
   >>> from outbreak_data import authenticate_user
   >>> authenticate_user.authenticate_new_user()

This opens a browser prompting you to enter your GISAID login
credentials. Once finished, create a new python file and paste the
following:

.. code:: python

   import pandas as pd
   from outbreak_data import outbreak_data

   df = outbreak_data.mutations_by_lineage(mutation='S:G593D')
   print(df.head())

::

     pangolin_lineage  lineage_count  mutation_count  proportion  proportion_ci_lower  proportion_ci_upper
   0            xbb.1          27879               1    0.000036             0.000004             0.000168

As we can see, this mutation has only been detected once, and in a
single sequence, belonging to the xbb.1 lineage. It would be good to
validate this SNV by seeing if it is found in other samples, but for
now, we can be fairly confident that this is a real mutation.
