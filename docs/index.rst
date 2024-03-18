.. Freyja documentation master file, created by
   sphinx-quickstart on Wed Feb 21 08:25:27 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Freyja Documentation
==================================
Freyja is a tool to recover relative lineage abundances from mixed SARS-CoV-2 samples from a sequencing dataset (BAM aligned to the Hu-1 reference). The method uses lineage-determining mutational "barcodes" derived from the UShER global phylogenetic tree as a basis set to solve the constrained (unit sum, non-negative) de-mixing problem.

Freyja is intended as a post-processing step after primer trimming and variant calling in iVar (Grubaugh and Gangavaparu et al., 2019). From measurements of SNV freqency and sequencing depth at each position in the genome, Freyja returns an estimate of the true lineage abundances in the sample.

To ensure reproducibility of results, we provide old (timestamped) barcodes and metadata in the separate Freyja-data repository. Barcode version can be checked using the freyja demix --version command.

.. toctree::
   :maxdepth: 2
   :caption: Freyja:
   
   installation
   usage

.. toctree::
   :maxdepth: 2
   :caption: Wiki:

   command_line_workflow
   cryptic_variants
   custom_plotting_tutorial
   lineage_barcode_extract
   read_analysis_tutorial
   terra_workflow
