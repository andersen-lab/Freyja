Freyja Documentation
==================================
Freyja is a tool to recover relative lineage abundances from mixed SARS-CoV-2 samples from a sequencing dataset (BAM aligned to the Hu-1 reference). The method uses lineage-determining mutational "barcodes" derived from the UShER global phylogenetic tree as a basis set to solve the constrained (unit sum, non-negative) de-mixing problem.

Freyja is intended as a post-processing step after primer trimming and variant calling in `iVar (Grubaugh and Gangavaparu et al., 2019) <https://github.com/andersen-lab/ivar>`_. From measurements of SNV freqency and sequencing depth at each position in the genome, Freyja returns an estimate of the true lineage abundances in the sample.

To ensure reproducibility of results, we provide old (timestamped) barcodes and metadata in the separate `Freyja-data <https://github.com/andersen-lab/Freyja-data>`_ repository. Barcode version can be checked using the ``freyja demix --version`` command.

.. toctree::
   :maxdepth: 2
   :caption: Usage:
   
   src/installation
   src/usage/demix
   src/usage/variants
   src/usage/barcode-build
   src/usage/get-lineage-def
   src/usage/update
   src/usage/boot
   src/usage/aggregate
   src/usage/plot
   src/usage/dash
   src/usage/relgrowthrate
   src/usage/extract
   src/usage/filter
   src/usage/covariants
   src/usage/plot-covariants

.. toctree::
   :maxdepth: 2
   :caption: Wiki:

   src/wiki/command_line_workflow
   src/wiki/cryptic_variants
   src/wiki/custom_plotting_tutorial
   src/wiki/freyja-measles
   src/wiki/custom-plots-with-R
   src/wiki/lineage_barcode_extract
   src/wiki/read_analysis_tutorial
   src/wiki/terra_workflow
   src/wiki/expanded_pathogens

