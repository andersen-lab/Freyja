# Freyja
[![freyja CI](https://github.com/andersen-lab/Freyja/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/andersen-lab/Freyja/actions/workflows/python-package-conda.yml) [![Anaconda-Server Badge](https://anaconda.org/bioconda/freyja/badges/version.svg)](https://anaconda.org/bioconda/freyja) ![docs](https://github.com/andersen-lab/Freyja/actions/workflows/update_docs.yml/badge.svg)

Detailed documentation, including installation, usage, and examples can be found [here](https://andersen-lab.github.io/Freyja/index.html#).

Freyja is a tool to recover relative lineage abundances from mixed SARS-CoV-2 samples from a sequencing dataset (BAM aligned to the Hu-1 reference). The method uses  lineage-determining mutational "barcodes" derived from the UShER global phylogenetic tree as a basis set to solve the constrained (unit sum, non-negative) de-mixing problem. 

Freyja is intended as a post-processing step after primer trimming and variant calling in [iVar (Grubaugh and Gangavaparu et al., 2019)](https://github.com/andersen-lab/ivar). From measurements of SNV freqency and sequencing depth at each position in the genome, Freyja returns an estimate of the true lineage abundances in the sample.

To ensure reproducibility of results, we provide old (timestamped) barcodes and metadata in the separate [Freyja-data](https://github.com/andersen-lab/Freyja-data) repository. Barcode version can be checked using the `freyja demix --version` command.


## Installation via conda
Freyja is entirely written in Python 3, but requires preprocessing by tools like iVar and [samtools](https://github.com/samtools/samtools) mpileup to generate the required input data. We recommend using python3.7, but Freyja has been tested on python versions up to 3.10.  First, create an environment for freyja
```
conda create -n freyja-env
```
then add the following channels
```
conda config --add channels bioconda
conda config --add channels conda-forge
```
and then install freyja
```
conda install freyja
```