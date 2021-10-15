# Depth-weighted De-Mixing (DWDM)
DWDM is a tool to recover relative lineage abundances from mixed SARS-CoV-2 samples from a sequencing dataset (BAM aligned to the Hu-1 reference). The method uses  lineage-determining mutational "barcodes" derived from the UShER global phylogenetic tree as a basis set to solve the constrained (unit sum, non-negative) de-mixing problem. 

DWDM is intended as a post-processing step after primer trimming and variant calling in [iVar (Grubaugh and Gangavaparu et al., 2019)](https://github.com/andersen-lab/ivar). From measurements of SNV freqency and sequencing depth at each position in the genome, DWDM returns an estimate of the true lineage abundances in the sample.   

## Installation
DWDM is entirely written in Python 3, but requires preprocessing by tools like iVar and samtools mpileup to generate the required input data. Successful installation of iVar (available via conda) should be sufficient to perform all required steps. 

### Dependencies
* [UShER](https://usher-wiki.readthedocs.io/en/latest/#)
* [cvxpy](https://www.cvxpy.org/)
* [tqdm](https://github.com/tqdm/tqdm)
* [numpy](https://numpy.org/)
* [pandas](https://pandas.pydata.org/)
