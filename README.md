# Freyja
Freyja is a tool to recover relative lineage abundances from mixed SARS-CoV-2 samples from a sequencing dataset (BAM aligned to the Hu-1 reference). The method uses  lineage-determining mutational "barcodes" derived from the UShER global phylogenetic tree as a basis set to solve the constrained (unit sum, non-negative) de-mixing problem. 

Freyja is intended as a post-processing step after primer trimming and variant calling in [iVar (Grubaugh and Gangavaparu et al., 2019)](https://github.com/andersen-lab/ivar). From measurements of SNV freqency and sequencing depth at each position in the genome, Freyja returns an estimate of the true lineage abundances in the sample.   

## Installation via conda
Freyja is entirely written in Python 3, but requires preprocessing by tools like iVar and [samtools](https://github.com/samtools/samtools) mpileup to generate the required input data. First, create an environment for freyja
```
conda create -n freyja-env
```
then add the following channels
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
and then install freyja
```
conda install freyja
```
### Dependencies
* [iVar](https://github.com/andersen-lab/ivar)
* [samtools](https://github.com/samtools/samtools)
* [UShER](https://usher-wiki.readthedocs.io/en/latest/#)
* [cvxpy](https://www.cvxpy.org/)
* [numpy](https://numpy.org/)
* [pandas](https://pandas.pydata.org/)

## Usage
After primer trimming in iVar, we get both variant call and sequencing depth information with the command:
```
freyja variants [bamfile] --variants [variant outfile name] --depths [depths outfile name]
```
which uses both samtools and iVar. 

We can then run Freyja on the output files using the commmand:
```
freyja demix [variants-file] [depth-file] --output [output-file]
```
This outputs to a tsv file that includes the lineages present, their corresponding abundances, and summarization by constellation. This method also includes a `--eps` option, which enables the user to define the minimum lineage abundance returned to the user (e.g. `--eps 0.0001`).  

---
### Additional options
By default, this method ships with an existing "data/usher_barcodes.csv" file for the barcodes, and the [outbreak.info](https://outbreak.info/) curated lineage metadata file for summarizing lineages by WHO designation. To update both of these we recommend running the command

```
freyja update
```
which downloads new versions of the curated lineage file as well as the UShER global phylogenetic [tree](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/), which is subsequently converted into barcodes and saved in "data/usher_barcodes.csv".

For rapid visualization of results, we also offer two utility methods for manipulating the "demixed" output files. The first is an aggregation method

```
freyja aggregate [directory-of-output-files] --output [aggregated-filename.tsv]
```
This resulting aggregated data can analyzed directly as a tsv file, or can be visualized using

```
freyja plot [aggregated-filename-tsv] --output [plot-filename(.pdf,.png,etc.)]
```
which provides a fractional abundance estimate for all aggregated samples. To modify the provide a lineage specific breakdown, the `--lineages` flag can be used. Example outputs:

|**Summarized** | **Lineage-Specific**|
|     :---:      |     :---:      |
|![Summarized](freyja/data/testSummary.png) | ![Lineage-Specific](freyja/data/test0.png)|

If users wish to include sample collection time information, this can be done using 

```
freyja plot [aggregated-filename-tsv] --output [plot-filename(.pdf,.png,etc.)] --times [times_metadata.csv(note csv!)] --interval [MS or D (month/day bins)]
```

When using the `--interval D` option, the `--windowsize NN` should also be specified, where `NN` is the width of the rolling average window. See `freyja/data/times_metadata.csv` for an example collection time metadata file. Example outputs:

|**Month Binning** | **Daily binning (with smoothing)|
|     :---:      |     :---:      |
|![Monthly](freyja/data/test2.png) | ![Daily-Smoothed](freyja/data/test.png)|