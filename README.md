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
freyja variants [bamfile] --variants [variant outfile name] --depths [depths outfile name] --ref [reference.fa]
```
which uses both samtools and iVar. Note that the reference should match the fasta file used for alignment.

We can then run Freyja on the output files using the commmand:
```
freyja demix [variants-file] [depth-file] --output [output-file]
```
This outputs to a tsv file that includes the lineages present, their corresponding abundances, and summarization by constellation. This method also includes a `--eps` option, which enables the user to define the minimum lineage abundance returned to the user (e.g. `--eps 0.0001`). A custom barcode file can be prvoided using the `--barcodes [path-to-barcode-file]` option.  An example output should have the format


|       | filename |
| ----------- | ----------- |
| summarized      | [('Delta', 0.65),  ('Other', 0.25),  ('Alpha', 0.1')] |
| lineages   | ['B.1.617.2' 'B.1.2' 'AY.6' 'Q.3']       |
| abundances   | "[0.5 0.25 0.15 0.1]"|
| resid   | 3.14159        |

Where ```summarized``` denotes a sum of all lineage abundances in a particular WHO designation (i.e. B.1.617.2 and AY.6 abundances are summed in the above example), otherwise they are grouped into "Other". The ```lineage``` array lists the identified lineages in descending order, and  ```abundances``` contains the corresponding abundances estimates. The value of ```resid``` corresponds to the residual of the weighted least absolute devation problem used to estimate lineage abundances. 

---
### Additional options
By default, this method ships with an existing "data/usher_barcodes.csv" file for the barcodes, and the [outbreak.info](https://outbreak.info/) curated lineage metadata file for summarizing lineages by WHO designation. To update both of these we recommend running the command

```
freyja update
```
which downloads new versions of the curated lineage file as well as the UShER global phylogenetic [tree](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/), which is subsequently converted into barcodes and saved in "data/usher_barcodes.csv".

We now provide a fast bootstrapping method for freyja, which can be run using the command

```
freyja boot [variants-file] [depth-file] --nt [number-of-cpus] --nb [number-of-bootstraps] --output_basename [base-name]
```
which results in two output files `base-name_lineages.csv` and `base-name_summarized.csv`, which contain the 0.05,0.25,0.5 (median),0.75, and 0.95 quantiles for each lineage and WHO designated VOI/VOC, respectively, as obtained via the bootstrap. We also provide the `--eps` and `--barcodes` options available in `freyja demix`. 

For rapid visualization of results, we also offer two utility methods for manipulating the "demixed" output files. The first is an aggregation method

```
freyja aggregate [directory-of-output-files] --output [aggregated-filename.tsv]
```
This resulting aggregated data can analyzed directly as a tsv file, or can be visualized using

```
freyja plot [aggregated-filename-tsv] --output [plot-filename(.pdf,.png,etc.)]
```
which provides a fractional abundance estimate for all aggregated samples. To modify the provide a lineage specific breakdown, the `--lineages` flag can be used. We now provide a `--colors [path-to-csv-of-hex-codes]` option so users can control the colors of the plot (see `freyja/data/colors.csv` for an example input file).   Example outputs:

|**Summarized** | **Lineage-Specific**|
|     :---:      |     :---:      |
|![Summarized](freyja/data/testSummary.png) | ![Lineage-Specific](freyja/data/test0.png)|

If users wish to include sample collection time information, this can be done using 

```
freyja plot [aggregated-filename-tsv] --output [plot-filename(.pdf,.png,etc.)] --times [times_metadata.csv(note csv!)] --interval [MS or D (month/day bins)]
```

When using the `--interval D` option, the `--windowsize NN` should also be specified, where `NN` is the width of the rolling average window. See `freyja/data/times_metadata.csv` for an example collection time metadata file. Example outputs:

|**Month binning** | **Daily binning (with smoothing)**|
|     :---:      |     :---:      |
|![Monthly](freyja/data/test2.png) | ![Daily-Smoothed](freyja/data/test.png)|
