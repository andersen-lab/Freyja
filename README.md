# Freyja
Freyja is a tool to recover relative lineage abundances from mixed SARS-CoV-2 samples from a sequencing dataset (BAM aligned to the Hu-1 reference). The method uses  lineage-determining mutational "barcodes" derived from the UShER global phylogenetic tree as a basis set to solve the constrained (unit sum, non-negative) de-mixing problem. 

Freyja is intended as a post-processing step after primer trimming and variant calling in [iVar (Grubaugh and Gangavaparu et al., 2019)](https://github.com/andersen-lab/ivar). From measurements of SNV freqency and sequencing depth at each position in the genome, Freyja returns an estimate of the true lineage abundances in the sample.   

## Installation
Freyja is entirely written in Python 3, but requires preprocessing by tools like iVar and [samtools](https://github.com/samtools/samtools) mpileup to generate the required input data. Successful installation of iVar (available via conda) should be sufficient to perform all required steps. 

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
samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f NC_045512_Hu-1.fasta [filename.trimmed.bam] | tee >(cut -f1-4 > [outpath.depth]) | ivar variants -p [outpath] -q 20 -r NC_045512_Hu-1.fa 
```

We can then run Freyja on the output files using the commmand:
```
freyja demix [variants-file] [depth-file] --output [output-file]
```
This outputs in a tsv file, which includes the lineages present, their corresponding abundances, and summarization by constellation. 

---
### Additional options
1. By default, this method will use the existing "data/usher_barcodes.csv" file for the barcodes, and the [outbreak.info](https://outbreak.info/) curated lineage metadata file for summarizing lineages by WHO designation. To update both of these we recommend running the command

```
freyja update

```
which downloads new versions of the UShER global phylogenetic [tree](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/) as well as the curated lineage file. 

Lineage defining mutation barcodes are extracted using 
```
matUtils extract -i public-latest.all.masked.pb.gz -C lineagePaths.txt
```
These are converted to a new barcode set by 
```
freyja barcode lineagePaths.txt
```
which saves the new barcodes as "usher_barcodes.csv". 

2. For summarizing of lineages by constellation, we pull directly from the [outbreak.info](https://outbreak.info/) curated lineage metadata file. To pull a new one, just run

```
wget -N https://raw.githubusercontent.com/outbreak-info/outbreak.info/master/web/src/assets/genomics/curated_lineages.json
```

---

Acknowledgements

