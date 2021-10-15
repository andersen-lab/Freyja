# Depth-Weighted De-Mixing (DWDM)
DWDM is a tool to recover relative lineage abundances from mixed SARS-CoV-2 samples from a sequencing dataset (BAM aligned to the Hu-1 reference). The method uses  lineage-determining mutational "barcodes" derived from the UShER global phylogenetic tree as a basis set to solve the constrained (unit sum, non-negative) de-mixing problem. 

DWDM is intended as a post-processing step after primer trimming and variant calling in [iVar (Grubaugh and Gangavaparu et al., 2019)](https://github.com/andersen-lab/ivar). From measurements of SNV freqency and sequencing depth at each position in the genome, DWDM returns an estimate of the true lineage abundances in the sample.   

## Installation
DWDM is entirely written in Python 3, but requires preprocessing by tools like iVar and [samtools](https://github.com/samtools/samtools) mpileup to generate the required input data. Successful installation of iVar (available via conda) should be sufficient to perform all required steps. 

### Dependencies
* [UShER](https://usher-wiki.readthedocs.io/en/latest/#)
* [cvxpy](https://www.cvxpy.org/)
* [tqdm](https://github.com/tqdm/tqdm)
* [numpy](https://numpy.org/)
* [pandas](https://pandas.pydata.org/)

## Usage
After primer trimming in iVar, we get both variant call and sequencing depth information with the command:
```
samtools mpileup -aa -A -d 600000 -B -Q 0 test.trimmed.bam | tee >(cut -f1-4 > test.depth) | ivar variants -p test -q 20 -r NC_045512_Hu-1.fa 
```

We can then run DWDM on the output files using the commmand:
```
python sample_deconv.py variant_tsvs/ depth_files/ output_result.tsv
```
This results in a tsv file, which includes the lineages present and their corresponding abundances. 

By default, this method will use the existing "usher_barcodes.csv" file for the barcodes. To make a new barcode library, download the latest global phylogenetic tree from UShER: http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/. 

Lineage defining mutation barcodes are extracted using 
```
matUtils extract -i public-latest.all.masked.pb.gz -C lineagePaths.txt
```
and these are converted to a new barcode set by 
```
python convert_paths2barcodes.py lineagePaths.txt
```
which saves the new barcodes as "usher_barcodes.csv". 


