# VariantAnnotator
The main goal of this pipeline is to annotate variants from a vcf file using a genbank reference file. The output is a csv file with the predicted coding changes along with the consequences of the change (synonymous/non-synonymous). 

## Authors
- William Shropshire (https://github.com/wshropshire)  
- Selvalakshmi (https://github.com/Selvalakshmi27)

## Dependencies
The R script was built using R version 4.2.0. The following are the required dependencies:
 + VariantAnnotation (1.42.1)
 + genbankr (1.24.0)
 + tidyverse (1.3.2)
 + GenomicFeatures (1.48.3)
 + randomcoloR (1.1.0.1)
 + gtools (3.9.2.2)
 + GenomeInfoDb (1.32.2)
 + dplyr (1.0.9)

## Usage
For linux 
```
Rscript VariantAnnotation.R [input.vcf] [reference.gb] [output.csv]
```
- Download the VariantAnnotation.R script from above
- input.vcf = The vcf file can be generated from tools such as snippy, samtools, freebayes etc
- reference.gb = The reference genome should be in genbank full format. 
- output.csv = Specify output file name in csv format. 

The example output file can be found above as coding-changes.csv
