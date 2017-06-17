# practical-microbiomess-2017
Practical session at the [2017 Microbiome Summer School](http://metagenomic.ca).

## Introduction

## GWAS with PLINK

### Data loading and pre-processing
Download the data using
```sh
wget TODO
```

Start by checking the files are intact and plink works, and get some basic statistics on your data.
```sh
plink --noweb --file simulated
```
From this command, PLINK understands it is going to find the genotype data in `simulated.ped` and SNP descriptions under `simulated.map`.

#### File formats
* `.ped`: TODO
* `.map`: TODO

#### Quality control

Now we're going to apply __quality control__ filters:
* SNPs with minor allele frequency (maf) lower than 1% will be removed.
* SNPs with missing data for more than 10% of individuals will be removed.
* SNPs that are not in Hardy-Weinberg equilibrium (p-value larger than 1e-6) will be removed.

We're also creating a binary ped file, called a `.bed` file, that will take up less space and speed up subsequent analyses.
```sh
plink --noweb --file simulated --maf 0.01 --hwe 1e-6 --geno 0.1 --make-bed --out mydata
```

The output should look like
```
TODO                           

```
We have hence created the following `mydata` files:
* `mydata.bed`
* `mydata.bim`
* `mydata.fam`
as well as `mydata.log` which contains the log that you also see on screen, and which can be quite useful to keep track of which parameters you used exactly.

#### File formats
* `.bed`: TODO
* `.bim`: TODO
* `.fam`: TODO

How many SNPs do we have left to analyze?

```sh
plink --noweb --bfile mydata
```
will give you the answer, which is TODO.

### GWAS
```sh
plink --noweb --bfile mydata --assoc --out assoc1
```
Results in `assoc1.assoc1`.

You can have a look at the contents of this file using
```sh
more assoc1.assoc
```

The results look like:
```
TODO
```
% where each row is a single SNP association result. The fields are:
%     Chromosome
%     SNP identifier
%     Code for allele 1 (the minor, rare allele based on the entire sample frequencies)
%     The frequency of this variant in cases
%     The frequency of this variant in controls
%     Code for the other allele
%     The chi-squared statistic for this test (1 df)
%     The asymptotic significance value for this test
%     The odds ratio for this test

Sort results by p-value

```sh
  sort --key=8 -r --general-numeric-sort assoc1.assoc | head -20
```

`--key=8` because CHISQ is the 8th column, `-r1` to reverse order (larger first), `--general-numeric-sort` to sort numerically and not alphabetically, `head -20` to only see the top 10 results.

The most significant SNP is rs1015896 with a p-value of 2.266e-18.

## Reproduce PLINK analysis in Python

## Epistasis detection

## Lasso

Recode SNPs for usage as a data matrix
```sh
  plink --noweb --bfile mydata --recodeA --out mydata1
```
Creates `mydata1.raw`

One line per individual. The first six columns are the usual
```
FID IID PAT MAT SEX PHENOTYPE
```
and then there is one column per SNP.
The first line is a line of header.



## Further resources
* Microbiome\ practical.ipynb is a Jupyter notebook containing all the data generation and processing steps.
