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
plink --noweb --bfile simulated --maf 0.01 --hwe 1e-6 --geno 0.1 --make-bed --out mydata
```

The output should look like
```
83534 (of 83534) markers to be included from [ simulated.map ]
89 individuals read from [ simulated.ped ]
89 individuals with nonmissing phenotypes
Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
Missing phenotype value is also -9
42 cases, 47 controls and 0 missing
89 males, 0 females, and 0 of unspecified sex
Before frequency and genotyping pruning, there are 83534 SNPs
89 founders and 0 non-founders found
0 markers to be excluded based on HWE test ( p <= 1e-06 )
      0 markers failed HWE test in cases
      0 markers failed HWE test in controls
Total genotyping rate in remaining individuals is 0.99441
859 SNPs failed missingness test ( GENO > 0.1 )
16994 SNPs failed frequency test ( MAF < 0.01 )
After frequency and genotyping pruning, there are 65803 SNPs
After filtering, 42 cases, 47 controls and 0 missing
After filtering, 89 males, 0 females, and 0 of unspecified sex
Writing pedigree information to [ mydata.fam ]
Writing map (extended format) information to [ mydata.bim ]
Writing genotype bitfile to [ mydata.bed ]
Using (default) SNP-major mode
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
CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
 1   rs6681049          1    1   0.2024   0.2234    2       0.1168       0.7326        0.882
 1   rs4074137          2    1  0.08333  0.07447    2      0.04811       0.8264         1.13
 1   rs1891905          4    1   0.4167   0.3936    2      0.09784       0.7544          1.1
 1   rs9729550          5    1   0.1548   0.1064    2       0.9227       0.3368        1.538
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

The most significant SNP is rs12045968.

## Reproduce PLINK analysis in Python

## Epistasis detection

## Lasso

Recode SNPs for usage as a data matrix
```sh
  plink --noweb --bfile mydata --recodeA
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
