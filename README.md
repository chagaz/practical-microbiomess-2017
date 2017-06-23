# practical-microbiomess-2017
Practical session on efficient biomarker discovery at the [2017 Microbiome Summer School](http://metagenomic.ca).

## Introduction
The goal of this practical is to manipulate quantitative GWAS data and start exploring how machine learning algorithms can be used to analyze this data. We will be working with the genotypes of 89 individuals from the [1000 Genomes Project](http://www.internationalgenome.org/data) (Han Chinese and Japanese ancestry), and a simulated quantitative phenotype. This phenotype can be imagined to represent the relative abundance of two microbial species in the gut of the host.

## GWAS with PLINK

One of the most well-known pieces of software for analyzing GWAS data is [PLINK](http://zzz.bwh.harvard.edu/plink/), developed maninly by [Shaun Purcell](http://zzz.bwh.harvard.edu/) at Harvard, MGH and the Broad Institute.

[PLINK 1.9](https://www.cog-genomics.org/plink/1.9/) is currently in beta version.

### Data loading and pre-processing
Uncompress the data using
```sh
tar zxvf simulated-gwas.tar.gz
```
The data is stored in PLINK format (a format that is commonly used to exchange genotype/phenotype data and that most GWAS software can manipulate).

#### File formats
* `.ped`: The samples data. Contains  as many lines as samples in the data and `6 + 2 x num_snps` columns.
The first 6 columns contain the following information: Family identifier (`FID`), individual identifier (`IID`), paternal identifier (`PAT`), maternal identifier (`MAT`), sex (`SEX`; male=1, female=2, unknown=other) and phenotype (`PHENOTYPE`). The following columns contain all bi-allelic SNP information. Each SNP is coded on 2 columns, each corresponding to one strand of DNA. The SNP can be encoded `A, T, C, G` or `1, 2` (corresponding to one or the other allele).
* `.map`: The markers data. Contains as many lines as SNPs, and 4 columns per SNP: chromosome, SNP identifier, genetic distance in morgans, and base-pair position.

PLINK formats are fully documented [here](http://zzz.bwh.harvard.edu/plink/data.shtml).


Start by checking the files are intact and plink works, and get some basic statistics on your data.
```sh
plink --noweb --file simulated
```
From this command, PLINK understands it is going to find the genotype data in `simulated.ped` and SNP descriptions under `simulated.map`.

#### Quality control

Now we're going to apply __quality control__ filters:
* SNPs with __minor allele frequency (MAF)__ lower than 1% will be removed. We focus on common variants for several reasons: the "common disease, common variant" hypothesis; the fact that rare variants are more likely to be technical artifacts; and, last but not least, because we have limited statistical power to detect the effect of rare SNPs.
* SNPs with __missing data__ for more than 10% of individuals will be removed.
* SNPs that are not in __Hardy-Weinberg equilibrium (HWE)__ (p-value larger than 1e-6) will be removed: departure from HWE is likely to be due to a genotyping error.

We're also creating a binary ped file, called a `.bed` file, that will take up less space and speed up subsequent analyses.
```sh
plink --file simulated --maf 0.01 --hwe 1e-6 --geno 0.1 --make-bed --out mydata
```

__Question__ How many SNPs passed quality control?

<details>
<summary>Click to see answer</summary>
Answer: 66 536 out of 83 534 (you can get this from the output of the PLINK command, either on screen or in `mydata.log`.
</details>


We have hence created the following `mydata` files:
* `mydata.bed`, a binary file that contains genotype information;
* `mydata.bim`, an extended `.map` file with to extra columns containing the allele names (in our case, always `1` and `2`);
* `mydata.fam`, the first 6 columns of the `.ped` file.

as well as `mydata.log` which contains the log that you also see on screen, and which can be quite useful to keep track of which parameters you used exactly.

### GWAS
Let us now use PLINK to test for statistical association between each SNP and the phenotype.
```sh
plink --noweb --bfile mydata --assoc --out assoc1
```
This creates a file called `assoc1.qassoc` (the `q` stands for "quantitative").

You can have a look at the contents of this file using
```sh
more assoc1.qassoc
```

The results look like:
```
CHR         SNP         BP    NMISS       BETA         SE         R2        T            P
  1   rs6681049          1       89   -0.04122    0.04808   0.008377  -0.8573       0.3936
  1   rs4074137          2       89   -0.01389    0.07988  0.0003473  -0.1738       0.8624
  1   rs1891905          4       89   -0.01533    0.04204   0.001527  -0.3647       0.7162
  1   rs9729550          5       89    0.03083    0.06171   0.002861   0.4996       0.6186
  1   rs3813196          6       89    0.09187     0.1362   0.005198   0.6743       0.5019
  1  rs12044597         11       89    0.02528     0.0501    0.00292   0.5047        0.615
[...]
```
For each SNP, PLINK fitted a regression between this SNP and the phenotype. It outputs among other things:
* BETA, the regression coefficient, which can be interpreted as _effect size_;
* R2, the correlation coefficient between the true phenotype and the phenotype predicted by the regression, which can be interpreted as the _proportion of (phentoypic) variance explained_ by the SNP;
* P, the _p-value_ of the [Wald test](https://www.wikiwand.com/en/Wald_test) that evaluates whether beta is significantly different from zero.

### GWAS results analysis in Python

We are going to use Python2.7 to analyze the output of PLINK.

Launch the interactive Python terminal with
```sh
  ipython2
```
__Tip:__ Once you have copied code, you can paste it in the ipython terminal by simply typing `%paste`. (This command is called a "magic" in ipython, and will keep code indentation.)

We're going to use a Python package called `pandas` to manipulate data. Pandas manipulates data in so-called "data frames". In your ipython terminal, type or copy-paste:

```python
# Load the pandas package
import pandas as pd
# Load the output of PLINK as a pandas dataframe called df
df = pd.read_csv('assoc1.qassoc', # File name
                 delim_whitespace=True) # Use any white space as delimiter
```

You can now visualize the object you have loaded with
```python
print df
```

Let us sort the SNPs by p-value (smallest first), and print the 10 SNPs with the smallest p-values.

```python
df_sorted = df.sort_values('P')
print df_sorted[:10]
```

What is our significance threshold? We need to account for __multiple hypothesis testing.__ If we call `r` the probability of getting a false positive when running our test, the probability of _not_ getting one is ``(1-r)``. Now if we run p tests, the probability of _not_ getting any false positive is `(1-r)**p`. Therefore the probability of getting at least one false positive is `1 - (1-r)**p`, which grows steadily with p and tends towards 1.

If you're curious, you can plot this function with matplotlib:
```python
import numpy as np
from matplotlib import pyplot as plt
# false positive probability
r = 0.01
# create a vector of 50 values, equally spaced between 1 and 500
x = np.arange(1, 500, 50)
# create a vector of (1- (1-r)**p) for all values of p stored in x
y = (1- (1-r)**x)
# plot y against x
plt.plot(x, y)
# visualize the plot
plt.show(block=False)
```

To compensate for this, we can use the __Bonferroni correction__, and divide the significance threshold (0.05) by the number of statistical tests we have run, that is to say, the number of SNPs.

```python
# df has as many lines as there are SNPs
num_snps = df.shape[0]
sig_thresh = 0.05/num_snps
print "Significance threshold = ", sig_thresh
```  
We can use this threshold to identify the SNPs that are significantly associated with the phenotype:
```python
sig_SNPs = df.SNP[df.P < sig_thresh]
print sig_SNPs
```
Hence we have found two SNPs that are significantly associated with the phenotype: rs1015896 and rs920160.

#### Manhattan plot
We can visualize our results with a Manhattan plot: we will plot all SNPs according to, along the x-axis, their genomic position (chromosomes will be situated next to each other) and, along the y-axis, their p-value. We will also color the SNPs by chromosome.

```python
import numpy as np
from matplotlib import pyplot as plt
plt.scatter(df.BP, # x = SNP position (in bp)
            -np.log10(df.P), # y = -log10 p-value (the higher the more significant)
            c=df.CHR) # color by chromosome

# Plot a line corresponding to -log10 the significance threshold
plt.plot([0, max(df.BP)], [-np.log10(sig_thresh), -np.log10(sig_thresh)], lw=2)

# Prettify
plt.xlabel("SNP", fontsize=14)
plt.ylabel("-log10 p-value", fontsize=14)
plt.xlim([0, max(df.BP)])

plt.show(block=False)

# If you want to save the figure:
# plt.savefig('manhattan.png', bbox_inches='tight')
```

What are the names of our significant SNPs?
```python
sig_SNPs = df.SNP[df.P < sig_thresh]
print sig_SNPs
```

#### Q-Q plot
Quantile-quantile (Q-Q) plots allow us to visualize the distribution of p-values, and whether it significantly deviates from the uniform. We expect the vast majority of SNPs _not_ to be associated with the phenotype, and hence their p-values to be uniformly distributed. A visible deviation from the uniform usually indicate that the analysis is confounded by population structure.

Let us plot a Q-Q plot in Python:
```python
import scipy.stats as ss
ss.probplot(df.P, dist="uniform", plot=plt)
plt.show(block=False)
```

Our Q-Q plot perfectly matches the diagonal line. This means there is no deviation from the uniform distribution, and very little chance for population structure confounding.

In practice, here this is due to how we simulated the phenotypes; it is quite unlikely that a cohort mixing Han Chinese and Japanese samples will not suffer from confounding due to population structure.  

The two main ways to correct for population structure are
* to use principal components (computed for example with [EIGENSTRAT](https://github.com/DReichLab/EIG/)) as covariates (`--covar` option of PLINK).
* to use linear mixed models, for example with [FastLMM](https://github.com/MicrosoftGenomics/FaST-LMM).


## Machine learning approaches

In what follows, we will work with an _additive_ encoding of the SNPs data:
* 0 for homozygous in the major allele;
* 1 for heterozygous;
* 2 for homozygous in the minor allele.

We can use PLINK to re-encode our SNPs.
```sh
  plink --noweb --bfile mydata --recodeA --out mydata1
```
This creates a file called `mydata1.raw` which contains one line per individual. The first six columns are the usual `FID IID PAT MAT SEX PHENOTYPE` fields,  and then there is one column per SNP. The first line is a line of header.

Excluding the first 6 columns and the header line, this is exactly the type of X matrix we usually represent data with in machine learning.

### Load the data in Python
We will be working with the `scikit-learn` machine learning suite for Python. This requires representing
* the phenotype as a 1-dimensional array `y` of length `num_samples`  
* the genotypes as a `num_snps x num_samples` 2-dimensional array `X`.

```python
import numpy as np
import pandas as pd
# Read the data file in a pandas data frame:
df = pd.read_csv('mydata1.raw',
                  delim_whitespace=True)

# Build X from the values of df, excluding its first 6  columns.
X = df.iloc[:, 6:].values
# Check the dimensions of X
print "X has shape", X.shape

# Build y from the 'PHENOTYPE' column of df
y = df.PHENOTYPE.values
print "y has shape", y.shape
```

PLINK has given to each SNP column a name of the format `rsXXX_A` where `rsXXX` is the SNP's name and `A` is `1` or `2` depending on which one is the minor allele in this data. To be able to map these column names to the names of SNPs we know, we need to remove the `_1` or `_2` suffix.

```python
new_col_names = list(df.columns[:6])
new_col_names.extend([col_name[:-2] for col_name in df.columns[6:]])
df.columns = new_col_names
```

### Percentage of variance explained by the significant SNPs
Before starting on more complicated models, we can use `scikit-learn` to determine the fraction of the phenotypic variance among our 89 individuals that our two significant  SNPs explain.

We first need to build a _linear model_ that uses only these two SNPs.
```python
# Identify the indices of the significant SNPs in X:
# they are the indices of the significant SNPs in df, minus 6 because we excluded the first 6 columns from df when building X.
sig_indices = [(df.columns.get_loc(snp)-6) for snp in sig_SNPs]
# Restrict X to only the 2 significant SNPs
X_sig = X[:, sig_indices]

# Create a linear model that uses the two significant SNPs
from sklearn import linear_model
model = linear_model.LinearRegression()
```

We then need to _fit_ this model to the data
```python
model.fit(X_sig, y)
```

Let us now use this model to make predictions.
```python
y_pred = model.predict(X_sig)
```

How close are we from the _true_ phenotype? We can first check this visually with a plot:
```python
plt.scatter(y, y_pred)
plt.xlabel("True phenotype", fontsize=14)
plt.ylabel("Predicted phenotype", fontsize=14)
plt.title("Phenotype predicted from 2 significant SNPs")
plt.show(block=False)
```

We can also quantify how well predictions match true values by the __proportion of variance explained__:
```python
from sklearn import metrics
print metrics.explained_variance_score(y, y_pred)
```

In this simulation, the two most significant SNPs explain about 76% of the phentoypic variance. To put this in perspective, the SNPs that have been associated with human height explain only about 5% of the phenotypic variance — although we known about 80% of height is inherited. This is typical of the _missing heritability_ problem.

We are now going to try to identify more SNPs to explain heritability.

### Lasso
For this purpose, we're going to use a Lasso. The Lasso finds a linear combination of the SNPs that best fits the phenotype, and has a built-in mechanism to create _sparse_ models, meaning that it will use at few SNPs as possible.

The Lasso tries to minimize a sum of two terms: the mean squared error between predicted phenotype and true phenotype, and a so-called regularization term. The regularization term is the sum of the absolute values of the coefficients of this linear regression, and minimizing it will encourage "useless" features to have a coefficient of zero.

The relative importance of both terms (error and regularization) is controlled by a _regularization parameter_, called `alpha` in `scikit-learn`. The greater `alpha`, the more important the regularization term, and the fewer SNPs will enter the model.

We will start by playing with several values of `alpha`, and then we'll see how to set it automatically.

```python
# Create a Lasso model with regularziation parameter 1.0
lasso = linear_model.Lasso(alpha=0.05)
# Fit the model to the data
lasso.fit(X, y)
# Indices of the SNPs with non-zero coefficients
selected_snps = np.where(lasso.coef_)[0]
# See how many SNPs have a non-zero coefficient
print len(selected_snps), "selected SNPs"
# Plot the coefficients
plt.scatter(range(lasso.coef_.shape[0]), # x-axis = SNPs
            lasso.coef_, # y-axis = SNP weight
            )
plt.show(block=False)
```
With `alpha=0.05`, the Lasso selected 3 SNPs.

__Question:__ What was the significance of each of these 3 SNPs on the Manhattan plot from PLINK?

<details>
<summary>Click to see answer</summary>
Answer: The Lasso built a model with 3 SNPs, two of which were the most significantly associated SNPs by PLINK. The third is not the third-most signficant SNP from the single-SNP analysis!
</details>

How well do these three SNPs explain the phenotype?

```python
# Fit a linear model with to the significant SNPs
model = linear_model.LinearRegression()
model.fit(X[:, selected_snps], y)

# Predict phenotype with this model
y_pred = model.predict(X[:, selected_snps])

# Plot predictions against true values
plt.scatter(y, y_pred)
plt.show(block=False)

# Percentage of variance explained
print metrics.explained_variance_score(y, y_pred)
```

We now explain 97% of the phenotypic variance with 3 SNPs!

__Question:__ What happens when `alpha=0.1`?

<details>
<summary>Click to see answer</summary>
Answer: We retrieve the linear model with the two most significant SNPs from earlier on.
</details>

__Question:__ What happens when `alpha=0.02`?

<details>
<summary>Click to see answer</summary>
Answer: We build a model that includes 6 SNPs and explains 98% of the variance.
</details>

#### Cross-validated Lasso
So, which value of `alpha` should we choose? One way to do this is by _cross-validation_: in a grid of possibilities, we will choose the one that gives the best performing model, in a cross-validation setting.

```python
# Use a variant of Lasso with inner 5-fold cross-validation to set the alpha parameter
lasso_cv = linear_model.LassoCV(cv=5)
# Fit the model to the data
lasso_cv.fit(X, y)

# What is the optimal alpha value?
print "Optimal alpha:", lasso_cv.alpha_

# Indices of the SNPs with non-zero coefficients
selected_snps = np.where(lasso_cv.coef_)[0]
# See how many SNPs have a non-zero coefficient
print len(selected_snps), "selected SNPs"

# Plot the coefficients
plt.scatter(range(lasso_cv.coef_.shape[0]), # x-axis = SNPs
           lasso_cv.coef_)
plt.show(block=False)

# Fit a linear model to the significant SNPs
model = linear_model.LinearRegression()
model.fit(X[:, selected_snps], y)
# Predict with this model
y_pred = model.predict(X[:, selected_snps])
# Percentage of variance explained
print metrics.explained_variance_score(y_pred, y)
```

__Question:__ What was the optimal value for `alpha`? How does it compare to the other values of `alpha` you have tested? In consequence, were you expecting to select more or fewer SNPs than before?

<details>
<summary>Click to see answer</summary>
Answer: The optimal value for alpha is 0.01. This is less than the values we tried before. With less regularization, the model will be less sparse: more SNPs will be selected. Indeed, we now select 32 SNPs. (Your answer may be different because of randomness in the cross-validation procedure.)
</details>

#### Overfitting
Wow! We're explaining almost all the phenotypic variance with only 32 SNPs! But could this be due to _overfitting_? We are evaluating our linear regression model on exactly the same data we used to build it. Would it really work that well on individuals we have not seen? To test for this, we will repeat the above experiments using only two thirds of our cohort for _discovery_ and the remaining third for _validation_.

__Caveat__ Using only two thirds of our samples for discovery reduces the statistical power of our analyses...

Let us separate our data `(X, y)` in a discovery (train) and validation sets:
```python
from sklearn import model_selection
X_disc, X_val, y_disc, y_val = model_selection.train_test_split(X, y, test_size=0.33,
                                random_state=42)
```

Let us fit our cross-validated Lasso again. Because there are only 59 samples in the discovery set, we will use a 3-fold cross-validation only.
```python
lasso_cv = linear_model.LassoCV(cv=3)
lasso_cv.fit(X_disc, y_disc)
```

__Question:__ What is now the optimized alpha value?

<details>
<summary>Click to see answer</summary><p>

Answer:
```python
print "Optimal alpha:", lasso_cv.alpha_
```
alpha=0.02. (Your answer may be different because of randomness in the cross-validation procedure.)

</p></details>

__Question:__ How many SNPs have we selected?

<details>
<summary>Click to see answer</summary><p>

Answer:
```python
# Indices of the SNPs with non-zero coefficients
selected_snps = np.where(lasso_cv.coef_)[0]
# See how many SNPs have a non-zero coefficient
print len(selected_snps), "selected SNPs"
```
18 SNPs. (Your answer may be different because of randomness in the cross-validation procedure.)

</p></details>

__Question:__ What percentage of the variance do we explain in the discovery set?
<details>
  <summary>Click to see answer</summary><p>

Answer:
```python
# Fit a linear model to the significant SNPs
model = linear_model.LinearRegression()
model.fit(X_disc[:, selected_snps], y_disc)
# Predict on the discovery set
y_pred = model.predict(X_disc[:, selected_snps])

# Percentage of variance explained
print metrics.explained_variance_score(y_pred, y_disc)
```
99.4% (Your answer may be different because of randomness in the cross-validation procedure.)

</p></details>

__Question:__ What percentage of the variance do we explain in the _validation_ set?
<details>
<summary>Click to see answer</summary><p>

Answer:
```python
# Use the model to predict on the validation set
y_pred = model.predict(X_val[:, selected_snps])

# Percentage of variance explained
print metrics.explained_variance_score(y_pred, y_val)
```
89.5% (Your answer may be different because of randomness in the cross-validation procedure.). This is significantly less than on the discovery set.

</p></details>

## Going further

Here are a few pointers towards questions you might want to explore:

### Stability
If you apply the above techniques (PLINK single-SNP GWAS, Lasso) to 90% of the samples, selected at random, and repeat this procedure multiple times, do you systematically select the same SNPs?

The Lasso is notoriously unstable in the presence of correlation between the variables, and it is quite likely that you will get different SNPs except maybe for the 2-3 top ones. What are the implications in terms of interpretability?

### Elastic Net
One way to stabilize the Lasso is to use the so-called Elastic Net, which mixes the sparsity effect of the Lasso with the grouped selection effect of an l2-regularizer (ridge regression).

It is implemented in `scikit-learn.linear_model` as `ElasticNet` (documented [here](http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html)). How does it affect the number of selected SNPs? The percentage of variance explained? The stability as observed above?

### Random forests
Random forests result in non-linear models, and have an intrinsic measure of feature importance. Applied to GWAS data, feature importance measures how important a SNP is towards explaining the phenotype, _in the context of all other SNPs_.

Unlike linear regressions, random forests are meant to deal with discrete variables and are therefore quite interesting for GWAS problems.

Unfortunately, random forests also suffer from the instability described above, and require large numbers of trees to be effective on data sets with large numbers of features.

In `scikit-learn`, random forests are implemented in `ensemble`, as [`RandomForestClassifier`](http://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html#sklearn.ensemble.RandomForestClassifier) for classification problems (case/control phenotypes) and [`RandomForestRegressor`](http://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestRegressor.html#sklearn.ensemble.RandomForestRegressor) for regression problems (quantitative phenotypes). Feature importance is given by the `feature_importances_` attribute.

The use of random forests for GWAS is an ongoing research topic.

### Epistasis
Phenotypic variation can also be explained by non-linear interaction effects between SNPs. (In fact, the phenotype for this tutorial was simulated with a multiplicative effect between `'rs1942455` and `rs1115764`.)

Epistasis detection is implemented in SNP with the `--epistasis` flag, but can be very long.

Many other software have been developed to address this problem, some of them using parallelization on clusters or GPUs. OMICtools maintains a list  [here](https://omictools.com/epistasis-detection-category). You can also read a recent review by [Niel et al., (2015)](http://journal.frontiersin.org/article/10.3389/fgene.2015.00285/full).

Epistasis detection (and definition...) are ongoing research topics.

### What about missing values?
In this tutorial, we've worked with full genotype data. In practice, there are always some individuals for which a given SNP has not been properly genotyped and hence is missing. The `X` matrice will therefore have missing data, and while for a single-SNP association you can just ignore, when analyzing SNP x, the individuals for which SNP x does not have a genotype, most machine learning algorithms will not be happy with missing values.

The typical way to address this is to impute the missing values using, for example, [IMPUTE2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html).

### I'm sure this Python business is nice, but I use R
I think you're missing out, but here are a few useful packages:
* `snpStats` for standard GWAS, Manhattan plots, Q-Q plots, etc.;
* `ggplot` for prettier plotting;
* `GWASTools` is also nice;
* `glmnet` for regularized linear regression.

### How was this data generated exactly?
The `Microbiome practical.ipynb` [Jupyter](http://jupyter.org) notebook in this repository contains all the data generation and processing steps.
