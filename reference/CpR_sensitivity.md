# Runs a sensitivity analysis for rates of accumulation of a given phylogenetic index

This function allows the evaluation of the sensitivity of the estimated
rates of accumulation of a given phylogenetic index (e.g.,
[`CpD()`](https://araujomat.github.io/treesliceR/reference/CpD.md),
[`CpE()`](https://araujomat.github.io/treesliceR/reference/CpE.md),
[`CpB()`](https://araujomat.github.io/treesliceR/reference/CpB.md),
[`CpB_RW()`](https://araujomat.github.io/treesliceR/reference/CpB_RW.md))
to the number of slices inputted by the user.

## Usage

``` r
CpR_sensitivity(tree, vec, mat, adj, rate, samp, comp, method, criterion, ncor)
```

## Arguments

- tree:

  phylo. An ultrametric phylogenetic tree in the "phylo" format.

- vec:

  numeric vector. A numeric vector containing a series of numbers of
  slices.

- mat:

  matrix. A presence/absence matrix containing all studied species and
  sites.

- adj:

  matrix. A square adjacency matrix containing the presence/absence
  information of all sites and their spatially adjacent ones.

- rate:

  character string. The desired cumulative phylogenetic rate to be
  assessed, which can be the phylogenetic diversity (CpD), phylogenetic
  endemism (CpE), phylogenetic B-diversity (CpB), or phylogenetic
  B-diversity range-weighted (CpB_RW). Default is NULL, but must be
  filled with "CpD", "CPE", "CpB_RW", or "CpB".

- samp:

  numeric. The number of assemblages, or sites, to be sampled to make
  the sensitivity analysis.

- comp:

  character string. The component of beta-diversity that the user wants
  to calculate the CpB. It can be either "sorensen", turnover" or
  "nestedness". This argument works only when "rate = CpB". Default is
  "sorensen".

- method:

  character string. The method for calculating the CpB-rate. It can be
  either "pairwise" or "multisite". This argument works only when the
  argument "rate" is set to run for "CpB" or "CpB_RW". Default is
  "multisite".

- criterion:

  character string. The method for cutting the tree. It can be either
  "my" (million years) or "PD" (accumulated phylogenetic diversity).
  Default is "my".

- ncor:

  numeric. A value indicating the number of cores the user wants to
  parallelize. Default is 0.

## Value

This function returns a data frame containing the sensitivity analysis
for a given rate of accumulation of a phylogenetic index. This outputted
data frame contains, for each row or assemblage, a column with the rate
value assessed for each inputted number of slices.

## Details

**Parallelization**

Users are advised to check the number of available cores within their
machines before running parallel programming.

**Plotting**

For plotting the sensitivity analysis output users can use
[`CpR_sensitivity_plot()`](https://araujomat.github.io/treesliceR/reference/CpR_sensitivity_plot.md).

## References

See the tutorial on how to use this function on our
[website](https://araujomat.github.io/treesliceR/articles/Passeriformes-diversification.html).

## See also

Other cumulative phylogenetic index rate analysis:
[`CpD()`](https://araujomat.github.io/treesliceR/reference/CpD.md),
[`CpE()`](https://araujomat.github.io/treesliceR/reference/CpE.md),
[`CpB()`](https://araujomat.github.io/treesliceR/reference/CpB.md),
[`CpB_RW()`](https://araujomat.github.io/treesliceR/reference/CpB_RW.md)

## Author

Matheus Lima de Araujo <matheusaraujolima@live.com>

## Examples

``` r
# Generate a random tree
tree <- ape::rcoal(20)

# Create a presence-absence matrix
mat <- matrix(sample(c(1,0), 20*10, replace = TRUE), ncol = 20, nrow = 10)
colnames(mat) <- tree$tip.label

# Calculate the CpD for 100 tree slices
CpD(tree, n = 100, mat = mat)
#> > The 100 number of pieces inputted equals to intervals of 0.0223240072262284 million of years.
#>         CpD       PD      pDO
#> 1  1.161595 7.742873 2.578982
#> 2  1.040543 6.839359 2.879008
#> 3  1.175988 7.836881 2.547418
#> 4  1.130367 7.545205 2.650229
#> 5  1.008798 6.337370 2.969607
#> 6  1.103703 7.355982 2.714256
#> 7  1.026489 6.723270 2.918426
#> 8  1.014933 6.547156 2.951657
#> 9  1.032840 6.761731 2.900481
#> 10 1.118490 7.448509 2.678372

# Create a vector of number of slices
vec <- c(25, 50, 75, 100, 125, 150)

# Calculate the sensitivity of the CpD
CpR_sensitivity(tree, vec, mat, rate = "CpD", samp = 5)
#> > The 25 number of pieces inputted equals to intervals of 0.0892960289049137 million of years.
#> > The 50 number of pieces inputted equals to intervals of 0.0446480144524569 million of years.
#> > The 75 number of pieces inputted equals to intervals of 0.0297653429683046 million of years.
#> > The 100 number of pieces inputted equals to intervals of 0.0223240072262284 million of years.
#> > The 125 number of pieces inputted equals to intervals of 0.0178592057809827 million of years.
#> > The 150 number of pieces inputted equals to intervals of 0.0148826714841523 million of years.
#>          25        50        75      100      125      150
#> 1 1.0156832 1.0733539 1.0934981 1.103703 1.109884 1.114029
#> 2 0.9599567 1.0128818 1.0312446 1.040543 1.046180 1.049955
#> 3 1.0389374 1.0988298 1.1197467 1.130367 1.136799 1.141109
#> 4 1.0289976 1.0876419 1.1081001 1.118490 1.124766 1.128979
#> 5 0.9323381 0.9825602 0.9999689 1.008798 1.014129 1.017700
```
