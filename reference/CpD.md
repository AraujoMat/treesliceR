# Calculates the rate of accumulation of phylogenetic diversity (CpD) over time slices

This function estimates the rates of accumulation of phylogenetic
diveristy (CpD) over time for inputted assemblages.

## Usage

``` r
CpD(tree, n, mat, criterion = "my", pDO = 5, ncor = 0)
```

## Arguments

- tree:

  phylo. An ultrametric phylogenetic tree in the "phylo" format.

- n:

  numeric. A numeric value indicating the number of temporal slices
  (method = 1) or the time interval in million years (or phylogenetic
  diversity) among the tree slices (method = 2). Default is 1.

- mat:

  matrix. A presence/absence matrix containing all studied species and
  sites.

- criterion:

  character string. The method for cutting the tree. It can be either
  "my" (million years) or "PD" (accumulated phylogenetic diversity).
  Default is "my".

- pDO:

  numeric. A value indicating the numeric proportion to define the
  temporal origin at which the phylogenetic diversity (PD) started to
  accumulate in a given assemblage. Default is 5%.

- ncor:

  numeric. A value indicating the number of cores the user wants to
  parallelize. Default is 0.

## Value

The function returns a data frame containing the assemblages' rates of
cumulative phylogenetic diversity (CpD), their total phylogenetic
diversity (PD), and their PD origin (pDO).

## Details

**Parallelization**

Users are advised to check the number of available cores within their
machines before running parallel programming.

## See also

Other cumulative phylogenetic rates analysis:
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
#> > The 100 number of pieces inputted equals to intervals of 0.025532819250821 million of years.
#>          CpD       PD      pDO
#> 1  0.8702400 7.891730 3.442421
#> 2  0.8197361 6.524105 3.654508
#> 3  0.8637876 7.826490 3.468135
#> 4  0.9577769 8.663085 3.127798
#> 5  0.9875099 9.186935 3.033622
#> 6  0.8235792 7.490258 3.637455
#> 7  0.8313366 6.608806 3.603513
#> 8  0.8421551 6.699907 3.557221
#> 9  0.7870920 6.283885 3.806076
#> 10 0.8422718 7.652187 3.556729
```
