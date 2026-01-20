# Calculates the rate of accumulation of phylogenetic endemism (CpE) over time slices

This function estimates the rates of accumulation of phylogenetic
endemism (CpE) over time for inputted assemblages.

## Usage

``` r
CpE(tree, n, mat, criterion = "my", pEO = 5, ncor = 0)
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

- pEO:

  numeric. A value indicating the numeric proportion to define the
  temporal origin at which the phylogenetic endemism (PE) started to
  accumulate in a given assemblage. Default is 5%.

- ncor:

  numeric. A value indicating the number of cores the user wants to
  parallelize. Default is 0.

## Value

The function returns a data frame containing the assemblages' rates of
cumulative phylogenetic endemism (CpE), their total phylogenetic
endemism (PE), and their PE origin (pEO).

## Details

**Parallelization**

Users are advised to check the number of available cores within their
machines before running parallel programming.

## See also

Other cumulative phylogenetic index analysis:
[`CpD()`](https://araujomat.github.io/treesliceR/reference/CpD.md),
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

# Calculate the CpE for 100 tree slices
CpE(tree, n = 100, mat = mat)
#> > The 100 number of pieces inputted equals to intervals of 0.0057614043538185 million of years.
#>         CpE        PE       pEO
#> 1  5.259451 0.3341360 0.5695902
#> 2  4.690477 0.4888579 0.6386840
#> 3  6.226204 0.4181528 0.4811491
#> 4  4.548108 0.2617525 0.6586766
#> 5  3.859096 0.3297783 0.7762782
#> 6  4.534381 0.4731461 0.6606707
#> 7  5.270881 0.3590355 0.5683552
#> 8  4.486627 0.4210215 0.6677026
#> 9  5.206627 0.3484977 0.5753690
#> 10 4.324126 0.3658145 0.6927948
```
