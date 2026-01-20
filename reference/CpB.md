# Calculates the rate of accumulation of phylogenetic B-diversity (CpB) over time slices

This function estimates the rates of accumulation of phylogenetic
B-diversity (CpB) over time for inputted assemblages.

## Usage

``` r
CpB(tree, n, mat, adj, comp, method, criterion, pBO, ncor)
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

- adj:

  matrix. A square adjacency matrix containing the presence/absence
  information of all sites and their spatially adjacent ones.

- comp:

  character string. The component of the phylogenetic beta-diversity to
  obtain the rates of accumulation. It can be either "sorensen",
  "turnover", or "nestedness". Default is "sorensen".

- method:

  character string. The method for calculating the phylogenetic
  beta-diversity. It can be either obtained through a "pairwise" or
  "multisite" approach. Default is "multisite".

- criterion:

  character string. The method for cutting the tree. It can be either
  "my" (million years) or "PD" (accumulated phylogenetic diversity).
  Default is "my".

- pBO:

  numeric. A value indicating the numeric proportion to define the
  temporal origin at which the phylogenetic B-diversity (PB) started to
  accumulate in a given assemblage. Default is 5%.

- ncor:

  numeric. A value indicating the number of cores the user wants to
  parallelize. Default is 0.

## Value

The function returns a data frame containing the assemblages' rates of
cumulative phylogenetic B-diversity (CpB), their total phylogenetic
B-diversity (PB), and their PB origin (pBO).

## Details

**Parallelization**

Users are advised to check the number of available cores within their
machines before running parallel programming.

## References

See the tutorial on how to use this function on our
[website](https://araujomat.github.io/treesliceR/articles/Passeriformes-diversification.html).

## See also

Other cumulative phylogenetic index analysis:
[`CpD()`](https://araujomat.github.io/treesliceR/reference/CpD.md),
[`CpE()`](https://araujomat.github.io/treesliceR/reference/CpE.md),
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

# Create a random adjacency matrix
adj <- matrix(sample(c(1,0), 10*10, replace = TRUE), ncol = 10, nrow = 10)

# Fill the diagonals with 1
diag(adj) <- 1

# Calculate their CpB (sorensen) for 100 tree slices
CpB(tree, n = 100, mat = mat, adj = adj, comp = "sorensen", method = "multisite")
#> > The 100 number of pieces inputted equals to intervals of 0.0421617378647582 million of years.
#>         CpB        PB       pBO
#> 1  5.680293 0.1252677 0.5273905
#> 2  5.116377 0.1527932 0.5855183
#> 3  5.264234 0.1420109 0.5690727
#> 4  5.703764 0.1173966 0.5252202
#> 5  5.319323 0.1897345 0.5631793
#> 6  5.303892 0.1840631 0.5648178
#> 7  5.280286 0.1929758 0.5673428
#> 8  5.276041 0.1442105 0.5677993
#> 9  5.702609 0.1510395 0.5253266
#> 10 5.836211 0.1127854 0.5133009
```
