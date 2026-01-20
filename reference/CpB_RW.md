# Calculates the range weighted rate of accumulation of phylogenetic B-diversity (CpB_RW) over time slices

This function estimates the range-weighted rates of accumulation of
phylogenetic B-diversity (CpB_RW) over time for inputted assemblages.

## Usage

``` r
CpB_RW(tree, n, mat, adj, method = "multisite", criterion = "my", pBO = 5, ncor = 0)
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
  temporal origin at which the range-weighted phylogenetic B-diversity
  (PB_RW) started to accumulate in a given assemblage. Default is 5%.

- ncor:

  numeric. A value indicating the number of cores the user wants to
  parallelize. Default is 0.

## Value

The function returns a data frame containing the assemblages' rates of
cumulative range-weighted phylogenetic B-diversity (CpB_RW), their total
range-weighted phylogenetic B-diversity (PB_RW), and their origin (pBO).

## Details

**Parallelization**

Users are advised to check the number of available cores within their
machines before running parallel programming.

## References

Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T.,
González-Orozco, C. E., Knerr, N., Thornhill, A. H., & Mishler, B. D.
(2016). Range-weighted metrics of species and phylogenetic turnover can
better resolve biogeographic transition zones. Methods in Ecology and
Evolution, 7(5), 580–588. https://doi.org/10.1111/2041-210x.12513

## See also

Other cumulative phylogenetic index analysis:
[`CpD()`](https://araujomat.github.io/treesliceR/reference/CpD.md),
[`CpE()`](https://araujomat.github.io/treesliceR/reference/CpE.md),
[`CpB()`](https://araujomat.github.io/treesliceR/reference/CpB.md)

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

# Calculate their CpB range weighted for 100 tree slices
CpB_RW(tree, n = 100, mat = mat, adj = adj, method = "multisite")
#> > The 100 number of pieces inputted equals to intervals of 0.0180041734100821 million of years.
#>      CpB_RW     PB_RW       pBO
#> 1  7.484801 0.3657205 0.4002421
#> 2  8.219977 0.2992264 0.3644453
#> 3  8.403562 0.2881230 0.3564836
#> 4  7.213133 0.3746163 0.4153164
#> 5  7.139844 0.3973057 0.4195795
#> 6  6.232263 0.3283471 0.4806813
#> 7  6.931592 0.3768518 0.4321853
#> 8  6.058951 0.3196225 0.4944308
#> 9  7.586139 0.3817362 0.3948955
#> 10 9.978632 0.1737816 0.3002147
```
