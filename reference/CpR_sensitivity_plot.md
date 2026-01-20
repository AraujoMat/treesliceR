# Make a plot of the sensitivity analysis evaluated under cumulative phylogenetic rate functions

This function generates a plot of the sensitivity analysis calculated
for a given cumulative phylogenetic rate (CpR) assessed through the
[`CpR_sensitivity()`](https://araujomat.github.io/treesliceR/reference/CpR_sensitivity.md)
function.

## Usage

``` r
CpR_sensitivity_plot(sst_output, rate = NULL, stc = "mean")
```

## Arguments

- sst_output:

  data frame. The outputted data frame from the
  [`CpR_sensitivity()`](https://araujomat.github.io/treesliceR/reference/CpR_sensitivity.md)
  function.

- rate:

  character string. The desired phylogenetic index rate to display. It
  can be filled with "CpD", "CPE", "CpB", or "CpB_RW".

- stc:

  character string. A statistical measure to summarize the phylogenetic
  rates and create the plot, which could be filled with "mean", "var",
  "median", "sd", "min", and "max". Default is "mean".

## Value

The function returns a ggplot graph.

## References

Users can use the `CpR_sensitivity_plot()` function for plotting
sensitivity analysis outputs.

## See also

CpR sensitivity analysis:
[`CpR_sensitivity()`](https://araujomat.github.io/treesliceR/reference/CpR_sensitivity.md).

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
#> > The 100 number of pieces inputted equals to intervals of 0.00821168057307371 million of years.
#>         CpD       PD       pDO
#> 1  3.223311 3.500258 0.9293960
#> 2  3.121716 2.971389 0.9596427
#> 3  3.268035 3.085609 0.9166769
#> 4  2.957965 3.172645 1.0127680
#> 5  3.152449 3.262462 0.9502874
#> 6  2.896482 2.898836 1.0342659
#> 7  3.158758 3.423315 0.9483892
#> 8  2.657978 2.310989 1.1270717
#> 9  3.444037 3.825147 0.8698317
#> 10 2.942772 3.105235 1.0179968

# Create a vector of number of slices
vec <- c(25, 50, 75, 100, 125, 150)

# Calculate the sensitivity of the CpD
Sens_out <- CpR_sensitivity(tree, vec, mat, rate = "CpD", samp = 5)
#> Warning: Removing the species in presence abscence matrix without any occurrence
#> Warning: Removing tips from phylogeny that are absent on species matrix
#> > The 25 number of pieces inputted equals to intervals of 0.0328467222922948 million of years.
#> Warning: Removing the species in presence abscence matrix without any occurrence
#> Warning: Removing tips from phylogeny that are absent on species matrix
#> > The 50 number of pieces inputted equals to intervals of 0.0164233611461474 million of years.
#> Warning: Removing the species in presence abscence matrix without any occurrence
#> Warning: Removing tips from phylogeny that are absent on species matrix
#> > The 75 number of pieces inputted equals to intervals of 0.0109489074307649 million of years.
#> Warning: Removing the species in presence abscence matrix without any occurrence
#> Warning: Removing tips from phylogeny that are absent on species matrix
#> > The 100 number of pieces inputted equals to intervals of 0.00821168057307371 million of years.
#> Warning: Removing the species in presence abscence matrix without any occurrence
#> Warning: Removing tips from phylogeny that are absent on species matrix
#> > The 125 number of pieces inputted equals to intervals of 0.00656934445845896 million of years.
#> Warning: Removing the species in presence abscence matrix without any occurrence
#> Warning: Removing tips from phylogeny that are absent on species matrix
#> > The 150 number of pieces inputted equals to intervals of 0.00547445371538247 million of years.

# Plot the sensitity analysis
CpR_sensitivity_plot(Sens_out, rate = "CpD", stc = "mean")

```
