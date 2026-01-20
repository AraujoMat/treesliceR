# Slices a temporal interval within a phylogenetic tree

This function slices a temporal interval located within an ultrametric
phylogenetic tree.

## Usage

``` r
squeeze_int(tree, from, to, invert = FALSE, criterion = "my", dropNodes = FALSE)
```

## Arguments

- tree:

  phylo. An ultrametric phylogenetic tree in the "phylo" format.

- from:

  numeric. A temporal threshold value that determines the time at which
  the interval should start.

- to:

  numeric. A temporal threshold value that determines the time at which
  the interval should end.

- invert:

  logical. A logical value indicating if the desired slice must be
  executed inside (invert = FALSE) or outside (invert = TRUE) the
  defined interval. Using the argument as TRUE will return a list
  containing a root and tip slices. Default is FALSE.

- criterion:

  character string. The method for cutting the tree. It can be either
  "my" (million years) or "PD" (accumulated phylogenetic diversity).
  Default is "my".

- dropNodes:

  logical. A logical value indicating whether the nodes that were sliced
  (void nodes, presenting no branch length) should be preserved in the
  node matrix. Default is FALSE.

## Value

The function returns an time-slice interval of a phylogenetic tree in
the "phylo" format.

## Details

**Slicing approach**

To return a given phylogenetic interval, this function simultaneously
applies simultaneously the same logic as
[`squeeze_tips()`](https://araujomat.github.io/treesliceR/reference/squeeze_tips.md)
and
[`squeeze_root()`](https://araujomat.github.io/treesliceR/reference/squeeze_root.md).
If "invert = TRUE", then the temporal interval set will be excluded from
the phylogeny, returning a list containing a tip and a root slices.

## References

See the tutorial on how to use this function on our
[website](https://araujomat.github.io/treesliceR/articles/Intro-treesliceR.html).

## See also

Other slicing methods:
[`squeeze_tips()`](https://araujomat.github.io/treesliceR/reference/squeeze_tips.md),
[`squeeze_root()`](https://araujomat.github.io/treesliceR/reference/squeeze_root.md),
[`phylo_pieces()`](https://araujomat.github.io/treesliceR/reference/phylo_pieces.md),
[`prune_tips()`](https://araujomat.github.io/treesliceR/reference/prune_tips.md)

## Author

Matheus Lima de Araujo <matheusaraujolima@live.com>

## Examples

``` r
# Generate a random tree
tree <- ape::rcoal(20)

# Slice an interval from 0.5 to 0.2 million years
tree <- squeeze_int(tree, from = 0.5, to = 0.2)

# Plot it
plot(tree)

```
