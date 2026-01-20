# Slices a phylogenetic tree following a 'tipward' orientation

This function slices a phylogenetic tree in a 'tipward' orientation,
starting from the tips and moving towards the root of the tree.

## Usage

``` r
squeeze_tips(tree, time, criterion = "my", dropNodes = FALSE)
```

## Arguments

- tree:

  phylo. An ultrametric phylogenetic tree in the "phylo" format.

- time:

  numeric. A value that determines the time, or accumulated phylogenetic
  diversity (PD), at which the tree should be cut.

- criterion:

  character string. The method for cutting the tree. It can be either
  "my" (million years) or "PD" (accumulated phylogenetic diversity).
  Default is "my".

- dropNodes:

  logical. A logical value indicating whether the nodes that were sliced
  (void nodes, presenting no branch length) should be preserved in the
  node matrix. Default is FALSE.

## Value

The function returns a time-slice of an inputted phylogenetic tree in
the "phylo" format, following a 'tipward' orientation.

## Details

**Slicing approach** The treesliceR package uses a simple approach for
cutting phylogenies, which reduces branch lengths in relation to an
inputted temporal threshold.

## References

See the tutorial on how to use this function on our
[website](https://araujomat.github.io/treesliceR/articles/Intro-treesliceR.html).

## See also

Other slicing methods:
[`squeeze_root()`](https://araujomat.github.io/treesliceR/reference/squeeze_root.md),
[`squeeze_int()`](https://araujomat.github.io/treesliceR/reference/squeeze_int.md),
[`phylo_pieces()`](https://araujomat.github.io/treesliceR/reference/phylo_pieces.md),
[`prune_tips()`](https://araujomat.github.io/treesliceR/reference/prune_tips.md).

## Author

Matheus Lima de Araujo <matheusaraujolima@live.com>

## Examples

``` r
# Generate a random tree
tree <- ape::rcoal(20)

# Slice "tipwardly" the phylogeny at 0.1 million years
tree <- squeeze_tips(tree, time = 0.1)

# Plot it
plot(tree)

```
