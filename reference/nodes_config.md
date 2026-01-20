# Capture general branching information for an inputted phylogenetic tree

This function captures general tree information, including branching,
node positions, and depths. It serves as a core function underlying all
algorithms for slicing phylogenies.

## Usage

``` r
nodes_config(tree)
```

## Arguments

- tree:

  phylo. An ultrametric phylogenetic tree in the "phylo" format.

## Value

The function returns a phylogenetic tree in the "phylo" format
containing three novel pieces of information (stored within "config",
"node_matrix", and "tree_depth").

## Details

This function captures node and edge information from an ultrametric
phylogenetic tree. The function provides a data frame containing
detailed branching information (within the internal "config" object), a
node matrix (within "node_matrix"), and the tree age (within
"tree_depth").

More specifically, the "config" object returns the following
information: NodeBegin: the node at which a given branch begins.
NodeEnd: the node at which a given branch ends. NodeLength: the branch
length of that nodes interval. YearBegin: the year at which a given node
begins. YearEnd: the year at which a given node ends.

## See also

Phylogenetic slicing methods:
[`squeeze_tips()`](https://araujomat.github.io/treesliceR/reference/squeeze_tips.md),[`squeeze_root()`](https://araujomat.github.io/treesliceR/reference/squeeze_root.md),[`squeeze_int()`](https://araujomat.github.io/treesliceR/reference/squeeze_int.md).

## Author

Matheus Lima de Araujo <matheusaraujolima@live.com>

## Examples

``` r
# Generate a random tree
tree <- ape::rcoal(20)

# Capture tree information
tree <- nodes_config(tree)

# Accessing these informations
tree$config # Nodes configurations
#>    NodeBegin NodeEnd  NodeLength YearBegin   YearEnd
#> 1         21       1 0.887949971 0.0000000 0.8879500
#> 22        21      22 0.116482812 0.0000000 0.1164828
#> 24        22      24 0.538438888 0.1164828 0.6549217
#> 27        24      27 0.098054404 0.6549217 0.7529761
#> 38        27      38 0.122722991 0.7529761 0.8756991
#> 2         38       2 0.012250875 0.8756991 0.8879500
#> 3         38       3 0.012250875 0.8756991 0.8879500
#> 4         27       4 0.134973866 0.7529761 0.8879500
#> 32        24      32 0.194947587 0.6549217 0.8498693
#> 5         32       5 0.038080683 0.8498693 0.8879500
#> 6         32       6 0.038080683 0.8498693 0.8879500
#> 23        22      23 0.368522309 0.1164828 0.4850051
#> 26        23      26 0.244753006 0.4850051 0.7297581
#> 7         26       7 0.158191844 0.7297581 0.8879500
#> 28        26      28 0.024654628 0.7297581 0.7544128
#> 33        28      33 0.102416707 0.7544128 0.8568295
#> 8         33       8 0.031120508 0.8568295 0.8879500
#> 9         33       9 0.031120508 0.8568295 0.8879500
#> 10        28      10 0.133537215 0.7544128 0.8879500
#> 25        23      25 0.221290601 0.4850051 0.7062957
#> 29        25      29 0.050021042 0.7062957 0.7563168
#> 31        29      31 0.087857646 0.7563168 0.8441744
#> 11        31      11 0.043775560 0.8441744 0.8879500
#> 34        31      34 0.015288546 0.8441744 0.8594630
#> 12        34      12 0.028487014 0.8594630 0.8879500
#> 13        34      13 0.028487014 0.8594630 0.8879500
#> 30        29      30 0.065628883 0.7563168 0.8219456
#> 14        30      14 0.066004323 0.8219456 0.8879500
#> 35        30      35 0.037932635 0.8219456 0.8598783
#> 15        35      15 0.028071688 0.8598783 0.8879500
#> 36        35      36 0.001048102 0.8598783 0.8609264
#> 16        36      16 0.027023586 0.8609264 0.8879500
#> 37        36      37 0.010240320 0.8609264 0.8711667
#> 17        37      17 0.016783266 0.8711667 0.8879500
#> 39        37      39 0.007903300 0.8711667 0.8790700
#> 18        39      18 0.008879966 0.8790700 0.8879500
#> 19        39      19 0.008879966 0.8790700 0.8879500
#> 20        25      20 0.181654248 0.7062957 0.8879500
tree$node_matrix # Node matrix
#> NULL
tree$tree_depth # Tree age
#> [1] 0.88795
```
