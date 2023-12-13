# Pruning the tips in a temporal threshold larger than supported by the phylogeny
test_that("Prune tips identify larger inputted thresholds than is supported by the phylogeny", {
  # Create a phylogentic tree
  phylo <- ape::rcoal(100)
  expect_error(prune_tips(phylo, time = 10^6))
})

# Using a quantile higher than its available by the phylogeny
test_that("Prune tips identify larger quantile thresholds than is supported by the phylogeny", {
  expect_error(prune_tips(phylo, time = 10, qtl = TRUE))
})
