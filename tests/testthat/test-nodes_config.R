# Testing the outputs
test_that("Nodes config returns information for all edge and nodes in the tree",{
  # Function to crate random phylogenies with a random number of tips
  generate_phylogenies <- function(num_trees) {
    # Initialize an empty list
    phylo_list <- vector("list", length = num_trees)
    # Create the trees and save them
    for (i in 1:num_trees) {
      size <- sample(15:100, 1)  # between 15 and 100 tips
      phylo <- ape::rcoal(size)
      phylo_list[[i]] <- phylo
    }
    # Return the final list
    return(phylo_list)
  }

  # Generate a list of 20 random phylogenies
  num_phylogenies <- 20
  phylos <- generate_phylogenies(num_phylogenies)

  # Adding nodes_config
  for(i in 1:length(phylos)){
    phylos[[i]] <- nodes_config(phylos[[i]])
  }

  for (i in 1:num_phylogenies) {
    expected_length <- length(phylos[[i]]$edge.length)
    expect_equal(expected_length, length(phylos[[i]]$config$NodeBegin))
    expect_equal(expected_length, length(phylos[[i]]$config$NodeEnd))
    expect_equal(expected_length, length(phylos[[i]]$config$NodeLength))
    expect_equal(expected_length, length(phylos[[i]]$config$YearBegin))
    expect_equal(expected_length, length(phylos[[i]]$config$YearEnd))
  }
})
