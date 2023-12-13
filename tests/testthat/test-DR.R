# Testing the output
test_that("DR returns the DR for all tips within the tree", {

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

  # Creating an empty output list
  DRvals <- list()

  # Calculating the DRs
  for(i in 1:length(phylos)){
    DRvals[[i]] <- DR(phylos[[i]])
  }

  # Test
  for(i in 1:20){
    expect_equal(length(phylos[[i]]$tip.label), nrow(DRvals[[i]]))
    }
  })
