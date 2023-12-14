# Checking if the interval size of the piece equals the difference between "from" and "to"
test_that("Squeeze int correctly slices the trees using criteria = my", {
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


  ## Doing the same for squeezes based on "my"
  # Initialize lists to store results
  squeeze_results <- vector("list", length = num_phylogenies)
  from <- runif(20, min = 0.21, max = 0.3)
  to <- runif(20, min = 0.1, max = 0.2)
  expected_v <- list()

  # Squeezing the tree and getting the size of the squeeze
  for (i in 1:num_phylogenies) {
    # Capturing an interval within the phylogeny
    squeeze_results[[i]] <- squeeze_int(phylos[[i]], from = from[[i]], to = to[[i]])
    # Capturing the new nodes config for each tree
    squeeze_results[[i]] <- nodes_config(squeeze_results[[i]])
    # Saving the interval size
    expected_v[[i]] <- max(squeeze_results[[i]]$config$YearEnd)
  }

  # Test
  for(i in 1:length(expected_v)){
    expect_equal(expected_v[[i]], (from[[i]] - to[[i]]))
  }
})
