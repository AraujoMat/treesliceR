# Testing if the size of the remaining slice correspond to the time inputted
test_that("Squeeze root correctly slices the trees using both criteria, my and pd", {
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
  time <- runif(20, min = 0.1, max = 0.3)
  expected_v <- list()

  # Squeezing the tree and getting the size of the squeeze
  for (i in 1:num_phylogenies) {
    # Squeezing the tree based on the sampled time
    squeeze_results[[i]] <- squeeze_root(phylos[[i]], time[[i]])
    # Capturing the new nodes_config
    squeeze_results[[i]] <- nodes_config(squeeze_results[[i]])
    # Checking the size of the slice
    expected_v[[i]] <- max(squeeze_results[[i]]$config$YearEnd)
  }

  # Test 1
  for(i in 1:length(expected_v)){
    expect_equal(expected_v[[i]], time[[i]])
  }

  ## Doing the same for squeezes based on "pd"
  # Initialize lists to store results
  squeeze_results <- vector("list", length = num_phylogenies)
  time <- runif(20, min = 0.1, max = 0.3)
  expected_v <- list()

  # Squeezing the tree and getting the size of the squeeze
  for (i in 1:num_phylogenies) {
    # Squeezing the tree based on the sampled time
    squeeze_results[[i]] <- squeeze_root(phylos[[i]], time[[i]], criteria = "pd")
    # Checking the amount of PD within removed from the slice
    expected_v[[i]] <- sum(phylos[[i]]$edge.length) - sum(squeeze_results[[i]]$edge.length)
  }

  # Testing if the pd removed correspond to the pd inputted
  # Test 2
  for(i in 1:length(squeeze_results)){
    expect_equal(expected_v[[i]], time[[i]])
  }
})




