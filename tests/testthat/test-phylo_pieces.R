# Testing for the size of the pieces
test_that("Phylo piece correctly create the inputted number of slicescreated slices of equal temporal width", {

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

  # Initialize lists to store results
  slice_results <- vector("list", length = num_phylogenies)
  number_slices <- sample(2:50, 20)

  # Squeezing the tree and getting the size of the squeeze
  for (i in 1:num_phylogenies) {
    # Slicing all the phylogenies into multiple pieces
    slice_results[[i]] <- phylo_pieces(phylos[[i]], n = number_slices[[i]])
  }

  # Test 1
  # Checking if the number of piece equals matches the input
  for(i in 1:length(slice_results)){
      expect_equal(length(slice_results[[i]]), number_slices[[i]])
    }

  # Checking if a single tree have all it's slices cut using equal temporal widths

  # Select one tree as example
  phylo <- slice_results[[1]]
  # Empty list
  expected_v <- NA
  # Capturing the size of the pieces
  for (i in 1:length(phylo)) {
    # Capturing the new nodes config for each tree
    phylo[[i]] <- nodes_config(phylo[[i]])
    # Saving the interval size
    expected_v[i] <- max(phylo[[i]]$config$YearEnd)
  }

  # Test2
  for(i in 1:length(expected_v)){
    expect_equal(all.equal(expected_v[1], expected_v[i]), TRUE)
  }
})
