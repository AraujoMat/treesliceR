# Evaluating the size of the outputs
test_that("r_phylo outputs a list with the same size as the inputted number of assemblages, and the correct inputted number of slices", {
  # Number of iterations
  n <- sample(10:20, 20, replace = TRUE)

  # Creating a presence-absence matrix
  mat_input <- c(0,1,0,0,0,0,1,1,0,1,1,0,1,1,0,1,0,1,0,1,1,0,1,0,1,0,0,1,0,1,
                 1,1,1,1,1,1,1,0,1,1,0,1,1,1,0,1,1,0,1,0,0,0,0,1,1,1,0,1,1,1,
                 1,0,1,1,1,0,0,0,0,1,0,0,1,0,1,0,0,1,1,0,0,1,1,0,0,0,0,1,0,0,
                 0,0,1,0,0,0,1,0,1,1,1,1,1,0,1,0,0,0,1,0,0,1,1,0,1,0,0,1,0,0,
                 1,1,1,1,0,1,1,1,0,1,1,0,1,0,1,0,0,1,1,0,1,1,0,0,1,1,0,0,1,0,
                 1,1,0,1,1,0,0,0,0,0,1,1,1,1,1,0,0,1,0,1,0,0,1,0,1,0,0,1,0,1,
                 0,1,1,0,1,0,1,0,1,0,1,1,1,1,1,1,0,1,1,1,0,0,1,0,0,0,0,0,1,1,
                 1,1,0,1,1,0,0,1,1,1,1,1,1,0,1,0,0,1,1,0,1,1,1,0,1,1,1,0,0,0,
                 0,1,1,0,1,0,1,0,1,1,0,1,1,0,0,1,0,0,0,1,1,1,0,1,1,1,0,1,0,0,
                 1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,1,0,1)
  # Transforming into a 10x30 matrix
  mat_input <- matrix(mat_input, nrow = 10, ncol = 30, byrow = TRUE)

  # Creating an adjacency matrix
  adj <- c(1,1,1,0,0,1,0,1,0,0,
           1,1,0,0,0,0,1,1,0,1,
           1,1,1,0,1,0,1,1,1,1,
           1,0,1,1,0,0,0,0,1,0,
           1,1,0,1,1,1,1,1,0,1,
           0,1,0,1,0,1,1,1,0,1,
           0,1,1,0,0,0,1,1,0,1,
           1,0,0,0,0,1,0,1,1,1,
           1,0,0,1,0,1,0,0,1,0,
           0,1,1,1,0,1,0,1,0,1)
  # Converting into a 10x10 matrix
  adj <- matrix(adj, nrow = 10, ncol = 10, byrow = TRUE)

  # Creating empty lists to store r_phylo outputs
  r_phylo_M1 <- list()
  r_phylo_M2 <- list()
  r_phylo_M3 <- list()
  r_phylo_M4 <- list()

  # And also lists to store phylogenies
  tree <- list()
  mat <- list()

  # Create 20 random assemblages
  for(i in 1:20){
    # Create a random phylogeny
    tree[[i]] <- ape::rcoal(30, br = c(min = 2, max = 10)) # Sets a minimum and maximum length for branch lengths
    # Create a random matrix
    mat[[i]] <- mat_input
    colnames(mat[[i]]) <- tree[[i]]$tip.label # Name its columns according to tip names
  }

  # Run the algorithm for different phylogentic indexes while suppressing some warnings
  # (related to tips no present within the matrix, and vice-versa)
  suppressWarnings({for(i in 1:20){
    r_phylo_M1[[i]] <- r_phylo(tree[[i]], n = n[i], mat = mat[[i]], index = "PD")
    r_phylo_M2[[i]] <- r_phylo(tree[[i]], n = n[i], mat = mat[[i]], index = "PE")
    r_phylo_M3[[i]] <- r_phylo(tree[[i]], n = n[i], mat = mat[[i]], adj = adj, index = "PB")
    r_phylo_M4[[i]] <- r_phylo(tree[[i]], n = n[i], mat = mat[[i]], adj = adj, index = "PB_RW")
  }})

  # Test 1
  for(i in 1:20){
    expect_equal(length(r_phylo_M1[[i]]), nrow(mat[[i]]))
    expect_equal(length(r_phylo_M2[[i]]), nrow(mat[[i]]))
    expect_equal(length(r_phylo_M3[[i]]), nrow(mat[[i]]))
    expect_equal(length(r_phylo_M4[[i]]), nrow(mat[[i]]))
  }

  # test 2
  for(i in 1:20){
    for(j in 1:length(r_phylo_M1[[1]])){  # i <- 1   j <- 1
      expect_equal(length(r_phylo_M1[[i]][[j]]), n[i])
      expect_equal(length(r_phylo_M2[[i]][[j]]), n[i])
      expect_equal(length(r_phylo_M3[[i]][[j]]), n[i])
      expect_equal(length(r_phylo_M4[[i]][[j]]), n[i])
    }
  }
})
