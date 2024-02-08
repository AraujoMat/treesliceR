# Evaluating the size of the outputs
test_that("r_phylo outputs a list with the same size as the inputted number of assemblages, and the correct inputted number of slices", {
  # Number of iterations
  n <- sample(5:15, 20, replace = TRUE)

  # Creating empty lists to store r_phylo outputs
  r_phylo_M1 <- list()
  r_phylo_M2 <- list()
  r_phylo_M3 <- list()
  r_phylo_M4 <- list()

  # And also lists to store phylogenies
  tree <- list()
  mat <- list()
  adj <- list()

  # Create 20 random assemblages
  for(i in 1:20){
    # Create a random phylogeny
    tree[[i]] <- ape::rcoal(30)
    # Create a random matrix
    mat[[i]] <- matrix(sample(c(1, 0), 30*15, replace = TRUE), ncol = 30, nrow = 15)
    colnames(mat[[i]]) <- tree[[i]]$tip.label # Name its columns according to tip names
    # Create a random adjacency matrix
    adj[[i]] <- matrix(sample(c(1,0), 15*15, replace = TRUE), ncol = 15, nrow = 15)
    # Fill the diagonals with 1
    diag(adj[[i]]) <- 1
  }

  # Run the algorithm for different phylogentic indexes while suppressing some warnings
  # (related to tips no present within the matrix, and vice-versa)
  suppressWarnings({for(i in 1:20){
    r_phylo_M1[[i]] <- r_phylo(tree[[i]], n = n[i], mat = mat[[i]], index = "PD")
    r_phylo_M2[[i]] <- r_phylo(tree[[i]], n = n[i], mat = mat[[i]], index = "PE")
    r_phylo_M3[[i]] <- r_phylo(tree[[i]], n = n[i], mat = mat[[i]], adj = adj[[i]], index = "PB")
    r_phylo_M4[[i]] <- r_phylo(tree[[i]], n = n[i], mat = mat[[i]], adj = adj[[i]], index = "PB_RW")
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
