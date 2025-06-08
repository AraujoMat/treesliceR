test_that("CpB_RW returns the expected data frame size", {
  # Number of iterations
  n <- sample(10:20, 20, replace = TRUE)

  # Creating empty lists to store CpB outputs
  CpB_RW_M1 <- list()

  # And also lists to store phylogenies
  tree <- list()
  mat <- list()
  adj <- list()

  # Create 20 random assemblages
  for(i in 1:20){
    # Create a random phylogeny
    tree[[i]] <- ape::rcoal(20, br = c(min = 2, max = 10)) # Sets a minimum and maximum length for branch lengths
    # Create a random matrix
    mat[[i]] <- matrix(sample(c(1, 0), 20*10, replace = TRUE), ncol = 20, nrow = 10)
    colnames(mat[[i]]) <- tree[[i]]$tip.label # Name its columns according to tip names
    # Create a random adjacency matrix
    adj[[i]] <- matrix(sample(c(1,0), 10*10, replace = TRUE), ncol = 10, nrow = 10)
    # Fill the diagonals with 1
    diag(adj[[i]]) <- 1
  }

  # Run the CpB_RW algorithm while suppressing some warnings
  # (related to tips no present within the matrix, and vice-versa)
  suppressWarnings({for(i in 1:20){
    CpB_RW_M1[[i]] <- CpB_RW(tree[[i]], n = n[i], mat = mat[[i]], adj = adj[[i]], criterion = "my")
  }})

  for(i in 1:20){
    expect_equal(dim(CpB_RW_M1[[i]]), c(nrow(adj[[i]]), 3))
  }
})
