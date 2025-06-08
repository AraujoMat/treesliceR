test_that("CpB returns the expected data frame size using any beta component", {

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

  # Creating empty lists to store CpB outputs
  CpB_M1 <- list()
  CpB_M2 <- list()

  # And also lists to store phylogenies
  tree <- list()
  mat <- list()

  # Create 20 random assemblages
  for(i in 1:20){
    # Create a random phylogeny
    tree[[i]] <- ape::rcoal(30, br = c(min = 3, max = 10)) # Sets a minimum and maximum length for branch lengths
    # Create a random matrix
    mat[[i]] <- mat_input
    colnames(mat[[i]]) <- tree[[i]]$tip.label # Name its columns according to tip names
  }

  # Run the CpB algorithm while suppressing some warnings
  # (related to tips no present within the matrix, and vice-versa)
  suppressWarnings({for(i in 1:20){
    CpB_M1[[i]] <- CpB(tree[[i]], n = n[i], mat = mat[[i]], adj = adj, comp = "sorensen", criterion = "my")
    CpB_M2[[i]] <- CpB(tree[[i]], n = n[i], mat = mat[[i]], adj = adj, comp = "turnover", criterion = "my")
  }})

  # Testing
  for(i in 1:20){
    expect_equal(dim(CpB_M1[[i]]), c(nrow(adj), 3))
    expect_equal(dim(CpB_M2[[i]]), c(nrow(adj), 3))
  }
})
