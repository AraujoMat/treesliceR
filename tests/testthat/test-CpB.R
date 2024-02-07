test_that("CpB returns the expected data frame size using any beta component", {

  # Number of iterations
  n <- sample(5:15, 20, replace = TRUE)

  # Creating empty lists to store CpB outputs
  CpB_M1 <- list()
  CpB_M2 <- list()
  CpB_M3 <- list()

  # And also lists to store phylogenies, species matrix, and adjacency matrix
  tree <- list()
  mat <- list()
  adj <- list()

  # Create 20 random assemblages
  for(i in 1:20){
    # Create a random phylogeny
    tree[[i]] <- ape::rcoal(20)
    # Create a random matrix
    mat[[i]] <- matrix(sample(c(1, 0), 20*10, replace = TRUE), ncol = 20, nrow = 10)
    colnames(mat[[i]]) <- tree[[i]]$tip.label # Name its columns according to tip names
    # Create a random adjacency matrix
    adj[[i]] <- matrix(sample(c(1,0), 10*10, replace = TRUE), ncol = 10, nrow = 10)
    # Fill the diagonals with 1
    diag(adj[[i]]) <- 1
  }

  # Run the CpB algorithm while suppressing some warnings
  # (related to tips no present within the matrix, and vice-versa)
  suppressWarnings({for(i in 1:20){
    CpB_M1[[i]] <- CpB(tree[[i]], n = n[i], mat = mat[[i]], adj = adj[[i]], comp = "sorensen", criterion = "my")
    CpB_M2[[i]] <- CpB(tree[[i]], n = n[i], mat = mat[[i]], adj = adj[[i]], comp = "turnover", criterion = "my")
    CpB_M3[[i]] <- CpB(tree[[i]], n = n[i], mat = mat[[i]], adj = adj[[i]], comp = "nestedness", criterion = "my")
  }})

  # Testing
  for(i in 1:20){
    expect_equal(dim(CpB_M1[[i]]), c(nrow(adj[[i]]), 3))
    expect_equal(dim(CpB_M2[[i]]), c(nrow(adj[[i]]), 3))
    expect_equal(dim(CpB_M3[[i]]), c(nrow(adj[[i]]), 3))
  }
})
