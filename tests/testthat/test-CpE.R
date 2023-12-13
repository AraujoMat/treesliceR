test_that("CpE returns the expected data frame size", {

  # Number of iterations
  n <- sample(10:20, 20, replace = TRUE)

  # Creating empty list to store CpE outputs
  CpE_M1 <- list()

  # And also lists to store phylogenies
  tree <- list()
  mat <- list()
?rcoal
  # Create 20 random assemblages
  for(i in 1:20){
    # Create a random phylogeny
    tree[[i]] <- ape::rcoal(20)
    # Create a random matrix
    mat[[i]] <- matrix(sample(c(1, 0), 20 * 10, replace = TRUE), ncol = 20, nrow = 10)
    colnames(mat[[i]]) <- tree[[i]]$tip.label # Name its columns according to tip names
  }

  # Run the CpE algorithm while suppressing some warnings
  # (related to tips no present within the matrix, and vice-versa)
  suppressWarnings({for(i in 1:20){
    CpE_M1[[i]] <- CpE(tree[[i]], n = n[i], mat = mat[[i]], criteria = "my")
  }})

  # Test
  for(i in 1:20){
    expect_equal(dim(CpE_M1[[i]]), c(nrow(mat[[i]]), 3))
  }
})
