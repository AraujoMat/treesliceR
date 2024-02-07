
test_that("CpD returns the expected data frame size", {

  # Number of iterations
  n <- sample(10:20, 20, replace = TRUE)

  # Creating empty lists to store CpD outputs
  CpD_M1 <- list()

  # And also lists to store phylogenies
  tree <- list()
  mat <- list()

  # Create 20 random assemblages
  for(i in 1:20){
    # Create a random phylogeny
    tree[[i]] <- ape::rcoal(20)
    # Create a random matrix
    mat[[i]] <- matrix(sample(c(1, 0), 20 * 10, replace = TRUE), ncol = 20, nrow = 10)
    colnames(mat[[i]]) <- tree[[i]]$tip.label # Name its columns according to tip names
  }

  # Run the CpD algorithm while suppressing some warnings
  # (related to tips no present within the matrix, and vice-versa)
  suppressWarnings({for(i in 1:20){
    CpD_M1[[i]] <- CpD(tree[[i]], n = n[i], mat = mat[[i]], criterion = "my")
  }})

  # Test
  for(i in 1:20){
    expect_equal(dim(CpD_M1[[i]]), c(nrow(mat[[i]]), 3))
  }
})
