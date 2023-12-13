test_that("CpB_RW returns the expected data frame size", {
  # Number of iterations
  n <- sample(10:20, 20, replace = TRUE)

  # Creating empty lists to store CpB outputs
  CpB_RW_M1 <- list()

  # And also lists to store phylogenies
  tree <- list()
  mat <- list()
  asb <- list()

  # Create 20 random assemblages
  for(i in 1:20){
    # Create a random phylogeny
    tree[[i]] <- ape::rcoal(20)
    # Create a random matrix
    mat[[i]] <- matrix(sample(c(1, 0), 20 * 10, replace = TRUE), ncol = 20, nrow = 10)
    colnames(mat[[i]]) <- tree[[i]]$tip.label # Name its columns according to tip names
    # Create an assemblage with its neighborhoods
    asb[[i]] <- list(mat[[i]][1:5, ], mat[[i]][6:10, ])
  }

  # Run the CpB_RW algorithm while suppressing some warnings
  # (related to tips no present within the matrix, and vice-versa)
  suppressWarnings({for(i in 1:20){
    CpB_RW_M1[[i]] <- CpB_RW(tree[[i]], n = n[i], mat = mat[[i]], asb = asb[[i]], criteria = "my")
  }})

  for(i in 1:20){
    expect_equal(dim(CpB_RW_M1[[i]]), c(length(asb[[i]]), 3))
  }
})
