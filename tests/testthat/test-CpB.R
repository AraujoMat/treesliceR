test_that("CpB returns the expected data frame size using any beta component", {

  # Number of iterations
  n <- sample(10:20, 20, replace = TRUE)

  # Creating empty lists to store CpB outputs
  CpB_M1 <- list()
  CpB_M2 <- list()
  CpB_M3 <- list()

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

  # Run the CpB algorithm while suppressing some warnings
  # (related to tips no present within the matrix, and vice-versa)
  suppressWarnings({for(i in 1:20){
    CpB_M1[[i]] <- CpB(tree[[i]], n = n[i], asb = asb[[i]], comp = "sorensen", criteria = "my")
    CpB_M2[[i]] <- CpB(tree[[i]], n = n[i], asb = asb[[i]], comp = "turnover", criteria = "my")
    CpB_M3[[i]] <- CpB(tree[[i]], n = n[i], asb = asb[[i]], comp = "nestedness", criteria = "my")
  }})

  # Testing
  for(i in 1:20){
    expect_equal(dim(CpB_M1[[i]]), c(length(asb[[i]]), 3))
    expect_equal(dim(CpB_M2[[i]]), c(length(asb[[i]]), 3))
    expect_equal(dim(CpB_M3[[i]]), c(length(asb[[i]]), 3))
  }
})
