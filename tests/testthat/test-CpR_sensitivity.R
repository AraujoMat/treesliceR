# Testing if the output returns the correct df dimension
test_that("CpR sensitivity returns the expected data frame size for all phylogenetic indexes approaches", {

  # Creating empty lists to store CpR_sensitivity outputs
  CpR_sensitivity_M1 <- list()
  CpR_sensitivity_M2 <- list()
  CpR_sensitivity_M3 <- list()
  CpR_sensitivity_M4 <- list()

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

  # Creating a vector containing the number of slices desired to run the sensitivity
  vec <- c(25, 50, 75, 100, 125)
  # and the number of samples to evaluate it
  samp <- 3

  # Run the CpR_sensitivity algorithm while suppressing some warnings
  # (related to tips no present within the matrix, and vice-versa)
  suppressWarnings({for(i in 1:20){
    CpR_sensitivity_M1[[i]] <- CpR_sensitivity(tree[[i]], vec = vec, mat = mat[[i]], samp = samp, rate = "CpD")
    CpR_sensitivity_M2[[i]] <- CpR_sensitivity(tree[[i]], vec = vec, mat = mat[[i]], samp = samp, rate = "CpE")
    CpR_sensitivity_M3[[i]] <- CpR_sensitivity(tree[[i]], vec = vec, mat = mat[[i]], adj = adj[[i]], samp = samp, rate = "CpB")
    CpR_sensitivity_M4[[i]] <- CpR_sensitivity(tree[[i]], vec = vec, mat = mat[[i]], adj = adj[[i]], samp = samp, rate = "CpB_RW")
  }})

  # Test
  for(i in 1:20){
    expect_equal(dim(CpR_sensitivity_M1[[i]]), c(samp, length(vec)))
    expect_equal(dim(CpR_sensitivity_M2[[i]]), c(samp, length(vec)))
    expect_equal(dim(CpR_sensitivity_M3[[i]]), c(samp, length(vec)))
    expect_equal(dim(CpR_sensitivity_M4[[i]]), c(samp, length(vec)))
  }
})

