
# Testing if all the saved outputs are a ggplot object
test_that("CpR_sensitivity_plot returns the correct class of outputs", {

  # Creating empty lists to store CpR_sensitivity outputs
  CpR_plot <- list()

  # And also lists to store phylogenies
  tree <- list()
  mat <- list()
  asb <- list()

  # Create 20 random assemblages
  for(i in 1:20){
    # Create a random phylogeny
    tree[[i]] <- ape::rcoal(20, br = c(min = 2, max = 10)) # Sets a minimum and maximum length for branch lengths
    # Create a random matrix
    mat[[i]] <- matrix(sample(c(1, 0), 20 * 10, replace = TRUE), ncol = 20, nrow = 10)
    colnames(mat[[i]]) <- tree[[i]]$tip.label # Name its columns according to tip names
  }

  # Creating a vector containing the number of slices desired to run the sensitivity
  vec <- c(25, 50, 75, 100, 125)
  # and the number of samples to evaluate it
  samp <- 10

  # Run the CpR_sensitivity algorithm for the CpD rate, and save its CpR_sensitivity_plot object,
  # while suppressing some warnings (related to tips no present within the matrix, and vice-versa)
  suppressWarnings({for(i in 1:20){
    CpR_obj <- CpR_sensitivity(tree[[i]], vec = vec, mat = mat[[i]], samp = samp, rate = "CpD")
    CpR_plot[[i]] <- CpR_sensitivity_plot(CpR_obj, rate = "CpD", stc = "mean")
  }})

  # Test
  for(i in 1:20){
    expect_s3_class(CpR_plot[[i]], "ggplot")
  }
})
