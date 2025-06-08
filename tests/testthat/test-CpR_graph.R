# Testing if all the saved outputs are a ggplot object
test_that("CpR_graph returns the correct class of outputs", {

  # Number of iterations
  n <- sample(10:20, 20, replace = TRUE)

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

  # Run the CpD rate function, and save its CpR_graph objects while suppressing some warnings
  # (related to tips no present within the matrix, and vice-versa)
  suppressWarnings({for(i in 1:20){
    CpR_obj <- CpD(tree[[i]], 20, mat = mat[[i]])
    CpR_plot[[i]] <- CpR_graph(CpR_obj, rate = "CpD")
  }})

  # Test
  for(i in 1:20){
    expect_s3_class(CpR_plot[[i]], "ggplot")
  }
})
