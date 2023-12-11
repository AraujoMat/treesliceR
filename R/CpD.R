#' Calculates the rate of accumulation of phylogenetic diversity (CpD) over time slices
#' @description
#' This function estimates the rates of accumulation of phylogenetic diveristy (CpD) over time for inputted assemblages.
#'
#' @usage CpD(tree, n, mat, criteria = "my", pDO = 5, ncor = 0)
#'
#' @param tree phylo. An ultrametric phylogenetic tree in the "phylo" format.
#' @param n numeric. A numeric value indicating the number of temporal slices (method = 1) or the time interval in million years (or phylogenetic diversity) among the tree slices (method = 2). Default is 1.
#' @param mat matrix. A presence/absence matrix containing all studied species and sites.
#' @param criteria character string. The method for cutting the tree. It can be either "my" (million years) or "PD" (accumulated phylogenetic diversity). Default is "my".
#' @param pDO numeric. A value indicating the numeric proportion to define the temporal origin at which the phylogenetic diversity (PD) started to accumulate in a given assemblage. Default is 5%.
#' @param ncor numeric. A value indicating the number of cores the user wants to parallelize. Default is 0.
#'
#' @return The function returns a data frame containing the assemblages' rates of cumulative phylogenetic diversity (CpD), their total phylogenetic diversity (PD), and their PD origin (pDO).
#'
#' @details
#'
#' \bold{Parallelization}
#'
#' Users are advised to check the number of available cores within their machines before running parallel programming.
#'
#' @seealso Other cumulative phylogenetic rates analysis: [CpE()], [CpB()], [CpB_RW()]
#'
#' @author Matheus Lima de Araujo <matheusaraujolima@live.com>
#'
#' @examples
#' # Generate a random tree
#' tree <- ape::rcoal(20)
#'
#' # Create a presence-absence matrix
#' mat <- matrix(sample(c(1,0), 20*10, replace = TRUE), ncol = 20, nrow = 10)
#' colnames(mat) <- tree$tip.label
#'
#' # Calculate the CpD for 100 tree slices
#' CpD(tree, n = 100, mat = mat)
#'
#' @export

CpD <- function(tree, n, mat, criteria = "my", pDO = 5, ncor = 0){

  ## Cleaning the phylogeny (if necessary) and cutting it into pieces ----------

  # If doesnt have a matrix, but a vector of species,
  # transform it into a species matrix
  if((is.matrix(mat) | is.data.frame(mat)) == FALSE){
    mat <- t(as.matrix(mat))
  }

  # Checking if there is species in the matrix without presence and removing them
  spps_pa <- colSums(mat)
  if(sum(spps_pa == 0) > 0){

    if(nrow(mat) == 1){
      # Removing those species
      mat <- t(as.matrix(mat[, -(which(spps_pa == 0))]))
      warning("Removing the species in presence abscence matrix without any occurrence")
    } else {
      # Removing those species
      mat <- mat[, -(which(spps_pa == 0))]
      warning("Removing the species in presence abscence matrix without any occurrence")
    }
  }

  # Dropping the lineage tips that arent in my spp matrix
  if(all(tree$tip.label %in% colnames(mat)) == FALSE){
    tree <- ape::keep.tip(tree, intersect(tree$tip.label, colnames(mat)))
    warning("Removing tips from phylogeny that are absent on species matrix")
  }

  # Cutting the phylogenetic tree into equal width slices
  branch_pieces <- phylo_pieces(tree, n, criteria = criteria,
                                timeSteps = TRUE, returnTree = TRUE)

  # Separating the time steps from the phylogenetic pieces
  age <- branch_pieces[[2]][length(branch_pieces[[2]]):1]
  tree <- branch_pieces[[3]]
  branch_pieces <- branch_pieces[[1]]
  i <- NULL

  ## Calculating assemblages CpD, PD and pDO -----------------------------------
  # The user wants to use more CPU cores?
  if(ncor > 0){
    # Register the number of desired clusters
    doParallel::registerDoParallel(ncor)
    # Loop and capture the values for each assemblage
    CPD <- foreach::foreach(i = 1:nrow(mat), .combine = rbind) %dopar% {

      # Which species are within this assemblage
      tips <- which(mat[i,] > 0) # (spp numbers follows same position on node matrix)  i <- 1000

      if(length(tips) == 0) {
        CPDrate <- c(NA, NA, NA)
        return(CPDrate)

      } else {
        # Obtaining the node matrix for those species
        if(length(tips) == 1){
          #  if there is a single spp, which nodes give origin to my species
          nodes <- as.numeric(names(which(rowSums(as.data.frame(tree$node_matrix[, tips])) > 0)))
        }

        if(length(tips) > 1){
          # Which nodes give origin to my species
          nodes <- as.numeric(names(which(rowSums(tree$node_matrix[, tips]) > 0)))
        }

        # Obtaining the tips and node positions
        positions <- which(tree$edge[,2] %in% c(tips, nodes))

        # Calculating the relative PE on each phylo slice
        CPD <- sapply(branch_pieces, function(x){     # x <- branch_pieces[[1800]]
          # Calculating the PE stored on each tree slice
          return(sum(x$edge.length[positions]))
        })

        # Saving the CPD
        CPD <- cumsum(CPD)/sum(CPD)  #  CPD[1809]

        # Fitting the non-linear model to obtain the CpB-rate:
        if(is.na(CPD[1]) == FALSE){ # Checking if there is any NA emerging due 0/0 fraction
          nonlinear_cum <- stats::nls(CPD ~ exp(-r*age),
                                      start = list(r = 0.2),
                                      control = stats::nls.control(maxiter = 1000))   # plot(CPE ~ age)
          # Saving the exponentital coefficient
          CPDrate <- stats::coef(nonlinear_cum)
        } else {
          # If it cant be calculated, return NA
          CPDrate <- NA
        }

        # Obtaining the PE
        PD <- sum(tree$edge.length[positions])

        # Creating a vec output containing the CpE-rate, pEO
        CPDrate <- c(CPDrate, PD, ((-1)*(log(pDO/100)/CPDrate)))
        return(CPDrate)
      }
    }
    # stop the clusters after running the algorithm
    doParallel::stopImplicitCluster()
  }

  # The user do not set clusters for being used
  if(ncor == 0){
    # Loop and capture the values for each assemblage
    CPD <- foreach::foreach(i = 1:nrow(mat), .combine = rbind) %do% {

      # Which species are within this assemblage
      tips <- which(mat[i,] > 0) # (spp numbers follows same position on node matrix)  i <- 1000

      if(length(tips) == 0) {
        CPDrate <- c(NA, NA, NA)
        return(CPDrate)

      } else {
        # Obtaining the node matrix for those species
        if(length(tips) == 1){
          #  if there is a single spp, which nodes give origin to my species
          nodes <- as.numeric(names(which(rowSums(as.data.frame(tree$node_matrix[, tips])) > 0)))
        }

        if(length(tips) > 1){
          # Which nodes give origin to my species
          nodes <- as.numeric(names(which(rowSums(tree$node_matrix[, tips]) > 0)))
        }

        # Obtaining the tips and node positions
        positions <- which(tree$edge[,2] %in% c(tips, nodes))

        # Calculating the relative PE on each phylo slice
        CPD <- sapply(branch_pieces, function(x){     # x <- branch_pieces[[1800]]
          # Calculating the PE stored on each tree slice
          return(sum(x$edge.length[positions]))
        })

        # Saving the CPD
        CPD <- cumsum(CPD)/sum(CPD)  #  CPD[1809]

        # Fitting the non-linear model to obtain the CpB-rate:
        if(is.na(CPD[1]) == FALSE){ # Checking if there is any NA emerging due 0/0 fraction
          nonlinear_cum <- stats::nls(CPD ~ exp(-r*age),
                                      start = list(r = 0.2),
                                      control = stats::nls.control(maxiter = 1000))   # plot(CPE ~ age)
          # Saving the exponentital coefficient
          CPDrate <- stats::coef(nonlinear_cum)
        } else {
          # If it cant be calculated, return NA
          CPDrate <- NA
        }

        # Obtaining the PE
        PD <- sum(tree$edge.length[positions])

        # Creating a vec output containing the CpE-rate, pEO
        CPDrate <- c(CPDrate, PD, ((-1)*(log(pDO/100)/CPDrate)))
        return(CPDrate)
      }
    }
  }

  # Renaming the columns and rows, and returning them as output
  if(nrow(mat) > 1){
    colnames(CPD) <- c("CpD", "PD", "pDO")
    rownames(CPD) <- 1:nrow(mat)
  } else {
    names(CPD) <- c("CpD", "PD", "pDO")
  }

  return(CPD)
}
