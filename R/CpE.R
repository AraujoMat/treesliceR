#' Calculates the rate of accumulation of phylogenetic endemism (CpE) over time slices
#' @description
#' This function estimates the rates of accumulation of phylogenetic endemism (CpE) over time for inputted assemblages.
#'
#' @usage CpE(tree, n, mat, criterion = "my", pEO = 5, ncor = 0)
#'
#' @param tree phylo. An ultrametric phylogenetic tree in the "phylo" format.
#' @param n numeric. A numeric value indicating the number of temporal slices (method = 1) or the time interval in million years (or phylogenetic diversity) among the tree slices (method = 2). Default is 1.
#' @param mat matrix. A presence/absence matrix containing all studied species and sites.
#' @param criterion character string. The method for cutting the tree. It can be either "my" (million years) or "PD" (accumulated phylogenetic diversity). Default is "my".
#' @param pEO numeric. A value indicating the numeric proportion to define the temporal origin at which the phylogenetic endemism (PE) started to accumulate in a given assemblage. Default is 5%.
#' @param ncor numeric. A value indicating the number of cores the user wants to parallelize. Default is 0.
#'
#' @return The function returns a data frame containing the assemblages' rates of cumulative phylogenetic endemism (CpE), their total phylogenetic endemism (PE), and their PE origin (pEO).
#'
#' @details
#'
#' \bold{Parallelization}
#'
#' Users are advised to check the number of available cores within their machines before running parallel programming.
#'
#' @seealso Other cumulative phylogenetic index analysis: [CpD()], [CpB()], [CpB_RW()]
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
#' # Calculate the CpE for 100 tree slices
#' CpE(tree, n = 100, mat = mat)
#'
#' @export

CpE <- function(tree, n, mat, criterion = "my", pEO = 5, ncor = 0){

  ## Cleaning the phylogeny (if necessary) and cutting it into pieces ----------

  # If a doesnt have a matrix, but a vector of species,
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
  branch_pieces <- phylo_pieces(tree, n, criterion = criterion,
                                timeSteps = TRUE, returnTree = TRUE)

  # Separating the time steps from the phylogenetic pieces
  age <- branch_pieces[[2]][length(branch_pieces[[2]]):1]
  tree <- branch_pieces[[3]]
  branch_pieces <- branch_pieces[[1]]
  i <- NULL

  ## Calculating the branch range size within each node branch -----------------

  # If there is only one site, the range of all nodes equals 1
  if(nrow(mat) == 1){
    r_sizes <- rep(1, length(tree$edge[,2]))
  } else {
    r_sizes <- sapply(tree$edge[,2], function(x){
      # If its a tip branch
      if(x < tree$edge[1, 1]){
        # There is only one species within the node, calculate and save its range
        deno <- sum(mat[, tree$tip.label[x]])
      } else {
        # If its a node branch
        # Which species share those nodes
        spps_node <- names(which(tree$node_matrix[as.character(x), ] > 0))  # node 1    x <- 689

        # Capturing the range size denominator from the node
        if(length(spps_node) == 1){
          # If there is only one species within the node
          deno <- sum(mat[, spps_node])
        } else { # If there is more species
          # Which is the shared range size of the species sharing this node
          deno <- sum(rowSums(mat[, spps_node]) > 0)
        }
      }
      return(deno)
    })
  }


  ## Calculating assemblages CpE, PE and pEO -----------------------------------
  # The user wants to use more CPU cores?
  if(ncor > 0){
    # Register the number of desired clusters
    doParallel::registerDoParallel(ncor)
    # Loop and capture the values for each assemblage
    CPE <- foreach::foreach(i = 1:nrow(mat), .combine = rbind) %dopar% {

      # Which species are within this assemblage
      tips <- which(mat[i,] > 0) # (spp numbers follows same position on node matrix)  i <- 1000

      if(length(tips) == 0) {
        CPErate <- c(NA, NA, NA)
        return(CPErate)

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
        CPE <- sapply(branch_pieces, function(x){     # x <- branch_pieces[[1800]]
          # Calculating the PE stored on each tree slice
          return(sum(x$edge.length[positions]/r_sizes[positions]))
        })

        # Saving the CPE
        CPE <- cumsum(CPE)/sum(CPE)  #  CPE[1809]

        # Fitting the non-linear model to obtain the CpB-rate:
        if(is.na(CPE[1]) == FALSE){ # Checking if there is any NA emerging due 0/0 fraction
          nonlinear_cum <- stats::nls(CPE ~ exp(-r*age),
                                      start = list(r = 0.2),
                                      control = stats::nls.control(maxiter = 1000))   # plot(CPE ~ age)
          # Saving the exponentital coefficient
          CPErate <- stats::coef(nonlinear_cum)
        } else {
          # If it cant be calculated, return NA
          CPErate <- NA
        }

        # Obtaining the PE
        PE <- sum(tree$edge.length[positions]/r_sizes[positions])

        # Creating a vec output containing the CpE-rate, pEO
        CPErate <- c(CPErate, PE, ((-1)*(log(pEO/100)/CPErate)))
        return(CPErate)
      }
    }
    # stop the clusters after running the algorithm
    doParallel::stopImplicitCluster()
  }

  # The user do not set clusters for being used
  if(ncor == 0){
    # Loop and capture the values for each assemblage
    CPE <- foreach::foreach(i = 1:nrow(mat), .combine = rbind) %do% {

      # Which species are within this assemblage
      tips <- which(mat[i,] > 0) # (spp numbers follows same position on node matrix)  i <- 1

      if(length(tips) == 0) {
        CPErate <- c(NA, NA, NA)
        return(CPErate)

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
        CPE <- sapply(branch_pieces, function(x){     # x <- branch_pieces[[1800]]
          # Calculating the PE stored on each tree slice
          return(sum(x$edge.length[positions]/r_sizes[positions]))
        })

        # Saving the CPE
        CPE <- cumsum(CPE)/sum(CPE)  #  CPE[1809]

        # Fitting the non-linear model to obtain the CpB-rate:
        if(is.na(CPE[1]) == FALSE){ # Checking if there is any NA emerging due 0/0 fraction
          nonlinear_cum <- stats::nls(CPE ~ exp(-r*age),
                                      start = list(r = 0.2),
                                      control = stats::nls.control(maxiter = 1000))   # plot(CPE ~ age)
          # Saving the exponentital coefficient
          CPErate <- stats::coef(nonlinear_cum)
        } else {
          # If it cant be calculated, return NA
          CPErate <- NA
        }

        # Obtaining the PE
        PE <- sum(tree$edge.length[positions]/r_sizes[positions])

        # Creating a vec output containing the CpE-rate, pEO
        CPErate <- c(CPErate, PE, ((-1)*(log(pEO/100)/CPErate)))
        return(CPErate)
      }
    }
  }

  # Renaming the columns and rows, and returning them as output
  if(nrow(mat) > 1){
    colnames(CPE) <- c("CpE", "PE", "pEO")
    rownames(CPE) <- 1:nrow(mat)
  } else {
    names(CPE) <- c("CpE", "PE", "pEO")
  }

  return(CPE)
}
