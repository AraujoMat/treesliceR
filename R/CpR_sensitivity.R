#' Runs a sensitivity analysis for rates of accumulation of a given phylogenetic index
#' @description
#' This function allows the evaluation of the sensitivity of the estimated rates of accumulation of a given phylogenetic index (e.g., [CpD()], [CpE()], [CpB()], [CpB_RW()]) to the number of slices inputted by the user.
#'
#' @usage CpR_sensitivity(tree, vec, mat, adj, rate, samp, comp, method, criterion, ncor)
#'
#' @param tree phylo. An ultrametric phylogenetic tree in the "phylo" format.
#' @param vec numeric vector. A numeric vector containing a series of numbers of slices.
#' @param mat matrix. A presence/absence matrix containing all studied species and sites.
#' @param adj matrix. A square adjacency matrix containing the presence/absence information of all sites and their spatially adjacent ones.
#' @param rate character string. The desired cumulative phylogenetic rate to be assessed, which can be the phylogenetic diversity (CpD), phylogenetic endemism (CpE), phylogenetic B-diversity (CpB), or phylogenetic B-diversity range-weighted (CpB_RW). Default is NULL, but must be filled with "CpD", "CPE", "CpB_RW", or "CpB".
#' @param samp numeric. The number of assemblages, or sites, to be sampled to make the sensitivity analysis.
#' @param comp character string. The component of beta-diversity that the user wants to calculate the CpB. It can be either "sorensen", turnover" or "nestedness". This argument works only when "rate = CpB". Default is "sorensen".
#' @param method character string. The method for calculating the CpB-rate. It can be either "pairwise" or "multisite". This argument works only when the argument "rate" is set to run for "CpB" or "CpB_RW". Default is "multisite".
#' @param criterion character string. The method for cutting the tree. It can be either "my" (million years) or "PD" (accumulated phylogenetic diversity). Default is "my".
#' @param ncor numeric. A value indicating the number of cores the user wants to parallelize. Default is 0.
#'
#' @return This function returns a data frame containing the sensitivity analysis for a given rate of accumulation of a phylogenetic index. This outputted data frame contains, for each row or assemblage, a column with the rate value assessed for each inputted number of slices.
#'
#' @details
#'
#' \bold{Parallelization}
#'
#' Users are advised to check the number of available cores within their machines before running parallel programming.
#'
#' \bold{Plotting}
#'
#' For plotting the sensitivity analysis output users can use [CpR_sensitivity_plot()].
#'
#' @seealso Other cumulative phylogenetic index rate analysis: [CpD()], [CpE()], [CpB()], [CpB_RW()]
#'
#' @author Matheus Lima de Araujo <matheusaraujolima@live.com>
#'
#' @references
#' See the tutorial on how to use this function on our [website](https://araujomat.github.io/treesliceR/articles/Passeriformes-diversification.html).
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
#' # Create a vector of number of slices
#' vec <- c(25, 50, 75, 100, 125, 150)
#'
#' # Calculate the sensitivity of the CpD
#' CpR_sensitivity(tree, vec, mat, rate = "CpD", samp = 5)
#'
#' @export


CpR_sensitivity <- function(tree, vec, mat = NULL, adj = NULL, rate = NULL, samp = 0,
                            comp = "sorensen", method = "multisite",
                            criterion = "my", ncor = 0){

  # The user provided the rate to be tested?
  if(is.null(rate) == TRUE){
    stop("The desired rate must be provided on rate argument")
  } else {
    # The user provided the number of sites to be sampled?
    if(samp == 0){
      stop("The number of sites for running the sensitivity analysis was not provided")
    } else {

      # Creating an empty "n"
      n <- NULL

      ## Cumulative Phylogenetic Diversity (CpD) -------------------------------

      if(rate == "CpD"){
        # If the samp inputed is bigger than the number of sites, return a error
        if(samp > nrow(mat)){
          stop("The number of site samples inputted in bigger than the sites available")
        } else {
          # Sampling random sites to run the algorithm
          sites <- sample(1:nrow(mat), samp)

          if(length(sites) == 1){
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpD(tree, n, mat[sites,], criterion = criterion, ncor = ncor)[1]
              return(output)
            }
            # Renaming the sensitivity data.frame object
            colnames(sensitivity) <- as.character(vec)
          } else {
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpD(tree, n, mat[sites,], criterion = criterion, ncor = ncor)[,1]
              return(output)
            }
            # Renaming the sensitivity obj
            colnames(sensitivity) <- as.character(vec)
            rownames(sensitivity) <- 1:nrow(sensitivity)
          }
        }
      }

      ## Cumulative Phylogenetic Endemism (CpE) --------------------------------

      if(rate == "CpE"){
        # If the samp inputed is bigger than the number of sites, return a error
        if(samp > nrow(mat)){
          stop("The number of site samples inputted in bigger than the sites available")
        } else {
          # Sampling random sites to run the algorithm
          sites <- sample(1:nrow(mat), samp)

          if(length(sites) == 1){
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpE(tree, n, mat[sites,], criterion = criterion, ncor = ncor)[1]
              return(output)
            }
            # Renaming the sensitivity data.frame object
            colnames(sensitivity) <- as.character(vec)

          } else {
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpE(tree, n, mat[sites,], criterion = criterion, ncor = ncor)[,1]
              return(output)
            }
            # Renaming the sensitivity obj
            colnames(sensitivity) <- as.character(vec)
            rownames(sensitivity) <- 1:nrow(sensitivity)
          }
        }
      }

      ## Cumulative Phylogenetic B-Diversity (CpB) -----------------------------

      if(rate == "CpB"){

        # If the samp inputed is bigger than the number of sites, return a error
        if(samp > nrow(mat)){
          stop("The number of site samples inputted in bigger than the sites available")
        } else {
          # Sampling random sites to run the algorithm
          sites <- sample(1:nrow(mat), samp)

          if(length(sites) == 1){
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpB(tree, n, mat, adj[sites, , drop = FALSE], method = method, comp = comp,
                            criterion = criterion, ncor = ncor)[1]
              return(output)
            }
            # Renaming the sensitivity data.frame object
            colnames(sensitivity) <- as.character(vec)

          } else {
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpB(tree, n, mat, adj[sites, , drop = FALSE], method = method, comp = comp,
                            criterion = criterion, ncor = ncor)[,1]
              return(output)
            }
            # Renaming the sensitivity obj
            colnames(sensitivity) <- as.character(vec)
            rownames(sensitivity) <- 1:nrow(sensitivity)
          }
        }
      }

      ## Cumulative Phylogenetic B-Diversity Range-Weighted (CpB RW) -----------

      if(rate == "CpB_RW"){

        # If the samp inputed is bigger than the number of sites, return a error
        if(samp > nrow(mat)){
          stop("The number of site samples inputted in bigger than the sites available")
        } else {
          # Sampling random sites to run the algorithm
          sites <- sample(1:nrow(mat), samp)

          if(length(sites) == 1){
            # Running the sensitivity algorithm
            sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
              output <- CpB_RW(tree, n, mat, adj[sites, , drop = FALSE], method = method,
                               criterion = criterion, ncor = ncor)[1]
              return(output)
            }
            # Renaming the sensitivity data.frame object
            colnames(sensitivity) <- as.character(vec)

            } else {
              # Running the sensitivity algorithm
              sensitivity <- foreach::foreach(n = vec, .combine = cbind) %do% {
                output <- CpB_RW(tree, n, mat, adj[sites, , drop = FALSE], method = method,
                                 criterion = criterion, ncor = ncor)[,1]
                return(output)
              }
              # Renaming the sensitivity obj
              colnames(sensitivity) <- as.character(vec)
              rownames(sensitivity) <- 1:nrow(sensitivity)
          }
        }
      }

      # Returning it
      return(sensitivity)
    }
  }
}
