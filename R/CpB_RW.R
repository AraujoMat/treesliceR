#' Calculates the range weighted rate of accumulation of phylogenetic B-diversity (CpB_RW) over time slices
#' @description
#' This function estimates the range-weighted rates of accumulation of phylogenetic B-diversity (CpB_RW) over time for inputted assemblages.
#'
#' @usage CpB_RW(tree, n, mat, adj, method = "multisite", criterion = "my", pBO = 5, ncor = 0)
#'
#' @param tree phylo. An ultrametric phylogenetic tree in the "phylo" format.
#' @param n numeric. A numeric value indicating the number of temporal slices (method = 1) or the time interval in million years (or phylogenetic diversity) among the tree slices (method = 2). Default is 1.
#' @param mat matrix. A presence/absence matrix containing all studied species and sites.
#' @param adj matrix. A square adjacency matrix containing the presence/absence information of all sites and their spatially adjacent ones.
#' @param method character string. The method for calculating the phylogenetic beta-diversity. It can be either obtained through a "pairwise" or "multisite" approach. Default is "multisite".
#' @param criterion character string. The method for cutting the tree. It can be either "my" (million years) or "PD" (accumulated phylogenetic diversity). Default is "my".
#' @param pBO numeric. A value indicating the numeric proportion to define the temporal origin at which the range-weighted phylogenetic B-diversity (PB_RW) started to accumulate in a given assemblage. Default is 5%.
#' @param ncor numeric. A value indicating the number of cores the user wants to parallelize. Default is 0.
#'
#' @return The function returns a data frame containing the assemblages' rates of cumulative range-weighted phylogenetic B-diversity (CpB_RW), their total range-weighted phylogenetic B-diversity (PB_RW), and their origin (pBO).
#'
#' @details
#'
#' \bold{Parallelization}
#'
#' Users are advised to check the number of available cores within their machines before running parallel programming.
#'
#' @seealso Other cumulative phylogenetic index analysis: [CpD()], [CpE()], [CpB()]
#'
#' @author Matheus Lima de Araujo <matheusaraujolima@live.com>
#'
#' @references
#' Laffan, S. W., Rosauer, D. F., Di Virgilio, G., Miller, J. T., González-Orozco, C. E., Knerr, N., Thornhill, A. H., & Mishler, B. D. (2016). Range-weighted metrics of species and phylogenetic turnover can better resolve biogeographic transition zones. Methods in Ecology and Evolution, 7(5), 580–588. https://doi.org/10.1111/2041-210x.12513
#'
#' @examples
#' # Generate a random tree
#' tree <- ape::rcoal(20)
#'
#' # Create a presence-absence matrix
#' mat <- matrix(sample(c(1,0), 20*10, replace = TRUE), ncol = 20, nrow = 10)
#' colnames(mat) <- tree$tip.label
#'
#' # Create a random adjacency matrix
#' adj <- matrix(sample(c(1,0), 10*10, replace = TRUE), ncol = 10, nrow = 10)
#'
#' # Fill the diagonals with 1
#' diag(adj) <- 1
#'
#' # Calculate their CpB range weighted for 100 tree slices
#' CpB_RW(tree, n = 100, mat = mat, adj = adj, method = "multisite")
#'
#' @export

CpB_RW <- function(tree, n, mat, adj, method = "multisite", criterion = "my", pBO = 5, ncor = 0){

  ## Cleaning the phylogeny (if necessary) and cutting it into pieces ----------

  # Checking if there is species in the matrix without presence and removing them
  spps_pa <- colSums(mat)
  if(sum(spps_pa == 0) > 0){
    # Removing those species
    mat <- mat[, -(which(spps_pa == 0))]
    warning("Removing the species in presence abscence matrix without any occurrence")
  }

  # Dropping the lineage tips that arent in my spp matrix
  if(all(tree$tip.label %in% colnames(mat)) == FALSE){
    tree <- ape::keep.tip(tree, intersect(tree$tip.label, colnames(mat)))
    warning("Removing tips from phylogeny that are absent on species matrix")
  }

  # Creating a list containing the focal and adjacent sites
  asb <- lapply(1:nrow(adj), function(x){
    # The matrix has only one site?
    if(nrow(adj) <= 1 & nrow(mat) <= 1){
      return(stop("At least one additional site must be provided to conduct B-diversity analysis."))
      # If it has more, list them
    } else if(all(which(adj[x,] == 1) %in% 1:nrow(mat)) == FALSE) {
      return(stop("Not all adjacent sites are included in the adjacency matrix."))
      # If everything is OK:
    } else {
      return(mat[which(adj[x,] == 1),, drop = FALSE])
    }
  })

  # Cutting the phylogenetic tree into equal width slices
  branch_pieces <- phylo_pieces(tree, n, criterion = criterion,
                                timeSteps = TRUE, returnTree = TRUE)

  # Separating the time steps from the phylogenetic pieces
  age <- branch_pieces[[2]][length(branch_pieces[[2]]):1]
  tree <- branch_pieces[[3]]
  branch_pieces <- branch_pieces[[1]]
  commu <- NULL

  ## Calculating the branch range size within each node branch -----------------
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


  ## Calculating the CpB_RW (using phylo-Sorensen) -----------------------------
  # The user wants to use more CPU cores?
  if(ncor > 0){

    # Register the number of desired clusters
    doParallel::registerDoParallel(ncor)

    # Calculating the CpB rate through the CpBrate algorithm
    CpBrate <- foreach::foreach(commu = asb, .combine = rbind) %dopar% {  # commu <- communities[[1]] asb <- communities commu <- communities[[8]]

      # Conditions, minimum of species, sites, etc
      if(nrow(commu) >= 2 & ncol(commu) >= 2 &
         all(commu == 1) == FALSE & all(commu == 0) == FALSE){

        # Combining the paired assemblages into a new species presence-absence matrix
        comb_commus <- t(apply(utils::combn(nrow(commu), 2), 2, function(x) colSums(commu[x,]) > 0))

        ## Find the species in each community
        species <- colnames(commu)

        ## Defining the nodes contained in each observed and paired matrix
        # OBS
        commus_nodes <- apply(commu, 1, function(x) { # x <- 1
          # Which species are present in my assemblage
          present_obs <- species[x > 0] # present_obs <- species[commu[1,] > 0]

          # Which is the species position in my node matrix
          spps <- which(colnames(tree$node_matrix) %in% present_obs)

          # if there is a single spp
          if(length(spps) == 1){
            # Which nodes give origin to my species
            nodes <- as.numeric(names(which(rowSums(as.data.frame(tree$node_matrix[,spps])) > 0)))

          } else {
            # Which nodes give origin to my species
            nodes <- as.numeric(names(which(rowSums(tree$node_matrix[,spps]) > 0)))
          }

          # Lines from the edge.length to preserve
          lines_prs <- which(tree$edge[, 1] %in% nodes &             # Which nodes are in my spliting side from matrix
                               tree$edge[, 2] %in% c(spps, nodes))   # Which nodes and tips are in my splitted side from matrix

        })
        # PWR
        p_commus_nodes <- lapply(1:nrow(comb_commus), function(x) { # x <- 1
          # Which species are present in my assemblage
          present_obs <- species[comb_commus[x,] > 0] # present_obs <- species[comb_commus[1,] > 0]

          # Which is the species position in my node matrix
          spps <- which(colnames(tree$node_matrix) %in% present_obs)


          # if there is a single spp
          if(length(spps) == 1){
            # Which nodes give origin to my species
            nodes <- as.numeric(names(which(rowSums(as.data.frame(tree$node_matrix[,spps])) > 0)))

          } else {
            # Which nodes give origin to my species
            nodes <- as.numeric(names(which(rowSums(tree$node_matrix[, spps]) > 0)))
          }

          # Lines from the edge.length to preserve
          lines_prs <- which(tree$edge[, 1] %in% nodes &             # Which nodes are in my spliting side from matrix
                               tree$edge[, 2] %in% c(spps, nodes))   # Which nodes and tips are in my splitted side from matrix

        })


        ## Calculating the PD contained on each community, but for the complete phylogenetic tree
        # For each commu
        pd_assemblages <- sapply(1:length(commus_nodes), function(x){
          return(sum(tree$edge.length[commus_nodes[[x]]]/r_sizes[commus_nodes[[x]]]))
        })
        # For pairwised comparisions among assemblages
        pw_pd_assemblages <- sapply(1:length(p_commus_nodes), function(x){
          return(sum(tree$edge.length[p_commus_nodes[[x]]]/r_sizes[p_commus_nodes[[x]]]))
        })


        ########################
        #### SORENSEN INDEX ####
        ########################

        if(method == "pairwise"){
          ## Capturing the community parameters a, b, c (CpB denominator)
          bc_sum <- apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, sum)
          # capturing "a"
          a <- pw_pd_assemblages - bc_sum
          # Sorensen
          pbd <- mean(bc_sum/((2*a) + bc_sum))
        }

        if(method == "multisite"){
          # Species nodes for all species within all assemblages (TOTAL)
          t_commus_nodes <- apply(t(as.matrix(colSums(commu) > 0)), 1, function(x) { # x <- 1
            # Which species are present in my assemblage
            present_obs <- species[x > 0] # present_obs <- species[commu[1,] > 0]

            # Which is the species position in my node matrix
            spps <- which(colnames(tree$node_matrix) %in% present_obs)
            # Which nodes give origin to my species
            nodes <- as.numeric(names(which(rowSums(tree$node_matrix[,spps]) > 0)))

            # Lines from the edge.length to preserve
            lines_prs <- which(tree$edge[, 1] %in% nodes &             # Which nodes are in my spliting side from matrix
                                 tree$edge[, 2] %in% c(spps, nodes))   # Which nodes and tips are in my splitted side from matrix

          })
          # Calculating the total PD for the complete tree
          t_pd_assemblages <- sum(tree$edge.length[t_commus_nodes]/r_sizes[t_commus_nodes])

          ## Capturing the community parameters a, b, c (CpB denominator)
          bc_sum <- sum(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)))
          a <- sum(sum(pd_assemblages) - t_pd_assemblages)   # t_pd_assemblages/5  if no one branch lenght were shared among assemblages, it is the value expected. In that case, a will equal 0.
          # Sorensen
          pbd <- bc_sum/((2*a) + bc_sum)
        }

        ## Calculating the PROPORTION of TURNOVER (sorensen) accumulation through time
        CpB <- sapply(branch_pieces, function(x){
          ## PDs for each phylogenetic tree slice
          # For each commu
          pd_piece <- sapply(1:length(commus_nodes), function(y){   # x <- branch_pieces[[1]]
            return(sum(x$edge.length[commus_nodes[[y]]]/r_sizes[commus_nodes[[y]]]))
          })
          # For pairwised comparisions among assemblages
          pd_pw_piece <- sapply(1:length(p_commus_nodes), function(y){
            return(sum(x$edge.length[p_commus_nodes[[y]]]/r_sizes[p_commus_nodes[[y]]]))
          })

          if(method == "pairwise"){
            # Calculating the bc sum within the slices (relatice (b + c) on slice) -> numerator sorensen
            r_bc_sum <- apply(pd_pw_piece - t(utils::combn(pd_piece, 2)), 1, sum)
            # Relative Sorensen
            r_pbd <- mean(r_bc_sum/((2*a) + bc_sum))
          }

          if(method == "multisite"){
            # Calculating the bc sum within the slices (relatice (b + c) on slice) -> numerator sorensen
            r_bc_sum <- sum(pd_pw_piece - t(utils::combn(pd_piece, 2)))
            # Relative Sorensen
            r_pbd <- r_bc_sum/((2*a) + bc_sum)
          }
          # Capturing the proportion of turnover in stored in this slice
          return(r_pbd/pbd)
        })

        ##################
        #### CpB rate ####
        ##################

        # Fitting the non-linear model to obtain the CpB-rate:
        if(is.na(CpB[1]) == FALSE){ # Checking if there is any NA emerging due 0/0 fraction
          CpB <- stats::nls(cumsum(CpB) ~ exp(-r*age),
                            start = list(r = 0.2),
                            control = stats::nls.control(maxiter = 1000))   # plot(cumsum(CpB) ~ age)
          # Mergin all data
          CpB <- c(stats::coef(CpB), pbd, ((-1)*(log(pBO/100)/stats::coef(CpB))))
          return(CpB)
        } else {
          CpB <- c(NA, NA, NA)
          return(CpB)
        }
      } else {
        CpB <- c(NA, NA, NA)
        return(CpB)
      }
    }

    # stop the clusters after running the algorithm
    doParallel::stopImplicitCluster()
  }

  # The user do not set clusters for being used
  if(ncor == 0){
    # Calculating the CpB rate through the CpBrate algorithm
    CpBrate <- foreach::foreach(commu = asb, .combine = rbind) %do% {  # commu <- communities[[1]] commu <- asb[[1]]

      # Conditions, minimum of species, sites, etc
      if(nrow(commu) >= 2 & ncol(commu) >= 2 &
         all(commu == 1) == FALSE & all(commu == 0) == FALSE){

        # Combining the paired assemblages into a new species presence-absence matrix
        comb_commus <- t(apply(utils::combn(nrow(commu), 2), 2, function(x) colSums(commu[x,]) > 0))

        ## Find the species in each community
        species <- colnames(commu)

        ## Defining the nodes contained in each observed and paired matrix
        # OBS
        commus_nodes <- apply(commu, 1, function(x) { # x <- 1
          # Which species are present in my assemblage
          present_obs <- species[x > 0] # present_obs <- species[commu[1,] > 0]

          # Which is the species position in my node matrix
          spps <- which(colnames(tree$node_matrix) %in% present_obs)

          # if there is a single spp
          if(length(spps) == 1){
            # Which nodes give origin to my species
            nodes <- as.numeric(names(which(rowSums(as.data.frame(tree$node_matrix[,spps])) > 0)))

          } else {
            # Which nodes give origin to my species
            nodes <- as.numeric(names(which(rowSums(tree$node_matrix[,spps]) > 0)))
          }

          # Lines from the edge.length to preserve
          lines_prs <- which(tree$edge[, 1] %in% nodes &             # Which nodes are in my spliting side from matrix
                               tree$edge[, 2] %in% c(spps, nodes))   # Which nodes and tips are in my splitted side from matrix

        })
        # PWR
        p_commus_nodes <- lapply(1:nrow(comb_commus), function(x) { # x <- 1
          # Which species are present in my assemblage
          present_obs <- species[comb_commus[x,] > 0] # present_obs <- species[comb_commus[1,] > 0]

          # Which is the species position in my node matrix
          spps <- which(colnames(tree$node_matrix) %in% present_obs)


          # if there is a single spp
          if(length(spps) == 1){
            # Which nodes give origin to my species
            nodes <- as.numeric(names(which(rowSums(as.data.frame(tree$node_matrix[,spps])) > 0)))

          } else {
            # Which nodes give origin to my species
            nodes <- as.numeric(names(which(rowSums(tree$node_matrix[, spps]) > 0)))
          }

          # Lines from the edge.length to preserve
          lines_prs <- which(tree$edge[, 1] %in% nodes &             # Which nodes are in my spliting side from matrix
                               tree$edge[, 2] %in% c(spps, nodes))   # Which nodes and tips are in my splitted side from matrix

        })


        ## Calculating the PD contained on each community, but for the complete phylogenetic tree
        # For each commu
        pd_assemblages <- sapply(1:length(commus_nodes), function(x){
          return(sum(tree$edge.length[commus_nodes[[x]]]/r_sizes[commus_nodes[[x]]]))
        })
        # For pairwised comparisions among assemblages
        pw_pd_assemblages <- sapply(1:length(p_commus_nodes), function(x){
          return(sum(tree$edge.length[p_commus_nodes[[x]]]/r_sizes[p_commus_nodes[[x]]]))
        })


        ########################
        #### SORENSEN INDEX ####
        ########################

        if(method == "pairwise"){
          ## Capturing the community parameters a, b, c (CpB denominator)
          bc_sum <- apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, sum)
          # capturing "a"
          a <- pw_pd_assemblages - bc_sum
          # Sorensen
          pbd <- mean(bc_sum/((2*a) + bc_sum))
        }

        if(method == "multisite"){
          # Species nodes for all species within all assemblages (TOTAL)
          t_commus_nodes <- apply(t(as.matrix(colSums(commu) > 0)), 1, function(x) { # x <- 1
            # Which species are present in my assemblage
            present_obs <- species[x > 0] # present_obs <- species[commu[1,] > 0]

            # Which is the species position in my node matrix
            spps <- which(colnames(tree$node_matrix) %in% present_obs)
            # Which nodes give origin to my species
            nodes <- as.numeric(names(which(rowSums(tree$node_matrix[,spps]) > 0)))

            # Lines from the edge.length to preserve
            lines_prs <- which(tree$edge[, 1] %in% nodes &             # Which nodes are in my spliting side from matrix
                                 tree$edge[, 2] %in% c(spps, nodes))   # Which nodes and tips are in my splitted side from matrix

          })
          # Calculating the total PD for the complete tree
          t_pd_assemblages <- sum(tree$edge.length[t_commus_nodes]/r_sizes[t_commus_nodes])

          ## Capturing the community parameters a, b, c (CpB denominator)
          bc_sum <- sum(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)))
          a <- sum(sum(pd_assemblages) - t_pd_assemblages)   # t_pd_assemblages/5  if no one branch lenght were shared among assemblages, it is the value expected. In that case, a will equal 0.
          # Sorensen
          pbd <- bc_sum/((2*a) + bc_sum)
        }

        ## Calculating the PROPORTION of TURNOVER (sorensen) accumulation through time
        CpB <- sapply(branch_pieces, function(x){
          ## PDs for each phylogenetic tree slice
          # For each commu
          pd_piece <- sapply(1:length(commus_nodes), function(y){   # x <- branch_pieces[[1]]
            return(sum(x$edge.length[commus_nodes[[y]]]/r_sizes[commus_nodes[[y]]]))
          })
          # For pairwised comparisions among assemblages
          pd_pw_piece <- sapply(1:length(p_commus_nodes), function(y){
            return(sum(x$edge.length[p_commus_nodes[[y]]]/r_sizes[p_commus_nodes[[y]]]))
          })

          if(method == "pairwise"){
            # Calculating the bc sum within the slices (relatice (b + c) on slice) -> numerator sorensen
            r_bc_sum <- apply(pd_pw_piece - t(utils::combn(pd_piece, 2)), 1, sum)
            # Relative Sorensen
            r_pbd <- mean(r_bc_sum/((2*a) + bc_sum))
          }

          if(method == "multisite"){
            # Calculating the bc sum within the slices (relatice (b + c) on slice) -> numerator sorensen
            r_bc_sum <- sum(pd_pw_piece - t(utils::combn(pd_piece, 2)))
            # Relative Sorensen
            r_pbd <- r_bc_sum/((2*a) + bc_sum)
          }
          # Capturing the proportion of turnover in stored in this slice
          return(r_pbd/pbd)
        })

        ##################
        #### CpB rate ####
        ##################

        # Fitting the non-linear model to obtain the CpB-rate:
        if(is.na(CpB[1]) == FALSE){ # Checking if there is any NA emerging due 0/0 fraction
          CpB <- stats::nls(cumsum(CpB) ~ exp(-r*age),
                            start = list(r = 0.2),
                            control = stats::nls.control(maxiter = 1000))   # plot(cumsum(CpB) ~ age)
          # Mergin all data
          CpB <- c(stats::coef(CpB), pbd, ((-1)*(log(pBO/100)/stats::coef(CpB))))
          return(CpB)
        } else {
          CpB <- c(NA, NA, NA)
          return(CpB)
        }
      } else {
        CpB <- c(NA, NA, NA)
        return(CpB)
      }
    }
  }

  # Renaming the columns and rows, and returning them as output
  if(length(asb) > 1){
    colnames(CpBrate) <- c("CpB_RW", "PB_RW", "pBO")
    rownames(CpBrate) <- 1:length(asb)
  } else {
    names(CpBrate) <- c("CpB_RW", "PB_RW", "pBO")
  }

  return(CpBrate)
}
