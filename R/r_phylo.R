#' Calculates the relative value of a phylogenetic index in a temporal sequence of phylogenetic slices.
#' @description
#' This function estimates the relative value of a phylogenetic index in a sequence of multiple phylogenetic slices cut from roots to tips.
#'
#' @usage r_phylo(tree, n, mat, adj, index = NULL, comp, method, criterion = "my", ncor = 0)
#'
#' @param tree phylo. An ultrametric phylogenetic tree in the "phylo" format.
#' @param n numeric. A numeric value indicating either the number of temporal slices (method = 1) or the time interval in million years (or phylogenetic diversity) among the tree slices (method = 2). Default is 1.
#' @param mat matrix. A presence/absence matrix containing all studied species and sites.
#' @param adj matrix. A square adjacency matrix containing the presence/absence information of all sites and their spatially adjacent ones.
#' @param index character string. The phylogenetic index to be calculated over the phylogenetic slices. It can be set as "PD" (phylogenetic diversity), "PE" (phylogenetic endemism), "PB" (phylogenetic B-diversity), or "PB_RW" (phylogenetic B-diversity range-weighted).
#' @param comp character string. The component of phylogenetic beta-diversity to obtain the relative value. It can be "sorensen", "turnover", or "nestedness". Default is "sorensen".
#' @param method character string. The method for calculating phylogenetic beta-diversity. It can be obtained through a "pairwise" or "multisite" approach. Default is "multisite".
#' @param criterion character string. The method for cutting the tree. It can be either "my" (million years) or "PD" (accumulated phylogenetic diversity). Default is "my".
#' @param ncor numeric. A value indicating the number of cores the user wants to parallel. Default is 0.
#'
#' @return The function returns a list where each object contains a vector (of length "n") with the relative phylogenetic index, from the phylogeny root to the tips, from the inputted assemblage.
#'
#' @details
#'
#' \bold{The "adj" argument}
#'
#' Must be filled only for phylogenetic B-diversity ("PB") or it's range weight version ("PB_RW", defined in "index").
#'
#' \bold{Parallelization}
#'
#' Users are advised to check the number of cores available within their machines before running in parallel programming.
#'
#' @seealso Other cumulative phylogenetic rates analysis: [CpD()], [CpE()], [CpB()], [CpB_RW()]
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
#' # Create a random adjacency matrix
#' adj <- matrix(sample(c(1,0), 10*10, replace = TRUE), ncol = 10, nrow = 10)
#'
#' # Fill the diagonals with 1
#' diag(adj) <- 1
#'
#' # Calculate the relative PD for 100 slices
#' rPD <- r_phylo(tree, n = 100, mat = mat, index = "PD")
#' # Plot the relative PD of the first assemblage
#' plot(rPD[[1]])
#'
#' # Calculate the relative PE for 100 slices
#' rPE <- r_phylo(tree, n = 100, mat = mat, index = "PE")
#' # Plot the relative PE of the first assemblage
#' plot(rPE[[1]])
#'
#' # Calculate the relative PB for 100 slices
#' rPB <-  r_phylo(tree, n = 100, mat = mat, adj = adj, index = "PB")
#' # Plot the relative PB of the first assemblage
#' plot(rPB[[1]])
#'
#' # Calculate the relative PB_RW for 100 slices
#' rPB_RW <- r_phylo(tree, n = 100, mat = mat, adj = adj, index = "PB_RW")
#' # Plot the relative PB_RW of the first assemblage
#' plot(rPB_RW[[1]])
#'
#' @export

r_phylo <- function(tree, n, mat, adj, index = NULL, comp = "sorensen",
                    method = "multisite", criterion = "my", ncor = 0){

  # Checking if the index was provided by the user:
  if(is.null(index) == TRUE){
    stop("The desired phylogenetic index must be provided on index argument")
  } else {
    # Relative PD on each phylogenetic slice -------------------------------------
    if(index == "PD"){

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
      tree <- branch_pieces[[3]]
      branch_pieces <- branch_pieces[[1]]
      i <- NULL

      if(ncor > 0){
        # Register the number of desired clusters
        doParallel::registerDoParallel(ncor)
        # Loop and capture the values for each assemblage
        CPD <- foreach::foreach(i = 1:nrow(mat)) %dopar% {

          # Which species are within this assemblage
          tips <- which(mat[i,] > 0) # (spp numbers follows same position on node matrix)  i <- 10

          if(length(tips) == 0) {
            CPD <- NA
            return(CPD)

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

            # returning CPD
            return(CPD)  # plot(CPD)
          }
        }
        # stop the clusters after running the algorithm
        doParallel::stopImplicitCluster()
      }

      # The user do not set clusters for being used
      if(ncor == 0){
        # Loop and capture the values for each assemblage
        CPD <- foreach::foreach(i = 1:nrow(mat)) %do% {

          # Which species are within this assemblage
          tips <- which(mat[i,] > 0) # (spp numbers follows same position on node matrix)  i <- 10

          if(length(tips) == 0) {
            CPD <- NA
            return(CPD)

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

            # returning CPD
            return(CPD)  # plot(CPD)
          }
        }
      }

      # Return
      return(CPD)
    }

    # Relative PE on each phylogenetic slice -------------------------------------
    if(index == "PE"){
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
        CPE <- foreach::foreach(i = 1:nrow(mat)) %dopar% {

          # Which species are within this assemblage
          tips <- which(mat[i,] > 0) # (spp numbers follows same position on node matrix)  i <- 10

          if(length(tips) == 0) {
            CPE <- NA
            return(CPE)

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

            # Returning the relative PE
            return(CPE)
          }
        }
        # stop the clusters after running the algorithm
        doParallel::stopImplicitCluster()
      }

      # The user do not set clusters for being used
      if(ncor == 0){
        # Loop and capture the values for each assemblage
        CPE <- foreach::foreach(i = 1:nrow(mat)) %do% {

          # Which species are within this assemblage
          tips <- which(mat[i,] > 0) # (spp numbers follows same position on node matrix)  i <- 10

          if(length(tips) == 0) {
            CPE <- NA
            return(CPE)

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

            # Returning the relative PE
            return(CPE)
          }
        }
      }

      # Return
      return(CPE)
    }

    # Relative PB on each phylogenetic slice -------------------------------------
    if(index == "PB"){

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
          return(mat[which(adj[x,] == 1), , drop = FALSE])
        }
      })

      # Capturing the tree configs and
      # Cutting the phylogenetic tree into equal width slices
      branch_pieces <- phylo_pieces(tree, n, criterion = criterion,
                                    timeSteps = TRUE, returnTree = TRUE)

      # Separating the time steps from the phylogenetic pieces
      tree <- branch_pieces[[3]]
      branch_pieces <- branch_pieces[[1]]
      commu <- NULL

      # The user wants to use more CPU cores?
      if(ncor > 0){

        # Register the number of desired clusters
        doParallel::registerDoParallel(ncor)

        # Calculating the CpB rate through the CpBrate algorithm
        CpB <- foreach::foreach(commu = asb) %dopar% {

          if(nrow(commu) >= 2 & ncol(commu) >= 2 &
             all(commu == 1) == FALSE & all(commu == 0) == FALSE){ # Minimum 2 neig and 2 spp

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


            ## Calculating the PD contained on each community, but for the complete phylogenetic tree (CpBnest denominator)
            # For each commu
            pd_assemblages <- sapply(1:length(commus_nodes), function(x){
              return(sum(tree$edge.length[commus_nodes[[x]]]))
            })
            # For pairwised comparisions among assemblages
            pw_pd_assemblages <- sapply(1:length(p_commus_nodes), function(x){
              return(sum(tree$edge.length[p_commus_nodes[[x]]]))
            })


            ##################
            #### SORENSEN ####
            ##################

            # If the component desired is the SORENSEN
            if(comp == "sorensen"){
              if(method == "pairwise"){
                ## Capturing the community parameters a, b, c (CpB denominator)
                bc_sum <- apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, sum) # minimun (b,c)
                # capturing "a"
                a <- pw_pd_assemblages - bc_sum
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
                t_pd_assemblages <- sum(tree$edge.length[t_commus_nodes])

                ## Capturing the community parameters a, b, c (CpB denominator)
                bc_sum <- sum(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2))) # minimun (b,c)
                a <- sum(sum(pd_assemblages) - t_pd_assemblages)   # t_pd_assemblages/5  if no one branch lenght were shared among assemblages, it is the value expected. In that case, a will equal 0.
              }
            }

            # If the component desired is the SORENSEN, calculate the relative PBD
            # stored in each phylogenetic slice following the imputed method
            if(comp == "sorensen"){
              ## Calculating the PROPORTION of NESTEDNESS accumulation over time
              CpB <- sapply(branch_pieces, function(x){
                ## PDs for each phylogenetic tree slice
                # For each commu
                pd_piece <- sapply(1:length(commus_nodes), function(y){   # x <- branch_pieces[[1]]
                  return(sum(x$edge.length[commus_nodes[[y]]]))
                })
                # For pairwised comparisions among assemblages
                pd_pw_piece <- sapply(1:length(p_commus_nodes), function(y){
                  return(sum(x$edge.length[p_commus_nodes[[y]]]))
                })


                if(method == "pairwise"){
                  # Calculating the bc sum within the slices (relatice (b + c) on slice) -> numerator sorensen
                  r_bc_sum <- apply(pd_pw_piece - t(utils::combn(pd_piece, 2)), 1, sum)
                  # Relative Sorensen
                  r_pbd <- r_bc_sum/((2*a) + bc_sum)
                }

                if(method == "multisite"){
                  # Calculating the bc sum within the slices (relatice (b + c) on slice) -> numerator sorensen
                  r_bc_sum <- sum(pd_pw_piece - t(utils::combn(pd_piece, 2)))
                  # Relative Sorensen
                  r_pbd <- r_bc_sum/((2*a) + bc_sum)
                }

                # Capturing the proportion of turnover in stored in this slice   (r_pbd/pbd)
                return(r_pbd)
              })
            }


            ####################
            ##### TURNOVER #####
            ####################

            # If the component desired is the TURNOVER, calculate the minimum
            # following the inputed method
            if(comp == "turnover"){
              if(method == "pairwise"){
                ## Capturing the community parameters min (b, c) - (CpB denominator)
                # Now the minimum within the community
                min_bc <- apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, min)
                ## Capturing the community parameters a, b, c (CpB denominator)
                bc_sum <- apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, sum) # minimun (b,c)
                # capturing "a"
                a <- pw_pd_assemblages - bc_sum
                # Turnover (Simpson)
                turn <- mean(min_bc)
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
                t_pd_assemblages <- sum(tree$edge.length[t_commus_nodes])

                ## Capturing the community parameters a, b, c (CpB denominator)
                a <- sum(sum(pd_assemblages) - t_pd_assemblages)   # t_pd_assemblages/5  if no one branch lenght were shared among assemblages, it is the value expected. In that case, a will equal 0.
                # Capturing the multisite community minimum (CpB denominator)
                turn <- sum(apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, min))
              }
            }

            # If the component desired is the TURNOVER, calculate the relative minimum
            # stored in each phylogenetic slice following the imputed method
            if(comp == "turnover"){
              ## Calculating the PROPORTION of TURNOVER stored accumulation over time
              CpB <- sapply(branch_pieces, function(x){
                ## PDs for each phylogenetic tree slice
                # For each commu
                piece <- sapply(1:length(commus_nodes), function(y){   # x <- branch_pieces[[1]]
                  return(sum(x$edge.length[commus_nodes[[y]]]))
                })
                # For pairwised comparisions among assemblages
                pw_piece <- sapply(1:length(p_commus_nodes), function(y){
                  return(sum(x$edge.length[p_commus_nodes[[y]]]))
                })

                if(method == "pairwise"){
                  r_turn <- mean(apply(pw_piece - t(utils::combn(piece, 2)), 1, min))/(a + turn)
                }

                if(method == "multisite"){
                  r_turn <- sum(apply(pw_piece - t(utils::combn(piece, 2)), 1, min))/(a + turn)
                }

                # Capturing the proportion of turnover in stored in this slice
                return(r_turn)
              })
            }


            ####################
            #### NESTEDNESS ####
            ####################

            # If the component desired is the NESTEDNESS, decompose its other components
            # to obtain the nestedness values following the imputed method
            if(comp == "nestedness"){
              if(method == "pairwise"){
                ## Capturing the community parameters a, b, c (CpB denominator)
                bc_sum <- apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, sum) # minimun (b,c)
                # capturing "a"
                a <- pw_pd_assemblages - bc_sum
                # Now the minimum within the community
                min_bc <- apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, min)
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
                t_pd_assemblages <- sum(tree$edge.length[t_commus_nodes])

                ## Capturing the community parameters a, b, c (CpB denominator)
                bc_sum <- sum(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2))) # minimun (b,c)
                a <- sum(sum(pd_assemblages) - t_pd_assemblages)   # t_pd_assemblages/5  if no one branch lenght were shared among assemblages, it is the value expected. In that case, a will equal 0.
                # Now the minimum within the community
                min_bc <- sum(apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, min))
              }
            }

            # If the component desired is the NESTEDNESS, calculate the relative nest
            # stored in each phylogenetic slice following the imputed method
            if(comp == "nestedness"){
              ## Calculating the PROPORTION of NESTEDNESS accumulation over time
              CpB <- sapply(branch_pieces, function(x){
                ## PDs for each phylogenetic tree slice
                # For each commu
                pd_piece <- sapply(1:length(commus_nodes), function(y){   # x <- branch_pieces[[1]]
                  return(sum(x$edge.length[commus_nodes[[y]]]))
                })
                # For pairwised comparisions among assemblages
                pd_pw_piece <- sapply(1:length(p_commus_nodes), function(y){
                  return(sum(x$edge.length[p_commus_nodes[[y]]]))
                })


                if(method == "pairwise"){
                  # Calculating the bc sum within the slices (relatice (b + c) on slice) -> numerator sorensen
                  r_bc_sum <- apply(pd_pw_piece - t(utils::combn(pd_piece, 2)), 1, sum)
                  # Calculating the minimum within the slices (relative min(b,c) on slice) -> numerator simpson
                  r_min_bc <- apply(pd_pw_piece - t(utils::combn(pd_piece, 2)), 1, min)

                  ## Mean relative value
                  # Relative Sorensen
                  r_pbd <- r_bc_sum/((2*a) + bc_sum)
                  # Relative Turnover (Simpson)
                  r_turn <- r_min_bc/(a + min_bc)
                  # Relative Nestedness
                  r_nest <- mean(r_pbd - r_turn)
                }

                if(method == "multisite"){
                  # Calculating the bc sum within the slices (relatice (b + c) on slice) -> numerator sorensen
                  r_bc_sum <- sum(pd_pw_piece - t(utils::combn(pd_piece, 2)))
                  # Calculating the minimum within the slices (relative min(b,c) on slice) -> numerator simpson
                  r_min_bc <- sum(apply(pd_pw_piece - t(utils::combn(pd_piece, 2)), 1, min))

                  # Relative Sorensen
                  r_pbd <- r_bc_sum/((2*a) + bc_sum)
                  # Relative Turnover (Simpson)
                  r_turn <- r_min_bc/(a + min_bc)
                  # Relative Nestedness
                  r_nest <- r_pbd - r_turn
                }
                # Return
                return(r_nest)
              })
            }

            # Returning
            return(CpB)
          }
        }

        # stop the clusters after running the algorithm
        doParallel::stopImplicitCluster()
      }

      # The user do not set clusters for being used
      if(ncor == 0){
        # Calculating the CpB rate through the CpBrate algorithm
        CpB <- foreach::foreach(commu = asb) %do% {

          if(nrow(commu) >= 2 & ncol(commu) >= 2 &
             all(commu == 1) == FALSE & all(commu == 0) == FALSE){ # Minimum 2 neig and 2 spp

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


            ## Calculating the PD contained on each community, but for the complete phylogenetic tree (CpBnest denominator)
            # For each commu
            pd_assemblages <- sapply(1:length(commus_nodes), function(x){
              return(sum(tree$edge.length[commus_nodes[[x]]]))
            })
            # For pairwised comparisions among assemblages
            pw_pd_assemblages <- sapply(1:length(p_commus_nodes), function(x){
              return(sum(tree$edge.length[p_commus_nodes[[x]]]))
            })


            ##################
            #### SORENSEN ####
            ##################

            # If the component desired is the SORENSEN
            if(comp == "sorensen"){
              if(method == "pairwise"){
                ## Capturing the community parameters a, b, c (CpB denominator)
                bc_sum <- apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, sum) # minimun (b,c)
                # capturing "a"
                a <- pw_pd_assemblages - bc_sum
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
                t_pd_assemblages <- sum(tree$edge.length[t_commus_nodes])

                ## Capturing the community parameters a, b, c (CpB denominator)
                bc_sum <- sum(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2))) # minimun (b,c)
                a <- sum(sum(pd_assemblages) - t_pd_assemblages)   # t_pd_assemblages/5  if no one branch lenght were shared among assemblages, it is the value expected. In that case, a will equal 0.
              }
            }

            # If the component desired is the SORENSEN, calculate the relative PBD
            # stored in each phylogenetic slice following the imputed method
            if(comp == "sorensen"){
              ## Calculating the PROPORTION of NESTEDNESS accumulation over time
              CpB <- sapply(branch_pieces, function(x){
                ## PDs for each phylogenetic tree slice
                # For each commu
                pd_piece <- sapply(1:length(commus_nodes), function(y){   # x <- branch_pieces[[1]]
                  return(sum(x$edge.length[commus_nodes[[y]]]))
                })
                # For pairwised comparisions among assemblages
                pd_pw_piece <- sapply(1:length(p_commus_nodes), function(y){
                  return(sum(x$edge.length[p_commus_nodes[[y]]]))
                })


                if(method == "pairwise"){
                  # Calculating the bc sum within the slices (relatice (b + c) on slice) -> numerator sorensen
                  r_bc_sum <- apply(pd_pw_piece - t(utils::combn(pd_piece, 2)), 1, sum)
                  # Relative Sorensen
                  r_pbd <- r_bc_sum/((2*a) + bc_sum)
                }

                if(method == "multisite"){
                  # Calculating the bc sum within the slices (relatice (b + c) on slice) -> numerator sorensen
                  r_bc_sum <- sum(pd_pw_piece - t(utils::combn(pd_piece, 2)))
                  # Relative Sorensen
                  r_pbd <- r_bc_sum/((2*a) + bc_sum)
                }

                # Capturing the proportion of turnover in stored in this slice   (r_pbd/pbd)
                return(r_pbd)
              })
            }


            ####################
            ##### TURNOVER #####
            ####################

            # If the component desired is the TURNOVER, calculate the minimum
            # following the inputed method
            if(comp == "turnover"){
              if(method == "pairwise"){
                ## Capturing the community parameters min (b, c) - (CpB denominator)
                # Now the minimum within the community
                min_bc <- apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, min)
                ## Capturing the community parameters a, b, c (CpB denominator)
                bc_sum <- apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, sum) # minimun (b,c)
                # capturing "a"
                a <- pw_pd_assemblages - bc_sum
                # Turnover (Simpson)
                turn <- mean(min_bc)
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
                t_pd_assemblages <- sum(tree$edge.length[t_commus_nodes])

                ## Capturing the community parameters a, b, c (CpB denominator)
                a <- sum(sum(pd_assemblages) - t_pd_assemblages)   # t_pd_assemblages/5  if no one branch lenght were shared among assemblages, it is the value expected. In that case, a will equal 0.
                # Capturing the multisite community minimum (CpB denominator)
                turn <- sum(apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, min))
              }
            }

            # If the component desired is the TURNOVER, calculate the relative minimum
            # stored in each phylogenetic slice following the imputed method
            if(comp == "turnover"){
              ## Calculating the PROPORTION of TURNOVER stored accumulation over time
              CpB <- sapply(branch_pieces, function(x){
                ## PDs for each phylogenetic tree slice
                # For each commu
                piece <- sapply(1:length(commus_nodes), function(y){   # x <- branch_pieces[[1]]
                  return(sum(x$edge.length[commus_nodes[[y]]]))
                })
                # For pairwised comparisions among assemblages
                pw_piece <- sapply(1:length(p_commus_nodes), function(y){
                  return(sum(x$edge.length[p_commus_nodes[[y]]]))
                })

                if(method == "pairwise"){
                  r_turn <- mean(apply(pw_piece - t(utils::combn(piece, 2)), 1, min))/(a + turn)
                }

                if(method == "multisite"){
                  r_turn <- sum(apply(pw_piece - t(utils::combn(piece, 2)), 1, min))/(a + turn)
                }

                # Capturing the proportion of turnover in stored in this slice
                return(r_turn)
              })
            }


            ####################
            #### NESTEDNESS ####
            ####################

            # If the component desired is the NESTEDNESS, decompose its other components
            # to obtain the nestedness values following the imputed method
            if(comp == "nestedness"){
              if(method == "pairwise"){
                ## Capturing the community parameters a, b, c (CpB denominator)
                bc_sum <- apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, sum) # minimun (b,c)
                # capturing "a"
                a <- pw_pd_assemblages - bc_sum
                # Now the minimum within the community
                min_bc <- apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, min)
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
                t_pd_assemblages <- sum(tree$edge.length[t_commus_nodes])

                ## Capturing the community parameters a, b, c (CpB denominator)
                bc_sum <- sum(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2))) # minimun (b,c)
                a <- sum(sum(pd_assemblages) - t_pd_assemblages)   # t_pd_assemblages/5  if no one branch lenght were shared among assemblages, it is the value expected. In that case, a will equal 0.
                # Now the minimum within the community
                min_bc <- sum(apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, min))
              }
            }

            # If the component desired is the NESTEDNESS, calculate the relative nest
            # stored in each phylogenetic slice following the imputed method
            if(comp == "nestedness"){
              ## Calculating the PROPORTION of NESTEDNESS accumulation over time
              CpB <- sapply(branch_pieces, function(x){
                ## PDs for each phylogenetic tree slice
                # For each commu
                pd_piece <- sapply(1:length(commus_nodes), function(y){   # x <- branch_pieces[[1]]
                  return(sum(x$edge.length[commus_nodes[[y]]]))
                })
                # For pairwised comparisions among assemblages
                pd_pw_piece <- sapply(1:length(p_commus_nodes), function(y){
                  return(sum(x$edge.length[p_commus_nodes[[y]]]))
                })


                if(method == "pairwise"){
                  # Calculating the bc sum within the slices (relatice (b + c) on slice) -> numerator sorensen
                  r_bc_sum <- apply(pd_pw_piece - t(utils::combn(pd_piece, 2)), 1, sum)
                  # Calculating the minimum within the slices (relative min(b,c) on slice) -> numerator simpson
                  r_min_bc <- apply(pd_pw_piece - t(utils::combn(pd_piece, 2)), 1, min)

                  ## Mean relative value
                  # Relative Sorensen
                  r_pbd <- r_bc_sum/((2*a) + bc_sum)
                  # Relative Turnover (Simpson)
                  r_turn <- r_min_bc/(a + min_bc)
                  # Relative Nestedness
                  r_nest <- mean(r_pbd - r_turn)
                }

                if(method == "multisite"){
                  # Calculating the bc sum within the slices (relatice (b + c) on slice) -> numerator sorensen
                  r_bc_sum <- sum(pd_pw_piece - t(utils::combn(pd_piece, 2)))
                  # Calculating the minimum within the slices (relative min(b,c) on slice) -> numerator simpson
                  r_min_bc <- sum(apply(pd_pw_piece - t(utils::combn(pd_piece, 2)), 1, min))

                  # Relative Sorensen
                  r_pbd <- r_bc_sum/((2*a) + bc_sum)
                  # Relative Turnover (Simpson)
                  r_turn <- r_min_bc/(a + min_bc)
                  # Relative Nestedness
                  r_nest <- r_pbd - r_turn
                }
                # Return
                return(r_nest)
              })
            }

            # Returning
            return(CpB)
          }
        }
      }

      # Returning
      return(CpB)
    }

    # Relative PB_RW on each phylogenetic slice ----------------------------------
    if(index == "PB_RW"){

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
          return(mat[which(adj[x,] == 1), , drop = FALSE])
        }
      })

      # Cutting the phylogenetic tree into equal width slices
      branch_pieces <- phylo_pieces(tree, n, criterion = criterion,
                                    timeSteps = TRUE, returnTree = TRUE)

      # Separating the time steps from the phylogenetic pieces
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
        CpB <- foreach::foreach(commu = asb) %dopar% {  # commu <- communities[[1]] asb <- communities commu <- communities[[8]]

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
              return(r_pbd)
            })

            # Return CpB
            return(CpB)

          } else {
            CpB <- NA
            return(CpB)
          }
        }

        # stop the clusters after running the algorithm
        doParallel::stopImplicitCluster()
      }

      # The user do not set clusters for being used
      if(ncor == 0){
        # Calculating the CpB rate through the CpBrate algorithm
        CpB <- foreach::foreach(commu = asb) %do% {  # commu <- communities[[1]] asb <- communities commu <- communities[[8]]

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
              return(r_pbd)
            })

            # Return CpB
            return(CpB)

          } else {
            CpB <- NA
            return(CpB)
          }
        }
      }

      # Return
      return(CpB)
    }
  }
}
