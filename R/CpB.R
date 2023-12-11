#' Calculates the rate of accumulation of phylogenetic B-diversity (CpB) over time slices
#' @description
#' This function estimates the rates of accumulation of phylogenetic B-diversity (CpB) over time for inputted assemblages.
#'
#' @usage CpB(tree, n, asb, comp, method, criteria, pBO, ncor)
#'
#' @param tree phylo. An ultrametric phylogenetic tree in the "phylo" format.
#' @param n numeric. A numeric value indicating the number of temporal slices (method = 1) or the time interval in million years (or phylogenetic diversity) among the tree slices (method = 2). Default is 1.
#' @param asb matrix, or list of matrices. A matrix (or list of matrices) containing a focal assemblage and its neighborhood assemblages (need at least two assemblages to run).
#' @param comp character string. The component of the phylogenetic beta-diversity to obtain the rates of accumulation. It can be either "sorensen", "turnover", or "nestedness". Default is "sorensen".
#' @param method character string. The method for calculating the phylogenetic beta-diversity. It can be either obtained through a "pairwise" or "multisite" approach. Default is "multisite".
#' @param criteria character string. The method for cutting the tree. It can be either "my" (million years) or "PD" (accumulated phylogenetic diversity). Default is "my".
#' @param pBO numeric. A value indicating the numeric proportion to define the temporal origin at which the phylogenetic B-diversity (PB) started to accumulate in a given assemblage. Default is 5%.
#' @param ncor numeric. A value indicating the number of cores the user wants to parallelize. Default is 0.
#'
#' @return The function returns a data frame containing the assemblages' rates of cumulative phylogenetic B-diversity (CpB), their total phylogenetic B-diversity (PB), and their PB origin (pBO).
#'
#' @details
#'
#' \bold{Parallelization}
#'
#' Users are advised to check the number of available cores within their machines before running parallel programming.
#'
#' @seealso Other cumulative phylogenetic index analysis: [CpD()], [CpE()], [CpB_RW()]
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
#' # And separate it into two assemblages with focal and neigs
#' asb <- list(mat[1:5,], mat[6:10,])
#'
#' # Calculate their CpB (sorensen) for 100 tree slices
#' CpB(tree, n = 100, asb = asb, comp = "sorensen", method = "multisite")
#'
#' @export

CpB <- function(tree, n, asb, comp = "sorensen",
                method = "multisite", criteria = "my", pBO = 5, ncor = 0){

  # Capturing the tree configs and
  # Cutting the phylogenetic tree into equal width slices
  branch_pieces <- phylo_pieces(tree, n, criteria = criteria,
                                timeSteps = TRUE, returnTree = TRUE)

  # Separating the time steps from the phylogenetic pieces
  age <- branch_pieces[[2]][length(branch_pieces[[2]]):1]
  tree <- branch_pieces[[3]]
  branch_pieces <- branch_pieces[[1]]
  commu <- NULL

  # Checking if a single matrix or a list of matrixes was provided
  if(sum(class(asb) != "list") >= 1){
    asb <- list(asb)  # Transforming the list into a single matrix
  }

  # The user wants to use more CPU cores?
  if(ncor > 0){

    # Register the number of desired clusters
    doParallel::registerDoParallel(ncor)

    # Calculating the CpB rate through the CpBrate algorithm
    CpBrate <- foreach::foreach(commu = asb, .combine = rbind) %dopar% {

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

            # Sorensen
            pbd <- bc_sum/((2*a) + bc_sum)
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

            # Sorensen
            pbd <- bc_sum/((2*a) + bc_sum)
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
            return(r_pbd/pbd)
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
            g_turn <- mean(min_bc/(a + min_bc))

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
            bc_sum <- sum(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2))) # minimun (b,c)
            a <- sum(sum(pd_assemblages) - t_pd_assemblages)   # t_pd_assemblages/5  if no one branch lenght were shared among assemblages, it is the value expected. In that case, a will equal 0.
            # Capturing the multisite community minimum (CpB denominator)
            turn <- sum(apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, min))

            # Turnover (Simpson)
            g_turn <- turn/(a + turn)
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
              r_turn <- mean(apply(pw_piece - t(utils::combn(piece, 2)), 1, min))/turn
            }

            if(method == "multisite"){
              r_turn <- sum(apply(pw_piece - t(utils::combn(piece, 2)), 1, min))/turn
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

            # Sorensen
            pbd <- bc_sum/((2*a) + bc_sum)
            # Turnover (Simpson)
            turn <- min_bc/(a + min_bc)
            # Nestedness
            nest <- mean(pbd - turn)
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

            # Sorensen
            pbd <- bc_sum/((2*a) + bc_sum)
            # Turnover (Simpson)
            turn <- min_bc/(a + min_bc)
            # Nestedness
            nest <- pbd - turn
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

            # Capturing the proportion of turnover in stored in this slice   (r_nest/nest)
            return(r_nest/nest)
          })
        }


        ####################
        ##### CpB rate #####
        ####################

        # Fitting the non-linear model to obtain the CpB-rate:
        if(is.na(CpB[1]) == FALSE){ # Checking if there is any NA emerging due 0/0 fraction
          nonlinear_cum <- stats::nls(cumsum(CpB) ~ exp(-r*age),
                                      start = list(r = 0.2),
                                      control = stats::nls.control(maxiter = 1000))   # plot(cumsum(CpB) ~ age)

          if(comp == "sorensen"){
            return(c(stats::coef(nonlinear_cum), pbd, ((-1)*(log(pBO/100)/stats::coef(nonlinear_cum)))))
          }
          if(comp == "turnover"){
            return(c(stats::coef(nonlinear_cum), g_turn, ((-1)*(log(pBO/100)/stats::coef(nonlinear_cum)))))
          }
          if(comp == "nestedness"){
            return(c(stats::coef(nonlinear_cum), nest, ((-1)*(log(pBO/100)/stats::coef(nonlinear_cum)))))
          }

        } else {
          return(c(NA, NA, NA))
        }
      } else {
        return(c(NA, NA, NA))
      }
    }

    # stop the clusters after running the algorithm
    doParallel::stopImplicitCluster()
  }

  # The user do not set clusters for being used
  if(ncor == 0){
    # Calculating the CpB rate through the CpBrate algorithm
    CpBrate <- foreach::foreach(commu = asb, .combine = rbind) %do% {

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

            # Sorensen
            pbd <- bc_sum/((2*a) + bc_sum)
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

            # Sorensen
            pbd <- bc_sum/((2*a) + bc_sum)
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
            return(r_pbd/pbd)
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
            g_turn <- mean(min_bc/(a + min_bc))

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
            bc_sum <- sum(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2))) # minimun (b,c)
            a <- sum(sum(pd_assemblages) - t_pd_assemblages)   # t_pd_assemblages/5  if no one branch lenght were shared among assemblages, it is the value expected. In that case, a will equal 0.
            # Capturing the multisite community minimum (CpB denominator)
            turn <- sum(apply(pw_pd_assemblages - t(utils::combn(pd_assemblages, 2)), 1, min))

            # Turnover (Simpson)
            g_turn <- turn/(a + turn)
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
              r_turn <- mean(apply(pw_piece - t(utils::combn(piece, 2)), 1, min))/turn
            }

            if(method == "multisite"){
              r_turn <- sum(apply(pw_piece - t(utils::combn(piece, 2)), 1, min))/turn
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

            # Sorensen
            pbd <- bc_sum/((2*a) + bc_sum)
            # Turnover (Simpson)
            turn <- min_bc/(a + min_bc)
            # Nestedness
            nest <- mean(pbd - turn)
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

            # Sorensen
            pbd <- bc_sum/((2*a) + bc_sum)
            # Turnover (Simpson)
            turn <- min_bc/(a + min_bc)
            # Nestedness
            nest <- pbd - turn
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

            # Capturing the proportion of turnover in stored in this slice   (r_nest/nest)
            return(r_nest/nest)
          })
        }


        ####################
        ##### CpB rate #####
        ####################

        # Fitting the non-linear model to obtain the CpB-rate:
        if(is.na(CpB[1]) == FALSE){ # Checking if there is any NA emerging due 0/0 fraction
          nonlinear_cum <- stats::nls(cumsum(CpB) ~ exp(-r*age),
                                      start = list(r = 0.2),
                                      control = stats::nls.control(maxiter = 1000))   # plot(cumsum(CpB) ~ age)

          if(comp == "sorensen"){
            return(c(stats::coef(nonlinear_cum), pbd, ((-1)*(log(pBO/100)/stats::coef(nonlinear_cum)))))
          }
          if(comp == "turnover"){
            return(c(stats::coef(nonlinear_cum), g_turn, ((-1)*(log(pBO/100)/stats::coef(nonlinear_cum)))))
          }
          if(comp == "nestedness"){
            return(c(stats::coef(nonlinear_cum), nest, ((-1)*(log(pBO/100)/stats::coef(nonlinear_cum)))))
          }

        } else {
          return(c(NA, NA, NA))
        }
      } else {
        return(c(NA, NA, NA))
      }
    }
  }


  # Renaming the columns
  # Renaming the columns and rows, and returning them as output
  if(length(asb) > 1){
    colnames(CpBrate) <- c("CpB", "PB", "pBO")
    rownames(CpBrate) <- 1:length(asb)
  } else {
    names(CpBrate) <- c("CpB", "PB", "pBO")
  }

  return(CpBrate)
}
