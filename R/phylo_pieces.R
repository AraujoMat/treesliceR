#' Slices a phylogenetic tree into multiple temporal slices
#' @description
#' This function slices a phylogenetic tree into multiple slices, spaced equally either in million years or intervals of phylogenetic diversity (PD).
#'
#' @usage phylo_pieces(tree, n, criteria, method, timeSteps, dropNodes, returnTree)
#'
#' @param tree phylo. An ultrametric phylogenetic tree in the "phylo" format.
#' @param n numeric. A numeric value indicating either the number of temporal slices (method = 1) or the time interval in million years (or phylogenetic diversity) among the tree slices (method = 2). Default is 1.
#' @param criteria character string. The method for slicing the tree. It can be either "my" (million years) or "PD" (accumulated phylogenetic diversity). Default is "my".
#' @param method numerical. A numerical value indicating the method to make the multiple slices. Setting "method = 1" will slice the phylogeny based on an "n" number of slices. If method = 2, the slices will be created based on a temporal interval. Default is 1.
#' @param timeSteps logical. A logical value indicating whether the vector containing the time-steps used for creating the multiple slices should be returned. If "timeSteps = TRUE", then it returns a list containing both the time vector and a list with the multiple tree slices. Default is FALSE.
#' @param dropNodes logical. A logical value indicating whether the nodes that were sliced (void nodes, presenting no branch length) should be preserved in the node matrix. Default is FALSE.
#' @param returnTree logical. A logical value indicating whether the original input tree should be returned with its slices in a list. Default is FALSE.
#'
#' @return The function returns a list containing multiple slices of a phylogenetic tree. The slices list is ordered from roots to tips. Thus, the first object within the outputted list is the root-slice, whereas the last is the tips-slice.
#'
#' @seealso Other slicing methods: [squeeze_root()], [squeeze_tips()], [squeeze_int()], [prune_tips()]
#'
#' @author Matheus Lima de Araujo <matheusaraujolima@live.com>
#'
#' @examples
#' # Generate a random tree
#' tree <- ape::rcoal(20)
#'
#' # Cuts a phylogeny into multiple temporal slices
#' tree <- phylo_pieces(tree, n = 3, criteria = "my", method = 1)
#'
#' # Plotting the three slices of our phylogeny
#' plot(tree[[1]])
#' plot(tree[[2]])
#' plot(tree[[3]])
#'
#' @export

phylo_pieces <- function(tree, n, criteria = "my", method = 1,
                         timeSteps = FALSE, dropNodes = FALSE, returnTree = FALSE){

  # Capturing the nodes and tips informations
  tree <- nodes_config(tree)

  # The choose criteria for making the slices are time
  if(criteria == "my"){

    if(method == 1){  # number of phylogenetic slices
      # If the number of slices provided were decimals...
      if(n < 1){
        stop("The n must be an integer under argument method = 1")
      } else {
        # Saving the number of slices into a new object
        nslices <- n

        # Capturing the time interval that it belongs
        n <- tree$tree_depth/nslices

        # Printing the threshold choose
        message(paste("> The", nslices, "number of pieces inputted equals to intervals of",
                  n, "million of years.\n"))
      }


    }

    if(method == 2){ # desired time interval among slices

      # if the imputed time is bigger than the tree time
      if(n > sum(tree$tree_depth)){
        stop("The time threshold inputted is older than the phylogeny")
      } else {

        ## Find the best time interval that equal the temporal width among slices
        # Finding the expected number of pieces under the imputed time-interval
        exp_pieces <- tree$tree_depth/n

        # Find the smallest rounding value from the expected
        f_piece <- tree$tree_depth/floor(exp_pieces)
        # Find the biggest rounding values from the expected
        cel_piece <- tree$tree_depth/ceiling(exp_pieces)

        # Check which of these are the closest value and select it
        if(abs(n - f_piece) < abs(n - cel_piece) & f_piece != tree$tree_depth){
          n <- f_piece # floor piece
          nslices <- floor(exp_pieces)
        } else {
          n <- cel_piece # ceiling piece
          nslices <- ceiling(exp_pieces)
        }

        # Printing the threshold choose
        message(paste("> The number of pieces was rounded to", tree$tree_depth/n,
                  ", with intervals of", n, "million of years.\n"))
      }
    }

    # Summing the cummulated values of the temporal slices (from roots -> tips)
    age <- cumsum(rep(n, nslices))

    ## Running the slicer algorithm
    cutted_tree <- lapply(age, function(j){ # j <- age[4]

      ## Cutting the phylogeny TIPWARDLY (TIPS -> ROOT):
      # Which nodes ends after j but are SMALLER or bigger than my thresholds?
      smaller <- which((tree$config$YearEnd[which(tree$config$YearEnd > j)] - j) >= tree$edge.length[which(tree$config$YearEnd > j)])
      bigger <- which((tree$config$YearEnd[which(tree$config$YearEnd > j)] - j) < tree$edge.length[which(tree$config$YearEnd > j)])

      # For those smaller, turn their edge.length 0 (ZERO)
      tree$edge.length[which(tree$config$YearEnd > j)][smaller] <- 0

      # Those values bigger than my threshold when subtracting j,
      # remove its remaining difference to reach the threshold
      tree$edge.length[which(tree$config$YearEnd > j)][bigger] <-
        tree$edge.length[which(tree$config$YearEnd > j)][bigger] - (tree$config$YearEnd[which(tree$config$YearEnd > j)][bigger] - j)

      if(j > n){  ## if j is smaller than n, make the slices specifically on nodes position

        ### Cutting the phylogeny ROOTWARDLY (ROOT -> TIPS):
        # The threshold occupy different nodes with different requirements for node slicing?
        # (ex: several nodes starting before and after a given threshold)
        if(length(unique(c(which(sort(unique(tree$config$YearEnd)) >= j)[1], which(sort(unique(tree$config$YearEnd)) >= (j-n))[1]))) > 1){

          ## Correcting the length of those nodes that ends inside the given interval,
          # but start before the threshold established in (j - n):
          tree$edge.length[tree$config$YearEnd >= (j-n) & tree$config$YearEnd <= (j) & tree$config$YearBegin < (j-n)] <-
            (tree$config$YearEnd[tree$config$YearEnd >= (j-n) & tree$config$YearEnd <= (j) & tree$config$YearBegin < (j-n)] - (j-n))

          ## All nodes ending before my threshold (j-n), turn into zero their lengths
          tree$edge.length[tree$config$YearEnd < (j-n)] <- 0

          ## All nodes originating befores my threshold (j - n) and ending after my threshold (j)
          # (an entire branch length within the interval), assign (n) to its edge.length
          tree$edge.length[tree$config$YearBegin < (j-n) & tree$config$YearEnd >= j] <- n
        }

        ## If all remaining nodes had its edges bigger than the interval:
        if(length(unique(c(which(sort(unique(tree$config$YearEnd)) >= j)[1], which(sort(unique(tree$config$YearEnd)) >= (j-n))[1]))) == 1){

          ## Turn into 0 those edges ending before my threshold
          tree$edge.length[which(tree$config$YearEnd < (j))] <- 0

          ## Those remaining nodes, turn into n their edge.lengths
          tree$edge.length[tree$edge.length > 0] <- n
        }
        return(tree)
      }
      ## If the j threshold is smaller than n, return the tree as it is
      if(j <= n){
        return(tree)
      }
    })
  }

  # The choose criteria for making the slices are PD
  if(criteria == "pd"){
    ## How many PD we have at each nodes spliting?
    # Separating only the nodes inital information
    nodes <- tree$config[!duplicated(tree$config$NodeBegin), c(1, 2, 4)]

    ## Making a matrix with these informations
    df <- data.frame(time = sort(unique(nodes$YearBegin)), # Time which init the node
                     nBranch = 2:(length(unique(nodes$YearBegin)) + 1), # Number of branches
                     timeLength = c(sort(unique(nodes$YearBegin)), max(tree$config$YearEnd))[-1] - # time length of the node
                       sort(unique(nodes$YearBegin)))

    # Calculating the cumulative TIME and PD per node split
    df$cumulativeTIME <- cumsum(df$timeLength)
    df$cumulativePD <- cumsum(df$nBranch * df$timeLength)


    if(method == 1){ # number of phylogenetic slices
      # If the number of slices provided were decimals...
      if(n < 1){
        stop("The n must be an integer under the argument method = 1")
      } else {
        # Saving the number of slices into a new object
        nslices <- n

        # Capturing the time interval that it belongs
        n <- sum(tree$edge.length)/nslices

        # Printing the threshold choose
        message(paste("> The", nslices, "number of pieces inputted equals to intervals of",
                  n, "phylogenetic diversity (PD).\n"))
      }
    }

    if(method == 2){ # desired PD interval among slices

      # if the imputed PD is bigger than the tree PD
      if(n > sum(tree$edge.length)){
        stop("The PD threshold inputted is older than the phylogeny")
      } else {

        ## Find the best PD interval that equal the temporal width among slices
        # Finding the expected number of pieces under the imputed PD-interval
        exp_pieces <- sum(tree$edge.length)/n

        # Find the smallest rounding value from the expected
        f_piece <- sum(tree$edge.length)/floor(exp_pieces)
        # Find the biggest rounding values from the expected
        cel_piece <- sum(tree$edge.length)/ceiling(exp_pieces)

        # Check which of these are the closest value and select it
        if(abs(n - f_piece) < abs(n - cel_piece) & f_piece != sum(tree$edge.length)){
          n <- f_piece # floor piece
          nslices <- floor(exp_pieces)
        } else {
          n <- cel_piece # ceiling piece
          nslices <- ceiling(exp_pieces)
        }

        # Printing the threshold choose
        message(paste("> The number of pieces was rounded to", sum(tree$edge.length)/n,
                  ", with intervals of", n, "phylogenetic diversity (PD).\n"))
      }
    }

    # Qual o PD dentro de cada intervalo desejado ?
    pd_slices <- sum(tree$edge.length)/nslices
    # Calculating the cumulative PD
    pd_slices <- cumsum(rep(pd_slices, nslices))

    # Retrieving the temporal slices of equal PD
    # Running the algorithm
    age <- sapply(pd_slices, function(z){
      # Separating the node information
      node_info <- df[which(df$cumulativePD >= z)[1],]
      # What is the PD difference within this node and the expected PD?
      rmv <- (node_info[, 5] - z)
      # Which is the temporal distance of this difference (based on the local node number)
      rmv <- rmv/node_info[, 2]
      # What time this expected PD are located in?
      time_pd <- node_info[, 4] - rmv
      return(time_pd)
    })

    # Configuring the time intervals into a data.frame
    pd_df <- data.frame(ID = 1:length(age), AGE_t1 = age, AGE_t0 = c(0, age[-length(age)]))


    ## Running the slicer algorithm
    cutted_tree <- lapply(pd_df[,2], function(k){ # k <- pd_df[2, 2]

      # Capturing the time of left-side slice
      l <- pd_df[which(pd_df[, 2] == k), 3]

      ## Cutting the phylogeny TIPWARDLY (TIPS -> ROOT):
      # Which nodes ends after j but are SMALLER or bigger than my thresholds?
      smaller <- which((tree$config$YearEnd[which(tree$config$YearEnd > k)] - k) >= tree$edge.length[which(tree$config$YearEnd > k)])
      bigger <- which((tree$config$YearEnd[which(tree$config$YearEnd > k)] - k) < tree$edge.length[which(tree$config$YearEnd > k)])

      # For those smaller, turn their edge.length 0 (ZERO)
      tree$edge.length[which(tree$config$YearEnd > k)][smaller] <- 0

      # Those values bigger than my threshold when subtracting k,
      # remove its remaining difference to reach the threshold
      tree$edge.length[which(tree$config$YearEnd > k)][bigger] <-
        tree$edge.length[which(tree$config$YearEnd > k)][bigger]-(tree$config$YearEnd[which(tree$config$YearEnd > k)][bigger] - k)

      if(l > 0){  ## if l is smaller than 0, make the slices specifically on nodes position

        ### Cutting the phylogeny ROOTWARDLY (ROOT -> TIPS):
        # The threshold occupy different nodes with different requirements for node slicing?
        # (ex: several nodes starting before and after a given threshold)
        if(length(unique(c(which(sort(unique(tree$config$YearEnd)) >= k)[1], which(sort(unique(tree$config$YearEnd)) >= l)[1]))) > 1){ # + de 1 TRUE

          ## Correcting the length of those nodes that ends inside the given interval,
          # but start before the threshold established in (k - l):
          tree$edge.length[tree$config$YearEnd >= l & tree$config$YearEnd <= (k) & tree$config$YearBegin < l] <-
            (tree$config$YearEnd[tree$config$YearEnd >= l & tree$config$YearEnd <= (k) & tree$config$YearBegin < l] - l)

          ## All nodes ending before my threshold (l), turn into zero their lengths
          tree$edge.length[tree$config$YearEnd < l] <- 0

          ## All nodes originating befores my threshold (l) and ending after my threshold (k)
          # (an entire branch length within the interval), assign (k - l) to its edge.length
          tree$edge.length[tree$config$YearBegin < l & tree$config$YearEnd >= k] <- k-l
        }

        ## If all remaining nodes had its edges bigger than the interval:
        if(length(unique(c(which(sort(unique(tree$config$YearEnd)) >= k)[1], which(sort(unique(tree$config$YearEnd)) >= l)[1]))) == 1){ # ! de 0, porque ele n?o foi cortado; o bra?o ?ntegro (sem recorte) que o valor de recorte

          ## Turn into 0 those edges ending before my threshold
          tree$edge.length[which(tree$config$YearEnd < (k))] <- 0

          ## Those remaining nodes, turn into k-l their edge.lengths
          tree$edge.length[tree$edge.length > 0] <- k-l
        }
        return(tree)
      }
      if(l == 0){
        return(tree)
      }
    })
  }

  ## If there are some void nodes/edges inside our tree pieces,
  # and the user want to remove them from each phylo piece
  if(dropNodes == TRUE){

    cutted_tree <- lapply(cutted_tree, function(tree_pieces){
      # Which are our void nodes with 0 length?
      rm_nd <- which(!(tree_pieces$edge[,2] %in% c(1:length(tree_pieces$tip.label))) & tree_pieces$edge.length == 0) # tree <- teste

      if(length(rm_nd) > 0){
        # Correcting their values
        for (i in 1:length(rm_nd)) { # i <- 1
          # which node comes after our node
          val <- tree_pieces$edge[rm_nd[i], 2]
          # all places with this value, need to be turned into its previous node
          tree_pieces$edge[which(tree_pieces$edge[,1] == val), 1] <- tree_pieces$edge[rm_nd[i], 1]
        }

        # Them, removing those edges and edgelenghts that are 0 and are node-to-node
        tree_pieces$edge <- tree_pieces$edge[-c(rm_nd), ]               # Removing from the edge matrix
        tree_pieces$edge.length <- tree_pieces$edge.length[-c(rm_nd)]   # Removing from the edge length vectors
        tree_pieces$Nnode <- length(unique(tree_pieces$edge[,1]))       # Adding the new number of nodes

        # Unique tree nodes
        oldnodes <- sort(unique(tree_pieces$edge[,1]))

        # Renaming the tree nodes
        for (i in 1:length(oldnodes)) { # i <- 1
          tree_pieces$edge[which(tree_pieces$edge[,1] == oldnodes[i]), 1] <- length(tree_pieces$tip.label) + i
          tree_pieces$edge[which(tree_pieces$edge[,2] == oldnodes[i]), 2] <- length(tree_pieces$tip.label) + i
        }
      }

      return(tree_pieces)
    })
  }

  # If the user want the time-steps used to make the slieces
  if(timeSteps == TRUE){

    # If the user want to return the original tree within the list
    if(returnTree == FALSE){
      # Return a list with all pieces and a time-step vector
      return(list(cutted_tree, age))
    } else {
      # Return a list with all pieces, a time-step vector, and the original tree
      return(list(cutted_tree, age, tree))
    }

  } else {

    # If the user want to return the original tree within the list
    if(returnTree == FALSE){
      # return only the phylogenetic pieces
      return(cutted_tree)
    } else {
      # Return a list with all pieces, a time-step vector, and the original tree
      return(list(cutted_tree, tree))
    }
  }
}
