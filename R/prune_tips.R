#' Collapse the tips of a phylogenetic tree based on a temporal threshold
#' @description
#' This function collapses, or prune, the tips of a phylogenetic tree based on an inputted temporal threshold.
#'
#' @usage prune_tips(tree, time, qtl = FALSE, method = 1)
#'
#' @param tree phylo. An ultrametric phylogenetic tree in the "phylo" format.
#' @param time numeric or numeric vector. A numeric value, or vector, containing the temporal threshold(s) in million years for pruning the tree tips. It must be within the interval of the tips ages present in the phylogenetic tree.
#' @param qtl logical. A logical value indicating whether the user wants to use the quantile values of the tip ages to prune the tree. Default is FALSE
#' @param method numerical. A numerical value indicating the method to prune the tips. If "method = 1", the function will prune the tips originating after an inputted temporal threshold. If method = 2, the function will prune the tips originating before a temporal threshold. Default is 1.
#'
#' @return The function returns a pruned tree in the "phylo" format if a single temporal threshold was inputted. Otherwise, if a vector of thresholds was inputted, it returns an list of pruned trees.
#'
#' @details
#' It uses the tip ages location in relation to an inputted time threshold to prune the phylogenetic tree. Setting "method = 2" makes slices based on the quantile distribution of tip ages available within the phylogenetic tree.
#'
#' @seealso Other slicing methods: [squeeze_root()],[squeeze_tips()],[squeeze_int()],[phylo_pieces()]
#'
#' @author Matheus Lima de Araujo <matheusaraujolima@live.com>
#'
#' @references
#' See the tutorial on how to use this function on our [website](https://araujomat.github.io/treesliceR/articles/Intro-treesliceR.html).
#'
#' @examples
#' # Generate a random tree
#' tree <- ape::rcoal(20)
#'
#' # Pruning the tips originating after 0.3 million years
#' tree1 <- prune_tips(tree, time = 0.1)
#'
#' # Plot it
#' plot(tree1)
#'
#'
#' # Pruning the tips based on quantiles of tip ages
#' tree2 <- prune_tips(tree, time = c(0.25, 0.75), qtl = TRUE)
#' plot(tree2[[1]])
#' plot(tree2[[2]])
#'
#' @export

prune_tips <- function(tree, time, qtl = FALSE, method = 1){

  # Tree configurations
  tree <- nodes_config(tree)

  # Tips position in matrix
  positions <- which(tree$config[,2] <= length(tree$tip.label))
  tip_depth <- tree$tree_depth - tree$config[positions, "YearBegin"]
  range_age <- range(tip_depth)

  if(method == 1){
    ## Is there only one time inputted?
    if(length(time) == 1){

      if(qtl == FALSE){
        # Is there an lineages depth bigger than the inputted time?
        if((sum(time < tip_depth) >= 1 & time >= 0) == TRUE){
          # Which linages to filter
          tips <- tree$tip.label[which(tip_depth >= time)]
          # Selecting only these tips
          return(ape::keep.tip(tree, tips))
        } else {
          stop(paste("The input time argument must be between", range_age[1], "and", range_age[2], sep = " "))
        }
      }

      if(qtl == TRUE){
        # Is the quantile inputted beyond the quantile intervals?
        if((time <= 1 & time >= 0) == TRUE){
          # Checking the tip depths
          thr <- stats::quantile(tip_depth, time)
          # Which linages to filter
          tips <- tree$tip.label[which(tip_depth >= thr)]
          # Selecting only these tips
          return(ape::keep.tip(tree, tips))
        } else {
          stop("The quantile value must fall within the range of 0 to 1")
        }
      }
    }

    ## Is there more than one time inputted?
    if(length(time) > 1){   # time = c(0.1, 0.2, 0.8)
      if(qtl == FALSE){
        # Is there an lineages depth bigger than the inputted time?
        if((sum(max(time) < tip_depth) >= 1 & min(time) >= 0) == TRUE){
          output <- lapply(time, function(x){
            # Which linages to filter
            tips <- tree$tip.label[which(tip_depth >= x)]
            # Selecting only these tips
            return(ape::keep.tip(tree, tips))
          })
          # Returning the output
          return(output)
        } else {
          stop(paste("The input time argument must be between", range_age[1], "and", range_age[2], sep = " "))
        }
      }

      if(qtl == TRUE){
        # Is the quantile inputted beyond the quantile intervals?
        if((max(time) <= 1 & min(time) >= 0) == TRUE){
          output <- lapply(time, function(x){
            # Checking the tip depths
            thr <- stats::quantile(tip_depth, x)
            # Which linages to filter
            tips <- tree$tip.label[which(tip_depth >= thr)]
            # Selecting only these tips
            return(ape::keep.tip(tree, tips))
          })
          # Returning the output
          return(output)
        } else {
          stop("The quantile value must fall within the range of 0 to 1")
        }
      }
    }
  }

  if(method == 2){
    ## Is there only one time inputted?
    if(length(time) == 1){

      if(qtl == FALSE){
        # Is there an lineages depth bigger than the inputted time?
        if((sum(time < tip_depth) >= 1 & time >= 0) == TRUE){
          # Which linages to filter
          tips <- tree$tip.label[which(tip_depth >= time)]
          # Selecting only these tips
          return(ape::drop.tip(tree, tips))
        } else {
          stop(paste("The input time argument must be between", range_age[1], "and", range_age[2], sep = " "))
        }
      }

      if(qtl == TRUE){
        # Is the quantile inputted beyond the quantile intervals?
        if((time <= 1 & time >= 0) == TRUE){
          # Checking the tip depths
          thr <- stats::quantile(tip_depth, time)
          # Which linages to filter
          tips <- tree$tip.label[which(tip_depth >= thr)]
          # Selecting only these tips
          return(ape::drop.tip(tree, tips))
        } else {
          stop("The quantile value must fall within the range of 0 to 1")
        }
      }
    }

    ## Is there more than one time inputted?
    if(length(time) > 1){   # time = c(0.1, 0.2, 0.8)
      if(qtl == FALSE){
        # Is there an lineages depth bigger than the inputted time?
        if((sum(max(time) < tip_depth) >= 1 & min(time) >= 0) == TRUE){
          output <- lapply(time, function(x){
            # Which linages to filter
            tips <- tree$tip.label[which(tip_depth >= x)]
            # Selecting only these tips
            return(ape::drop.tip(tree, tips))
          })
          # Returning the output
          return(output)
        } else {
          stop(paste("The input time argument must be between", range_age[1], "and", range_age[2], sep = " "))
        }
      }

      if(qtl == TRUE){
        # Is the quantile inputted beyond the quantile intervals?
        if((max(time) <= 1 & min(time) >= 0) == TRUE){
          output <- lapply(time, function(x){
            # Checking the tip depths
            thr <- stats::quantile(tip_depth, x)
            # Which linages to filter
            tips <- tree$tip.label[which(tip_depth >= thr)]
            # Selecting only these tips
            return(ape::drop.tip(tree, tips))
          })
          # Returning the output
          return(output)
        } else {
          stop("The quantile value must fall within the range of 0 to 1")
        }
      }
    }
  }
}
