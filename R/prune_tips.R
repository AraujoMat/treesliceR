# Description:
# Collapses, or prune, the tips of a phylogeny based on a given temporal threshold.

# Arguments:
# tree = a phylogenetic tree.
# time = a numeric or vector of temporal thrsholds in millions of years (numeric)
# or quantile values (set argument qtl = TRUE). If a vector is inputted, it returns
# a list of trees based on the vector values set.
# method = could be 1 or 2. Setting "method = 1" defines makes the algorithm drops
# (or prune) the tips originating after a given threshold set. Otherwise, if the
# user wants to drop those tips branches originating before a given time threshold,
# set "method = 2".

# Errors:
# There is no tip with branching event deeper than the value inputted
# Quantile value doesnt exist, need to be within probabilities between 0 and 1.

#' Title
#'
#' @param tree
#' @param time
#' @param qtl
#' @param method
#'
#' @return
#' @export
#'
#' @examples

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
