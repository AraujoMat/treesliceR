#' Slices a temporal interval within a phylogenetic tree
#' @description
#' This function slices a temporal interval located within an ultrametric phylogenetic tree.
#'
#' @usage squeeze_int(tree, from, to, invert = FALSE, criteria = "my", dropNodes = FALSE)
#'
#' @param tree phylo. An ultrametric phylogenetic tree in the "phylo" format.
#' @param from numeric. A temporal threshold value that determines the time at which the interval should start.
#' @param to numeric. A temporal threshold value that determines the time at which the interval should end.
#' @param invert logical. A logical value indicating if the desired slice must be executed inside (invert = FALSE) or outside (invert = TRUE) the defined interval. Using the argument as TRUE will return a list containing a root and tip slices. Default is FALSE.
#' @param criteria character string. The method for cutting the tree. It can be either "my" (million years) or "PD" (accumulated phylogenetic diversity). Default is "my".
#' @param dropNodes logical. A logical value indicating whether the nodes that were sliced (void nodes, presenting no branch length) should be preserved in the node matrix. Default is FALSE.
#'
#' @return The function returns an time-slice interval of a phylogenetic tree in the "phylo" format.
#'
#' @details
#'
#' \bold{Slicing approach}
#'
#' To return a given phylogenetic interval, this function simultaneously applies simultaneously the same logic as [squeeze_tips()] and [squeeze_root()]. If "invert = TRUE", then the temporal interval set will be excluded from the phylogeny, returning a list containing a tip and a root slices.
#'
#' @seealso Other slicing methods: [squeeze_tips()], [squeeze_root()], [phylo_pieces()], [prune_tips()]
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
#' # Slice an interval from 0.5 to 0.2 million years
#' tree <- squeeze_int(tree, from = 0.5, to = 0.2)
#'
#' # Plot it
#' plot(tree)
#'
#' @export

squeeze_int <- function(tree, from, to, invert = FALSE, criteria = "my", dropNodes = FALSE){

  # The used want a phylogenetic interval?
  if(invert == FALSE){
    # The order for making the slices are important depending on the criteria used
    if(criteria == "my"){
      if(from > to){
        tree <- squeeze_root(tree, from, criteria = criteria, dropNodes = dropNodes)
        tree <- squeeze_tips(tree, to, criteria = criteria, dropNodes = dropNodes)
      } else {
        stop("The thresholds set in arguments [from] and [to] are incompatible")
      }

    }
    if(criteria == "pd"){
      if(from < to){
        tree <- squeeze_tips(tree, to, criteria = criteria, dropNodes = dropNodes)
        tree <- squeeze_root(tree, from, criteria = criteria, dropNodes = dropNodes)
      } else {
        stop("The thresholds set in arguments [from] and [to] are incompatible")
      }

    }
    return(tree)
  }

  # or he wants to remove a phylogenetic interval?
  if(invert == TRUE){

    if(criteria == "my"){
      if(from < to){
        stop("The thresholds set in arguments [from] and [to] are incompatible")
      }
    } else {
      tree1 <- squeeze_root(tree, to, criteria = criteria, dropNodes = dropNodes)
      tree2 <- squeeze_tips(tree, from, criteria = criteria, dropNodes = dropNodes)
    }

    if(criteria == "pd"){
      if(from > to){
        stop("The thresholds set in arguments [from] and [to] are incompatible")
      }
    } else {
      tree1 <- squeeze_root(tree, to, criteria = criteria, dropNodes = dropNodes)
      tree2 <- squeeze_tips(tree, from, criteria = criteria, dropNodes = dropNodes)
    }

    return(list(tree1, tree2))
  }
}
