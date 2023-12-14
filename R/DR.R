#' Calculate the tip diversification rates (DR) for a phylogenetic tree
#' @description
#' This function computes the tip diversification rates Jetz et al., 2012, or DR, for an inputted ultrametric phylogenetic tree.
#'
#' @usage DR(tree)
#'
#' @param tree phylo. An ultrametric phylogenetic tree in the "phylo" format.
#'
#' @return The function returns a data frame containing the tip diversification rates for all species within the inputted phylogenetic tree.
#'
#' @author Matheus Lima de Araujo <matheusaraujolima@live.com>
#'
#' @references
#' See the tutorial on how to use this function on our [website](https://araujomat.github.io/treesliceR/articles/Passeriformes-diversification.html).
#' Jetz, Walter, et al. "The global diversity of birds in space and time." Nature 491.7424 (2012): 444-448. <doi:10.1038/nature11631>
#'
#' @examples
#' # Generate a random tree
#' tree <- ape::rcoal(20)
#'
#' # Computing the tip-DR
#' DR(tree)
#'
#' @export

DR <- function(tree){

  # Capturing nodes information
  df <- as.data.frame(nodes_config(tree)$node_matrix)

  drs <- sapply(1:ncol(df), function(x){  #  x <- 8

    # Which nodes give origin to this specie
    nodes <- as.numeric(row.names(df[which(df[, x] == 1),]))

    # Which lines from the edge matrix it occupies?
    lines <- which(tree$edge[, 1] %in% nodes &             # Which nodes are in my spliting side from matrix
                     tree$edge[, 2] %in% c(x, nodes))

    # Save the branch lengths into a vector
    vec <- tree$edge.length[lines]

    # If there is an node of length 0, it is the species node conserved from the slice,
    # which could not be accounted for the DR calculation. Thus, must be removed.
    if(length(which(vec == 0)) > 0){
      vec <- vec[-c(which(vec == 0))]
    }

    # Calculate the DR
    dr_spp <- (sum(vec * sapply(length(vec):1, function(z){return(1/(2^(z-1)))})))^-1
    return(dr_spp)
  })

  # Saving it into a data.frame
  drs <- data.frame(Species = tree$tip.label, DR = drs)
  # Return the DRs of all tree species
  return(drs)
}
