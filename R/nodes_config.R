#' Capture general branching information for an inputted phylogenetic tree
#' @description
#' This function captures general tree information, including branching, node positions, and depths. It serves as a core function underlying all algorithms for slicing phylogenies.
#'
#' @usage nodes_config(tree)
#'
#' @param tree phylo. An ultrametric phylogenetic tree in the "phylo" format.
#'
#' @details
#' This function captures node and edge information from an ultrametric phylogenetic tree. The function provides a data frame containing detailed branching information (within the internal "config" object), a node matrix (within "node_matrix"), and the tree age (within "tree_depth").
#'
#' More specifically, the "config" object returns the following information:
#' NodeBegin: the node at which a given branch begins.
#' NodeEnd: the node at which a given branch ends.
#' NodeLength: the branch length of that nodes interval.
#' YearBegin: the year at which a given node begins.
#' YearEnd: the year at which a given node ends.
#'
#' @return The function returns a phylogenetic tree in the "phylo" format containing three novel pieces of information (stored within "config", "node_matrix", and "tree_depth").
#'
#' @seealso Phylogenetic slicing methods: [squeeze_tips()],[squeeze_root()],[squeeze_int()].
#'
#' @author Matheus Lima de Araujo <matheusaraujolima@live.com>
#'
#' @examples
#' # Generate a random tree
#' tree <- ape::rcoal(20)
#'
#' # Capture tree information
#' tree <- nodes_config(tree)
#'
#' # Accessing these informations
#' tree$config # Nodes configurations
#' tree$node_matrix # Node matrix
#' tree$tree_depth # Tree age
#'
#' @export

nodes_config <- function(tree){

  # if the tree is ultrametric, make the evaluations
  if(ape::is.ultrametric(tree) == TRUE){
    ## Capturing nodes distances and configurations
    matrix_nodes <- ape::dist.nodes(tree)
    # Putting it into my tree
    tree$config <- data.frame(NodeBegin = tree$edge[, 1],
                              NodeEnd = tree$edge[, 2],
                              NodeLength = tree$edge.length,
                              YearBegin = matrix_nodes[tree$edge[, 2], tree$edge[1, 1]] - tree$edge.length,
                              YearEnd = matrix_nodes[tree$edge[, 2], tree$edge[1, 1]])

    ## Creating a node path matrix
    node_mat <- matrix(0, ncol = length(tree$tip.label), nrow = tree$Nnode)
    # Naming its columns as species and its rows as nodes
    row.names(node_mat) <- unique(tree$edge[,1])
    colnames(node_mat) <- tree$tip.label

    # Creating a presence-abscense nodes-matrix
    paths <- ape::nodepath(tree)
    for(i in 1:length(tree$tip.label)){
      node_mat[which(row.names(node_mat) %in% paths[[i]]), i] <- 1
    }
    # Adding it to the tree
    tree$node_matrix <- node_mat

    # Calculating the total tree depth
    tree$tree_depth <- matrix_nodes[tree$edge[1, 1], 1] # distance from root to the species 1 (Just for ultrametric trees)

    return(tree)
  } else {

    stop("The inputted phylogenetic tree is not ultrametric.")
  }
}
