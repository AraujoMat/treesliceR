% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nodes_config.R
\name{nodes_config}
\alias{nodes_config}
\title{Capture general branching information for an inputted phylogenetic tree}
\usage{
nodes_config(tree)
}
\arguments{
\item{tree}{phylo. An ultrametric phylogenetic tree in the "phylo" format.}
}
\value{
The function returns a phylogenetic tree in the "phylo" format containing three novel pieces of information (stored within "config", "node_matrix", and "tree_depth").
}
\description{
This function captures general tree information, including branching, node positions, and depths. It serves as a core function underlying all algorithms for slicing phylogenies.
}
\details{
This function captures node and edge information from an ultrametric phylogenetic tree. The function provides a data frame containing detailed branching information (within the internal "config" object), a node matrix (within "node_matrix"), and the tree age (within "tree_depth").

More specifically, the "config" object returns the following information:
NodeBegin: the node at which a given branch begins.
NodeEnd: the node at which a given branch ends.
NodeLength: the branch length of that nodes interval.
YearBegin: the year at which a given node begins.
YearEnd: the year at which a given node ends.
}
\examples{
# Generate a random tree
tree <- ape::rcoal(20)

# Capture tree information
tree <- nodes_config(tree)

# Accessing these informations
tree$config # Nodes configurations
tree$node_matrix # Node matrix
tree$tree_depth # Tree age

}
\seealso{
Phylogenetic slicing methods: \code{\link[=squeeze_tips]{squeeze_tips()}},\code{\link[=squeeze_root]{squeeze_root()}},\code{\link[=squeeze_int]{squeeze_int()}}.
}
\author{
Matheus Lima de Araujo \href{mailto:matheusaraujolima@live.com}{matheusaraujolima@live.com}
}
