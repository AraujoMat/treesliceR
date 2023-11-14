#' Australia grid map
#'
#' This is a gridded shapefile of Australia, with grid cells at a 0.5-degree resolution.
#'
#' @format Shapefile
"AU_grid"

#' Presence-absence matrix of Australian passeriformes
#'
#' This is a presence-absence matrix for passeriformes obtained from Australian grid-cells. Each column represents a specific passeriforme species, while each row corresponds to an assemblage (or grid cell) within the Australian grid-map. All passerine range maps were sourced from the BirdLife platform.
#'
#' @format Matrix
#' @source <https://www.birdlife.org/>
"pass_mat"

#' Presence-absence matrices for focal and adjacent assemblages
#'
#' This is a list of presence-absence matrices for passeriformes from Australia. Each object within the list represents a passerine presence-absence matrix contaning the assemblage of a focal cell and its respective neighborhood assemblages.
#'
#' @format List
"pass_asb"

#' Phylogenetic tree of passerines from Australia
#'
#' This is a list contaning 100 phylogenetic trees of passerines from Australia. Each object within the list corresponds to a sampled phylogeny of class "phylo". The phylogenies are subsets obtained from the global bird phylogeny (https://birdtree.org/) made available by Jetz et al. (2012) in their published research, "The global diversity of birds in space and time".
#'
#' @format List of phylogenetic trees of format phylo
#' @source <https://doi.org/10.1038/nature11631>
"pass_trees"