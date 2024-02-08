# treesliceR 1.0.1

* The functions for calculating beta-diversity indexes have undergone changes in their arguments. Previously, they required a list containing each focal assemblage along with their respective adjacent cells separately, filled in the `asb` argument. Now, for simplicity, these functions need to be filled with a species presence-absence matrix (i.e., sites in rows and species in columns) in the `mat` argument and an adjacency matrix containing the focal assemblages and their respective adjacent cells in the `adj` argument, being the `asb` argument removed.

* The internal object `pass_asb` has been removed. This object previously contained a list of where each object had a matrix of a focal cell and its adjacent neighborhoods for Australian passeriformes. A spatial adjacency matrix of the Australian grid assemblages has been added and stored within the `AU_adj` internal object.

* The argument `criteria` in cutting functions has been changed to its singular form, `criterion`.

# treesliceR 1.0.0

* Initial CRAN submission.
