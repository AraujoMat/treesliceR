#' Make a plot of the sensitivity analysis evaluated under cumulative phylogenetic rate functions
#' @description
#' This function generates a plot of the sensitivity analysis calculated for a given cumulative phylogenetic rate (CpR) assessed through the [CpR_sensitivity()] function.
#'
#' @usage CpR_sensitivity_plot(sst_output, rate = NULL, stc = "mean")
#'
#' @param sst_output data frame. The outputted data frame from the [CpR_sensitivity()] function.
#' @param stc character string. A statistical measure to summarize the phylogenetic rates and create the plot, which could be filled with "mean", "var", "median", "sd", "min", and "max". Default is "mean".
#' @param rate character string. The desired phylogenetic index rate to display. It can be filled with "CpD", "CPE", "CpB", or "CpB_RW".
#'
#' @return The function returns a ggplot graph.
#'
#' @seealso CpR sensitivity analysis: [CpR_sensitivity()].
#'
#' @author Matheus Lima de Araujo <matheusaraujolima@live.com>
#'
#' @references
#' Users can use the [CpR_sensitivity_plot()] function for plotting sensitivity analysis outputs.
#'
#' @examples
#' # Generate a random tree
#' tree <- ape::rcoal(20)
#'
#' # Create a presence-absence matrix
#' mat <- matrix(sample(c(1,0), 20*10, replace = TRUE), ncol = 20, nrow = 10)
#' colnames(mat) <- tree$tip.label
#'
#' # Calculate the CpD for 100 tree slices
#' CpD(tree, n = 100, mat = mat)
#'
#' # Create a vector of number of slices
#' vec <- c(25, 50, 75, 100, 125, 150)
#'
#' # Calculate the sensitivity of the CpD
#' Sens_out <- CpR_sensitivity(tree, vec, mat, rate = "CpD", samp = 5)
#'
#' # Plot the sensitity analysis
#' CpR_sensitivity_plot(Sens_out, rate = "CpD", stc = "mean")
#'
#' @export


CpR_sensitivity_plot <- function(sst_output, rate = NULL, stc = "mean"){

  if(is.null(rate) == TRUE){
    stop("The user must inform the rate inputted to display it graphically")
  } else {
    if(rate == "CpD"){
      lab_R <- bquote(.(stc) ~  CpD["rate"])
    }
    if(rate == "CpE"){
      lab_R <- bquote(.(stc) ~  CpE["rate"])
    }
    if(rate == "CpB"){
      lab_R <- bquote(.(stc) ~  CpB["rate"])
    }
    if(rate == "CpB_RW"){
      lab_R <- bquote(.(stc) ~  CpB["rate"])
    }
  }

  # Removing NA's from the outputs
  if(sum(is.na(sst_output)) > 0){
    # Finding the rows contaning NA's and removing them
    sst_output <- sst_output[-(which(rowSums(is.na(sst_output)) > 0)),]
  }

  # Creating a df with summarizing statistic informations
  val <- apply(sst_output, 2, stc)
  vec <- as.numeric(colnames(sst_output))
  df <- data.frame(val, vec)

  # Plotting those informations
  return(ggplot2::ggplot(df, ggplot2::aes(x = vec)) +
           ggplot2::geom_line(ggplot2::aes(y = val), colour = "grey50", size = 1.2) +
           ggplot2::geom_point(ggplot2::aes(y = val), col = "black", fill = "#000000", size = 2, shape = 21) +
           ggplot2::labs(x = "Number of phylogenetic slices", y = lab_R) +
           ggplot2::scale_x_continuous(breaks = vec,
                                       labels = vec) +
           ggplot2::theme_classic() +
           ggplot2::theme(axis.text.y = ggplot2::element_text(colour = "black", size=9),
                          axis.text.x = ggplot2::element_text(colour = "black", size=9),
                          axis.title = ggplot2::element_text(size=14),
                          axis.ticks = ggplot2::element_line(colour = "black")))
}
