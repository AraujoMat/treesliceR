#' Make a line graph or a grid map of estimated rates of accumulation of a given phylogenetic index.
#' @description
#' This function creates a line graph, or a grid map, depicting the estimated rates of accumulation of a given phylogenetic index (e.g., phylogenetic diversity, endemism, etc.), obtained from functions such as [CpD()], [CpE()], [CpB()], or [CpB_RW()].
#'
#' @usage CpR_graph(data, rate = NULL, map = NULL, pal = NULL, qtl = FALSE)
#'
#' @param data data frame. The outputted data frame from a CpR-rate function.
#' @param rate character string. The desired cumulative phylogenetic rate to plot, which can be phylogenetic diversity (CpD), phylogenetic endemism (CpE), phylogenetic B-diversity (CpB), or phylogenetic B-diversity range-weighted (CpB_RW). Default is NULL, but must be filled with "CpD", "CPE", "CpB", or "CpB_RW".
#' @param map spatial data. A grid map containing the assemblages at which the phylogenetic rates were assessed. Default is NULL.
#' @param pal character vector. A vector containing a color palette. If none provided, a default color palette will be used.
#' @param qtl boolean. Should the color palette be displayed according to CpR-rates quantiles? It can be either TRUE or FALSE. Default is FALSE.
#'
#' @return The function returns a ggplot graph.
#'
#' @seealso Other cumulative phylogenetic rate analysis: [CpD()], [CpE()], [CpB()], [CpB_RW()], [CpR_sensitivity()]
#' Other plotting: [CpR_sensitivity_plot()].
#'
#' @author
#' Matheus Lima de Araujo <matheusaraujolima@live.com>
#'
#' @references
#' See the tutorial on how to use this function on our [website](https://araujomat.github.io/treesliceR/articles/Passeriformes-diversification.html).
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
#' CpD_DF <- CpD(tree, n = 100, mat = mat)
#'
#' # Plot it using the CpR_graph
#' CpR_graph(CpD_DF, rate = "CpD")
#'
#' @export

CpR_graph <- function(data, rate = NULL, map = NULL, pal = NULL, qtl = FALSE){

  if(is.null(rate) == TRUE){
    stop("The user must inform the rate inputted to display it graphically")
  } else {

    if(rate == "CpD"){
      lab1 <- "CpD (%)"
      lab2 <- "CpD[rate]"
      lab3 <- expression(CpD["rate"])
      pOR <- "pDO"
    }
    if(rate == "CpE"){
      lab1 <- "CpE (%)"
      lab2 <- "CpE[rate]"
      lab3 <- expression(CpE["rate"])
      pOR <- "pEO"
    }
    if(rate == "CpB"){
      lab1 <- "CpB (%)"
      lab2 <- "CpB[rate]"
      lab3 <- expression(CpB["rate"])
      pOR <- "pBO"
    }
    if(rate == "CpB_RW"){
      lab1 <- "CpB (%)"
      lab2 <- "CpB[rate]"
      lab3 <- expression(CpB["rate"])
      pOR <- "pBO"
    }
  }

  ## Configuring the colour pallete
  if(is.null(pal) == TRUE){
    # Saving the default pallete on "pal" object
    pal <- grDevices::colorRampPalette(c("#DD0000", "#FE6230", "#FEF600", "#00BC00","#009BFE", "#000083", "#30009B"))

    # Checking if it need to be discretized
    if(qtl == TRUE){
      # Creating a vector with dicretized colors
      vec <- stats::quantile(data[,1], probs = seq(0, 1, (1/100)), na.rm = TRUE)
      # Creating a label to plot it
      names(vec)[-c(1, 21, 41, 61, 81, 101)] <- ""
    }
  } else {

    # Creating the pallete
    pal <- grDevices::colorRampPalette(pal)

    # Checking if it need to be discretized
    if(qtl == TRUE){
      # Creating a vector with dicretized colors
      vec <- stats::quantile(data[,1], probs = seq(0, 1, (1/100)), na.rm = TRUE)
      # Creating a label to plot it
      names(vec)[-c(1, 21, 41, 61, 81, 101)] <- ""
    }
  }

  # If the user doesnt want to map their CpB-rate
  if(is.null(map) == TRUE){

    # Removing NA's
    if(sum(is.na(data)) > 0){
      # Finding the rows contaning NA's and removing them
      data <- data[-(which(rowSums(is.na(data)) > 0)),]
    }

    # Create a data.frame containing time-slices with predicted CpB and its respective rate
    n_rows <- nrow(data)*1000
    id <- numeric(n_rows)
    CpR <- numeric(n_rows)
    pO <- numeric(n_rows)
    age <- numeric(n_rows)
    pred <- numeric(n_rows)
    # Combining them
    df <- data.frame(id, CpR, pO, age, pred)

    # Saving the initial iterator
    j <- 1
    # looping it
    for (i in 1:nrow(data)) {
      df[j:(j+999), 1] <- rep(i, 1000)
      df[j:(j+999), 2] <- rep(data[, 1][i], 1000)
      df[j:(j+999), 3] <- rep(data[, 3][i], 1000)
      df[j:(j+999), 4] <- cumsum(rep(max(data[, 3], na.rm = T)/1000, 1000))
      df[j:(j+999), 5] <- exp(-(data[, 1][i])*cumsum(rep(max(data[, 3], na.rm = T)/1000, 1000)))  # i <- 2
      j <- j + 1000
    }

    ## Creating the pBO histogram
    hist_pBo <- ggplot2::ggplot() +
      ggplot2::geom_histogram(ggplot2::aes(x = data[,3]), fill = "white", colour = "black", bins = 23, size = 1) +
      ggplot2::labs(x = pOR, y = "Frequency") +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                     axis.text = ggplot2::element_text(size = 8))


    # The user want discrete quantile intervals of colour palletes in their plot
    if(qtl == TRUE){

      # Plot
      return(ggplot2::ggplot(data = df, ggplot2::aes(age, pred, group = id)) +
               ggplot2::geom_line(ggplot2::aes(color = CpR), size = 1, alpha = 0.7) +
               ggplot2::binned_scale(name = NULL,
                            aesthetics = "color",
                            scale_name = "custom",
                            palette = function (x) c(pal(100)),
                            guide = "colorsteps",
                            breaks = c(vec),
                            labels = names(vec)) +
               ggplot2::ylab(lab1) + ggplot2::xlab("Years before present (My)") +
               ggplot2::theme_classic() +
               ggplot2::theme(legend.position = c(.38, 0.83),
                              legend.direction = "horizontal",
                              legend.key.size = ggplot2::unit(0.7, "cm"),
                              legend.title = ggplot2::element_blank(),
                              legend.text = ggplot2::element_text(size = 6),   # insert some breaks here?
                              legend.background = ggplot2::element_blank(),
                              axis.text = ggplot2::element_text(size = 12),
                              axis.title = ggplot2::element_text(size = 16)) +
               ggplot2::annotate("text", x = (max(data[,3], na.rm = T)*0.38), y = 0.95, label = lab2, size = 5,
                                 parse=TRUE, colour = "black") +
               ggplot2::annotation_custom(ggplot2::ggplotGrob(hist_pBo), xmin = (max(data[,3], na.rm = T)*0.53),
                                 xmax = (max(data[,3], na.rm = T)*0.99), ymin = 0.53, ymax = 0.99))
    } else {

      # Plot
      return(ggplot2::ggplot(data = df, ggplot2::aes(age, pred, group = id)) +
               ggplot2::geom_line(ggplot2::aes(color = CpR), size = 1, alpha = 0.7) +
               ggplot2::scale_colour_gradientn(colours = pal(100)) +
               ggplot2::ylab(lab1) + ggplot2::xlab("Years before present (My)") +
               ggplot2::theme_classic() +
               ggplot2::theme(legend.position = c(.38, 0.83),
                     legend.direction = "horizontal",
                     legend.key.size = ggplot2::unit(0.7, "cm"),
                     legend.title = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size = 6),
                     legend.background = ggplot2::element_blank(),
                     axis.text = ggplot2::element_text(size = 12),
                     axis.title = ggplot2::element_text(size = 16)) +
               ggplot2::annotate("text", x = (max(data[,3], na.rm = T)*0.38), y = 0.95, label = lab2, size = 5,
                        parse=TRUE, colour = "black") +
               ggplot2::annotation_custom(ggplot2::ggplotGrob(hist_pBo), xmin = (max(data[,3], na.rm = T)*0.53),
                                 xmax = (max(data[,3], na.rm = T)*0.99), ymin = 0.53, ymax = 0.99))
    }
  } else {

    ## If a grid map is provided in the argument "map"

    # Turn into a sf spatial object
    map <- sf::st_as_sf(map)

    # Adding the CpB-rate information on the grid map cells
    map$CpR <- data[,1]

    if(qtl == TRUE){

      # Return the provided grid map with reescaled to quantile gradient based on "pal" colours
      return(ggplot2::ggplot() +
               ggplot2::geom_sf(data = map, ggplot2::aes(fill = CpR), colour = "grey20", alpha = 1) +
               ggplot2::binned_scale(name = lab3,
                            aesthetics = "fill",
                            scale_name = "custom",
                            palette = function (x) c(pal(100)),
                            guide = "colorsteps",
                            breaks = c(vec),
                            labels = names(vec)) +
               ggplot2::theme_void() + ggplot2::labs(x = NULL, y = NULL) +
               ggplot2::theme(legend.position = c(.15, .14),
                              legend.title = ggplot2::element_text(size = 13),
                              legend.key.size = ggplot2::unit(0.75, "cm"),
                              legend.direction = "horizontal"))

    } else {

      # Return the provided grid map without reescaling into quantiles the response variable
      return(ggplot2::ggplot() +
               ggplot2::geom_sf(data = map, ggplot2::aes(fill = CpR), colour = "grey20", alpha = 1) +
               ggplot2::scale_fill_gradientn(name = lab3, colours = pal(100), na.value = "white") +
               ggplot2::theme_void() + ggplot2::labs(x = NULL, y = NULL) +
               ggplot2::theme(legend.position = c(.15, .14),
                              legend.title = ggplot2::element_text(size = 13),
                              legend.key.size = ggplot2::unit(0.75, "cm"),
                              legend.direction = "horizontal"))
    }
  }
}
