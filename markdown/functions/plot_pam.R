#' 
#' @title plot pam
#'
#' @description creates the neccessary data and plots the pam results
#'
#' @param df is the dataframe normally ts.clean
#' @param n is the model to run
#' @param y0 the y axis start
#' @param y1 the y axis end
#' @param x0 the x axis start
#' @param x0 the x axis end
#'
#' @return plots species optima and tolerances
#'
#' @keywords tolerances
#' @export
#' @examples
#' species.tol(ts.clean, 1, 60, 240, 0, 38)


col.ramp <- function(x, n, cols, b) { # x is clustered data, n is the number of breaks, cols = colours, b = where the breaks occur or the number of breaks
  rbPal <- colorRampPalette(cols)
  env$col <- rbPal(n)[as.numeric(cut(env$SWLI,breaks = b))]
  sil.info <- as.data.frame(x$pamobject$silinfo$widths)
  sil.info$id  <- 1:nrow(sil.info)
  merg <- merge(as.data.frame(sil.info), as.data.frame(env), by='row.names', all=TRUE)
  pam.col <- merg[order(merg$id), ]
  clus.info <- as.data.frame(x$pamobject$clustering)
  clus.info$id  <- 1:nrow(clus.info)
  merg <- merge(as.data.frame(clus.info), as.data.frame(env), by='row.names', all=TRUE)
  clus.col <- merg[order(merg$id), ]
  return(list(pam.col = pam.col, clus.col = clus.col))
}


plot.pamk <- function(k, p2 = TRUE, p3 = TRUE, boxplot = TRUE) {
  # to create colours for plotting
  col_ramp <- col.ramp(cluster, n = 6, cols = c('#d73027','#fc8d59','#fee090','#e0f3f8','#91bfdb','#4575b4'), b = c(50, 100, 125, 150, 175, 200, 270) )
  # to plot the clusters
  str(si <- silhouette(cluster$pamobject))
  row.names(si) <- col_ramp$pam.col$SWLI
  plot(si, nmax = 180, col = col_ramp$pam.col$col, cex.names = 0.4)
  if (p2 == TRUE) plot(col_ramp$clus.col$SWLI, col_ramp$clus.col$`x$pamobject$clustering`, col = col_ramp$clus.col$col, xlab = "SWLI", ylab = "Cluster", las=1)
  if (p3 == TRUE) plot(col_ramp$clus.col$SWLI, col_ramp$clus.col$`x$pamobject$clustering`, col = colr, xlab = "SWLI", ylab = "Cluster", las=1)
  if (boxplot == TRUE) {
    boxplot(col_ramp$clus.col$SWLI ~ col_ramp$clus.col$`x$pamobject$clustering` , ylab="SWLI" , boxwex=0.4 , main="", las = 1, bty = "n", xlab = "Cluster", axes=F)
    axis(side=2, lwd= 1, las = 2, font=1, tck= -0.03)
    axis(side=2,  at=seq(100, 350, by=10), lwd=1, las = 2, tck= -0.01, labels = NA)
    axis(side=1, at=seq(1,3,1),lwd= 1, las = 1, font=1, tck= -0.03)
  }
  clus.df <- data.frame("cluster" = 1:cluster$nc)
  clus.df$size <- cluster$pamobjec$clusinfo[,1]
  clus.df$av.widths <- round(cluster$pamobjec$silinfo$clus.avg.widths, digits = 2)
  clus.df$max_diss <- round(cluster$pamobjec$clusinfo[, 2], digits = 2)
  clus.df$av_diss <- round(cluster$pamobjec$clusinfo[, 3], digits = 2)
  row.names(col_ramp$clus.col) <- col_ramp$clus.col$Row.names
  clus.4 <- col_ramp$clus.col[cluster$pamobjec$medoids, ]
  clus.df$med_swli <- clus.4$SWLI
  clus.df$min_swli <- tapply(col_ramp$clus.col$SWLI, col_ramp$clus.col$`x$pamobject$clustering`, min)
  clus.df$max_swli <- tapply(col_ramp$clus.col$SWLI, col_ramp$clus.col$`x$pamobject$clustering`, max)
  assign(paste0("clusters", cluster$nc), clus.df)
}
