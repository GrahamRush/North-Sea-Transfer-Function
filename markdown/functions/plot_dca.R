#' 
#' @title plot the dca
#'
#' @description plots the dca results
#'
#' @param origin is logical to plot by origin of 0
#' @param scale is the scale by which to plot the site by swli
#'
#' @return plots the dca for site and species
#'
#' @keywords tolerances
#' @export
#' @examples
#' plot_dca(origin = TRUE, scale = 50000, gradients = FALSE)

plot_dca <- function (origin = TRUE, scale = 50000, gradients = F) {
  df.split (training.set, ts.sp)    # split the data
  env$colour <- training_sets$cols$northsea[env$Site]

  axis1 <- c(paste("Axis 1 (Eigenvalue: ", round (ns.dca$evals.decorana[1], digits = 2), ")", sep = ""))    # create axis labels
  axis2 <- c(paste("Axis 2 (Eigenvalue: ", round (ns.dca$evals.decorana[2], digits = 2), ")", sep = ""))    # create axis labels
  
  # to add species include the following 
  m.2 <- apply(spec, 2, function(x) max(x > 10, na.rm=TRUE)) # calculate which species have a max abundace > 10 % 
  m3 <- as.vector(which(m.2 ==1)) # position of species if > 10 %
  # calculate names if needed
  nam <- colnames (spec)
  nam4 <- as.character()
  for (i in 1:ncol(spec)) {
  nam2 <- unlist(strsplit(nam[i], "[.]"))
  nam3 <- paste(nam2[1], nam2[2]) 
  nam4[i] <- nam3
  }
  if (origin == FALSE) {
    plot(ns.dca, type = "n", lwd = 2, xlab = "", ylab = "", las = 1)    #  create empty plot
    axis(1)    # add axis
    axis(2, las = 2)
    mtext(side = 1, axis1, line = 3, cex =1)
    mtext(side = 2, axis2, line = 2.5, cex =1)
    t1 <- c(paste(reg, "DCA", sep = " "))    # create main title
    mtext(t1, font = 2, side=3, outer=TRUE, line=-3)    # add main title
    points(ns.dca, "sites", cex = env$SWLI * (env$SWLI * 2) / scale, bg = env$colour, pch=21)    # add sites
    legend(x = "topleft", name, col=colr, pch = 20, pt.cex=2) # set to appropriate
    if (gradients == TRUE) {
      swli.grad <- envfit(ns.dca, env$SWLI, permu = 999)    # create environemtal gradients
      plot(swli.grad, font = 2, col = "blue", labels = "SWLI")    #  add environemtal gradients
    }
    #  re-plot with species
    plot(ns.dca, type = "n", lwd = 2, xlab = "", ylab = "", las = 1)    #  create empty plot
    axis(1)    # add axis
    axis(2, las = 2)
    mtext(side = 1, axis1, line = 3, cex =1)
    mtext(side = 2, axis2, line = 2.5, cex =1)
    points(ns.dca, "sites", cex = env$SWLI * (env$SWLI * 2) / scale, bg=env$colour, pch=21)    # add sites
    if (gradients == TRUE) {
      swli.grad <- envfit(ns.dca, env$SWLI, permu = 999)    # create environemtal gradients
      plot(swli.grad, font = 2, col = "blue", labels = "SWLI")    #  add environemtal gradients
    }
    points(ns.dca, "spec", pch=4, lwd= 2, scaling = 3, col = 2)    # add species
    text (ns.dca, "spec", cex = 0.7, pos = 3)   # add text
  }
  if (origin == TRUE) {    # plot samples to an origin of 0
    plot.new()
    par(fig=c(0,0.5,0,1), new=TRUE, xpd = TRUE) # , mar = c(0, 2, 4, 1))
    plot(x = c(ns.dca$rproj[,1]), y = c(ns.dca$rproj[,2]), type = "n", lwd = 2, yaxt = "n", xaxt="n", xlab = "", ylab = "", bty = "n", las = 1)    # empty plot
    par(lwd =1, font = 1)     # add axis
    axis(1)
    axis(2, las = 2)
    mtext(side = 1, axis1, line = 3, cex =1)
    mtext(side = 2, axis2, line = 2.5, cex =1)
    points(x = c(ns.dca$rproj[,1]), y = c(ns.dca$rproj[,2]), cex= env$SWLI * (env$SWLI * 2) / scale, bg = env$colour, pch = 21)
    # plot the names
    par(fig=c(0.5,1,0,1), new=TRUE, xpd = TRUE) # , mar = c(0, 2, 4, 1))
    plot(x = c(ns.dca$rproj[,1]), y = c(ns.dca$rproj[,2]), type = "n", lwd = 2, yaxt = "n", xaxt="n", xlab = "", ylab = "", bty = "n", las = 1)
    axis(1)     # add axis
    axis(2, las = 2)
    mtext(side = 1, axis1, line = 3, cex =1)
    mtext(side = 2, axis2, line = 2.5, cex =1)
    points(x = c(ns.dca$rproj[,1]), y = c(ns.dca$rproj[,2]), cex= env$SWLI * (env$SWLI * 2) / scale, bg = env$colour, pch = 21)    # add sites scaled to swli value
    text (x = c(ns.dca$cproj[m3,1]), y = c(ns.dca$cproj[m3,2]), cex = 0.7, labels = c(nam4[c(m3)]), pos = 3)    # add species text
    points (x = c(ns.dca$cproj[m3,1]), y = c(ns.dca$cproj[m3,2]), cex = 0.7, col = 2, pch = 4)     # add species as X
    legend(x = "topleft", name, col=colr, pch = 20, pt.cex=2)    # add a legen and set to appropriate
  }
}
