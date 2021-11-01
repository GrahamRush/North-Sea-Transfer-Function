#' 
#'
#' @title plot the residuals
#'
#' @description plot the residuals and predictions following cross-validation
#'
#' @param zoom is whether to zoom to the data or retain consistency
#' @param wa.mod describes which model to use. 2 = wa = classical
#' 
#' @return files using the names supplied above. 1. training set with data. 2.
#'
#' @keywords compare
#' @export
#' @examples
#' plot_residuals(wa.mod = 2, zoom = TRUE)
#' 

plot_residuals <- function(wa.mod = 2, zoom = TRUE) {
  par(mfcol=c(2,2), xpd=FALSE)    # set parameters
  tf.cols <- colr[env$Site]    # set colours
  #  plot first run
  # wapls.mod <- 2 # adjust this to fit the best choice of components
  if (zoom == FALSE) {
    plot (env$SWLI, cv$predicted[ ,wa.mod], xlim = c(40, 300), ylim = c(40, 300), pch = 16, cex = 0.7, abline(0, 1, lty = 2), xlab = "Elevation SWLI", ylab = "Predicted SWLI", col = tf.cols, las = 1) #  set sub heading and model number cv$predicted [,!]
    plot (env$SWLI, cv$residuals.cv[ ,wa.mod], ylim = c(-100, 100), abline(0, 0, lty = 2), pch = 16, cex = 0.7, ylab = "Residuals (SWLI)", xlab = "Elevation (SWLI)", col = tf.cols, las = 1, main = "WA Classical Tolerance downweighting")
  }
  if (zoom == TRUE) {
    plot (env$SWLI, cv$predicted[ ,wa.mod], xlim = c(min(env$SWLI), max(env$SWLI)), ylim = c(min(env$SWLI), max(env$SWLI)), pch = 16, cex = 0.7, abline(0, 1, lty = 2), xlab = "Elevation SWLI", ylab="Predicted SWLI", col = tf.cols, las = 1) #  set sub heading and model number cv$predicted [,!]
    plot (env$SWLI, cv$residuals.cv[ ,wa.mod], ylim = c(-100, 100), xlim = c(min(env$SWLI), max(env$SWLI)), abline(0, 0, lty = 2), pch = 16, cex = 0.7, ylab="Residuals (SWLI)", xlab = "Elevation (SWLI)", col=tf.cols, las=1, main = "WA Classical Tolerance downweighting")
  }
  if (zoom == FALSE) {
    plot (env$SWLI, cv.loo$predicted[ ,wapls.mod], xlim = c(40, 300), ylim = c(40, 300), pch= 16, cex = 0.7, abline(0, 1, lty = 2), xlab = "Elevation SWLI", ylab = "Predicted SWLI", col = tf.cols, las=1) #  set sub heading and model number cv$predicted [,!]
    plot (env$SWLI, cv.loo$residuals.cv[ ,wapls.mod], ylim = c(-100, 100), xlim = c(40, 300), abline(0, 0, lty = 2), pch = 16, cex = 0.7, ylab= "Residuals (SWLI)", xlab = "Elevation (SWLI)", col = tf.cols, las = 1, main = main1)
  }
  abline(h = sd(env$SWLI), lty = 3, col ="grey")
  abline(h = -2*sd(env$SWLI), lty = 4, col ="grey")
  abline(h = -sd(env$SWLI), lty = 3, col ="grey")
  abline(h = 2*sd(env$SWLI), lty = 4, col ="grey")
  if (zoom == TRUE) {
    plot (env$SWLI , cv.loo$predicted[ ,wapls.mod], xlim = c(min(env$SWLI), max(env$SWLI)), ylim = c(min(env$SWLI), max(env$SWLI)), pch = 16, cex = 0.7, abline(0, 1, lty = 2), xlab = "Elevation SWLI", ylab="Predicted SWLI", col=tf.cols, las=1) #  set sub heading and model number cv$predicted [,!]
    plot (env$SWLI , cv.loo$residuals.cv [ ,wapls.mod], ylim = c(-100, 100), abline(0, 0, lty = 2), pch = 16, cex = 0.7, ylab = "Residuals (SWLI)", xlab = "Elevation (SWLI)", col = tf.cols, las=1, main = main1)
  }
  abline(h = sd(env$SWLI), lty = 3, col ="grey")
  abline(h = -2*sd(env$SWLI), lty = 4, col ="grey")
  abline(h = -sd(env$SWLI), lty = 3, col ="grey")
  abline(h = 2*sd(env$SWLI), lty = 4, col ="grey")
  mtext(reg, font = 2, side = 3, outer = TRUE, line = -3) # add main title
}


#' 
#'
#' @title plot the residuals
#'
#' @description plot the residuals and predictions following locally weighted transfer functions
#'
#' @param zoom is whether to zoom to the data or retain consistency
#' @param model describes which model to use. 2 = wa classical, or wapls componment 2
#' 
#' @return files using the names supplied above. 1. training set with data. 2.
#'
#' @keywords compare
#' @export
#' @examples
#' plot_residuals(wa.mod = 2, zoom = TRUE)
#' 

plot_residuals_lw <- function(model = 2, zoom = TRUE) {
  par(mfcol=c(2,2), xpd=FALSE)    # set parameters
  tf.cols <- colr[env$Site]    # set colours
  #  plot first run
  if (zoom == FALSE) {
    plot (env$SWLI, cv.lw$predicted[ ,wapls.mod], xlim = c(40, 300), ylim = c(40, 300), pch= 16, cex = 0.7, abline(0, 1, lty = 2), xlab = "Elevation SWLI", ylab = "Predicted SWLI", col = tf.cols, las=1) #  set sub heading and model number cv$predicted [,!]
    plot (env$SWLI, cv.lwresiduals.cv[ ,wapls.mod], ylim = c(-100, 100), xlim = c(40, 300), abline(0, 0, lty = 2), pch = 16, cex = 0.7, ylab= "Residuals (SWLI)", xlab = "Elevation (SWLI)", col = tf.cols, las = 1, main = main1)
  }
  abline(h = sd(env$SWLI), lty = 3, col ="grey")
  abline(h = -2*sd(env$SWLI), lty = 4, col ="grey")
  abline(h = -sd(env$SWLI), lty = 3, col ="grey")
  abline(h = 2*sd(env$SWLI), lty = 4, col ="grey")
  if (zoom == TRUE) {
    plot (env$SWLI , cv.lw$predicted[ ,wapls.mod], xlim = c(min(env$SWLI), max(env$SWLI)), ylim = c(min(env$SWLI), max(env$SWLI)), pch = 16, cex = 0.7, abline(0, 1, lty = 2), xlab = "Elevation SWLI", ylab="Predicted SWLI", col=tf.cols, las=1) #  set sub heading and model number cv$predicted [,!]
    plot (env$SWLI , cv.lw$residuals.cv [ ,wapls.mod], ylim = c(-100, 100), abline(0, 0, lty = 2), pch = 16, cex = 0.7, ylab = "Residuals (SWLI)", xlab = "Elevation (SWLI)", col = tf.cols, las=1, main = main1)
  }
  abline(h = sd(env$SWLI), lty = 3, col ="grey")
  abline(h = -2*sd(env$SWLI), lty = 4, col ="grey")
  abline(h = -sd(env$SWLI), lty = 3, col ="grey")
  abline(h = 2*sd(env$SWLI), lty = 4, col ="grey")
  mtext(reg, font = 2, side = 3, outer = TRUE, line = -3) # add main title
}

#' 
#'
#' @title plot the loso
#'
#' @description plot the loso and bootstrapped cross-validations
#'
#' @param zoom is whether to zoom to the data or retain consistency
#' @param text describes whether to label points or include a legend
#' 
#' @return files using the names supplied above. 1. training set with data. 2.
#'
#' @keywords compare
#' @export
#' @examples
#' plot_loso(zoom = TRUE, text = FALSE)
#' 


plot_loso <- function(zoom = TRUE, text = FALSE) {
  par(mfrow=c(1,1))
  #  set colours
  loso_site <- colr[compare$Site]    # create colours
  if (zoom == TRUE) {
    plot(compare$LOSO , compare$LOO, xlim = c(0, max(compare[ ,1:2])), ylim = c(0, max(compare[ ,1:2])), pch = 16, cex = 2, col= loso_site, ylab="Bootstrapped cross validation (RMSEP)", xlab = "Leave one site out cross validation (RMSEP)", las = 1)    # plot loso
    if (text == TRUE) {
      text (compare$LOSO , compare$LOO, labels = name, pos = 1,font =2, cex = 0.8)    # add text to plot
    }
    if (text == FALSE) {
    legend(0, max(compare[,1:2]), name, col = colr, pch = 20, pt.cex=2)    #  add legend
    }
  }
  if (zoom == FALSE) {
    plot(compare$LOSO , compare$LOO, xlim = c(0, 65), ylim = c(0, 65), pch = 16, cex = 2, col= loso_site, ylab="Bootstrapped cross validation (RMSEP)", xlab = "Leave one site out cross validation (RMSEP)", las = 1)    # plot loso
    if (text == TRUE) {
      text (compare$LOSO , compare$LOO, labels = name, pos = 1,font =2, cex = 0.8)    # add text to plot
    }
    if (text == FALSE) {
      legend(0, 65, name, col= colr, pch = 20, pt.cex=2)    #  add legend
    }
  }
  
  abline(0,1, lty=2)    # add 1:1 line
  op <- par(font=2)
  t1 <- c(paste(reg, " WAPLS", sep = " "))
  mtext(t1, font = 2, side=3, outer=TRUE, line=-3) # add main title
}
