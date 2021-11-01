#' A clean function for the core
#'
#' @title calculate mat
#'
#' @description calculates the mat of samples
#'
#' @param z is the number of rows in training set -1
#' @param d is the dissimilarity coefficient, nominally set ot square-chord
#'
#' @return file with mat results. good is better than the 5th percentile, close is better than the  20th percentile
#'
#' @keywords mat
#' @export
#' @examples
#' run.mat(185)

run.mat <- function(z = 30, d = "sq.chord", plot = TRUE) { # z is number is the samples -1
  fit <- MAT(spec, env$SWLI, k = (nrow(spec)-1), dist.method = d)
  perc.20 <- quantile(fit$dist.n, c(.20))
  perc.5 <- quantile(fit$dist.n, c(.05))
  fit <- MAT(spec, env$SWLI, k = z, dist.method = d)
  pred.core.mat <- predict(fit, core$spec)
  mat <- ifelse(pred.core.mat$diagnostics$minD < perc.5, "good", ifelse(pred.core.mat$diagnostics$minD > perc.5 & pred.core.mat$diagnostics$minD < perc.20, "close", "poor" ))
  mat.fac <- as.factor(mat)
  sizes <- ordered(mat.fac, levels = c("good", "close", "poor"))
  my.pallete <- c("green", "orange", "red")
  cols.mat <- my.pallete[sizes]
  if (plot == TRUE) {
    x.m <- ifelse (max(pred.core.mat$diagnostics$minD) > perc.20, max(pred.core.mat$diagnostics$minD), perc.20)
    plot(pred.core.mat$diagnostics$minD, core$env$Depth..od., xlim = c(0, x.m), type = "n" )
    abline(v = perc.20, lty = 2)
    abline(v = perc.5, lty = 2)
    lines(pred.core.mat$diagnostics$minD, core$env$Depth..od.)
    points(pred.core.mat$diagnostics$minD, core$env$Depth..od.,  col = cols.mat, pch = 16)
  }
  return(list(mat = mat, cols_mat = cols.mat, value = pred.core.mat))
}
