#' test of goodness of fit
#'
#' @title goodness of fit 
#'
#' @description calculates the goodness of fit following Telford and Kemp (2015)
#'
#' @param z is the number of rows in training set -1
#' @param plot is logical to plot the reuslts
#'
#' @return file with mat results. good is better than the 5th percentile, close is better than the  20th percentile
#'
#' @keywords mat
#' @export
#' @examples
#' run.mat(185)

goodness_of_fit <- function(plot = TRUE) {
  res.length <- analogue::residLen(spec, env$SWLI, core$spec, method = c("cca"))
  g0f.perc <- quantile(res.length$train,   c(.90, 0.95))
  gof <- ifelse(res.length$passive < g0f.perc[1], "good", ifelse(res.length$passive > g0f.perc[1] & res.length$passive < g0f.perc[2], "close", "poor"))
  gof.fac <- as.factor(gof)
  sizes <- ordered(gof.fac, levels = c("good", "close", "poor"))
  my.pallete <- c("green","orange","red")
  cols.gof <- my.pallete[sizes]
  if (plot == TRUE) {
    plot(res.length$passive, core$env$Depth..od., type = "n")
    abline(v= res.length$train[1], lty = 2)
    abline(v= res.length$train[2], lty = 2)
    points(res.length$passive, core$env$Depth..od., col = cols.gof, pch = 16)
    lines(res.length$passive, core$env$Depth..od.)
  }
  return(gof)
}
