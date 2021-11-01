#' Save the output
#'
#' @title Save Output
#'
#' @description Saves the important output from the transfer function
#'
#' @param region is the region
#' @param method is the model used
#' @return files using the names supplied above. 1. training set with data. 2.
#'
#' @keywords compare
#' @export
#' @examples
#' saved("northsea", "WAPLS-C2")
#'
saved <- function (region, method)
{
  a1 <- c(paste(method, "-loo", sep = ""))
  a2 <- c(paste(method, "-loso", sep = ""))
  b <- c("NA")
  c <- c("NA")
  d <- c("NA")
  e <- c("NA")
  f <- c("NA")
  g <- c("NA")
  model1 <- data.frame(a1, b, c, d, e, f, g)
  model2 <- data.frame(a2, b, c, d, e, f, g)
  names(model1) <- c("RMSE", "R2", "Avg.Bias", "Max.Bias", "Skill", "delta.RMSE", "p")
  names(model2) <- c("RMSE", "R2", "Avg.Bias", "Max.Bias", "Skill", "delta.RMSE", "p")
  results <- rbind(test.wa, model1, test.wapls.loo, model2, test.wapls.loso)
  return(results)
}
