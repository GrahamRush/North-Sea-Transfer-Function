#' A clean fucntion
#'
#' @title Clean
#'
#' @description comparison between LOO and LOSO cross-validation
#'
#' @param comp is the model or component
#' @param a-e are the sites e.g "Alnmouth"
#' @return compare a dataframe with site comparisons
#'
#' @keywords compare
#' @export
#' @examples
#' loso(2, "Alnmouth","Brancaster","Cowpen","Kjelst","Rantum","Sonderho","Thornham","Welwick","Ythan")

loso <- function (comp, a=NULL,b=NULL,c=NULL,d=NULL,e=NULL,f=NULL,g=NULL,h=NULL,i=NULL) {
  looloso <- sapply (unique(training.set$Site), function(n) {
    N <- training.set$Site==n
    c(LOO=sqrt(mean((cv.loo$predict[N,comp] - training.set$SWLI[N])^2)),LOSO=sqrt(mean((cv.loso$predict[N,comp] - training.set$SWLI[N])^2)))
  })
  ## turn matrix into data frame
  compare <- as.data.frame(t(looloso), row.names = c(a,b,c,d,e,f,g,h,i))
  compare$Site <- as.factor(c(a,b,c,d,e,f,g,h,i))
  print(compare)
}
