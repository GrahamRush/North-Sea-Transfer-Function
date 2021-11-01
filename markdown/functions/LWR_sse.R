#' @title LWR sse
#'
#' @description function to run locally weighted TF and extract sample specific errors
#'
#' @param component is the wapls component or wa model
#' @param tf is the transfer function. WA or WAPLS
#' @param n is the number of closest analogues. Usually between 30 and 50
#'
#' @return list with the fit and sample specific errro for each fossil sample
#'
#' @keywords sse
#' @export
#' @examples
#' lw_tf(component = 1, tf = WAPLS, n = 50)

lw_sse <- function (component = 1, tf = WAPLS, n = 50) {
  fit.lw <- LWR(spec, env$SWLI, FUN = tf, dist.method = "sq.chord", k = n, tolDW = T)    # run the locally weighted TF
  p.lw <- predict(fit.lw, spl.c1$spec, sse = F, nboot = 100)    # run locally weighted TF to create the names
  lw.see <- vector()    # create empty vector
  lw.swli <- vector()    # create empty vector
  lw.name <- matrix(ncol = n, nrow = nrow(spl.c1$spec))
  for (i in 1:nrow(spl.c1$spec)) {    # loop to run on all rows of the core
    analogues <- p.lw$match.name[i,]    # extract the k analogues
    spec.lw <- spec[rownames(spec) %in% analogues, ]    # extract the species data for the associated analogues
    spec.lw.zero <- spec.lw[, colSums(spec.lw != 0) > 0]    # remove species with zero abundances
    env.lw <- env[rownames(env) %in% keep_rows, ]    # extract the enviro data for the associated analogues
    fit.wapls <- tf(spec.lw.zero, env.lw$SWLI)    # run the TF for the new data
    core.wapls <- predict(fit.wapls, spl.c1$spec[i,], sse=T, nboot=1000)    # predict the sample from the TF
    lw.see[i] <- core.wapls$SEP.boot[, component]
    lw.swli[i] <- core.wapls$fit[, component]
    lw.name[i,] <- analogues
  }
  return(list(fit = lw.swli, sse = lw.see))
}


 # function to run locally weighted TF and extract sample specific errors
  lw_tf <- function (component = 1, tf = tf, n = 50) {
    fit.lw <- LWR(spec, env$SWLI, FUN = tf, dist.method = "sq.chord", k = n)    # run the locally weighted TF
    p.lw <- predict(fit.lw, core.spec, sse = F, nboot = 100)    # run locally weighted TF to create the names
    lw.see <- vector()    # create empty vector
    lw.swli <- vector()    # create empty vector
    lw.name <- matrix(ncol = n, nrow = nrow(core.spec))
    lw.swli.mod <- matrix(ncol = m, nrow = nrow(core.spec))
  lw.tf <-  for (i in 1:nrow(core.spec)) {    # loop to run on all rows of the core
      analogues <- p.lw$match.name[i,]    # extract the k analogues
      spec.lw <- spec[rownames(spec) %in% analogues, ]    # extract the species data for the associated analogues
      env.lw <- env[rownames(env) %in% analogues, ]    # extract the enviro data for the associated analogues
      spec.lw.zero <- spec.lw[, colSums(spec.lw != 0) > 0]    # remove species with zero abundances
      fit.wapls <- WAPLS(spec.lw.zero, env.lw$SWLI)    # run the TF for the new data
      core.wapls <- predict(fit.wapls, core.spec[i,], sse=T, nboot=1000)    # predict the sample from the TF
      lw.see[i] <- core.wapls$SEP.boot[, component]
      lw.swli[i] <- core.wapls$fit[, component]
      lw.name[i,] <- analogues
      lw.swli.mod[i,] <- env.lw$SWLI
    }
    return(list(fit = lw.swli, sse = lw.see, analogues = lw.name, analogue.swli = lw.swli.mod))
  }
  
  
analogues_sensitivty <- function(k){
  analogue_sensitivty <- list()
  for(i in c(k)) { #Loop through the numbers of ID's instead of the ID's
    dat <-  run_lw_tf(n_analogues = i)
    name <- paste('k = ',i,sep='')
    analogues_sensitivty[[name]] <- dat
  }
  return(analogue_sensitivty) #Return the list of dataframes.
}

