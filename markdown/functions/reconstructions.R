#' @title core reconstructions
#'
#' @description function to reconstruct sea level for the core samples
#'
#' @param core is the core data
#' @param tf is the transfer function
#' @param samples determines whether to manually set the swli
#' @param b is the value for above samples
#' @param species is the final row of species data if the core needs splitting
#' @param mhhw is the palaeo mhhw value 
#' @param mtl is the palaeo mtl value 
#' @param mod is the name of the model used 
#'
#' @return data frame with sea-level predictions
#'
#' @keywords sse
#' @export
#' @examples
#' reconstruct_core(core = core, tf = fit.wapls, samples = 18:21, b = 230, components = 2, species = 15, mhhw = 1.95, mtl = 0.3, mod = "wapls.w", plot = TRUE)

reconstruct_core <- function (core, tf, samples = NULL, b = NULL, components, change.sp = FALSE, species, mhhw, mtl, mod, plot = FALSE) {
  core.pred <- predict(tf, core$spec, sse=TRUE, nboot=1000)
  if (change.sp == TRUE) {
   core.pred$fit [samples,1:4] <- c(b)    # add the swli for the peat layer as min e.g 230.
   }
  reconstruction <- ( ( (core.pred$fit[,components] - 100) * (mhhw - mtl) ) / 100) + mtl     # convert swli into (sea level m od)
  sea.level.pred <- (core$env$Depth..od. / 100) - ( ( ( (core.pred$fit[,components] - 100) * (mhhw - mtl) ) / 100) + mtl )
  range <- (mhhw - mtl)/100*core.pred$SEP.boot[,components]     # [.,components] corresponds to the number of components chosen
  range <- round (range, digits = 2)
  if (plot == TRUE) {
    y <- sea.level.pred
    x <- core$env$Depth..od.
    plot(x, y, xlab = "sea level (m)", ylab = "Depth (m)", pch = 16)
    arrows(x0 = core$env$Depth..od., y0 = sea.level.pred + range, y1 = sea.level.pred - range, angle = 90, code = 3, length = 0.05)
  }
  sl <- data.frame("s1" = round(sea.level.pred, 2), "s2" = round(range, 2), "s3" = round(core.pred$fit[,components], 2), "s4" =  round(reconstruction, 2), "s5" = round(core.pred$SEP.boot[,2], 2))
  s1 <- c(paste(mod, ".rsl", sep = ""))
  s2 <- c(paste(mod, ".SSE.m", sep = ""))
  s3 <- c(paste(mod, ".sl.swli", sep = ""))
  s4 <- c(paste(mod, ".sl.m", sep = ""))
  s5 <- c(paste(mod, ".SSE.swli", sep = ""))
  sl2 <- setNames (sl, c(s1, s2, s3, s4, s5))
  return(sl2)
}


#' @title LWR sse
#'
#' @description function to run locally weighted TF and extract sample specific errors
#'
#' @param component is the wapls component or wa model
#' @param tf is the transfer function. WA or WAPLS
#' @param n_analogues is the number of closest analogues. Usually between 30 and 50
#' @param samples determines whether to manually set the swli
#' @param b is the value for above samples
#' @param species is the final row of species data if the core needs splitting
#' @param mhhw is the palaeo mhhw value 
#' @param mtl is the palaeo mtl value 
#' @param mod is the name of the model used 
#'
#' @return list with the fit and sample specific error for each fossil sample
#'
#' @keywords sse
#' @export
#' @examples
#' run_lw_tf(n_analogues = 30, components = 2, tf = WAPLS, samples = 18:21, b = 230, species = 15, mhhw = 1.95, mtl = 0.3, mod = "lwr.w", method = "sq.chord")

# function to predict sea level for Locally weighted transfer fucntions
run_lw_tf <- function (n_analogues = 30, components = 2, tf = WAPLS, samples = 18:21, b = 230, species = 15, mhhw = 1.95, mtl = 0.3, mod = "lwr.w", plot = TRUE, method = "sq.chord") {
  
  lw_tf <- function (cpt = components, trans.f = tf, m = n_analogues) {    # function to run locally weighted TF and extract sample specific errors
    fit.lw <- LWR(spec, env$SWLI, FUN = WAPLS, dist.method = method, k = m, tolDW = T)    # run the locally weighted TF
    p.lw <- predict(fit.lw, core$spec, sse = F, nboot = 100)    # run locally weighted TF to create the names
    lw.see <- vector()    # create empty vector
    lw.swli <- vector()    # create empty vector
    lw.name <- matrix(ncol = m, nrow = nrow(core$spec))
    lw.site.mod <- data.frame(matrix(ncol = m, nrow = nrow(core$spec)))
    names <- as.character()
    for (i in 1:m) {
      names[i] <- paste("analogue", i, sep = "")
    }
    colnames(lw.site.mod) <- c(names)
    rownames(lw.site.mod) <- c(rownames(core$spec))
    lw.swli.mod <- lw.site.mod
    lw.site.mod[] <- lapply(lw.site.mod, as.character)
    lw.tf <-  for (i in 1:nrow(core$spec)) {    # loop to run on all rows of the core
        analogues <- p.lw$match.name[i,]    # extract the k analogues
        spec.lw <- spec[rownames(spec) %in% analogues, ]    # extract the species data for the associated analogues
        env.lw <- env[rownames(env) %in% analogues, ]    # extract the enviro data for the associated analogues
        spec.lw.zero <- spec.lw[, colSums(spec.lw != 0) > 0]    # remove species with zero abundances
        fit.wapls <- WAPLS(spec.lw.zero, env.lw$SWLI)    # run the TF for the new data
        core.pred <- predict(fit.wapls, core$spec[i,], sse=T, nboot=1000)    # predict the sample from the TF
        lw.see[i] <- core.pred$SEP.boot[, cpt]
        lw.swli[i] <- core.pred$fit[, cpt]
        lw.name[i,] <- analogues
        lw.swli.mod[i,] <- env.lw$SWLI
        lw.site.mod[i,] <- as.factor(env.lw$Site)
      }
    return(list(fit = lw.swli, sse = lw.see, analogues = lw.name, analogue.swli = lw.swli.mod, analogue.site = lw.site.mod))
  }

  lw.se <- lw_tf(m = n_analogues)    # run the function
  
  # create a table showing the number of analogues per site
  # sites <- as.data.frame(table(lw.se$analogue.site)) # Levels: Alnmouth Brancaster Cowpen Thornham Welwick Ythan
  # rownames(sites) <- w.names
  # sites
  
  lw.se$fit[samples] <- c(b)    # change the last column if applying b
  reconstruction <- (((lw.se$fit - 100) * (mhhw - mtl)) / 100) + mtl    # convert SWLI to sea level 
  sea.level.pred <- (core$env$Depth..od. / 100) - ((((lw.se$fit - 100) * (mhhw - mtl)) / 100) + mtl)    # convert to sea level
     # number in [,!] is the chosen transfer function model. The later numbers are the tidal datums to convert swli back to elevation 1.95 = HAT, 0.3 = MTL
  range <- round((mhhw - mtl) / 100 * lw.se$sse, digits = 2) # [.,1] corresponds to the number of components chosen
  sl <- data.frame("s1" = round(sea.level.pred, 2), "s2" = round(range, 2), "s3" = round(lw.se$fit, 2), "s4" =  round(reconstruction, 2), "s5" = round(lw.se$sse, 2))
  s1 <- c(paste(mod, n_analogues, "rsl", sep = "-"))
  s2 <- c(paste(mod, n_analogues, "SSE-m", sep = "-"))
  s3 <- c(paste(mod, n_analogues, "sl-swli", sep = "-"))
  s4 <- c(paste(mod, n_analogues, "sl-m", sep = "-"))
  s5 <- c(paste(mod, n_analogues, "SSE-swli", sep = "-"))
  sl.lw <- setNames (sl, c(s1, s2, s3, s4, s5))
  
  if (plot == TRUE) {
    y <- sea.level.pred
    x <- core$env$Depth..od.
    plot(x, y, xlab = "sea level (m)", ylab = "Depth (m")
    arrows(x0 = core$env$Depth..od., y0 = sea.level.pred + range, y1 = sea.level.pred - range, angle = 90, code = 3, length = 0.05)
  }
  
  return(list(sea_level = sl.lw, analysis = lw.se))
}


#' @title LWR sesnisitivity test
#'
#' @description function to run locally weighted TF and test the sesnitivity to k closest analogues
#'
#' @param k is a vector of number of closest analogues. Usually between 30 and 50
#' @param plot is logical to plot the results 
#'
#' @return list with the fit and sample specific error for each fossil sample
#'
#' @keywords sse
#' @export
#' @examples
#' run_lw_tf(n_analogues = 30, components = 2, tf = WAPLS, samples = 18:21, b = 230, species = 15, mhhw = 1.95, mtl = 0.3, mod = "lwr.w")

analogue_sensitivity_test <- function (k = c(30, 40, 50), plot = TRUE) {
  analogue_sensitivity <- function (m = k) {
    analogues_sensitivty <- list()
    for(i in m) {    # Loop through the numbers of ID's instead of the ID's
      dat <-  run_lw_tf(n_analogues = i)
      name <- paste('k = ',i,sep='')
      analogues_sensitivty[[name]] <- dat$sea_level
    }
    return(analogues_sensitivty)   # Return the list of dataframes.
  }
  analogues_sensitivty_test <- analogue_sensitivity(m = k)
  if (plot == TRUE) {
    analogues_sensitivity_sse <- lapply(analogues_sensitivty_test, '[', 1:17, 2)
    analogues_sensitivity_sse <- as.data.frame(matrix(unlist(analogues_sensitivity_sse),nrow=length(analogues_sensitivity_sse),byrow=TRUE))
    
    analogues_sensitivity_sse_mf <- lapply(analogues_sensitivty_test, '[', 1:10, 2) # extract the mudflat samples
    analogues_sensitivity_sse_mf <- as.data.frame(matrix(unlist(analogues_sensitivity_sse_mf),nrow=length(analogues_sensitivity_sse_mf),byrow=TRUE))
    
    analogues_sensitivity_sse_sm <- lapply(analogues_sensitivty_test, '[', 11:17, 2) # extract the salt marsh samples
    analogues_sensitivity_sse_sm <- as.data.frame(matrix(unlist(analogues_sensitivity_sse_sm),nrow=length(analogues_sensitivity_sse_sm),byrow=TRUE))
    
    # extract the means of each k test
    all_mean <- apply(analogues_sensitivity_sse, 1, mean)
    mf_mean <- apply(analogues_sensitivity_sse_mf, 1, mean)
    sm_mean <- apply(analogues_sensitivity_sse_sm, 1, mean)
    
    # plot the results
    plot(x = k, y = all_mean, type = "p", pch = 16, ylim = c(0.2, 0.4), ylab = "Mean sample specific error (m)", xlab = "Number of closest analogues")
    points(x = k, y = mf_mean, pch = 17)
    points(x = k, y = sm_mean, pch = 15)
    text(x = rep(46, 3), y = c(0.25, 0.32, 0.38), labels = c("Salt-marsh samples", "All samples", "Tidal-flat samples"), pos = 4)
  }
  return(analogues_sensitivty_test) #Return the list of dataframes.
}
  
