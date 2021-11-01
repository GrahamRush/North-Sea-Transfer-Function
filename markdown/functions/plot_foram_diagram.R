#' @title foram diagram
#'
#' @description function to produce a foram diagram locally weighted TF and extract sample specific errors
#'
#' @param data is the training set
#' @param min.species is the minimum % for which to show the species number of closest analogues. Usually between 30 and 50
#' @param clusters is the desired number of clusters
#' @param moethod is the methdo by whihc ot cluster
#' @param plot_zones add zones to the plot
#' @param zones if plot_zones is TRUE then set the swli value for the zones
#' @param save_plot logical to save the plot
#' 
#' @return list with the fit and sample specific errro for each fossil sample
#'
#' @keywords sse
#' @export
#' @examples
#' plot_foram_diagram(data = training.set, min.species = 10, clusters = 3, method = "euclidean", plot_zones = TRUE,  zones = c(200, 150, 100), save_plot = TRUE)
#' 

  
plot_foram_diagram <- function (data = training.set, min.species = 10, clusters = 3, method = "euclidean", plot_zones = TRUE,  zones = c(200, 150, 100), save_plot = FALSE) {

  ## function to order the data by optima and create labels
  species.tol.clu.rev <- function (df, max.species) {     # m = the min percentage to retain
    
    tf.wa <- function (x) {
      spec <- subset (x, select = -c(Site, Total.dead, Elevation, SWLI, HoF))
      env <- subset (x, select = c(Site, Total.dead, Elevation, SWLI, HoF))
      fit.wa <- WA(spec,env$SWLI, tolDW = T)
      cv <- crossval(fit.wa, cv.method="boot", nboot = 1000)      # to cross validate the transfer function with leave-one out
      return(cv)
    }
    
    cv <- tf.wa(df)    # run the nested function
    df.split(training.set, ts.sp)
    env$colour <- colr[env$Site]

    ## now set up the data for plotting
    sp.coeficients <- data.frame(cv$coefficients)
    
    ## add rows and columns to work on in the plot
    sp.coeficients$names <- rownames(sp.coeficients)    # add names
    sp.coeficients <- sp.coeficients[order(sp.coeficients$Optima),]     # order by value
    sp.coeficients$code <- c(1:nrow(sp.coeficients))    # give a code to refer to later
    ords <- as.character(sp.coeficients$names)
    spec1 <- spec[, ords]    # return ordered dataframe 
    ## Remove taxa where total abundanace is less than 10%
    s.max <- apply(spec1, 2, function(x) max(x > max.species, na.rm=TRUE)) # calculate which species have a max abundace > n %
    spec2 <- spec1[, which(s.max == 1)] # Subset data where 1 means species < 10 %
    ## calculate names
    nam <- colnames (spec2)
    nam4 <- as.character()    # create empty vector
    for (i in 1:ncol(spec2)) {     # loop to remove "." from species names 
      nam2 <- unlist(strsplit(nam[i], "[.]"))
      nam3 <- paste(nam2[1], nam2[2]) 
      nam4[i] <- nam3
    }
    nam5 <- c("SWLI", nam4)
    
    return(list(tols = spec2,
                labs = nam5))
  }
  
  species_tolerance <- species.tol.clu.rev(df = data, max.species = min.species)    # run nested function
  

  ## function to create dataframe to order and colour plotting
  spec.ab <- function (clusts, dist.m) { 
    df.split(training.set, ts.sp)
    env$colour <- colr[env$Site]
    dis <- vegdist(spec, dist.m) # Set this e.g. "bray" is bray curtis
    # Test the suggested number of clusters
    clusters <- pamk(dis, krange = clusts, criterion = "multiasw", critout = FALSE, metric = dist.m) # krange sets the number of clusters to be tested. critout prints all of the results
    spec4 <- cbind(species_tolerance$tols, env) # join the datasets
    # arrange and join cluster data
    pam.name <- row.names(clusters$pamobject$silinfo$widths)
    pam.ord <- as.data.frame(clusters$pamobject$silinfo$widths)
    pam.ord$sample <- pam.name
    spec4$sample <- row.names(spec4)
    spec6 <- join(pam.ord, spec4) # join to set by pam order
    spec6 <- spec6[order(spec6$SWLI),] # order by SWLI
    spec6$order <- c(1:nrow(spec6)) # add code for plotting order
    
    spec7 <- cbind(spec6$SWLI / 2, spec6[,5:(3 + length(species_tolerance$labs))])    # join with SWLI for ploting and create a swli equivalent to 100
    spec7$pam1 <- as.numeric(if_else(spec6$cluster == 1, 10, 0))
    spec7$pam2 <- as.numeric(if_else(spec6$cluster == 2, 10, 0))
    if (clusts > 2) {
    spec7$pam3 <- as.numeric(if_else(spec6$cluster == 3, 10, 0))
    }
    if (clusts > 3) {
    spec7$pam4 <- as.numeric(if_else(spec6$cluster == 4, 10, 0))
    }
    if (clusts > 4) {
    spec7$pam4 <- as.numeric(if_else(spec6$cluster == 5, 10, 0))
    }
    
    return(list (spec6 = spec6, spec7 = spec7))
  }
  
  species_order <- spec.ab(clusts = clusters, dist.m = method)    # run nested function


  ## fucntion to make the plot
  plot_foram <- function (zones,  z, save) {
  
    plot.new()
    par(fig=c(0,1,0,1), new=TRUE, xpd = TRUE) # , mar = c(0, 2, 4, 1))
    
    swli.line <- c("black", rep(NA, times = ncol(species_order$spec7) - 1))    # create swli colours
  
    if (save == TRUE) {
      file_name <- paste("results/", reg, "_foram_diagram.pdf", sep = "")
      pdf(file = file_name, width = 26/2.54, height = 19/2.54, useDingbats = FALSE, pointsize = 7)
    }
    
    for.plot <- strat.plot(species_order$spec7[, 2:ncol(species_order$spec7)], yvar = species_order$spec6$SWLI, plot.line = FALSE, plot.poly = FALSE, plot.bar = TRUE, cex.xlabel = 0.4, col.bar = species_order$spec6$colour, lwd.bar = 1, sep.bar = TRUE, scale.minmax  =  F, scale.percent = T, xSpace = 0.01, x.pc.inc = 10, x.pc.lab = TRUE, x.pc.omit0 = TRUE, las = 2, srt.xlabel = 45, y.axis = TRUE, x.names = species_tolerance$labs[2:length(species_tolerance$labs)])    # plot against SWLI
    
    if (zones == TRUE) {    # add zones
      addZone(for.plot, upper = z, lty = 2)
    }
    
    legend(x = "topleft", legend = name, col=colr, pch = 15)    # set to appropriate
  
    if (save == TRUE) {
      dev.off()
    }
    
  } 
 
  plot_foram(zones = plot_zones,  z = zones, save = save_plot)    # run nested function
  
}





