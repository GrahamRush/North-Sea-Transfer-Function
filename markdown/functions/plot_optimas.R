
#' 
#' @title site-species tolerances
#'
#' @description creates the neccessary data and plots individual site-species optima
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

species.tol <- function (df, n, y0, y1, x0, x1) {
  tf.wa <- function (x) {
    spec <- subset (x, select = -c(Site, Total.dead, Elevation, SWLI, HoF))
    env <- subset (x, select = c(Site, Total.dead, Elevation, SWLI, HoF))
    fit.wa <- WA(spec,env$SWLI, tolDW = T)
    ### to cross validate the transfer function with leave-one out
    cv <<- crossval(fit.wa, cv.method="loo")
  }
  tf.wa(df)
  #tf.wa.v2(n)
  ns.ts <- training.set [, colSums(spec != 0) > 0]
  ## now set up the data for plotting
  sp.coeficients <<- data.frame(cv$coefficients)
  ## add rows and columns to work on in the plot
  sp.coeficients$names <- rownames(sp.coeficients)
  sp.coeficients <- sp.coeficients[order(sp.coeficients$Optima),]
  n <- nrow(sp.coeficients)
  sp.coeficients$code <- c(1:n)
  ##prepare data for plotting
  y.ns <- sp.coeficients$Optima
  se.y.ns <- sp.coeficients$Tolerances
  x.ns <- sp.coeficients$code
  sp.coeficients.ns <- sp.coeficients
  sp.coeficients.ns$Original <- sp.coeficients.ns$Optima
  ns.tols <<- sp.coeficients.ns[,3:5]
  plot(x.ns,y.ns, pch=3,lwd=2, yaxt = "n", xaxt="n", xlab = "", ylab = "", ylim = c(y0, y1), xlim = c(x0, x1), bty = "n")
  # mtext(side = 3, "WA Inv TOL", font=2, line = -1, cex =1.5, lwd =2)
  axis(2, lwd = 2, font = 2, labels=c(seq(from=y0,to=y1,by=20)), at=c(seq(from=y0,to=y1,by=20)),las=1)
  mtext(side = 2, "SWLI", font=2, line = 3, cex =1,lwd =2)
  ## To plot on top of grey bar
  arrows(x.ns, y.ns+se.y.ns, x.ns, y.ns-se.y.ns, code=3, angle=90, length=0, lwd = 4, col = "#D2D2D2")
  # add text
  text(seq(x0, x1, by=1), y=65, adj = c(1,0), labels = sp.coeficients.ns$names, srt = 90, xpd = TRUE, cex = 0.7)
}


#'
#' @title site-species tolerances
#'
#' @description creates the neccessary data and plots individual site-species optima
#'
#' @param site is the site in ""
#' @param d is the colour to plot
#' @param n is the model of choice i.e3 is wa.inv.tol
#' @param y0 the y axis start
#' @param y1 the y axis end
#' @param x0 the x axis start
#' @param x0 the x axis end
#' @param ex is logical to return table
#'
#' @return adds site optima to the species tolerances plot
#'
#' @keywords tolerances
#' @export
#' @examples
#' ind.site.tol("Alnmouth", "#B0E2FF", 1, 0, 38, 60, 240)

ind.site.tol <- function (site, d, n, x0, x1, y0, y1, ex) { # 's' is the site e.g. "Alnmouth". 'd' is the colour to plot
  specie <- training.set[training.set$Site %in% c(site),]
  y <- specie [, colSums(specie != 0) > 0]
  tf.wa <- function (x) {
    spec <<- subset (x, select = -c(Site, Total.dead, Elevation, SWLI, HoF))
    env <<- subset (x, select = c(Site, Total.dead, Elevation, SWLI, HoF))
    fit.wa <- WA(spec,env$SWLI, tolDW = T)
    ### to cross validate the transfer function with leave-one out
    cv.loo <<- crossval(fit.wa, cv.method="loo")
  }
  tf.wa(y)
  #tf.wa.v2(n) # make sure that you clear the loso f
  sp.coeficients <- data.frame(cv.loo$coefficients)
  ## add rows and columns to work on in the plot
  sp.coeficients$names <- rownames(sp.coeficients)
  sp.coeficients <- sp.coeficients[order(sp.coeficients$Optima),]
  sp.coeficients <- merge(sp.coeficients,ns.tols, all=FALSE)
  ## prepare data for plotting
  y <- sp.coeficients$Optima
  se.y <- sp.coeficients$Tolerances
  x <- sp.coeficients$code
  par(new=T)
  plot(x, y, pch=3, lwd=2,col= d, ylim = c(y0, y1), xlim = c(x0, x1), ylab ="", xlab = "", axes = F)
  out <<- subset(sp.coeficients, Optima < 60 | Optima > 260)
  sp.coeficients$diff <- sp.coeficients$Optima - sp.coeficients$Original
  if (ex==FALSE) return()
  if (ex==TRUE) return(sp.coeficients)
}


#'
#' @title WAPLS site-species tolerances
#'
#' @description creates the neccessary data and plots individual site-species optima
#'
#' @param site is the site in ""
#' @param d is the colour to plot
#' @param n is the model of choice i.e 3 is wa.inv.tol
#' @param y0 the y axis start
#' @param y1 the y axis end
#' @param x0 the x axis start
#' @param x0 the x axis end
#' @param ex is logical to return table
#'
#' @return adds site optima to the species tolerances plot
#'
#' @keywords tolerances
#' @export
#' @examples
#' ind.site.tol.wapls("Alnmouth", "#B0E2FF", 1, 0, 38, 60, 240)

ind.site.tol.wapls <- function (site, d, n, x0, x1, y0, y1, ex) { # 's' is the site e.g. "Alnmouth". 'd' is the colour to plot
  specie <- training.set[training.set$Site %in% c(site),]
  y <- specie [, colSums(specie != 0) > 0]
  tf.wa <- function (x) {
    spec <<- subset (x, select = -c(Site, Total.dead, Elevation, SWLI, HoF))
    env <<- subset (x, select = c(Site, Total.dead, Elevation, SWLI, HoF))
    fit.wa <- WAPLS(spec,env$SWLI)
    ### to cross validate the transfer function with leave-one out
    cv.loo <<- crossval(fit.wa, cv.method="loo")
    }
  tf.wa(y)
  #tf.wa.v2(n) # make sure that you clear the loso f
  sp.coeficients <- data.frame(optima = cv.loo$coefficients[,n])
  ## add rows and columns to work on in the plot
  sp.coeficients$names <- rownames(sp.coeficients)
  sp.coeficients
  sp.coeficients <- sp.coeficients[order(sp.coeficients$optima),]
  sp.coeficients <- merge(sp.coeficients,ns.tols, all=FALSE)
  ## prepare data for plotting
  y <- sp.coeficients$optima
  x <- sp.coeficients$code
  par(new=T)
  plot(x, y, pch=3, lwd=2, col= d, ylim = c(y0, y1), xlim = c(x0, x1), ylab ="", xlab = "", axes = F)
  sp.coeficients$diff <- sp.coeficients$optima - sp.coeficients$Original
  diffs <<- sp.coeficients
  if (ex==FALSE) return()
  if (ex==TRUE) return(sp.coeficients)
}


#'
#' @title plot all of the site-species tolerances
#'
#' @description runs the above functions as defined by the user
#'
#' @param model is the site in ""
#' @param d is the colour to plot
#' @param model is the model of choice i.e 3 is wa.inv.tol
#' @param y_0 the y axis start
#' @param y_1 the y axis end
#' @param x_0 the x axis start
#' @param x_1 the x axis end
#' @param wa logical whether to plot WA
#' @param wapls logical whether to plot WA
#' @param component the wapls component to plot
#' @param wa logical whether to plot individual sit updates
#' 
#' @return plots the species optima with various updates
#'
#' @keywords tolerances
#' @export
#' @examples
#' plot_species_tolerances(model = 2, y_0 = 60, y_1 = 280, x_1 = ncol(spec), wa = TRUE, wapls = TRUE, component = 2, individual = TRUE)



### Species tolerances
# plot the site specific tolerances and optima or use "plot.species.tols.revison.R" to plot as regions
plot_species_tolerances <- function (model = 2, y_0 = 60, y_1 = 280, x_1 = ncol(spec), wa = TRUE, wapls = TRUE, component = 2, individual = FALSE) {
  par(xpd=TRUE)
  if (wa == TRUE) { # plot the site optimas using wa
    df <- df.split(training.set, ts.sp)
    species.tol(df = training.set, n = model, y0 = y_0, y1 = y_1, x0 = 1, x1 = x_1)
    # plot the site specific optima, choosing the reuired sites
    for (i in 1:9) {
      if (training_sets$names$northsea[i] %in% name) {
        ind.site.tol(training_sets$names$northsea[i], d = training_sets$cols$northsea[i], n = model, x0 = 1,x1 = x_1, y0 = y_0, y1 = y_1, ex = TRUE)
      }
    }
    mtext("WA Classical", font = 2, side=3, outer=TRUE, line=-3)
    # add a legend
    legend(1, y_1, name, col=colr,pch = 20, pt.cex=2)
  }
  if (wapls == TRUE) { # plot the site optimas using wapls
    df <- df.split(training.set, ts.sp)
    species.tol(df = training.set, n = model, y0 = y_0, y1 = y_1, x0 = 1, x1 = x_1)
    # plot the site specific optima, choosing the reuired sites
    for (i in 1:9) {
      if (training_sets$names$northsea[i] %in% name) {
        ind.site.tol.wapls(training_sets$names$northsea[i], d = training_sets$cols$northsea[i], n = component, x0 = 1, x1 = x_1, y0 = y_0, y1 = y_1, ex = TRUE)
      }
    }
    mtext("WAPLS Component 2", font = 2, side=3, outer=TRUE, line=-3)
    # add a legend
    legend(1, y_1, name, col=colr,pch = 20, pt.cex=2)
  }
  if (individual == TRUE) { # plot the updated differences in individual sites for 1,2 and 3 components
    for (i in 1:9) {
      if (training_sets$names$northsea[i] %in% name) {
        species.tol(training.set, n = model, y0 = -200, y1 = 400, x0 = 1, x1 = x_1)
        ind.site.tol(training_sets$names$northsea[i], d = training_sets$cols$northsea[1],n = model, x0 = 1, x1 = x_1, -200, 400, ex = T)
        ind.site.tol.wapls(training_sets$names$northsea[i], d = training_sets$cols$northsea[3], n = component, x0 = 1, x1 = x_1, y0 = -200, y1 = 400, ex = T)
        ind.site.tol.wapls(training_sets$names$northsea[i], d = training_sets$cols$northsea[6], n = component + 1, x0 = 1, x1 = x_1, y0 = -200, y1 = 400, ex = T)
        mtext(training_sets$names$northsea[i], font = 2, side=3, outer=TRUE, line=-3)
        legend (1, 400, c("WA-Classical", "WAPLS-C2", "WAPLS-C3"), col= c(training_sets$cols$northsea[1],training_sets$cols$northsea[3], training_sets$cols$northsea[6]),pch = 3, pt.cex=1)
      }
    }
  }
}
