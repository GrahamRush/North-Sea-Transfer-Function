pam_hockey_stick <- function () {
  plot(x=c(1:10), y = clusters$crit, type = "n", lwd = 1, axes = F, xlab = "", ylab = "", bty = "n")
  axis(1, at=seq(1,10, by=1), lwd=1)
  axis(2, las =2)
  mtext(side = 1, "Number of clusters", line = 3, cex =1)
  mtext(side = 2, "Average sillhouette width", line = 2.5, cex =1)
  lines(x=c(1:10), y = clusters$crit)
  points(x=c(1:10), y = clusters$crit, pch = 16)
}  
