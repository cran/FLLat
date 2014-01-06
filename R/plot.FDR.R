plot.FDR <- function(x,xlab="Threshold",ylab="FDR",...) {

  plot(x=x$thresh.values,y=x$FDRs,xlab=xlab,ylab=ylab,...)

}
