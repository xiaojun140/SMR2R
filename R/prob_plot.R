#' Plot a figure for probe
#'
#' @param probeinfobuf probeinfobuf
#' @param k k
#' @param colplot colplot
#' @param x.min x.min
#' @param x.max x.max
#' @param y.min y.min
#' @param y.max y.max
#' @param pchbuf pchbuf
#' @param heidi heidi
#'
#' @export plot_probe
plot_probe = function(probeinfobuf, k, colplot, x.min, x.max, y.min, y.max,pchbuf,heidi) {
  # parameters for plot
  genemove = 0.01; txt=1.1;  cex =1.3; lab=1.1; axis=1; top_cex=1.2;
  xcenter = as.numeric(probeinfobuf[k,3])
  pvalbuf = as.numeric(probeinfobuf[k,8])
  strbuf = probeinfobuf[k,1]
  par(new=TRUE)
  if(heidi==TRUE) {
    plot(xcenter, pvalbuf, ylim=c(y.min,y.max),  xlim=c(x.min,x.max),cex.axis=axis,
         xlab="", ylab="", col=colplot, bg=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
  } else {
    plot(xcenter, pvalbuf, ylim=c(y.min,y.max),  xlim=c(x.min,x.max),cex.axis=axis,
         xlab="", ylab="", col=colplot, bty="n", pch=pchbuf, cex=1, axes=F)
  }
}
