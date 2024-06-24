#' The script of plotting SMREffectPlot
#'
#' @param data outcome data of SMR
#' @param trait_name name of outcome
#' @param cisWindow cisWindow
#' @param transWindow transWindow
#' @param pointsize pointsize
#' @examples
#' # see "https://yanglab.westlake.edu.cn/software/smr/#SMRlocusplot19"
#' @export SMREffectPlot
#'
SMREffectPlot = function(data=SMRData, trait_name="",cisWindow=2000, transWindow=5000, pointsize=20){
  # parameters for plot
  genemove = 0.01; txt=1.1;  cex =1.3; lab=1.1; axis=1; top_cex=1.2;
  # parameters for plot
  pch_top = 24; pch_cis = 21; pch_trans = 22
  col_top = "red"; col_cis = "Navy"; col_trans = "green"
  cex_coeff=3/4 * pointsize/15;

  # Extract the probe for plot
  snbuf = which(as.character(data$eQTL[,1])==data$probeID)
  if(length(snbuf)==0) {
    print(paste("ERROR: no eQTL infomation found for probe",data$probeID,"!",sep=""));
    quit();
  }
  plotData = data$eQTL[snbuf,]
  idx=which(is.na(plotData[,5]))
  if(length(idx)>0) plotData=plotData[-idx,]

  # SNPs in common
  snpbuf = Reduce(intersect, list(as.character(plotData[,2]), data$GWAS[,1]))
  plotData = plotData[match(snpbuf, as.character(plotData[,2])),]
  plotGWAS = data$GWAS[match(snpbuf, as.character(data$GWAS[,1])),]
  # Effect size
  snplist = as.character(plotData[,2])
  bZX = as.numeric(as.character(plotData[,3]));
  seZX = as.numeric(as.character(plotData[,4]));
  snpCorr=as.numeric(as.character(plotData[,5]));
  bZY = as.numeric(as.character(plotGWAS[,2]));
  seZY = as.numeric(as.character(plotGWAS[,3]));
  # Limit
  xmin =  min(bZX - seZX, na.rm=T)
  xmax =  max(bZX + seZX, na.rm=T)
  ymin =  min(bZY - seZY, na.rm=T)
  ymax =  max(bZY + seZY, na.rm=T)

  if(xmin>0) xmin = -xmax/2
  if(xmax<0) xmax = -xmin/2
  if(ymin>0) ymin = -ymax/2
  if(ymax<0) ymax = -ymin/2

  # Plots
  par(mar=c(5,6.5,5,1), xpd=FALSE)
  layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(4.5,1), heights=c(1,1))

  # Start to plot
  nsnps = dim(plotData)[1]
  # Split data to cis- and trans-
  idx=which(data$SMR[,1]==data$probeID);
  if(length(idx)!=1)
  {
    print("ERROR: plot file is not correct!");
    quit();
  }
  if(is.na(data$SMR[idx,8]))
  {
    print(paste("ERROR: no SMR reslult for probe",data$probeID,"!",sep=""));
    quit();
  }
  probeChr = as.numeric(as.character(data$SMR[idx,2]))
  probeBP = as.numeric(as.character(data$SMR[idx,3]))
  GeneID =as.character(data$SMR[idx,4])
  idx=match(snplist,data$SNP[,1]);
  snpChr = as.numeric(as.character(data$SNP[idx,2]))
  snpBP = as.numeric(as.character(data$SNP[idx,3]))
  cisIndx = which(probeChr==snpChr & abs(snpBP-probeBP)<cisWindow*1000);
  ncis = length(cisIndx)
  transIndx = which(probeChr!=snpChr | (probeChr==snpChr & abs(snpBP-probeBP)>transWindow*1000));
  ntrans = length(transIndx)
  # Plot the cis-eQTL
  snplist_tmp = snplist[cisIndx]
  maxsnpCorr = snpCorr[cisIndx]
  bZX_tmp = bZX[cisIndx]; seZX_tmp = seZX[cisIndx]; zZX_tmp = bZX_tmp/seZX_tmp;
  bZY_tmp = bZY[cisIndx]; seZY_tmp = seZY[cisIndx]; zZY_tmp = bZY_tmp/seZY_tmp;
  maxid = which.max(zZX_tmp^2)
  maxsnp = snplist[maxid]
  maxsnpCorr = maxsnpCorr^2;
  for( k in 1 : ncis ) {
    # effect sizes
    colbuf = rgb(0, 0, 128/255, maxsnpCorr[k])
    colcir = rgb(0, 0, 1-maxsnpCorr[k]);
    cex = 1
    plot(bZX_tmp[k], bZY_tmp[k], pch=pch_cis, col=colcir, bg=colbuf,
         bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
         cex=cex, xlab="", ylab="", xaxt="n", yaxt="n")
    par(new=TRUE)
    plot(bZX_tmp[k], bZY_tmp[k], pch=20, col=colcir, bg=colbuf,
         bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
         cex=0.1, xlab="", ylab="", xaxt="n", yaxt="n")
    par(new=TRUE)
  }

  # standard error
  # cis eQTL
  for( k in 1 : ncis ) {
    colcir = rgb(0, 0, 1-maxsnpCorr[k]);
    segments(bZX_tmp[k]-seZX_tmp[k], bZY_tmp[k], bZX_tmp[k]+seZX_tmp[k], bZY_tmp[k],
             col=colcir, lwd=0.5+maxsnpCorr[k])
    segments(bZX_tmp[k], bZY_tmp[k]-seZY_tmp[k], bZX_tmp[k], bZY_tmp[k]+seZY_tmp[k],
             col=colcir, lwd=0.5+maxsnpCorr[k])
  }

  # line
  colline = rgb(244/255,164/255,96/255,1)
  bXY = bZY_tmp[maxid]/bZX_tmp[maxid]
  abline(0, bXY, col=colline, lwd=2, lty=2)

  # plot effect size of the top SNP
  colbuf = "white"
  colcir = col_top
  cex=2.3
  par(new=TRUE)
  plot(bZX_tmp[maxid], bZY_tmp[maxid], pch=pch_top, col=colcir, bg=colbuf,
       bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
       cex=cex, xlab="", ylab="", xaxt="n", yaxt="n")
  colbuf = col_top
  colcir = col_top
  cex = 1
  par(new=TRUE)
  plot(bZX_tmp[maxid], bZY_tmp[maxid], pch=pch_top, col=colcir, bg=colbuf,
       bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
       cex=cex, xlab="", ylab="", xaxt="n", yaxt="n")

  # se of the top SNP
  colcir = rgb(1,0,0)
  segments(bZX_tmp[maxid]-seZX_tmp[maxid], bZY_tmp[maxid],
           bZX_tmp[maxid]+seZX_tmp[maxid], bZY_tmp[maxid],
           col=colcir, lwd=1.5)
  segments(bZX_tmp[maxid], bZY_tmp[maxid]-seZY_tmp[maxid],
           bZX_tmp[maxid], bZY_tmp[maxid]+seZY_tmp[maxid],
           col=colcir, lwd=1.5)

  # Plot the trans-eQTLs
  if(ntrans>0) {
    snplist_tmp = snplist[transIndx]
    bZX_tmp = bZX[transIndx]; seZX_tmp = seZX[transIndx]; zZX_tmp = bZX_tmp/seZX_tmp;
    bZY_tmp = bZY[transIndx]; seZY_tmp = seZY[transIndx]; zZY_tmp = bZY_tmp/seZY_tmp;
    par(new=TRUE)
    for( k in 1 : ntrans ) {
      # effect sizes
      colbuf = col_trans;
      colcir = col_trans;
      cex = 1
      plot(bZX_tmp[k], bZY_tmp[k], pch=pch_cis, col=colcir, bg=colbuf,
           bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
           cex=cex, xlab="", ylab="", xaxt="n", yaxt="n")
      par(new=TRUE)
      plot(bZX_tmp[k], bZY_tmp[k], pch=20, col=colcir, bg=colbuf,
           bty="n", xlim=c(xmin, xmax), ylim=c(ymin, ymax),
           cex=0.1, xlab="", ylab="", xaxt="n", yaxt="n")
      par(new=TRUE)
    }
    # standard error
    # trans-eQTL
    for( k in 1 : ntrans ) {
      colcir = col_trans
      segments(bZX_tmp[k]-seZX_tmp[k], bZY_tmp[k], bZX_tmp[k]+seZX_tmp[k], bZY_tmp[k],
               col=colcir, lwd=1)
      segments(bZX_tmp[k], bZY_tmp[k]-seZY_tmp[k], bZX_tmp[k], bZY_tmp[k]+seZY_tmp[k],
               col=colcir, lwd=1)
    }
  }

  # plot the axis
  # x axis
  devbuf = (xmax - xmin)/5
  if(xmax!=0 & xmin!=0) {
    numbuf = min(abs(xmin), abs(xmax))
    if( devbuf > numbuf ) devbuf = numbuf
  }
  numbuf = as.numeric()
  if( xmin < 0 ) numbuf = c(numbuf, -seq(0, abs(xmin), devbuf))
  if( xmax > 0 ) numbuf = c(numbuf, seq(0, xmax, devbuf))
  axis(1, at=numbuf, labels=round(numbuf,2), las=1, cex.axis=axis)
  xmid = (xmax+xmin)/2
  mtext("eQTL effect sizes", side=1, at=xmid, line=3, cex=cex_coeff)
  # y axis
  devbuf = (ymax - ymin)/5
  if(ymax!=0 & ymin!=0) {
    numbuf = min(abs(ymin), abs(ymax))
    if( devbuf > numbuf ) devbuf = numbuf
  }
  numbuf = as.numeric()
  if( ymin < 0 ) numbuf = c(numbuf, -seq(0, abs(ymin), devbuf))
  if( ymax > 0 ) numbuf = c(numbuf, seq(0,ymax,devbuf))
  axis(2, at=numbuf, labels=round(numbuf,3), las=1, cex.axis=axis)
  ymid = (ymax + ymin)/2
  mtext("GWAS effect sizes", side=2, at=ymid, line=4.5, cex=cex_coeff)

  mainstr1 = trait_name
  mainstr2 = substitute(paste(probeid, " (", italic(gene), ")", sep=""),
                        list(probeid=as.character(data$probeID),
                             gene=as.character(GeneID)))
  mtext(mainstr1, side=3, at=xmin, adj=0, line=2.5, cex=cex_coeff)
  mtext(mainstr2, side=3, at=xmin, adj=0, line=0.5, cex=cex_coeff)
  # Plot legend
  lstr = c("top cis-eQTL", "cis-eQTL")
  col_led = c(col_top, col_cis); pch_led = c(pch_top, pch_cis)
  if(ntrans>0) {
    lstr=c(lstr, "trans-eQTL"); col_led = c(col_led, col_trans)
    pch_led = c(pch_led, pch_trans)
  }

  if(bXY>0) {
    legend("topleft", lstr, bty="n", border="white", pch=pch_led, col=col_led, pt.bg=col_led, cex=axis)
  } else {
    legend("topright", lstr, bty="n", border="white", pch=pch_led, col=col_led, pt.bg=col_led, cex=axis)
  }

  # add the scale bar
  par(mar=c(5,1,5,4.5))
  pal=colorRampPalette(c(rgb(0, 0, 1),  rgb(0, 0, 0)))
  breaks <- seq(min(snpCorr^2), max(snpCorr^2),length.out=100)
  imageScale(snpCorr^2, col=pal(length(breaks)-1), breaks=breaks, horiz=FALSE, axis.pos=4, yaxt="n")
  dvd = (max(snpCorr^2) - min(snpCorr^2))/5
  pos = seq(min(snpCorr^2), max(snpCorr^2), dvd)
  axis(4, at=pos, label=sprintf("%.2f", pos), las=2)
  mtext(expression(italic(r)^2), side=4, line=3.5, las=2 )
}

