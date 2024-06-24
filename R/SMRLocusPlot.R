
#' The script of plotting SMR Locus Plot
#'
#' @param data outcome data of SMR
#' @param probeNEARBY probe
#' @param smr_thresh smr_thresh
#' @param smr_thresh_plot smr_thresh_plot
#' @param heidi_thresh thresh of heidi
#' @param plotWindow plotWindow
#' @param pointsize pointsize
#' @param max_anno_probe max_anno_probe
#' @param anno_selfdef anno_selfdef
#' @examples
#' # see "https://yanglab.westlake.edu.cn/software/smr/#SMRlocusplot19"
#'
#' @export SMRLocusPlot

SMRLocusPlot = function(data=SMRData, probeNEARBY=NULL,smr_thresh=NULL, smr_thresh_plot=NULL, heidi_thresh=NULL, plotWindow=NULL,pointsize=20,max_anno_probe=16,anno_selfdef=TRUE)
{
  genemove = 0.01; txt=1.1;  cex =1.3; lab=1.1; axis=1; top_cex=1.2;
  cex_coeff=3/4 * pointsize/15;
  if(length(smr_thresh)==0){
    print("ERROR: please specify the threshold of SMR test!");
    quit();
  }
  if(length(heidi_thresh)==0){
    print("ERROR: please specify the threshold of HEIDI test!");
    quit();
  }
  if(length(plotWindow)==0){
    print("ERROR: please specify the plot window size!");
    quit();
  }
  if(length(which(is.na(data$SMR[,3])))>0)
  {
    print("ERROR: Some probes' physical positon is missing!");
    quit();
  }
  idx=match(data$probeID,data$SMR[,1]);
  if(length(idx)==0){
    print("ERROR: Plot file is not generated correctly, can't find target probe!");
    quit();
  }
  if(length(smr_thresh_plot)==0){
    smr_thresh_plot=smr_thresh;
  }
  cis_start=data$SMR[idx,3]-plotWindow*1000;
  if(cis_start<0) cis_start=0
  cis_end=data$SMR[idx,3]+plotWindow*1000;
  idx=which(data$SMR[,3]>=cis_start & data$SMR[,3]<=cis_end)
  data$SMR=data$SMR[idx,]
  idx=match(data$GWAS[,1],data$SNP[,1])
  tmpsnpbp=data$SNP[idx,3]
  idx=which(tmpsnpbp>=cis_start &tmpsnpbp<=cis_end)
  data$GWAS=data$GWAS[idx,]
  idx=match(data$eQTL[,2],data$SNP[,1])
  tmpsnpbp=data$SNP[idx,3]
  idx=which(tmpsnpbp>=cis_start &tmpsnpbp<=cis_end)
  data$eQTL=data$eQTL[idx,]

  if(!is.null(data$Gene))
  {
    idx=which(data$Gene[,2]>=cis_start & data$Gene[,3]<=cis_end )
    data$Gene=data$Gene[idx,]
  }

  #start to plot
  smrindx = which(data$SMR[,8] <= smr_thresh_plot)
  #heidiindx = which((data$SMR[,8] <= smr_thresh_plot) & (data$SMR[,9] >= heidi_thresh_plot))
  smrprobes = NULL; heidiprobes = NULL;
  if(length(smrindx)>0) { smrprobes =  as.character(data$SMR[smrindx,1]) }
  #if(length(heidiindx)>0) { heidiprobes = as.character(data$SMR[heidiindx,1]) }

  smrindx_bonferr = which(data$SMR[,8] <= smr_thresh)
  heidiindx_strengent = which((data$SMR[,9] >= heidi_thresh))
  smrprobes_red = NA; heidiprobes_solid = NA;
  if(length(smrindx_bonferr)>0) { smrprobes_red =  as.character(data$SMR[smrindx_bonferr,1]) }
  if(length(heidiindx_strengent)>0) { heidiprobes_solid = as.character(data$SMR[heidiindx_strengent,1]) }

  if(length(probeNEARBY)>0)
  {
    idx=match(probeNEARBY,data$SMR[,1])
    idxx=which(is.na(idx))
    if(length(idxx)>0)
    {
      for(ii in 1:length(idxx)) {
        print(paste("WARNING: cann't find probe ",probeNEARBY[idxx[ii]], " in plot region.",sep=""))
      }
      probeNEARBY=probeNEARBY[-idxx]
    }

  }
  probePLOT=smrprobes #draw the eQTL of all the probes that passed smr_thresh_plot
  probePLOT=unique(c(data$probeID,probePLOT,probeNEARBY)) # draw the target probe anyway
  nprobePLOT = length(probePLOT)

  idx=which(is.na(data$GWAS[,2]) | is.na(data$GWAS[,3]))
  if(length(idx)>0) data$GWAS=data$GWAS[-idx,]
  pZY=-log10(pchisq((data$GWAS[,2]/data$GWAS[,3])^2,1,lower.tail=F))

  idx=match(data$probeID,data$SMR[,1]);
  if(length(idx)>0){
    chrPLOT = data$SMR[idx,2]
  }else{
    print("ERROR: Plot file is not generated correctly, please report this bug!");
    quit();
  }
  idx=which(is.na(data$SMR[,8]) )
  if(length(idx)>0) {
    probeINFO=data$SMR[-idx,];
  }else{
    probeINFO=data$SMR;
  }
  idx=which(is.na(probeINFO[,5]) | is.na(probeINFO[,6]));
  idx2=which(is.na(probeINFO[,3]));
  if(length(intersect(idx,idx2))>0)
  {
    print("ERROR: Some probes' physical positon is missing!");
    quit();
  }
  probeINFO[idx,5]=probeINFO[idx,3]-7500;
  probeINFO[idx,6]=probeINFO[idx,3]+7500;
  probeINFO[,8]=-log10(probeINFO[,8]);
  probeINFO[,3]=probeINFO[,3]/1e6;
  probeINFO[,5]=probeINFO[,5]/1e6;
  probeINFO[,6]=probeINFO[,6]/1e6;
  pXY=probeINFO[,8];
  yMAX = ceiling(max(c(pZY, pXY), na.rm=T)) + 1;
  if(is.null(data$Gene))
  {
    glist=cbind(probeINFO[,2],probeINFO[,5:6],as.character(probeINFO[,4]),probeINFO[,7]);
  } else {
    glist=data$Gene;
    glist[,2]=glist[,2]/1e6;
    glist[,3]=glist[,3]/1e6;
  }
  colnames(glist)=c("CHR", "GENESTART",  "GENEEND",   "GENE", "ORIENTATION");
  idx=which(is.na(glist[,2]) | is.na(glist[,3]));
  if(length(idx>0)) glist=glist[-idx,];
  generow = GeneRowNum(glist);
  num_row = max(as.numeric(generow$ROW));
  offset_map = ceiling(yMAX);
  offset_probe = yMAX / 2.5;
  num_probe = nprobePLOT
  offset_eqtl = ceiling(yMAX / 2.5) + 0.5;
  dev_axis = 0.1*yMAX;
  if(dev_axis<1.5) dev_axis = 1.5;
  yaxis.min = -offset_map - offset_eqtl*num_probe - dev_axis*(num_probe+1);
  yaxis.max = yMAX + ceiling(offset_probe) + 1;
  # scales of x-axis
  idx=match(data$GWAS[,1],data$SNP[,1]);
  gwasBP = as.numeric(data$SNP[idx,3])/1e6;
  #min.pos = min(gwasBP);
  #max.pos = max(gwasBP);
  min.pos = cis_start/1e6
  max.pos = cis_end/1e6
  start = min(as.numeric(glist[,2]));
  end = max(as.numeric(glist[,3]));
  bp = c(min.pos, max.pos, start, end);
  xmin = min(bp, na.rm=T) - 0.001;  xmax = max(bp, na.rm=T) +0.001;
  xmax=xmax+(xmax-xmin)*0.1 #extend
  ylab = expression(paste("-", log[10], "(", italic(P), " GWAS or SMR)", sep=""));
  xlab = paste("Chromosome", chrPLOT, "Mb");
  # plot GWAS p value
  par(mar=c(5,5,3,2), xpd=TRUE)
  plot(gwasBP, pZY, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max),
       ylab="", xlab=xlab, cex.lab=lab, cex.axis=axis,cex=0.6,
       xlim=c(xmin, xmax), pch=20, col="gray68");

  # y1 axis
  devbuf1 = yMAX/4
  axis(2, at=seq(0,yMAX,devbuf1), labels=round(seq(0,yMAX,devbuf1),0), las=1, cex.axis=axis);
  mtext(ylab, side=2, line=3, at=(yMAX*2/3), cex=cex_coeff);
  eqtl.lab = expression(paste("-", log[10], "(", italic(P), " eQTL)", sep=""));
  axis.start = 0; axis.down = offset_eqtl + dev_axis;
  for( k in 1 : nprobePLOT ) {
    axis.start = axis.start - axis.down
    eqtlinfobuf = data$eQTL[which(data$eQTL[,1]==probePLOT[k]),]
    if(dim(eqtlinfobuf)[1]==0) next;
    pvalbuf=-log10(pchisq((eqtlinfobuf[,3]/eqtlinfobuf[,4])^2,1,lower.tail=F));
    pvalbuf[which(is.infinite(pvalbuf))]=1e-300;
    if(length(which(smrprobes_red==probePLOT[k]))==0) {
      col_eqtl = "navy"
    } else col_eqtl = "maroon"
    eqtl.min = 0; eqtl.max = ceiling(max(pvalbuf))
    eqtl.max =ceiling(eqtl.max *1.25) #extend
    pvalbuf = pvalbuf/eqtl.max * offset_eqtl + axis.start
    idx=match(eqtlinfobuf[,2],data$SNP[,1]);
    eqtlbp = as.numeric(data$SNP[idx,3])/1e6;
    probegene = unique(as.character(data$SMR[which(data$SMR[,1]==probePLOT[k]),4]))
    par(new=TRUE)
    pchbuf = 4;
    #if(k%%2==0) pchbuf = 20;
    plot(eqtlbp, pvalbuf, yaxt="n", bty="n", ylim=c(yaxis.min,yaxis.max), xaxt="n",
         ylab="", xlab="", cex=0.8, pch=pchbuf, col=col_eqtl, xlim=c(xmin, xmax))
    # annotate the eQTLs
    text(xmin, axis.start+offset_eqtl-dev_axis/2 , label=substitute(paste(probeid, " (",italic(geneid), ")", sep=""),list(probeid=probePLOT[k], geneid=probegene)),col="black", cex=1, adj=0)
    # axis
    devbuf1 = offset_eqtl/3; devbuf2 = eqtl.max/3
    axis(2, at=seq(axis.start,(axis.start+offset_eqtl),devbuf1),
         labels=round(seq(0,eqtl.max,devbuf2),0),
         las=1, cex.axis=axis)
    # add separator line
    segments(xmin, axis.start+offset_eqtl+dev_axis/2, xmax, axis.start+offset_eqtl+dev_axis/2,
             col="dim grey", lty="24", lwd=1)
  }
  #ypos = (axis.start - dev_axis)/2
  ypos = (axis.start - dev_axis)*2/3
  mtext(eqtl.lab, side=2, at=ypos, line=3, cex=cex_coeff)

  # plot p value of bTG
  # all the probes
  num_gene = dim(generow)[1]
  dist = offset_map/num_row
  for( k in 1 : num_row ) {
    generowbuf = generow[which(as.numeric(generow[,5])==k),]
    xstart = as.numeric(generowbuf[,3])
    xend = as.numeric(generowbuf[,4])
    snbuf = which(xend-xstart< 1e-3)
    if(length(snbuf)>0) {
      xstart[snbuf] = xstart[snbuf] - 0.0025
      xend[snbuf] = xend[snbuf] + 0.0025
    }
    xcenter = (xstart+xend)/2
    xcenter = spread.labs(xcenter, mindiff=0.01, maxiter=1000, min = xmin, max = xmax)
    num_genebuf = dim(generowbuf)[1]
    for( l in 1 : num_genebuf ) {
      ofs=0.3
      if(l%%2==0) ofs=-0.8
      m = num_row - k
      ypos = m*dist + yaxis.min
      code = 1
      if(generowbuf[l,2]=="+") code = 2;
      arrows(xstart[l], ypos, xend[l], ypos, code=code, length=0.07, ylim=c(yaxis.min,yaxis.max),
             col=colors()[75], lwd=1)
      movebuf = as.numeric(generowbuf[l,6])*genemove
      text(xcenter[l]+movebuf, ypos,label=substitute(italic(genename), list(genename=as.character(generowbuf[l,1]))), pos=3, offset=ofs, col="black", cex=0.9)
    }
  }

  # plot the probes
  probeINFO=probeINFO[order(probeINFO[,8],decreasing = TRUE),];
  nprobeINFO=dim(probeINFO)[1];
  if(nprobeINFO>max_anno_probe){
    probeINFO=probeINFO[c(1:max_anno_probe),]
    nprobeINFO=dim(probeINFO)[1];
  }
  if(anno_selfdef) probeINFO=probeINFO[order(probeINFO[2],probeINFO[3]),] ####20170217
  xcenter = as.numeric(probeINFO[,3])
  xcbuf = xcenter
  ####20170217####
  if(anno_selfdef)
  {
    reginlength=(xmax-(xmax-xmin)*0.15)-xmin
    leftspot=xmin+reginlength/20
    rightspot=(xmax-(xmax-xmin)*0.15)-reginlength/20
    itvl=(rightspot-leftspot)/dim(probeINFO)[1]
    if(dim(probeINFO)[1]==1) {
      xcenter=as.numeric(probeINFO[,3])
    } else {
      xcenter=leftspot+itvl/2
      for( k in 2:dim(probeINFO)[1]) xcenter=c(xcenter,leftspot+k*itvl)
    }

  } else {
    xcenter = spread.labs(xcenter[1:nprobeINFO], mindiff=0.08, maxiter=1000, min = xmin, max = xmax-1)
  }
  # adjust the line position

  adjflag = rep(0, nprobeINFO)
  if(nprobeINFO>1) {
    dbuf = c(0, xcbuf[1:(nprobeINFO-1)])
    mflag = as.numeric(abs(xcbuf[1:(nprobeINFO)] - dbuf) < 0.01)
    adjflag = as.numeric( mflag | c(mflag[2:nprobeINFO],0) )
  }

  for( k in 1 : nprobeINFO)  {
    hitflag=FALSE
    if(length(which(heidiprobes_solid==probeINFO[k,1]))>0 & length(which(smrprobes_red==probeINFO[k,1]))>0) {
      hitflag=TRUE
      colplot = "maroon"; colfont=2; pchbuf=23;
    } else if(length(which(smrprobes_red==probeINFO[k,1]))>0) {
      colplot = "maroon"; colfont=2; pchbuf=5
      #} else if (length(which(heidiprobes_solid==probeINFO[k,1]))>0) {
      #hitflag=TRUE
      # colplot = "navy"; colfont=1; pchbuf=23
    } else {
      colplot = "navy"; colfont=1; pchbuf=5
    }
    if( as.numeric(probeINFO[k,8]) < 0 ) {
      colplot = "black"; colfont=1;
    }
    # plot p value of bxy
    plot_probe(probeINFO, k, colplot, xmin, xmax, yaxis.min, yaxis.max,pchbuf,hitflag)
    # annotate the probes
    if(k<=max_anno_probe)
    {
      ypos = 1.02*yMAX
      strbuf =
        text(xcenter[k], ypos,
             labels=substitute(paste(probeid, " (", italic(genename), ")", sep=""),
                               list(probeid=as.character(probeINFO[k,1]),
                                    genename=as.character(probeINFO[k,4]))),
             ylim=c(yaxis.min, yaxis.max),
             srt=30, col=colplot, font=colfont, cex=1, adj=0)
      # plot the lines
      # 1st step
      xstart = xcbuf[k]
      ystart = as.numeric(probeINFO[k,8]); yend = yMAX*(1-1/20);
      if( nprobeINFO > 1 ) {
        if(adjflag[k]==1) {
          xstart = (xcbuf[k] + xcenter[k])/2
          segments(xcbuf[k], ystart, xstart, ystart, col=colplot, lwd=axis, lty=3)
        }
      }
      segments(xstart, ystart, xstart, yend, col=colplot, lwd=axis, lty=3)
      # 2nd step
      xend = xcenter[k]; ystart = yMAX*(1-1/20); yend = yMAX*1.01;
      segments(xstart, ystart, xend, yend, col=colplot, lwd=axis, lty=3)
    }
  }
  # plot the threshold
  # SMR threshold
  ybuf = -log10(as.numeric(smr_thresh)); dev_anno = yMAX/9;
  strbuf = paste("pSMR = ",smr_thresh, sep="")
  segments(xmin, ybuf, xmax, ybuf, col="maroon", lty=2, lwd=1);
  text(xmax, ybuf+dev_anno, labels=strbuf, adj=1, col="maroon", cex=axis,font=3);

}

