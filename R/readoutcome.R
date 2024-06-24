#' Read the plot file
#'
#' @param plotfile The text file for plot
#'
#' @return SMRData
#' @export ReadSMRData
ReadSMRData = function(plotfile){
  SMRData = list();
  key=c("$probe","$SNP","$GWAS","$eQTL");
  skiplines=0;
  keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
  if(keywords[1]!=key[1])
  {
    print("ERROR: plot file is not correct!");
    quit();
  }
  nprobes=as.numeric(keywords[2]);
  SMRData$probeID=keywords[3];


  skiplines=skiplines+1;
  SMRData$SMR=read.table(plotfile, header=F, nrows=nprobes, skip=skiplines);
  skiplines=skiplines+nprobes;
  keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
  if(keywords[1]!=key[2])
  {
    print("ERROR: plot file is not correct!");
    quit();
  }
  nrs=as.numeric(keywords[2]);
  skiplines=skiplines+1;
  SMRData$SNP=read.table(plotfile, header=F, nrows=nrs, skip=skiplines);
  skiplines=skiplines+nrs;
  keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
  if(keywords[1]!=key[3])
  {
    print("ERROR: plot file is not correct!");
    quit();
  }
  ngwas=as.numeric(keywords[2]);
  skiplines=skiplines+1;
  SMRData$GWAS=read.table(plotfile, header=F, nrows=ngwas, skip=skiplines);
  skiplines=skiplines+ngwas;
  keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
  if(keywords[1]!=key[4])
  {
    print("ERROR: plot file is not correct!");
    quit();
  }
  neqtl=as.numeric(keywords[2]);
  skiplines=skiplines+1;

  keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
  prbname=keywords[1];
  neqtlsnp=as.numeric(keywords[2]);
  skiplines=skiplines+1;
  SMRData$eQTL=read.table(plotfile, header=F, nrows=neqtlsnp, skip=skiplines);
  SMRData$eQTL=cbind(prbname,SMRData$eQTL)
  skiplines=skiplines+neqtlsnp;
  if(neqtl>1)
  {
    for(i in 2:neqtl)
    {
      keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
      prbname=keywords[1];
      neqtlsnp=as.numeric(keywords[2]);
      skiplines=skiplines+1;
      raweQTLtmp=read.table(plotfile, header=F, nrows=neqtlsnp, skip=skiplines);
      raweQTLtmp=cbind(prbname,raweQTLtmp);
      SMRData$eQTL=rbind(SMRData$eQTL,raweQTLtmp);
      skiplines=skiplines+neqtlsnp;
    }
  }

  keywords=scan(plotfile, what="", nlines=1, skip=skiplines);
  if(length(keywords)>0)
  {
    if(keywords[1]!="$Gene")
    {
      print("ERROR: plot file is not correct!");
      quit();
    }
    ngenes=as.numeric(keywords[2]);
    skiplines=skiplines+1;
    SMRData$Gene=read.table(plotfile, header=F, nrows=ngenes, skip=skiplines);
  }
  return(SMRData)
}
