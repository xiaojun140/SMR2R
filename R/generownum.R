#' Title
#'
#' @param GENELIST Here is a function for plotting within the package, no need to run it.
#'
#' @export GeneRowNum
#'

GeneRowNum = function(GENELIST) {
  BP_THRESH = 0.03; MAX_ROW = 5
  # get the start and end position
  GENELIST = GENELIST[!duplicated(GENELIST$GENE),]
  START1 = as.numeric(GENELIST$GENESTART); END1 = as.numeric(GENELIST$GENEEND)
  STRLENGTH = nchar(as.character(GENELIST$GENE))
  MIDPOINT = (START1 + END1)/2
  START2 = MIDPOINT-STRLENGTH/250; END2 = MIDPOINT+STRLENGTH/250
  START = cbind(START1, START2); END = cbind(END1, END2);
  START = apply(START, 1, min); END = apply(END, 1, max)
  GENELIST = data.frame(GENELIST, START, END)
  GENELIST = GENELIST[order(as.numeric(GENELIST$END)),]
  START = as.numeric(GENELIST$START); END = as.numeric(GENELIST$END)
  # get the row index for each gene
  NBUF = dim(GENELIST)[1]
  ROWINDX = rep(1, NBUF)
  ROWEND = as.numeric(rep(0, MAX_ROW))
  MOVEFLAG = as.numeric(rep(0, NBUF))
  if(NBUF>1) {
    for( k in 2 : NBUF ) {
      ITERFLAG=FALSE
      if(START[k] < END[k-1]) {
        INDXBUF=ROWINDX[k-1]+1
      } else INDXBUF = 1
      if(INDXBUF>MAX_ROW) INDXBUF=1;
      REPTIME=0
      repeat{
        if( ROWEND[INDXBUF] > START[k] ) {
          ITERFLAG=FALSE
          INDXBUF=INDXBUF+1
          if(INDXBUF>MAX_ROW) INDXBUF = 1
        } else {
          ITERFLAG=TRUE
        }
        if(ITERFLAG) break;
        REPTIME = REPTIME+1
        if(REPTIME==MAX_ROW) break;
      }
      ROWINDX[k]=INDXBUF;

      if( (abs(ROWEND[ROWINDX[k]]-START[k]) < BP_THRESH)
          | ((ROWEND[ROWINDX[k]]-START[k])>0) ) {
        MOVEFLAG[k] = 1
        SNBUF = tail(which(ROWINDX[c(1:k)]==ROWINDX[k]), n=2)[1]
        MOVEFLAG[SNBUF] = MOVEFLAG[SNBUF] - 1
      }
      if(ROWEND[ROWINDX[k]]<END[k]) {
        ROWEND[ROWINDX[k]] = END[k]  }
    }
  }
  GENEROW = data.frame(as.character(GENELIST$GENE),
                       as.character(GENELIST$ORIENTATION),
                       as.numeric(GENELIST$GENESTART),
                       as.numeric(GENELIST$GENEEND),
                       ROWINDX, MOVEFLAG)
  colnames(GENEROW) = c("GENE", "ORIENTATION", "START", "END", "ROW", "MOVEFLAG")
  return(GENEROW)
}
