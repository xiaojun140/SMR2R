#' Title Standardizing outcome format
#'
#' @param file The path or dataframe of outcome.
#' @param sep The delimiter of the input data.
#' @param SNP The column name of SNP.
#' @param A1 The column name of alt.
#' @param A2 The column name of ref.
#' @param freq The column name of af_alt.
#' @param b The column name of beta.
#' @param se The column name of sebeta.
#' @param p The column name of pval.
#' @param n The column name of ncase+ncontrol.
#' @param output The name of output file.
#'
#' @return A formatted file of outcome in current path.
#' @export format_SMR
#'
format_SMR <- function(file,sep="\t",SNP,A1,A2,freq,b,se,p,n=NULL,output="mygwas"){
  if (is.data.frame(file)) {
    df <- file
  } else if (is.character(file)) {
    if (file.exists(file)) {
      cat("Reading the outcome file \n")
      df <- data.table::fread(file,sep = sep)
    } else {
      stop("File not found.\n")
    }
  } else {
    stop("Invalid input. Please provide either a data frame or a file path.\n")
  }

  if (is.null(file) || is.null(SNP) || is.null(A1) || is.null(A2) || is.null(freq) || is.null(b) || is.null(se) || is.null(p)) {
    stop("Error: You must provide necessary parameter.\n")
  }

  df <- as.data.frame(df)
  df <- dplyr::select(df,SNP=!!rlang::enquo(SNP),A1=!!rlang::enquo(A1),A2=!!rlang::enquo(A2),
                  freq=!!rlang::enquo(freq),b=!!rlang::enquo(b),se=!!rlang::enquo(se),p=!!rlang::enquo(p))
  df <- dplyr ::filter(df,SNP!="" ,SNP!=".",!is.na(SNP))
  df <- tidyr::as_tibble(df) |> tidyr::separate_rows(SNP, sep = ",")
  df=df[!duplicated(df$SNP),]
  df=df[df$SNP!="",]

if(is.null(n)){df$n="NA"}else{df <-dplyr::select(df,n=!!rlang::enquo(n))}

  data.table::fwrite(df,file = stringr::str_c(output,".ma"),sep = "\t",row.names = F,col.names = T,quote = F)
  return("DONE")
}
