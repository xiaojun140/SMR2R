#' The function of making genes' qlt file
#'
#' @param qlt_path qlt data's path
#' @param genes Gene list file, usually a gene symbol ID file or gene Ensembl id file
#' @param query Gene-associated SNP threshold, usually 5e-8
#' @param out_file output file path
#'
#' @export get_gene_qlt
#' @details
#' The demo code detailed in the test file
#'
get_gene_qlt <- function(qlt_path,#qlt路径
                   genes,
                   query=5.0e-8,
                   out_file#输出路径
                   ){
  # .\smr-1.3.1-win.exe --beqtl-summary ../SMR/cis-eQTLs-full_eQTLGen --genes os.list --query 5.0e-8 --out os_eqtl --make-besd
   command <- paste0(smrfile(),
                    " --beqtl-summary ",
                    qlt_path,
                    " --genes ",
                    genes,
                    " --query ",
                    base::as.character(query),
                    " --out ",
                    paste0(getwd(),'/',out_file),
                    " --make-besd")

  system(command)
}
