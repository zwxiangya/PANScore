#'#' This is some description of this function.
#'@title  Calculate the PANoptosis related risk score
#' @param expr log2 transformed TPM matrix derived from  normalized gene expression matrix from bulk samples; rows are genes, cols are samples.
#'
#' @author Wei Zhang
#'
#' @return A dataframe containing two columns: continuous numerical variables with Riskscore, the ID of the samples.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("exprset.hovon")
#' RCDscore <- PANScore.calculate(exprset  = exprset.hovon)
#'}


PANScore.calculate = function(exprset){

  exprset.t = exprset %>% t() %>% as.data.frame()
  data("df")
  df = df[,2:3]
  colnames(df) = c("Gene","Coeff")
  rownames(df) =df$Gene
  comsa = intersect(colnames(exprset.t), df$Gene)
  df = df[comsa, ]

  risk.score <- apply(exprset.t[,comsa],1,function(x) {x %*% as.numeric(df$Coeff)})
  RS.table = data.frame(ID= rownames(exprset.t), Riskscore =risk.score )
  return(RS.table)
}
