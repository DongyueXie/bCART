#' @title Calculate Accuracy of binary classification
#' @param label True labels
#' @param pred predicted labels
#' @export

accuracy=function(label,pred){
  sum(diag(table(label,pred)))/length(label)
}
