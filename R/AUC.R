#' @title Calculate AUC of binary classification
#' @param label true labels
#' @param pred predicted lables
#' @return AUC
#' @importFrom pROC roc
#' @importFrom pROC auc
#' @export

AUC=function(label,pred){
  roc_obj=roc(label,pred,quiet=TRUE)
  return(auc(roc_obj))
}
