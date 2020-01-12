#' @title Calculate AUC of binary classification
#' @param label true labels
#' @param pred predicted lables
#' @return AUC
#' @importFrom pROC roc
#' @export

AUC=function(label,pred){
  roc_obj=roc(label,pred)
  return(auc(roc_obj))
}
