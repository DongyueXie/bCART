#' @title AUC
#' @export

AUC=function(label,pred){
  roc_obj=pROC::roc(label,pred)
  return(auc(roc_obj))
}
