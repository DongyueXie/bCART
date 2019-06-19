#' @title Accuracy
#' @param label True labels
#' @export

accuracy=function(label,pred){
  sum(diag(table(label,pred)))/length(label)
}
