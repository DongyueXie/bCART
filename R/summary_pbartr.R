#' @title Summary pBART fit
#' @export

summary_pbartr = function(lbart.obj,y,y.test){
  fit = lbart.obj
  message(sprintf('Training accuracy %.3f',accuracy(fit$yhat,y)))
  message(sprintf('Training AUC %.3f',AUC(fit$yhat,y)))
  message(sprintf('Testing accuracy %.3f',accuracy(fit$ypred,y.test)))
  message(sprintf('Testing AUC %.3f',AUC(fit$ypred,y.test)))
  ar = apply(fit$tree_proposal_accept/fit$tree_proposal_total, 2, mean)
  message(sprintf('Acceptance ratio: grow %.3f, prune %.3f, change %.3f',ar[1],ar[2],ar[3]))
  message(sprintf('Tree leaf count:'))
  print(table(fit$tree_leaf_count))
}
