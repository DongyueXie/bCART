#' @title Summary of BART fit
#' @export


summary.bartr = function(lbart.obj,y,Ey,Ey.test){
  fit = lbart.obj
  message(sprintf('Training error %.3f',rmse(fit$yhat.train.mean,y)))
  message(sprintf('In-sample estimaiton error %.3f',rmse(fit$yhat.train.mean,Ey)))
  message(sprintf('Out-sample estimaiton error %.3f',rmse(fit$yhat.test.mean,Ey.test)))
  ar = apply(fit$tree_proposal_accept/fit$tree_proposal_total, 2, mean)
  message(sprintf('Acceptance ratio: grow %.3f, prune %.3f, change %.3f',ar[1],ar[2],ar[3]))
  message(sprintf('Tree leaf count:'))
  print(table(fit$tree_leaf_count))
  plot(fit$sigma,ylab = 'sigma',main='sigma draws')
}
