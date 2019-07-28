#' @title Draw of y training sample only
#' @return Draws of yhat
#' @export

yhat.draw.train=function(btree_obj,Rj,tau,sigma2){
  t_data=btree_obj$t_data
  t_R<-lapply(t_data,function(x) Rj[x])

  mean.draw<-lapply(t_R,function(x){
  pmean=length(x)*mean(x)/sigma2/(length(x)/sigma2+1/tau^2)
  pvar=1/(length(x)/sigma2+1/tau^2)
  draw.mu=rnorm(1,pmean,sqrt(pvar))
  return(draw.mu)
  })

  yhat=c()
  for (dd in 1:length(mean.draw)) {
    yhat[t_data[[dd]]]=mean.draw[[dd]]
  }
  return(yhat)
}
