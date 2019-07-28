#' @title Draw of y when reject tree movement
#' @return Draws of yhat and ypred
#' @export

yhat.draw2=function(btree_obj,x.test,Rj,tau,sigma2){
  t_data=btree_obj$t_data

  t_R<-lapply(t_data,function(x) Rj[x])
  mean.draw<-lapply(t_R,function(x){

    pmean=length(x)*mean(x)/sigma2/(length(x)/sigma2+1/tau^2)
    pvar=1/(length(x)/sigma2+1/tau^2)
    draw.mu=rnorm(1,pmean,sqrt(pvar))

    return(draw.mu)
  })


  if(is.null(btree_obj$t_test_data)){
    t_idx = apply(x.test,1,function(x){find_terminal_idx(x,btree_obj)})
  }else{
    t_idx = btree_obj$t_test_data
  }



  yhat=c()
  ypred=c()
  for (dd in 1:length(mean.draw)) {
    yhat[t_data[[dd]]]=mean.draw[[dd]]
    ypred[which(t_idx==dd)]=mean.draw[[dd]]
  }
  return(list(yhat=yhat,ypred=ypred))



}
