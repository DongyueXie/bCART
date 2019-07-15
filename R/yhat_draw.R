#' @title Draw of y
#' @return Draws of yhat and ypred
#' @export

yhat.draw=function(btree_obj,x.test,Rj,tau,sigma2,draw.test=T){
  t_data=btree_obj$t_data

  t_R<-lapply(t_data,function(x) Rj[x])
  mean.draw<-lapply(t_R,function(x){

    pmean=length(x)*mean(x)/sigma2/(length(x)/sigma2+1/tau^2)
    pvar=1/(length(x)/sigma2+1/tau^2)
    draw.mu=rnorm(1,pmean,sqrt(pvar))

    return(draw.mu)
  })

  #x.test.list = lapply(seq_len(nrow(x.test)),function(z){x.test[z,]})
  #t_idx=mclapply(x.test.list,function(x){find_terminal_idx(x,btree_obj)},mc.cores=6)
  #t_idx = unlist(t_idx)

  if(draw.test){
    t_idx = apply(x.test,1,function(x){find_terminal_idx(x,btree_obj)})
    yhat=c()
    ypred=c()
    for (dd in 1:length(mean.draw)) {
      yhat[t_data[[dd]]]=mean.draw[[dd]]
      ypred[which(t_idx==dd)]=mean.draw[[dd]]
    }
    return(list(yhat=yhat,ypred=ypred))
  }else{
    yhat=c()
    for (dd in 1:length(mean.draw)) {
      yhat[t_data[[dd]]]=mean.draw[[dd]]
    }
    return(list(yhat=yhat))
  }


}
