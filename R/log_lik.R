#' @title Evaluate log lik of a tree
#' @export

log_lik=function(t_data,Rj,Tmin,sigma2,tau){
  t_R<-lapply(t_data, function(x) Rj[x])
  n_points<-lapply(t_R,length)
  n_points<-unlist(n_points)

  # each terminal node should contain at least 5 points
  if(min(n_points)<Tmin){
    # print("not valid tree")
    return(-Inf)
  }else{
    p_R<-lapply(t_R,llik_leave,sigma2=sigma2,tau=tau)
    p_R<-unlist(p_R)
    return(sum(p_R))
  }
}
