#' @title Generate sparse random orojection matrix
#' @param p,k matrix dimension p*k
#' @param nr number of matrices
#' @param lb prop of 0's in matrix
#' @export

Gen_SRP=function(p,k,nr,lb=1/3){
  mats=array(dim = c(p,k,nr))
  pp=p^2
  for (ir in 1:nr) {
    #generate one rotation matrix
    Q=sample(c(-1,0,1),p*k,replace = T,prob=c(lb/2,1-lb,lb/2))
    mats[,,ir]=matrix(Q,nrow=p,ncol=k)
  }
  return(mats)
}
