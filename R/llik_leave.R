#' @title Evalute log lik of a leave
#' @importFrom mvtnorm dmvnorm
#' @export

llik_leave<-function(y,sigma2,tau){
  ny = length(y)
  llik = dmvnorm(y,rep(0,ny),sigma2*diag(ny)+tau^2,log=T)
  llik
}
