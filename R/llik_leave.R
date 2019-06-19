#' @title Evalute log lik of a leave
#' @export

llik_leave<-function(y,sigma2,tau){
  ny<-length(y)
  # llik= -ny/2*log(2*pi*sigma2)+
  #   0.5*(log(sigma2)-log(sigma2+ny*tau^2))-
  #   1/(2*sigma2)*(sum((y-mean(y))^2)-mean(y)^2*ny^2/(ny+sigma2/tau^2)+ny*mean(y)^2)
  llik = dmvnorm(y,rep(0,ny),sigma2*diag(ny)+tau^2,log=T)
  return(llik)
}
