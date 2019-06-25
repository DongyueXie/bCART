#' @title Generate Random Rotation matrix
#' @export

Rmat = function(p){
  Q=gramschmidt(matrix(rnorm(p^2),nrow=p))$Q
  if(det(Q)<0){idx=sample(1:p,2);Q[c(idx[1],idx[2]),]=Q[c(idx[2],idx[1]),]}
  return(Q)
}
