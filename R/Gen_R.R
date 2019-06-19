#' @title Generate an array of random rotation matrices
#' @param p p*p matrix
#' @param nr number of matrices
#' @export
#'
Gen_R=function(p,nr){
  mats=array(dim = c(p,p,nr))
  for (ir in 1:nr) {
    mats[,,ir]=Rmat(p)
  }
  return(mats)
}
