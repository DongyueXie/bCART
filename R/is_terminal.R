#' @title Check if a node is a leaf node
#' @export

is.terminal<-function(kid_pos,t_pos){
  return((kid_pos %in% t_pos))
}
