#' @title Check if a node is terminal
#' @export

is.terminal<-function(kid_pos,t_pos){
  return((kid_pos %in% t_pos))
}
