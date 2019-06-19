#' @title number of internal nodes that have two children terminal nodes.
#' @export

w2=function(btree_obj){
  x1=btree_obj$t_pos/2
  in.node=floor(x1)
  sum(duplicated(in.node))
}
