#' @title Update a tree after tree moves
#' @export
#'

update_tree<-function(kid_pos,kid_data,btree_obj,X){
  #why add is.terminal? I don't think it can be terminla node since we selected from internal nodes.
  #i see; this function itself is a recursive function
  if(is.terminal(kid_pos,btree_obj$t_pos)){
    kid_idx<-which(btree_obj$t_pos==kid_pos)
    btree_obj$t_data[[kid_idx]]<-kid_data
    return(btree_obj)
  }else{
    kid_idx<-which(btree_obj$s_pos==kid_pos)
    btree_obj$s_data[[kid_idx]]<-kid_data
    btree_obj$s_obs[[kid_idx]]=length(kid_data)
    left_kid<-kid_pos*2
    right_kid<-kid_pos*2+1
    X_sub<-X[kid_data,]
    after_split<-new_split(X_sub,kid_data,btree_obj$s_dir[,kid_idx],btree_obj$s_rule[kid_idx])
    left_data<-after_split$left_data
    right_data<-after_split$right_data

    # update the left branch
    left_update=update_tree(left_kid,left_data,btree_obj,X)

    # then we update the right branch
    right_update=update_tree(right_kid,right_data,left_update,X)

    return(right_update)
  }

}
