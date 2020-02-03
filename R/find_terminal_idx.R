#' @title Find which leaf node a test sample belongs to
#' @export
#'

find_terminal_idx<-function(X_test,btree_obj){

  #start with the top node
  flag_pos=1

  while(!is.terminal(flag_pos,btree_obj$t_pos)){
    split_idx=which(btree_obj$s_pos == flag_pos)
    if(is.null(dim(btree_obj$s_dir))){
      # check if rule is bart or oblque. BART split supports factor x.
      if(sum(btree_obj$s_dir)==1){
        split_proj=X_test[which(btree_obj$s_dir==1)]
      }else{
        split_proj = sum(X_test*(btree_obj$s_dir))
      }
    }else{
      s_dir = btree_obj$s_dir[,split_idx]
      if(sum(s_dir)==1){
        split_proj = X_test[which( s_dir== 1)]
      }else{
        split_proj = sum(X_test*s_dir)
      }

    }
    split_rule=btree_obj$s_rule[split_idx]

    if(is.character(split_rule)){
      if(split_proj == split_rule){
        flag_pos=flag_pos*2
      }else{
        flag_pos=flag_pos*2+1
      }
    }else{
      if(split_proj<=split_rule){
        flag_pos=flag_pos*2
      }else{
        flag_pos=flag_pos*2+1
      }
    }
  }

  t_idx=which(btree_obj$t_pos==flag_pos)

  return(t_idx)

}
