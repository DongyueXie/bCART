#' @title Grow tree proposal
#' @export


grow_tree=function(btree_obj,X,Rj,Tmin,rule){

  # find the split node, var and rule
  # We only want to find large terminal node(more than Tmin samples)
  t_length= lapply(btree_obj$t_data,length)
  t_length=unlist(t_length)
  t_large=which(t_length>=(2*Tmin))
  if(length(t_large)>0){
    s_new<-sample(t_large,1)
  }else{
    s_new<-sample(1:length(t_length),1)
  }
  pos_new=btree_obj$t_pos[s_new]
  depth_new=btree_obj$t_depth[s_new]
  # split the data
  sub_data<-btree_obj$t_data[[s_new]]
  X_sub=X[sub_data,]
  Rj_sub=Rj[sub_data]
  after_select<-select_rule(X_sub,Rj_sub,Tmin,rule)
  dir_new<-after_select$dir_new
  #proj=after_select$proj
  rule_new<-after_select$rule_new
  after_split<-new_split(X_sub,sub_data,dir_new,rule_new)
  # record old tree terminal data index
  t_data_old = btree_obj$t_data[s_new]
  # remove this new split node from the terminal node record
  btree_obj$t_pos=btree_obj$t_pos[-s_new]
  #btree_obj$s_obs=btree_obj$s_obs[-s_new]
  old_t_depth=btree_obj$t_depth[s_new]
  btree_obj$t_depth=btree_obj$t_depth[-s_new]
  btree_obj$t_data=btree_obj$t_data[-s_new]

  # update its information in split node record
  btree_obj$s_pos<-c(btree_obj$s_pos,pos_new)
  btree_obj$s_obs=c(btree_obj$s_obs,length(sub_data))
  btree_obj$s_dir<-cbind(btree_obj$s_dir,dir_new)
  btree_obj$s_rule<-c(btree_obj$s_rule,rule_new)
  btree_obj$s_depth<-c(btree_obj$s_depth,depth_new)
  btree_obj$s_data[[length(btree_obj$s_pos)]]<-sub_data

  # push the new terminal nodes into record
  k=length(btree_obj$t_pos)
  btree_obj$t_pos=c(btree_obj$t_pos,2*pos_new,2*pos_new+1)
  btree_obj$t_depth=c(btree_obj$t_depth,old_t_depth+1,old_t_depth+1)
  btree_obj$t_data[[(k+1)]]=after_split$left_data
  btree_obj$t_data[[(k+2)]]=after_split$right_data
  t_data_new = list(after_split$left_data,after_split$right_data)
  return(list(btree_obj=btree_obj,d=old_t_depth,t_data_old = t_data_old,t_data_new = t_data_new))
}
