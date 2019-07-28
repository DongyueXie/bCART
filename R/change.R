#' @title Change proposal
#' @export

change_tree=function(btree_obj,X,Tmin,rule){

  parent_pos<-floor(btree_obj$t_pos/2)
  parent_freq<-data.frame(table(parent_pos))
  parent_2_pos<-parent_freq[which(parent_freq$Freq>1),1]
  # factor to numeric
  parent_2_pos<-as.numeric(as.character(parent_2_pos))

  #select the death split node
  idx_change<-sample(1:length(parent_2_pos),1)
  pos_change=parent_2_pos[idx_change]

  s_change=which(btree_obj$s_pos==pos_change)


  # repartition the data on this split node
  sub_data<-btree_obj$s_data[[s_change]]
  X_sub=X[sub_data,]
  s_change_dir = btree_obj$s_dir[,s_change]
  after_select<-select_rule(X_sub,Tmin,rule)
  dir_new<-after_select$dir_new
  #proj=after_select$proj
  rule_new<-after_select$rule_new

  # update the dir and rule
  btree_obj$s_dir[,s_change]=dir_new
  btree_obj$s_rule[s_change]=rule_new

  # update the kids
  after_update<-update_tree(pos_change,sub_data,btree_obj,X)

  #check
  left_kid = pos_change*2
  right_kid = pos_change*2+1

  t_data_old = btree_obj$t_data[match(c(left_kid,right_kid),btree_obj$t_pos)]
  t_data_new = after_update$t_data[match(c(left_kid,right_kid),after_update$t_pos)]
  return(list(btree_obj=after_update,t_data_old=t_data_old,t_data_new=t_data_new))
}
