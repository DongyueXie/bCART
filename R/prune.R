#' @title Prune tree proposal in bCART
#' @export


prune_tree<-function(btree_obj){

  # we find the proper death node from the terminal node
  # (we only want the split node with two terminal node kids)
  parent_pos<-floor(btree_obj$t_pos/2)
  parent_freq<-data.frame(table(parent_pos))
  parent_2_pos<-parent_freq[which(parent_freq$Freq>1),1]
  # factor to numeric
  parent_2_pos<-as.numeric(as.character(parent_2_pos))

  #select the death split node
  idx_death<-sample(1:length(parent_2_pos),1)
  pos_death=parent_2_pos[idx_death]

  # the kids of the death node
  death_left<-2*pos_death
  death_left<-which(btree_obj$t_pos==death_left)
  death_right<-2*pos_death+1
  death_right=which(btree_obj$t_pos==death_right)

  # remove the death node from split record
  s_death=which(btree_obj$s_pos==pos_death)
  btree_obj$s_pos=btree_obj$s_pos[-s_death]
  btree_obj$s_obs=btree_obj$s_obs[-s_death]
  btree_obj$s_dir=btree_obj$s_dir[,-s_death,drop=F]
  btree_obj$s_rule=btree_obj$s_rule[-s_death]
  old_s_depth=btree_obj$s_depth[s_death]
  btree_obj$s_depth=btree_obj$s_depth[-s_death]
  sub_data=btree_obj$s_data[[s_death]]
  btree_obj$s_data=btree_obj$s_data[-s_death]

  # remove the death terminal node
  btree_obj$t_pos=btree_obj$t_pos[-c(death_left,death_right)]
  t_data_old = btree_obj$t_data[c(death_left,death_right)]
  btree_obj$t_data=btree_obj$t_data[-c(death_left,death_right)]
  btree_obj$t_depth=btree_obj$t_depth[-c(death_left,death_right)]

  k=length(btree_obj$t_pos)
  # push the new terminal node into the series
  btree_obj$t_pos=c(btree_obj$t_pos,pos_death)
  btree_obj$t_depth=c(btree_obj$t_depth,old_s_depth)
  btree_obj$t_data[[k+1]]=sub_data
  t_data_new = list(sub_data)
  return(list(btree_obj=btree_obj,d=old_s_depth,t_data_old = t_data_old, t_data_new = t_data_new))
}
