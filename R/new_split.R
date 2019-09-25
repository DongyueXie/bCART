#' @title Genrate new splitted data index
#' @return left and right data index
#' @export

new_split=function(X_sub,sub_data,dir_new,rule_new){
  proj=X_sub[,which(dir_new==1),drop=F]

  #print(proj)
  #print(rule_new)

  if(is.character(rule_new)){
    left_idx = which(proj==rule_new)
    right_idx = which(proj!=rule_new)
  }else{
    left_idx = which(proj<=rule_new)
    right_idx = which(proj>rule_new)
  }

  # if(is.matrix(X_sub)){
  #   left_idx<-which(proj<=rule_new)
  #   right_idx<-which(proj>rule_new)
  # }else{
  #   left_idx<-which(X_sub<=rule_new)
  #   right_idx<-which(X_sub>rule_new)
  # }
  #data in the left child(index)
  left_data = sub_data[left_idx]
  left_data = na.omit(left_data)

  #data in the right child(index)
  right_data = sub_data[right_idx]
  right_data = na.omit(right_data)

  return(list(left_data=left_data,right_data=right_data))
}
