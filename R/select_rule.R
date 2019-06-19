#' @title Select splitting vairable and value
#' @export

select_rule=function(X,y,Tmin){
  p_var=ncol(X)
  v=rmultinom(1,1,rep(1/p_var,p_var))
  proj=X%*%v
  n_rule=length(proj)
  if(n_rule>=Tmin*2){
    rule_order=sort.int(proj,index.return=TRUE)$ix
    rule_order=rule_order[-c(1:Tmin,(n_rule-Tmin+2):n_rule)]
    if(length(rule_order)==1){
      rule_new=rule_order
    }else{
      rule_new=sample(rule_order,1)
    }
  }else{
    rule_new=sample(1:n_rule,1)
  }
  rule_new=proj[rule_new]
  return(list(dir_new=v,rule_new=rule_new,proj=proj))
}
