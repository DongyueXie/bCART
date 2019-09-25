#' @title select splitting variable and value
#' @export

select_rule=function(X,Tmin,rule){
  p_var=ncol(X)

  if(missing(rule)){
    rule='bart'
  }

  if(rule == 'bart'){
    unique_check = apply(X,2,unique)
    if(is.list(unique_check)){
      avail_var = which(unlist(lapply(unique_check,length))!=1)
      p_var_avail = length(avail_var)
      v=rmultinom(1,1,rep(1/p_var_avail,p_var_avail))
      split_var = avail_var[which(v==1)]
      proj = X[,split_var]
      v = rep(0,p_var)
      v[split_var] = 1
    }else{
      v=rmultinom(1,1,rep(1/p_var,p_var))
      proj = X[,which(v==1)]
    }

    if(is.factor(proj)){
      avail_split_value = levels(proj)[which(table(proj)>=Tmin)]
      if(length(avail_split_value)==0){
        rule_new = as.character(sample(proj,1))
      }else if(length(avail_split_value)==1){
        rule_new = avail_split_value
      }else{
        rule_new = sample(avail_split_value,1)
      }
    }else{
      table_proj = table(proj)
      n_split_value = length(table_proj)
      split_value = as.numeric(names(table_proj))
      cusum1 = cumsum(table_proj)
      cusum2 = cumsum(rev(table_proj))
      avail_split_value = intersect(split_value[cusum1>=Tmin],rev(split_value)[-1][cusum2[-1]>=Tmin])
      if(length(avail_split_value)==0){
        rule_new = sample(split_value,1)
      }else if(length(avail_split_value)==1){
        rule_new = avail_split_value
      }else{
        rule_new = sample(avail_split_value,1)
      }
    }
  }else if(rule=='hyperplane'){
    v.idx=sample(1:nrow(X),2)
    v=X[v.idx[1],]-X[v.idx[2],]
    proj=X%*%v
    eps = runif(1)
    rule_new = eps*proj[v.idx[1]]+(1-eps)*proj[v.idx[2]]
  }else{
    if(rule=='grp'){
      v = rnorm(p_var)
      v = v/sqrt(sum(v^2))
      proj=X%*%v
    }
    if(rule=='sgrp'){
      a = 1
      b = p_var
      gama=0
      while(sum(gama)==0){
        pp=rbeta(1,a,b)
        gama = rbinom(p_var,1,pp)
      }
      v = rnorm(p_var)*gama
      v = v/sqrt(sum(v^2))
      proj=X%*%v
    }
    ##  version not considering factors for oblique split
    n_rule=length(proj)
    if(n_rule>=Tmin*2){
      rule_order=sort.int(proj,index.return=TRUE)$ix
      if(n_rule==2){
        rule_order=rule_order[1]
      }else{
        rule_order=rule_order[-c(1:(Tmin-1),(n_rule-Tmin+1):n_rule)]
      }
      if(length(rule_order)==1){
        rule_new=rule_order
      }else{
        rule_new=sample(rule_order,1)
      }
    }else{
      rule_new=sample(1:n_rule,1)
    }
    rule_new=proj[rule_new]
  }

  # v: one-hot vector indicate which variable to split
  # rule_new: split value
  return(list(dir_new=v,rule_new=rule_new,proj=proj))
}
