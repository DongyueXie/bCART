#' @title select splitting variable and value
#' @export

select_rule=function(X,Tmin,rule){
  p_var=ncol(X)

  if(missing(rule)){
    rule='bart'
  }

  if(rule=='hyperplane'){
    v.idx=sample(1:nrow(X),2)
    v=X[v.idx[1],]-X[v.idx[2],]
    proj=X%*%v
    eps = runif(1)
    rule_new = eps*proj[v.idx[1]]+(1-eps)*proj[v.idx[2]]
  }else{
    if(rule=='bart'){
      v=rmultinom(1,1,rep(1/p_var,p_var))
    }
    if(rule=='grp'){
      v = rnorm(p_var)
      v = v/sqrt(sum(v^2))
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
    }

    proj=X%*%v
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

  return(list(dir_new=v,rule_new=rule_new,proj=proj))
}
