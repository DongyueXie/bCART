#' @title BART model train workhourse
#' @param rule grp: Gaussian random projection; sgrp: sparse Gaussian random projection; bart: originla bart; hyperplane: connect two points
#' @export
#'

BARTr_train = function(X,Rj,treej,p_modify,Tmin,rule,sigma2,tau,base,power){


  move=which(rmultinom(1,1,p_modify)==1)

  # when we have no split node, only grow tree
  if(length(treej$s_pos)<1){move=1}

  if(move==1){
    #grow
    grown_tree=grow_tree(treej,X,Tmin,rule)
    new_treej=grown_tree$btree_obj
    #calculate acceptance probablity
    lik_ratio = exp(log_lik(grown_tree$t_data_new,Rj,Tmin,sigma2,tau)
                    - log_lik(grown_tree$t_data_old,Rj,Tmin,sigma2,tau))
    trans_ratio=p_modify[2]/p_modify[1]*length(treej$t_pos)/w2(new_treej)
    #use new p_split
    #prior_ratio = (1-r^(-grown_tree$d-1))^2*r^(-grown_tree$d)/(1-r^(-grown_tree$d))
    prior_ratio=base*(1-base/(2+grown_tree$d)^power)^2/((1+grown_tree$d)^power-base)
    alpha=lik_ratio*trans_ratio*prior_ratio
    #print(sprintf('lik_ratio %.3f,trans_ratio %.3f,prior.ratio %.3f,alpha %.3f',lik_ratio,trans_ratio,prior_ratio,alpha))
    #print(sprintf('loglik_new %.3f, old %.3f',log_lik(grown_tree$t_data_new,X,Rj,Tmin,sigma_draw[i]^2,V),log_lik(grown_tree$t_data_old,X,Rj,Tmin,sigma_draw[i]^2,V)))
  }else if(move==2){
    #prune
    pruned_tree=prune_tree(treej)
    new_treej=pruned_tree$btree_obj

    lik_ratio = exp(log_lik(pruned_tree$t_data_new,Rj,Tmin,sigma2,tau)
                    - log_lik(pruned_tree$t_data_old,Rj,Tmin,sigma2,tau))

    trans_ratio=p_modify[1]/p_modify[2]*w2(new_treej)/(length(new_treej$t_pos))
    prior_ratio=((1+pruned_tree$d)^power-base)/(base*(1-base/(2+pruned_tree$d)^power)^2)
    #prior_ratio = (1-r^(-pruned_tree$d))/((1-r^(-pruned_tree$d-1))^2*r^(-pruned_tree$d))
    alpha=lik_ratio*trans_ratio*prior_ratio

  }else{
    # change(simple)
    changed_tree=change_tree(treej,X,Tmin,rule)
    new_treej = changed_tree$btree_obj
    lik_ratio = exp(log_lik(changed_tree$t_data_new,Rj,Tmin,sigma2,tau)
                    - log_lik(changed_tree$t_data_old,Rj,Tmin,sigma2,tau))
    alpha = lik_ratio
  }

  return(list(alpha=alpha,move=move,new_treej=new_treej))

}



