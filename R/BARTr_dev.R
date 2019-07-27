#' @title Fit a BART model, developer version
#' @param rule grp: Gaussian random projection; sgrp: sparse Gaussian random projection; bart: originla bart; hyperplane: connect two points
#' @export


BARTr_dev=function(X,y,x.test,sigdf=3, sigquant=.90,
               k=2.0, lambda=NA, sigest=NA,sigmaf=NA,
               power=2.0, base=.95,w=rep(1,length(y)),
               ntree=50,ndpost=700,nskip=300,Tmin=5,printevery=100,p_modify=c(2.5, 2.5, 4)/9,
               save_trees=F,rule='bart'){

  n=nrow(X)
  p=ncol(X)
  nt=nrow(x.test)

  #center
  fmean=mean(y)
  y.train = y-fmean
  #priors: nu,lambda,tau
  nu=sigdf
  if(is.na(lambda)) {
    if(is.na(sigest)) {
      if(p < n) {
        df = data.frame(X,y.train)
        lmf = lm(y.train~.,df)
        sigest = summary(lmf)$sigma
      } else {
        sigest = sd(y.train)
      }
    }
    qchi = qchisq(1.0-sigquant,nu)
    lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
  } else {
    sigest=sqrt(lambda)
  }

  if(is.na(sigmaf)) {
    tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
  } else {
    tau = sigmaf/sqrt(ntree)
  }

  #a list of ntree empty lists(trees).
  treelist=vector(ntree,mode='list')
  #give each tree a list speicifying the parameters
  treelist=lapply(treelist, function(x){
    x=list(s_pos=NULL,s_dir=NULL,s_rule=NULL,s_data=NULL,s_depth=NULL,s_obs=NULL,
           t_pos=1,t_data=list(1:n),t_depth=0)
  })

  #statistics to save
  tree_history=list()#a list of ndpost lists and each of ndpost lists is a list of ntree lists.

  tree_proposal_total=matrix(rep(0,ntree*length(p_modify)),nrow = ntree,ncol = length(p_modify))
  tree_proposal_accept=matrix(rep(0,ntree*length(p_modify)),nrow=ntree,ncol=length(p_modify))#record proposal acceptance number

  total_iter=nskip+ndpost

  sigma_draw=c(sigest)

  yhat.train=matrix(nrow=ndpost,ncol=n)
  yhat.test=matrix(nrow=ndpost,ncol=nt)
  #initilize single terminal node trees

  yhat.train.j=matrix(rnorm(ntree*n,0,sqrt(1/(n/sigest^2+1/tau^2))),nrow=ntree,ncol=n)
  yhat.test.j=matrix(rep(0,ntree*nt),nrow=ntree,ncol=nt)

  #####run bart for 100 iters then switch to 'rule'
  split_rule = rule
  #####

  for (i in 1:(total_iter)) {
    if(i%%printevery==0){print(sprintf("done %d (out of %d)",i,total_iter))};
    if(save_trees){tree_history[[i]]=treelist}
    #propose modification to each tree

    if(i<=100){rule = 'bart'}else{rule = split_rule}

    for (j in 1:ntree) {
      Rj=y.train-colSums(yhat.train.j[-j,,drop=F])

      # propose a modification
      move=which(rmultinom(1,1,p_modify)==1)

      # when we have no split node, only grow tree
      if(length(treelist[[j]]$s_pos)<1){move=1}
      tree_proposal_total[j,move]=tree_proposal_total[j,move]+1
      if(move==1){
        #grow

        #record time
        t1 = Sys.time()

        grown_tree=grow_tree(treelist[[j]],X,Rj,Tmin,rule)

        t2 = Sys.time()

        #print(paste('Grow tree takes', t2-t1,'second'))

        new_treej=grown_tree$btree_obj
        #calculate acceptance probablity
        t1 = Sys.time()
        lik_ratio = exp(log_lik(grown_tree$t_data_new,Rj,Tmin,sigma_draw[i]^2,tau)
                        - log_lik(grown_tree$t_data_old,Rj,Tmin,sigma_draw[i]^2,tau))
        t2 = Sys.time()

        #print(paste('Evaluate likelihood after grow takes', t2-t1,'second'))

        t1 = Sys.time()
        trans_ratio=p_modify[2]/p_modify[1]*length(treelist[[j]]$t_pos)/w2(new_treej)
        t2 = Sys.time()
        #print(paste('Evaluate trans ratio takes', t2-t1,'second'))
        #use new p_split
        #prior_ratio = (1-r^(-grown_tree$d-1))^2*r^(-grown_tree$d)/(1-r^(-grown_tree$d))
        prior_ratio=base*(1-base/(2+grown_tree$d)^power)^2/((1+grown_tree$d)^power-base)
        alpha=lik_ratio*trans_ratio*prior_ratio
        #print(sprintf('lik_ratio %.3f,trans_ratio %.3f,prior.ratio %.3f,alpha %.3f',lik_ratio,trans_ratio,prior_ratio,alpha))
        #print(sprintf('loglik_new %.3f, old %.3f',log_lik(grown_tree$t_data_new,X,Rj,Tmin,sigma_draw[i]^2,V),log_lik(grown_tree$t_data_old,X,Rj,Tmin,sigma_draw[i]^2,V)))
      }else if(move==2){
        #prune
        t1=Sys.time()
        pruned_tree=prune_tree(treelist[[j]])
        t2 = Sys.time()

        #print(paste('Prune tree takes', t2-t1,'second'))

        new_treej=pruned_tree$btree_obj

        t1 = Sys.time()
        lik_ratio = exp(log_lik(pruned_tree$t_data_new,Rj,Tmin,sigma_draw[i]^2,tau)
                        - log_lik(pruned_tree$t_data_old,Rj,Tmin,sigma_draw[i]^2,tau))
        t2=Sys.time()
        #print(paste('Evaluate likelihood after prune takes', t2-t1,'second'))

        t1 = Sys.time()
        trans_ratio=p_modify[1]/p_modify[2]*w2(new_treej)/(length(new_treej$t_pos))
        t2=Sys.time()
        #print(paste('Evaluate trans ratio takes', t2-t1,'second'))

        prior_ratio=((1+pruned_tree$d)^power-base)/(base*(1-base/(2+pruned_tree$d)^power)^2)
        #prior_ratio = (1-r^(-pruned_tree$d))/((1-r^(-pruned_tree$d-1))^2*r^(-pruned_tree$d))
        alpha=lik_ratio*trans_ratio*prior_ratio

      }else{
        # change(simple)
        t1 = Sys.time()
        changed_tree=change_tree(treelist[[j]],X,Rj,Tmin,rule)
        t2 = Sys.time()
        #print(paste('Change tree takes', t2-t1,'second'))

        new_treej = changed_tree$btree_obj

        t1 = Sys.time()
        lik_ratio = exp(log_lik(changed_tree$t_data_new,Rj,Tmin,sigma_draw[i]^2,tau)
                        - log_lik(changed_tree$t_data_old,Rj,Tmin,sigma_draw[i]^2,tau))
        t2 = Sys.time()
        #print(paste('Evaluate likelihood after change takes', t2-t1,'second'))


        alpha = lik_ratio
      }

      A=runif(1)

      if(is.nan(alpha)){
        alpha=0
      }
      #if a tree has a leaf node with no obs, discard it.
      #this can happen

      #if((0%in%as.numeric(unlist(lapply(new_treej$t_data, length))))){
      #  alpha=0
      #}

      #
      #

      t1 = Sys.time()

      if(A<alpha){
        # we accept the new tree
        #accept=accept+1
        tree_proposal_accept[j,move]=tree_proposal_accept[j,move]+1
        treelist[[j]]=new_treej
        #hat=yhat.draw(new_treej,x.test,Rj,tau,sigma_draw[i]^2)
        #hat=yhat.draw.linear(new_treej,X,x.test,Rj)
        if(i<=nskip){
          hat=yhat.draw(new_treej,x.test,Rj,tau,1,draw.test=F)
          yhat.train.j[j,] = hat$yhat
        }else{
          hat=yhat.draw(new_treej,x.test,Rj,tau,1)
          yhat.train.j[j,] = hat$yhat
          yhat.test.j[j,] = hat$ypred
        }


      }else{
        if(i<=nskip){
          hat=yhat.draw(treelist[[j]],x.test,Rj,tau,1,draw.test=F)
          yhat.train.j[j,] = hat$yhat
        }else{
          hat=yhat.draw(treelist[[j]],x.test,Rj,tau,1)
          yhat.train.j[j,] = hat$yhat
          yhat.test.j[j,] = hat$ypred
        }
      }

      t2 = Sys.time()
      #print(paste('Draw yhat and ypred takes', t2-t1,'second'))


    }

    if(i>nskip){
      yhat.train[i-nskip,]=colSums(yhat.train.j)
      yhat.test[i-nskip,]=colSums(yhat.test.j)
      #draw sigma
      res=y.train-yhat.train[i-nskip,]
    }else{
      res=y.train-colSums(yhat.train.j)
    }


    sigma_draw[i+1]=sqrt((nu*lambda + sum(res^2))/rchisq(1,n+nu))


  }

  yhat.train=yhat.train+fmean
  yhat.test=yhat.test+fmean
  sigma_draw=sigma_draw

  tree_leaf_count=as.numeric(unlist(lapply(treelist,function(x){length(x$t_data)})))

  return(list(yhat.train=yhat.train,yhat.test=yhat.test,
              yhat.train.mean=colSums(yhat.train)/nrow(yhat.train),
              yhat.test.mean=colSums(yhat.test)/nrow(yhat.test),sigma=sigma_draw,
              tree_history=tree_history,
              tree_proposal_total=tree_proposal_total,tree_proposal_accept=tree_proposal_accept,
              tree_leaf_count=tree_leaf_count))
}
