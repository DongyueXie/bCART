#' @title Rotation BART for classification
#' @param rotate 'rr'='random rotation'','rraug'='random rotation+augmentation','srp'='sparse random projection'
#' @param srpk the number of cols of sparse projection matrix
#' @export


RotpBART=function(X,y,x.test,cutoff=0.5,
                  k=2.0, binaryOffset=NULL,
                  power=2.0, base=.95,w=rep(1,length(y)),
                  ntree=50,ndpost=700,nskip=300,Tmin=5,printevery=100,p_modify=c(2.5, 2.5, 4)/9,
                  save_trees=F,rotate = 'rr',srpk=2*ncol(X)){

  n=nrow(X)
  p=ncol(X)
  nt=nrow(x.test)

  #center
  fmean=mean(y)

  tau = 3/(k*sqrt(ntree))
  if(length(binaryOffset)==0){binaryOffset=qnorm(fmean)}

  #a list of ntree empty lists(trees).
  treelist=vector(ntree,mode='list')
  #give each tree a list speicifying the parameters
  treelist=lapply(treelist, function(x){
    x=list(s_pos=NULL,s_dir=NULL,s_rule=NULL,s_data=NULL,s_depth=NULL,s_obs=NULL,
           t_pos=1,t_data=list(1:n),t_depth=0)
  })

  #a list of ndpost lists and each of ndpost lists is a list of ntree lists.
  tree_history=list()

  #record proposal acceptance number
  tree_proposal_total=matrix(rep(0,ntree*length(p_modify)),nrow = ntree,ncol = length(p_modify))
  tree_proposal_accept=matrix(rep(0,ntree*length(p_modify)),nrow=ntree,ncol=length(p_modify))

  total_iter=nskip+ndpost

  yhat.train=matrix(nrow=ndpost,ncol=n)
  yhat.test=matrix(nrow=ndpost,ncol=nt)
  #initilize single terminal node trees

  yhat.train.j=matrix(rnorm(ntree*n,0,sqrt(1/(n+1/tau^2))),nrow=ntree,ncol=n)
  yhat.test.j=matrix(rep(0,ntree*nt),nrow=ntree,ncol=nt)

  #rotate each data set
  if(rotate=='rr'){
    RotationMat = Gen_R(p,ntree)
    x.train.r=array(NA,dim=c(n,p,ntree))
    x.train.r[]=apply(RotationMat, 3, function(r){X%*%r})
    x.test.r=array(NA,dim=c(nt,p,ntree))
    x.test.r[]=apply(RotationMat,3,function(r){x.test%*%r})
  }
  if(rotate=='rraug'){
    RotationMat = Gen_R(p,ntree)
    x.train.r=array(NA,dim=c(n,2*p,ntree))
    x.test.r=array(NA,dim=c(nt,2*p,ntree))
    x.train.r[]=apply(RotationMat, 3, function(r){cbind(X,X%*%r)})
    x.test.r[]=apply(RotationMat,3,function(r){cbind(x.test,x.test%*%r)})
  }
  if(rotate=='srp'){
    RotationMat = Gen_SRP(p,srpk,ntree)
    x.train.r=array(NA,dim=c(n,srpk,ntree))
    x.test.r=array(NA,dim=c(nt,srpk,ntree))
    x.train.r[]=apply(RotationMat, 3, function(r){X%*%r})
    x.test.r[]=apply(RotationMat,3,function(r){x.test%*%r})
  }

  X = x.train.r
  x.test = x.test.r

  for (i in 1:(total_iter)) {
    if(i%%printevery==0){print(sprintf("done %d (out of %d)",i,total_iter))};
    if(save_trees){tree_history[[i]]=treelist}

    bm.f = colSums(yhat.train.j)
    y.train = c()
    for(ni in 1:n){
      if(y[ni]==1){
        y.train[ni] = rtnorm(bm.f[ni],-binaryOffset,1)
      }else{
        y.train[ni] = -rtnorm(-bm.f[ni],binaryOffset,1)
      }
    }

    #propose modification to each tree
    for (j in 1:ntree) {
      Rj=y.train-colSums(yhat.train.j[-j,,drop=F])

      # propose a modification
      move=which(rmultinom(1,1,p_modify)==1)

      # when we have no split node, only grow tree
      if(length(treelist[[j]]$s_pos)<1){move=1}
      tree_proposal_total[j,move]=tree_proposal_total[j,move]+1
      if(move==1){
        #grow
        grown_tree=grow_tree(treelist[[j]],X[,,j],Rj,Tmin)
        new_treej=grown_tree$btree_obj
        #calculate acceptance probablity
        lik_ratio = exp(log_lik(grown_tree$t_data_new,Rj,Tmin,1,tau)
                        - log_lik(grown_tree$t_data_old,Rj,Tmin,1,tau))
        trans_ratio=p_modify[2]/p_modify[1]*length(treelist[[j]]$t_pos)/w2(new_treej)
        #use new p_split
        #prior_ratio = (1-r^(-grown_tree$d-1))^2*r^(-grown_tree$d)/(1-r^(-grown_tree$d))
        prior_ratio=base*(1-base/(2+grown_tree$d)^power)^2/((1+grown_tree$d)^power-base)
        alpha=lik_ratio*trans_ratio*prior_ratio
        #print(sprintf('lik_ratio %.3f,trans_ratio %.3f,prior.ratio %.3f,alpha %.3f',lik_ratio,trans_ratio,prior_ratio,alpha))
        #print(sprintf('loglik_new %.3f, old %.3f',log_lik(grown_tree$t_data_new,X,Rj,Tmin,1,V),log_lik(grown_tree$t_data_old,X,Rj,Tmin,1,V)))
      }else if(move==2){
        #prune
        pruned_tree=prune_tree(treelist[[j]])
        new_treej=pruned_tree$btree_obj

        lik_ratio = exp(log_lik(pruned_tree$t_data_new,Rj,Tmin,1,tau)
                        - log_lik(pruned_tree$t_data_old,Rj,Tmin,1,tau))

        trans_ratio=p_modify[1]/p_modify[2]*w2(new_treej)/(length(new_treej$t_pos))
        prior_ratio=((1+pruned_tree$d)^power-base)/(base*(1-base/(2+pruned_tree$d)^power)^2)
        #prior_ratio = (1-r^(-pruned_tree$d))/((1-r^(-pruned_tree$d-1))^2*r^(-pruned_tree$d))
        alpha=lik_ratio*trans_ratio*prior_ratio

      }else{
        # change(simple)
        changed_tree=change_tree(treelist[[j]],X[,,j],Rj,Tmin)
        new_treej = changed_tree$btree_obj
        lik_ratio = exp(log_lik(changed_tree$t_data_new,Rj,Tmin,1,tau)
                        - log_lik(changed_tree$t_data_old,Rj,Tmin,1,tau))
        alpha = lik_ratio
      }

      A=runif(1)

      if(is.nan(alpha)){
        alpha=0
      }
      #if a tree has a leaf node with no obs, discard it.
      #this can happen
      if((0%in%as.numeric(unlist(lapply(new_treej$t_data, length))))){
        alpha=0
      }
      #
      #
      if(A<alpha){
        # we accept the new tree
        #accept=accept+1
        tree_proposal_accept[j,move]=tree_proposal_accept[j,move]+1
        treelist[[j]]=new_treej
        hat=yhat.draw(new_treej,x.test[,,j],Rj,tau,1)
        yhat.train.j[j,] = hat$yhat
        yhat.test.j[j,] = hat$ypred

      }else{
        hat=yhat.draw(treelist[[j]],x.test,Rj,tau,1)
        yhat.train.j[j,] = hat$yhat
        yhat.test.j[j,] = hat$ypred
      }


    }

    if(i>nskip){
      yhat.train[i-nskip,]=colSums(yhat.train.j)
      yhat.test[i-nskip,]=colSums(yhat.test.j)
      #draw sigma
      res=y.train-yhat.train[i-nskip,]
    }else{
      res=y.train-colSums(yhat.train.j)
    }

  }

  tree_leaf_count=as.numeric(unlist(lapply(treelist,function(x){length(x$t_data)})))

  phat.train = pnorm(yhat.train+binaryOffset)
  phat.test = pnorm(yhat.test+binaryOffset)

  phat.train.mean=colSums(phat.train)/nrow(phat.train)
  phat.test.mean=colSums(phat.test)/nrow(phat.test)

  yhat = 1*((phat.train.mean)>=cutoff)
  ypred = 1*((phat.test.mean)>=cutoff)

  return(list(yhat=yhat,ypred=ypred,
              yhat.train=yhat.train,yhat.test=yhat.test,
              phat.train.mean=phat.train.mean,
              phat.test.mean=phat.test.mean,
              tree_history=tree_history,
              tree_proposal_total=tree_proposal_total,tree_proposal_accept=tree_proposal_accept,
              tree_leaf_count=tree_leaf_count))
}
