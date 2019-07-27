#' @title Fit a bayesian additive decision stump model
#' @param rule grp: Gaussian random projection; sgrp: sparse Gaussian random projection; bart: originla bart; hyperplane: connect two points
#' @export


BADS=function(X,y,x.test,sigdf=3, sigquant=.90,k=2,
               lambda=NA, sigest=NA,sigmaf=NA,
               ntree=50,ndpost=200,nskip=100,Tmin=2,printevery=100,
               save_trees=F,rule='bart',pre_train=T){

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

  total_iter=nskip+ndpost

  sigma_draw=c(sigest)

  yhat.train=matrix(nrow=ndpost,ncol=n)
  yhat.test=matrix(nrow=ndpost,ncol=nt)
  #initilize single terminal node trees

  yhat.train.j=matrix(rnorm(ntree*n,0,sqrt(1/(n/sigest^2+1/tau^2))),nrow=ntree,ncol=n)
  yhat.test.j=matrix(rep(0,ntree*nt),nrow=ntree,ncol=nt)

  # initialize decision stumps
  for(j in 1:ntree){
    Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
    grown_tree=grow_tree(treelist[[j]],X,Rj,Tmin,rule)
    new_treej=grown_tree$btree_obj
    treelist[[j]]=new_treej
    hat=yhat.draw(new_treej,x.test,Rj,tau,sigest^2,draw.test=F)
    yhat.train.j[j,] = hat$yhat
  }

  #####run bart for 100 iters then switch to 'rule'
  split_rule = rule
  #####

  for (i in 1:(total_iter)) {
    if(i%%printevery==0){print(sprintf("done %d (out of %d)",i,total_iter))};
    if(save_trees){tree_history[[i]]=treelist}
    #propose modification to each tree

    if(pre_train){
      if(i<=100){rule = 'bart'}else{rule = split_rule}
    }


    if(i<=nskip){

      #######
      for (j in 1:ntree) {
        Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
        changed_tree=change_tree(treelist[[j]],X,Rj,Tmin,rule)
        new_treej = changed_tree$btree_obj
        lik_ratio = exp(log_lik(changed_tree$t_data_new,Rj,Tmin,sigma_draw[i]^2,tau)
                        - log_lik(changed_tree$t_data_old,Rj,Tmin,sigma_draw[i]^2,tau))
        alpha = lik_ratio

        A=runif(1)

        if(is.nan(alpha)){
          alpha=0
        }

        if(A<alpha){
          # we accept the new tree
          treelist[[j]]=new_treej
          #hat=yhat.draw(new_treej,x.test,Rj,tau,sigma_draw[i]^2)
          #hat=yhat.draw.linear(new_treej,X,x.test,Rj)
          hat=yhat.draw(new_treej,x.test,Rj,tau,sigma_draw[i]^2,draw.test=F)
          yhat.train.j[j,] = hat$yhat


        }else{
          hat=yhat.draw(treelist[[j]],x.test,Rj,tau,sigma_draw[i]^2,draw.test=F)
          yhat.train.j[j,] = hat$yhat
        }
      }
      res=y.train-colSums(yhat.train.j)

      ######
    }else{

      for (j in 1:ntree) {
        Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
        changed_tree=change_tree(treelist[[j]],X,Rj,Tmin,rule)
        new_treej = changed_tree$btree_obj
        lik_ratio = exp(log_lik(changed_tree$t_data_new,Rj,Tmin,sigma_draw[i]^2,tau)
                        - log_lik(changed_tree$t_data_old,Rj,Tmin,sigma_draw[i]^2,tau))
        alpha = lik_ratio

        A=runif(1)

        if(is.nan(alpha)){
          alpha=0
        }

        if(A<alpha){
          # we accept the new tree
          treelist[[j]]=new_treej
          #hat=yhat.draw(new_treej,x.test,Rj,tau,sigma_draw[i]^2)
          #hat=yhat.draw.linear(new_treej,X,x.test,Rj)
          hat=yhat.draw(new_treej,x.test,Rj,tau,sigma_draw[i]^2)
          yhat.train.j[j,] = hat$yhat
          yhat.test.j[j,] = hat$ypred

        }else{
          hat=yhat.draw(treelist[[j]],x.test,Rj,tau,sigma_draw[i]^2)
          yhat.train.j[j,] = hat$yhat
          yhat.test.j[j,] = hat$ypred
        }


      }


        yhat.train[i-nskip,]=colSums(yhat.train.j)
        yhat.test[i-nskip,]=colSums(yhat.test.j)
        #draw sigma
        res=y.train-yhat.train[i-nskip,]


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
              tree_leaf_count=tree_leaf_count))
}
