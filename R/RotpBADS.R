#' @title Rotation BADS for classification
#' @param rotate 'rr'='random rotation'','rraug'='random rotation+augmentation','srp'='sparse random projection'
#' @param srpk the number of cols of sparse projection matrix
#' @export


RotpBADS=function(X,y,x.test,cutoff=0.5,
                  k=2.0, binaryOffset=NULL,
                  ntree=50,ndpost=700,nskip=300,Tmin=2,printevery=100,
                  save_trees=F,rotate = 'rr',srpk=2*ncol(X),pre_train=T){

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
           t_pos=1,t_data=list(1:n),t_depth=0,RotMat=diag(p))
  })

  #a list of ndpost lists and each of ndpost lists is a list of ntree lists.
  tree_history=list()

  total_iter=nskip+ndpost

  yhat.train=matrix(nrow=ndpost,ncol=n)
  yhat.test=matrix(nrow=ndpost,ncol=nt)
  #initilize single terminal node trees

  yhat.train.j=matrix(rnorm(ntree*n,0,sqrt(1/(n+1/tau^2))),nrow=ntree,ncol=n)
  yhat.test.j=matrix(rep(0,ntree*nt),nrow=ntree,ncol=nt)

  # initialize decision stumps

  bm.f = colSums(yhat.train.j)
  y.train = c()
  for(ni in 1:n){
    if(y[ni]==1){
      y.train[ni] = rtnorm(bm.f[ni],-binaryOffset,1)
    }else{
      y.train[ni] = -rtnorm(-bm.f[ni],binaryOffset,1)
    }
  }

  for(j in 1:ntree){
    Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
    grown_tree=grow_tree(treelist[[j]],X,Rj,Tmin)
    new_treej=grown_tree$btree_obj
    treelist[[j]]=new_treej
    hat=yhat.draw(new_treej,x.test,Rj,tau,1,draw.test=F)
    yhat.train.j[j,] = hat$yhat
  }


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

      RotMat_j = Rmat(p)
      X_j = X%*%RotMat_j
      changed_tree=change_tree(treelist[[j]],X_j,Rj,Tmin)
      new_treej = changed_tree$btree_obj
      lik_ratio = exp(log_lik(changed_tree$t_data_new,Rj,Tmin,1,tau)
                        - log_lik(changed_tree$t_data_old,Rj,Tmin,1,tau))
      alpha = lik_ratio


      A=runif(1)

      if(is.nan(alpha)){
        alpha=0
      }

      if(A<alpha){
        new_treej$RotMat = RotMat_j
        treelist[[j]]=new_treej

        #hat=yhat.draw(new_treej,x.test,Rj,tau,sigma_draw[i]^2)
        #hat=yhat.draw.linear(new_treej,X,x.test,Rj)
        if(i<=nskip){
          hat=yhat.draw(new_treej,x.test%*%RotMat_j,Rj,tau,1,draw.test=F)
          yhat.train.j[j,] = hat$yhat
        }else{
          hat=yhat.draw(new_treej,x.test%*%RotMat_j,Rj,tau,1)
          yhat.train.j[j,] = hat$yhat
          yhat.test.j[j,] = hat$ypred
        }


      }else{
        if(i<=nskip){
          hat=yhat.draw(treelist[[j]],x.test%*%treelist[[j]]$RotMat,Rj,tau,1,draw.test=F)
          yhat.train.j[j,] = hat$yhat
        }else{
          hat=yhat.draw(treelist[[j]],x.test%*%treelist[[j]]$RotMat,Rj,tau,1)
          yhat.train.j[j,] = hat$yhat
          yhat.test.j[j,] = hat$ypred
        }
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
              tree_leaf_count=tree_leaf_count))
}