#' @title Rotation BADS for classification
#' @param rotate 'rr'='random rotation'','rraug'='random rotation+augmentation'
#' @param Others see ?pBADS
#' @return See ?pBADS
#' @importFrom BART rtnorm
#' @export


RotpBADS=function(X,y,x.test,cutoff=0.5,
                  k=2.0, binaryOffset=NULL,
                  ntree=50,ndpost=700,nskip=300,Tmin=2,printevery=100,
                  save_trees=F,rule='bart',rotate = 'rr'){

  n=nrow(X)
  p=ncol(X)
  nt=nrow(x.test)

  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  if(is.factor(y)){
    y = as.numeric.factor(y)
  }

  fmean=mean(y)

  #whcih y = 1; y=0
  y1.idx = which(y==1)
  y0.idx = which(y==0)


  tau = 3/(k*sqrt(ntree))
  if(length(binaryOffset)==0){binaryOffset=qnorm(fmean)}

  #a list of ntree empty lists(trees).
  treelist=vector(ntree,mode='list')
  #give each tree a list speicifying the parameters


  treelist=lapply(treelist, function(x){
    x=list(s_pos=NULL,s_dir=NULL,s_rule=NULL,s_data=NULL,s_depth=NULL,s_obs=NULL,
           t_pos=1,t_data=list(1:n),t_depth=0,RotMat=diag(p),t_test_data=NULL)
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
  y.train[y1.idx] = rtnorm(length(y1.idx),bm.f[y1.idx],1,-binaryOffset)
  y.train[y0.idx] = -rtnorm(length(y0.idx),-bm.f[y0.idx],1,binaryOffset)

  for(j in 1:ntree){
    Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
    if(rotate=='rr'){
      treelist[[j]]=grow_tree(treelist[[j]],X,Tmin)$btree_obj
    }
    if(rotate=='rraug'){
      treelist[[j]]=grow_tree(treelist[[j]],cbind(X,X%*%Rmat(p)),Tmin)$btree_obj
    }

    yhat.train.j[j,]=yhat.draw.train(treelist[[j]],Rj,tau,1)
  }


  for (i in 1:(total_iter)) {
    if(i%%printevery==0){print(sprintf("done %d (out of %d)",i,total_iter))};
    if(save_trees){tree_history[[i]]=treelist}

    bm.f = colSums(yhat.train.j)
    y.train = c()
    y.train[y1.idx] = rtnorm(length(y1.idx),bm.f[y1.idx],1,-binaryOffset)
    y.train[y0.idx] = -rtnorm(length(y0.idx),-bm.f[y0.idx],1,binaryOffset)

    #propose modification to each tree
    for (j in 1:ntree) {
      Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
      sig2 = 1

      RotMat_j = Rmat(p)

      if(rotate=='rr'){
        X_j = X%*%RotMat_j
      }
      if(rotate=='rraug'){
        X_j = cbind(X,X%*%RotMat_j)
      }

      changed_tree=change_tree(treelist[[j]],X_j,Tmin)
      new_treej = changed_tree$btree_obj
      alpha = exp(log_lik(changed_tree$t_data_new,Rj,Tmin,sig2,tau)
                      - log_lik(changed_tree$t_data_old,Rj,Tmin,sig2,tau))

      A=runif(1)

      if(is.nan(alpha)){
        alpha=0
      }

      if(A<alpha){
        new_treej$RotMat = RotMat_j

        if(i<=nskip){
          hat=yhat.draw.train(new_treej,Rj,tau,sig2)
          yhat.train.j[j,] = hat
        }else{
          hat=yhat.draw.bads(new_treej,x.test,Rj,tau,sig2,rotate)
          yhat.train.j[j,] = hat$yhat
          yhat.test.j[j,] = hat$ypred
          new_treej$t_test_data = hat$t_idx
        }

        treelist[[j]]=new_treej

      }else{
        if(i<=nskip){
          hat=yhat.draw.train(treelist[[j]],Rj,tau,sig2)
          yhat.train.j[j,] = hat
        }else{
          hat=yhat.draw2.bads(treelist[[j]],x.test,Rj,tau,sig2,rotate)
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
