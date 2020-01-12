#' @title Rotation BADS for classification
#' @param rotate 'rr'='random rotation'','rraug'='random rotation+augmentation'
#' @param srpk the number of cols of sparse projection matrix
#' @param Others see ?BADS
#' @return See ?BADS
#' @export


RotBADS=function(X,y,x.test,sigdf=3, sigquant=.90,k=2,
                 lambda=NA, sigest=NA,sigmaf=NA,
                 ntree=50,ndpost=200,nskip=100,Tmin=2,printevery=100,
                 save_trees=F,rotate='rr',rule='bart'){

  n=nrow(X)
  p=ncol(X)
  nt=nrow(x.test)

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
           t_pos=1,t_data=list(1:n),t_depth=0,RotMat=diag(p),t_test_data=NULL)
  })

  #a list of ndpost lists and each of ndpost lists is a list of ntree lists.
  tree_history=list()

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
    treelist[[j]]=grow_tree(treelist[[j]],X,Tmin,rule)$btree_obj
    yhat.train.j[j,] = yhat.draw.train(treelist[[j]],Rj,tau,sigest^2)
  }


  for (i in 1:(total_iter)) {
    if(i%%printevery==0){print(sprintf("done %d (out of %d)",i,total_iter))};
    if(save_trees){tree_history[[i]]=treelist}

    #propose modification to each tree
    for (j in 1:ntree) {
      Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
      sig2 = sigma_draw[i]^2

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

    sigma_draw[i+1]=sqrt((nu*lambda + sum(res^2))/rchisq(1,n+nu))
  }

  tree_leaf_count=as.numeric(unlist(lapply(treelist,function(x){length(x$t_data)})))

  yhat.train=yhat.train+fmean
  yhat.test=yhat.test+fmean
  sigma_draw=sigma_draw

  return(list(yhat.train=yhat.train,yhat.test=yhat.test,
              yhat.train.mean=colSums(yhat.train)/nrow(yhat.train),
              yhat.test.mean=colSums(yhat.test)/nrow(yhat.test),sigma=sigma_draw,
              tree_history=tree_history,
              tree_leaf_count=tree_leaf_count))
}
