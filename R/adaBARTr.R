#'@title Fit an adaptive bayesian additive regression tree(adaBART) model for classification
#'@param all see ?pBARTr
#'@export


adaBARTr=function(X,y,x.test,cutoff=0.5,
                k=2.0, binaryOffset=NULL,
                power=2.0, base=.95,p_split='CGM',r=2,
                ntree=50,ndpost=700,nskip=300,Tmin=2,printevery=100,p_modify=c(0.5, 0.5, 0),
                save_trees=F,rule='bart',pre_train=F,n_pre_train=100){



  n=nrow(X)
  p=ncol(X)
  nt=nrow(x.test)


  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  if(is.factor(y)){
    y = as.numeric.factor(y)
  }



  #whcih y = 1; y=0
  y1.idx = which(y==1)
  y0.idx = which(y==0)

  #priors: nu,lambda,tau
  tau = 3/(k*sqrt(ntree))
  if(length(binaryOffset)==0){binaryOffset=qnorm(mean(y))}

  #a list of ntree empty lists(trees).
  treelist=vector(ntree,mode='list')
  #give each tree a list specifying the parameters
  treelist=lapply(treelist, function(x){
    x=list(s_pos=NULL,s_dir=NULL,s_rule=NULL,s_data=NULL,s_depth=NULL,s_obs=NULL,
           t_pos=1,t_data=list(1:n),t_depth=0,t_test_data=NULL)
  })

  #statistics to save
  tree_history=list()#a list of ndpost lists and each of ndpost lists is a list of ntree lists.

  n_move = length(p_modify)

  tree_proposal_total=matrix(rep(0,ntree*n_move),nrow = ntree,ncol = n_move)
  tree_proposal_accept=matrix(rep(0,ntree*n_move),nrow=ntree,ncol=n_move)

  total_iter=nskip+ndpost

  yhat.train=matrix(nrow=ndpost,ncol=n)
  yhat.test=matrix(nrow=ndpost,ncol=nt)
  #initilize single terminal node trees

  yhat.train.j=matrix(rnorm(ntree*n,0,sqrt(1/(n+1/tau^2))),nrow=ntree,ncol=n)
  yhat.test.j=matrix(rep(0,ntree*nt),nrow=ntree,ncol=nt)

  #####run bart for 100 iters then switch to 'rule'
  split_rule = rule
  #####

  # initialize weights for each data points, and err, and alpha value, as in standard asaboost algorithm.
  wgt = rep(1,n)
  err.i = rep(0,ndpost)
  alpha.i = rep(0,ndpost)

  for (i in 1:(total_iter)) {
    if(i%%printevery==0){print(sprintf("done %d (out of %d)",i,total_iter))};
    if(save_trees){tree_history[[i]]=treelist}

    if(pre_train){
      if(i<=n_pre_train){rule = 'bart'}else{rule = split_rule}
    }

    bm.f = colSums(yhat.train.j)
    y.train = c()

    y.train[y1.idx] = rtnorm(length(y1.idx),bm.f[y1.idx],sqrt(wgt[y1.idx]),-binaryOffset)
    y.train[y0.idx] = -rtnorm(length(y0.idx),-bm.f[y0.idx],sqrt(wgt[y0.idx]),binaryOffset)

    #propose modification to each tree
    for (j in 1:ntree) {
      Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
      sig2 = 1

      BART_draw = BARTr_train(X,Rj,treelist[[j]],p_modify,Tmin,
                              rule,sig2,tau,base,power,p_split,r)
      alpha = BART_draw$alpha
      new_treej = BART_draw$new_treej
      move = BART_draw$move

      tree_proposal_total[j,move]=tree_proposal_total[j,move]+1

      A=runif(1)

      if(is.nan(alpha)){
        alpha=0
      }

      if(A<alpha){
        # we accept the new tree
        tree_proposal_accept[j,move]=tree_proposal_accept[j,move]+1

        if(i<=nskip){
          hat=yhat.draw.train(new_treej,Rj,tau,sig2)
          yhat.train.j[j,] = hat
        }else{
          hat=yhat.draw(new_treej,x.test,Rj,tau,sig2)
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
          hat=yhat.draw2(treelist[[j]],x.test,Rj,tau,sig2)
          yhat.train.j[j,] = hat$yhat
          yhat.test.j[j,] = hat$ypred
        }
      }


    }

    if(i>nskip){
      yhat.train[i-nskip,]=1*(pnorm(colSums(yhat.train.j)+binaryOffset)>=cutoff)
      yhat.test[i-nskip,]=1*(pnorm(colSums(yhat.test.j)+binaryOffset)>=cutoff)
      #update weights
      err.i[i-nskip] = sum(wgt*(yhat.train[i-nskip,]!=y))/sum(wgt)
      alpha.i[i-nskip] = log((1-err.i[i-nskip])/err.i[i-nskip])
      wgt = wgt * exp(alpha.i[i-nskip]*(yhat.train[i-nskip,]!=y))
      wgt = wgt/sum(wgt)*n
    }

  }

  tree_leaf_count=as.numeric(unlist(lapply(treelist,function(x){length(x$t_data)})))

  yhat.train.weighted = diag(alpha.i)%*%(2*yhat.train-1)
  yhat.test.weighted = diag(alpha.i)%*%(2*yhat.test-1)

  yhat.train.mean=sign(colSums(yhat.train.weighted))
  yhat.test.mean=sign(colSums(yhat.test.weighted))

  yhat = ((yhat.train.mean)+1)/2
  ypred = ((yhat.test.mean)+1)/2



  return(list(yhat=yhat,ypred=ypred,
              yhat.train=yhat.train,yhat.test=yhat.test,
              sample_weight = wgt,
              classifier_err = err.i,
              tree_history=tree_history,
              tree_proposal_total=tree_proposal_total,
              tree_proposal_accept=tree_proposal_accept,
              tree_leaf_count=tree_leaf_count))
}
