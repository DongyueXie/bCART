#' @title Fit a bayesian additive decision stump model for regression
#' @description Bayesian additive decision stump(BADS) is a Bayesian sum of two-leaf-node trees model.
#' @param X samples by features matrix
#' @param y response
#' @param x.test test samples by feature matrix
#' @param sigdf,sigquant,k,lambda,sigest,sigmaf see ?BART::wbart
#' @param ntree number of decison stumps
#' @param nskip,ndpost number of burn-in and posterior draws
#' @param Tmin minimum number of samples in a leaf node allowed
#' @param printevery print progress for every 'printevery' iterations
#' @param save_trees whether save all the trees from each iteration as a list
#' @param rule The splitting rule of a node. Choices are: 1. "grp": Gaussian random projection, randomly draw a length p vector from standard normal as the linear combination coefficients of p variables; 2. sgrp: sparse Gaussian random projection, which generates sparse linear combination coefficients; 3. bart: originla bart splits, which are axis-aligned splits; 4. hyperplane: randomly connect two points from the node as the partiton of node space.
#' @param pre_train whether pre-train the model using 'bart' rule before switching to another splitting rule.
#' @param n_pre_train number of iterations of pre-train
#' @return BADS returns a list of the following elements.
#' \item{yhat.train}{A matrix with ndpost rows and nrow(X) columns.}
#' \item{yhat.test}{A matrix with ndpost rows and nrow(x.test) columns.}
#' \item{yhat.train.mean}{Posterior mean of MCMC draws of traning data fits}
#' \item{yhat.test.mean}{Posterior mean of MCMC draws of testing data fits}
#' \item{sigma}{draws of random error vairaince, length = nskip+ndpost}
#' \item{tree_history}{If save_trees = TRUE, then a list of all trees}
#' @author Dongyue Xie: \email{dongyxie@gmail.com}
#' @references Chipman, H., George, E., and McCulloch R. (2010) Bayesian Additive Regression Trees. The Annals of Applied Statistics, 4,1, 266-298 <doi:10.1214/09-AOAS285>.
#' @export


BADS=function(X,y,x.test,sigdf=3, sigquant=.90,k=2,
               lambda=NA, sigest=NA,sigmaf=NA,
               ntree=50,ndpost=200,nskip=100,Tmin=2,printevery=100,
               save_trees=F,rule='bart',pre_train=T,n_pre_train=100){

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
           t_pos=1,t_data=list(1:n),t_depth=0,t_test_data=NULL)
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
    treelist[[j]]=grow_tree(treelist[[j]],X,Tmin,rule)$btree_obj
    yhat.train.j[j,] = yhat.draw.train(treelist[[j]],Rj,tau,sigest^2)
  }

  #####run bart for 100 iters then switch to 'rule'
  split_rule = rule
  #####

  for (i in 1:(total_iter)) {
    if(i%%printevery==0){print(sprintf("done %d (out of %d)",i,total_iter))};
    if(save_trees){tree_history[[i]]=treelist}
    #propose modification to each tree

    if(pre_train){
      if(i<=n_pre_train){rule = 'bart'}else{rule = split_rule}
    }


      for (j in 1:ntree) {
        Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
        sig2 = sigma_draw[i]^2
        changed_tree=change_tree(treelist[[j]],X,Tmin,rule)
        new_treej = changed_tree$btree_obj
        alpha = exp(log_lik(changed_tree$t_data_new,Rj,Tmin,sig2,tau)
                        - log_lik(changed_tree$t_data_old,Rj,Tmin,sig2,tau))

        A=runif(1)

        if(is.nan(alpha)){
          alpha=0
        }

        if(A<alpha){
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
        yhat.train[i-nskip,]=colSums(yhat.train.j)
        yhat.test[i-nskip,]=colSums(yhat.test.j)
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
              tree_leaf_count=tree_leaf_count))
}
