#' @title Fit a bayesian additive regression tree(BART) model for classification
#' @param X training samples by features matrix
#' @param y response, factor with level 0,1
#' @param x.test testing samples by feature matrix
#' @param cutoff label = 1 if p>cutoff; else label = 0.
#' @param k,power,base see ?BART::pbart
#' @param binaryOffset The model is P(Y=1 | x) = F(f(x) + binaryOffset).
#' @param p_split choice of 'CGM', 'RS'. 'CGM' splits an internal node with probability base*(1+d)^(-power); 'RS': r^(-d), 2 <= r <= n. d is the depth of an internal node.
#' @param r 'RS' splits an internal node with probability r^(-d), 2 <= r <= n. d is the depth of an internal node.
#' @param ntree number of trees
#' @param nskip,ndpost number of burn-in and posterior draws
#' @param Tmin minimum number of samples in a leaf node allowed
#' @param printevery print progress for every 'printevery' iterations
#' @param p_modify proportion of three moves: grow, prune, change.
#' @param save_trees whether save all the trees from each iteration as a list
#' @param rule The splitting rule of an internal node. Choices are: 1. "grp": Gaussian random projection, randomly draw a length p vector from standard normal as the linear combination coefficients of p variables; 2. sgrp: sparse Gaussian random projection, which generates sparse linear combination coefficients; 3. bart: originla bart splits, which are axis-aligned splits; 4. hyperplane: randomly connect two points from the node as the partiton of node space.
#' @param pre_train whether pre-train the BART model using 'bart' rule before switching to another splitting rule.
#' @param n_pre_train number of iterations of pre-train
#' @return BARTr returns a list of the following elements.
#' \item{yhat.train}{A matrix with ndpost rows and nrow(X) columns.}
#' \item{yhat.test}{A matrix with ndpost rows and nrow(x.test) columns.}
#' \item{yhat.train.mean}{Posterior mean of MCMC draws of traning data fits}
#' \item{yhat.test.mean}{Posterior mean of MCMC draws of testing data fits}
#' \item{sigma}{draws of random error vairaince, length = nskip+ndpost}
#' \item{tree_history}{If save_trees = TRUE, then a list of all trees}
#' \item{tree_proposal_total}{A ntree by length(p_modify) matrix, the (i,j)th entry is the total number of jth proposal of ith tree}
#' \item{tree_proposal_accept}{A ntree by length(p_modify) matrix, the (i,j)th entry is the total number of accepted jth proposal of ith tree}
#' \item{tree_leaf_count}{Number of leaf nodes in each tree}
#' @author Dongyue Xie: \email{dongyxie@gmail.com}
#' @references Chipman, H., George, E., and McCulloch R. (2010) Bayesian Additive Regression Trees. The Annals of Applied Statistics, 4,1, 266-298 <doi:10.1214/09-AOAS285>.
#' @importFrom BART rtnorm
#' @export


pBARTr=function(X,y,x.test,cutoff=0.5,
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

  fmean=mean(y)

  #whcih y = 1; y=0
  y1.idx = which(y==1)
  y0.idx = which(y==0)

  #priors: nu,lambda,tau
  tau = 3/(k*sqrt(ntree))
  if(length(binaryOffset)==0){binaryOffset=qnorm(fmean)}

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

  for (i in 1:(total_iter)) {
    if(i%%printevery==0){print(sprintf("done %d (out of %d)",i,total_iter))};
    if(save_trees){tree_history[[i]]=treelist}

    if(pre_train){
      if(i<=n_pre_train){rule = 'bart'}else{rule = split_rule}
    }

    bm.f = colSums(yhat.train.j)
    y.train = c()

    y.train[y1.idx] = rtnorm(length(y1.idx),bm.f[y1.idx],1,-binaryOffset)
    y.train[y0.idx] = -rtnorm(length(y0.idx),-bm.f[y0.idx],1,binaryOffset)

    #propose modification to each tree
    for (j in 1:ntree) {
      Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
      sig2 = 1

      #t1 = Sys.time()

      BART_draw = BARTr_train(X,Rj,treelist[[j]],p_modify,Tmin,
                              rule,sig2,tau,base,power,p_split,r)

      #t2 = Sys.time()

      #print('---------')
      #if((t2-t1)>0.01){
      #  print(BART_draw$move)
      #}
      #print(t2-t1)

      alpha = BART_draw$alpha
      new_treej = BART_draw$new_treej
      move = BART_draw$move

      tree_proposal_total[j,move]=tree_proposal_total[j,move]+1

      A=runif(1)

      if(is.nan(alpha)){
        alpha=0
      }

      t3 = Sys.time()

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

        t4 = Sys.time()

        #print(t4-t3)
        #print('aaaaaaaaa')


      }else{
        if(i<=nskip){
          hat=yhat.draw.train(treelist[[j]],Rj,tau,sig2)
          yhat.train.j[j,] = hat
        }else{
          hat=yhat.draw2(treelist[[j]],x.test,Rj,tau,sig2)
          yhat.train.j[j,] = hat$yhat
          yhat.test.j[j,] = hat$ypred
        }

        t4 = Sys.time()

        #print(t4-t3)
        #print('rrrrrrrrr')
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
