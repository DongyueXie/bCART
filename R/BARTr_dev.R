#' @title Fit a BART model
#' @param rule grp: Gaussian random projection; sgrp: sparse Gaussian random projection; bart: originla bart; hyperplane: connect two points
#' @export


BARTr_dev=function(X,y,x.test,sigdf=3, sigquant=.90,
               k=2.0, lambda=NA, sigest=NA,sigmaf=NA,
               power=2.0, base=.95,w=rep(1,length(y)),
               ntree=50,ndpost=300,nskip=100,Tmin=2,printevery=100,p_modify=c(2.5, 2.5, 4)/9,
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
  # s_pos: internal node index
  # s_dir: internal node projection
  # s_rule: internal node splitting value
  # s_data: internal node data point index
  # s_depth: internal node depth
  # s_obs: internal node number of observations
  # t_pos: leaf nodei index
  # t_data: terminal node data point index
  # t_depth: terminal node depth
  # t_test_data: terminal node testing data point index
  treelist=lapply(treelist, function(x){
    x=list(s_pos=NULL,s_dir=NULL,s_rule=NULL,s_data=NULL,s_depth=NULL,s_obs=NULL,
           t_pos=1,t_data=list(1:n),t_depth=0,t_test_data=NULL)
  })

  #statistics to save
  tree_history=list()#a list of ndpost lists and each of ndpost lists is a list of ntree lists.

  n_move = length(p_modify)

  tree_proposal_total=matrix(rep(0,ntree*n_move),nrow = ntree,ncol = n_move)
  tree_proposal_accept=matrix(rep(0,ntree*n_move),nrow=ntree,ncol=n_move)

  sigma_draw=c(sigest)

  yhat.train=matrix(nrow=ndpost,ncol=n)
  yhat.test=matrix(nrow=ndpost,ncol=nt)
  #initilize single terminal node trees

  yhat.train.j=matrix(rnorm(ntree*n,0,sqrt(1/(n/sigest^2+1/tau^2))),nrow=ntree,ncol=n)
  yhat.test.j=matrix(rep(0,ntree*nt),nrow=ntree,ncol=nt)

  sigma_draw_all = c()

  ############# pre-train standard BART ################

  if(pre_train){
    print(sprintf('Pre_training BART, total %d iterations',n_pre_train))
    for(i in 1:n_pre_train){
      for(j in 1:ntree){
        Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
        sig2 = sigma_draw[i]^2
        BART_draw = BARTr_train(X,Rj,treelist[[j]],p_modify,Tmin,
                                rule='bart',sig2,tau,base,power)
        Rj=y.train-colSums(yhat.train.j[-j,,drop=F])

        BART_draw = BARTr_train(X,Rj,treelist[[j]],p_modify,Tmin,
                                rule,sig2,tau,base,power)



        alpha = BART_draw$alpha
        new_treej = BART_draw$new_treej

        A=runif(1)

        if(is.nan(alpha)){
          alpha=0
        }
        if(A<alpha){
           hat=yhat.draw.train(new_treej,Rj,tau,sig2)
           yhat.train.j[j,] = hat
           treelist[[j]]=new_treej
        }else{
          hat=yhat.draw.train(treelist[[j]],Rj,tau,sig2)
          yhat.train.j[j,] = hat
        }
      }
      res=y.train-colSums(yhat.train.j)
      sigma_draw[i+1]=sqrt((nu*lambda + sum(res^2))/rchisq(1,n+nu))
    }
    print('Pre-training done')
    sigma_draw_all = sigma_draw
    sigma_draw = sigma_draw[length(sigma_draw)]
  }

################### pre-train done #################

################### start burn-in period ###########

  for(i in 1:nskip){
    if(i%%printevery==0){print(sprintf("Burn-in: done %d (out of %d)",i,nskip))};
    for(j in 1:ntree){
      Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
      sig2 = sigma_draw[i]^2

      BART_draw = BARTr_train(X,Rj,treelist[[j]],p_modify,Tmin,
                              rule,sig2,tau,base,power)



      alpha = BART_draw$alpha
      new_treej = BART_draw$new_treej
      A=runif(1)
      if(is.nan(alpha)){
        alpha=0
      }

      if(A<alpha){
       hat=yhat.draw.train(new_treej,Rj,tau,sig2)
       yhat.train.j[j,] = hat
       treelist[[j]]=new_treej
      }else{
       hat=yhat.draw.train(treelist[[j]],Rj,tau,sig2)
       yhat.train.j[j,] = hat
      }
    }
    res=y.train-colSums(yhat.train.j)
    sigma_draw[i+1]=sqrt((nu*lambda + sum(res^2))/rchisq(1,n+nu))
  }
  sigma_draw_all = c(sigma_draw_all,sigma_draw)
  sigma_draw = sigma_draw[length(sigma_draw)]


################## burn in done #################

################# start MCMC draws #############

  for (i in 1:(ndpost)) {
    if(i%%printevery==0){print(sprintf("MCMC draws: done %d (out of %d)",i,ndpost))};
    if(save_trees){tree_history[[i]]=treelist}
    #propose modification to each tree

    for (j in 1:ntree) {
      Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
      sig2 = sigma_draw[i]^2

      BART_draw = BARTr_train(X,Rj,treelist[[j]],p_modify,Tmin,
                              rule,sig2,tau,base,power)



      alpha = BART_draw$alpha
      new_treej = BART_draw$new_treej
      move = BART_draw$move

      tree_proposal_total[j,move]=tree_proposal_total[j,move]+1

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
      if(A<alpha){
        # we accept the new tree
        tree_proposal_accept[j,move]=tree_proposal_accept[j,move]+1
        hat=yhat.draw(new_treej,x.test,Rj,tau,sig2)
        yhat.train.j[j,] = hat$yhat
        yhat.test.j[j,] = hat$ypred
        new_treej$t_test_data = hat$t_idx


        treelist[[j]]=new_treej


      }else{
        hat=yhat.draw2(treelist[[j]],x.test,Rj,tau,sig2)
        yhat.train.j[j,] = hat$yhat
        yhat.test.j[j,] = hat$ypred
      }


    }


      yhat.train[i,]=colSums(yhat.train.j)
      yhat.test[i,]=colSums(yhat.test.j)
      #draw sigma
      res=y.train-yhat.train[i,]



    sigma_draw[i+1]=sqrt((nu*lambda + sum(res^2))/rchisq(1,n+nu))

  }


  yhat.train=yhat.train+fmean
  yhat.test=yhat.test+fmean
  sigma_draw=c(sigma_draw_all,sigma_draw)

  tree_leaf_count=as.numeric(unlist(lapply(treelist,function(x){length(x$t_data)})))

  return(list(yhat.train=yhat.train,yhat.test=yhat.test,
              yhat.train.mean=colSums(yhat.train)/nrow(yhat.train),
              yhat.test.mean=colSums(yhat.test)/nrow(yhat.test),sigma=sigma_draw,
              tree_history=tree_history,
              tree_proposal_total=tree_proposal_total,tree_proposal_accept=tree_proposal_accept,
              tree_leaf_count=tree_leaf_count))
}
