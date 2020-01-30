#' @title Rotation BART for classification
#' @param rotate 'rr'='random rotation'','rraug'='random rotation+augmentation','srp'='sparse random projection'
#' @param srpk the number of cols of sparse projection matrix
#' @param Others See ?BARTr
#' @return See ?BARTr
#' @importFrom BART rtnorm
#' @export


RotpBART=function(X,y,x.test,cutoff=0.5,
                  k=2.0, binaryOffset=NULL,
                  power=2.0, base=.95,w=rep(1,length(y)),
                  ntree=50,ndpost=700,nskip=300,Tmin=2,printevery=100,p_modify=c(0.5,0.5,0),
                  save_trees=F,rule='bart',p_split='CGM',rotate = 'rr',srpk=2*ncol(X)){

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
    y.train[y1.idx] = rtnorm(length(y1.idx),bm.f[y1.idx],1,-binaryOffset)
    y.train[y0.idx] = -rtnorm(length(y0.idx),-bm.f[y0.idx],1,binaryOffset)
    # for(ni in 1:n){
    #   if(y[ni]==1){
    #     y.train[ni] = rtnorm(bm.f[ni],-binaryOffset,1)
    #   }else{
    #     y.train[ni] = -rtnorm(-bm.f[ni],binaryOffset,1)
    #   }
    # }

    #propose modification to each tree
    for (j in 1:ntree) {
      Rj=y.train-colSums(yhat.train.j[-j,,drop=F])
      sig2 = 1

      BART_draw = BARTr_train(X[,,j],Rj,treelist[[j]],p_modify,Tmin,
                              rule,sig2,tau,base,power,p_split,r)


      alpha = BART_draw$alpha
      new_treej = BART_draw$new_treej
      move = BART_draw$move

      A=runif(1)

      tree_proposal_total[j,move]=tree_proposal_total[j,move]+1

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

        #hat=yhat.draw(new_treej,x.test,Rj,tau,sigma_draw[i]^2)
        #hat=yhat.draw.linear(new_treej,X,x.test,Rj)
        if(i<=nskip){
          hat=yhat.draw.train(new_treej,Rj,tau,sig2)
          yhat.train.j[j,] = hat
        }else{
          hat=yhat.draw(new_treej,x.test[,,j],Rj,tau,sig2)
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
          hat=yhat.draw2(treelist[[j]],x.test[,,j],Rj,tau,sig2)
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
              tree_proposal_total=tree_proposal_total,tree_proposal_accept=tree_proposal_accept,
              tree_leaf_count=tree_leaf_count))
}
