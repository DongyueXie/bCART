#'@title Run reps of pbart on real data set
#'@export

pbart_reps = function(x,y,ntrain = floor(nrow(x)*0.5),nmaxtest=1000,nreps=20,seed=12345,
                      ntrees=c(50,200),ks=c(0.5,2,5),psplit='CGM',powers = c(2,3,5,10,100),k_folds=5){

  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")

  if(is.factor(y)){
    y = as.numeric(y)
    y = (y - min(y))/(max(y)-min(y))
    }

  cv_results = c()
  cv_best_params = c()

  if(psplit=='CGM'){
    powers = c(powers,ntrain)
  }
  if(psplit=='RS'){
    powers = 1/c(powers,ntrain)
  }

  for(rp in 1:nreps){
    train.idx = sample(1:nrow(x),ntrain)
    xtrain = x[train.idx,]
    ytrain = y[train.idx]
    if((nrow(x)-ntrain)>nmaxtest){
      test.idx = sample((1:nrow(x))[-train.idx],nmaxtest)
      xtest = x[test.idx,]
      ytest = y[test.idx]
    }else{
      xtest = x[-train.idx,]
      ytest = y[-train.idx]
    }

    cv_fit = pbartCV(xtrain,ytrain,xtest,ytest,ntrees,ks,powers,k_folds)

    cv_results[rp] = cv_fit$Test_AUC
    cv_best_params = rbind(cv_best_params,cv_fit$best_params)

  }

  return(list(cv_results=cv_results,cv_best_params=cv_best_params))
}
