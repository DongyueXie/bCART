#'@title cross-validated pbart
#'@return predicted label and test AUC: ypred,Test_AUC,best_params
#'@export

pbartCV = function(x,y,x.test,y.test,ntrees=c(50,200),ks=c(0.5,2,5),powers = c(2,3,5,10,nrow(x)),
                   k_folds=5,verbose=F){
  n_combns = length(ntrees)*length(ks)*length(powers)
  results = matrix(nrow=n_combns,ncol=4)
  #create folds
  folds = lapply(split(sample(1:nrow(x)),1:k_folds),sort)
  #folds = createFolds(1:nrow(x),k_folds)

  y_folds = lapply(folds,function(z){y[z]})

  while(sum(unlist(lapply(y_folds,function(z){length(unique(z))}))==1)>0){
    folds = lapply(split(sample(1:nrow(x)),1:k_folds),sort)
    #folds = createFolds(1:nrow(x),k_folds)
    y_folds = lapply(folds,function(z){y[z]})
    }

  idx = 0
  for(it in (ntrees)){
    for(jk in (ks)){
      for(kp in (powers)){
        idx = idx + 1
        if(verbose & idx%%5==0){print(idx)}


          cv_AUC = lapply(folds,function(z){
            x_train = x[-z,]
            y_train = y[-z]
            x_test = x[z,]
            y_test = y[z]
            prob.test = quiet(pbart(x_train,y_train,x_test,ntree=it,k=jk,power=kp,nkeeptrain=0,nkeeptreedraws=0))$prob.test.mean
            Test_AUC = AUC(y_test,as.numeric(prob.test>0.5))
            Test_AUC
          })

        cv_AUC = mean(unlist(cv_AUC))
        results[idx,] = c(cv_AUC,it,jk,kp)
      }
    }
  }

  best_params = results[which.max(results[,1]),]

  prob.test = quiet(pbart(x,y,x.test,ntree=best_params[2],k=best_params[3],power=best_params[4],nkeeptrain=0,nkeeptreedraws=0))$prob.test.mean

  ypred = as.numeric(prob.test>0.5)
  Test_AUC = AUC(y.test,ypred)

  return(list(ypred=ypred,Test_AUC=Test_AUC,best_params=best_params))

}

