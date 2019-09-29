#'@title cross-validated bart
#'@return predicted label and test AUC: ypred,Test_AUC,best_params
#'@export

pbartCV = function(x,y,x.test,y.test,ntrees=c(50,200),ks=c(0.5,2,5),powers = c(2,3,5,10,nrow(x)),
                   k_folds=5,verbose=F,isparallel=T){
  n_combns = length(ntrees)*length(ks)*length(powers)
  results = matrix(nrow=n_combns,ncol=4)
  #create folds
  folds = createFolds(1:nrow(x),k_folds)

  idx = 0
  for(it in (ntrees)){
    for(jk in (ks)){
      for(kp in (powers)){
        idx = idx + 1
        if(verbose & idx%%5==0){print(idx)}

        if(isparallel){
          cv_AUC = mclapply(folds,function(z){
            x_train = x[-z,]
            y_train = y[-z]
            x_test = x[z,]
            y_test = y[z]
            fit = quiet(pbart(x_train,y_train,x_test,ntree=it,k=jk,power=kp))
            Test_AUC = AUC(y_test,as.numeric(colSums(pnorm(fit$yhat.test))/nrow(fit$yhat.test)>0.5))
            Test_AUC
          },mc.cores = k_folds)
        }else{
          cv_AUC = lapply(folds,function(z){
            x_train = x[-z,]
            y_train = y[-z]
            x_test = x[z,]
            y_test = y[z]
            fit = quiet(pbart(x_train,y_train,x_test,ntree=it,k=jk,power=kp))
            Test_AUC = AUC(y_test,as.numeric(colSums(pnorm(fit$yhat.test))/nrow(fit$yhat.test)>0.5))
            Test_AUC
          })
        }
        cv_AUC = mean(unlist(cv_AUC))
        results[idx,] = c(cv_AUC,it,jk,kp)
      }
    }
  }

  best_params = results[which.max(results[,1]),]

  fit = quiet(pbart(x,y,x.test,ntree=best_params[2],k=best_params[3],power=best_params[4]))

  ypred = as.numeric(colSums(pnorm(fit$yhat.test))/nrow(fit$yhat.test)>0.5)
  Test_AUC = AUC(y.test,ypred)

  return(list(ypred=ypred,Test_AUC=Test_AUC,best_params=best_params))

}

