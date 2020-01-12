#'@title Evaluate liklihood ratio of new tree and old tree for GROW
#'@param old_t_data a list of old tree leaf data index
#'@param new_t_data a list of new tree leaf data index
#'@export

lik_ratio_grow = function(old_t_data,new_t_data,Rj,Tmin,sigma2,tau){
  tau2 = tau^2
  Rl = Rj[unlist(old_t_data)]
  nl = length(Rl)

  RlL = lapply(new_t_data, function(x) Rj[x])
  RlR = unlist(RlL[2])
  nlR = length(RlR)
  RlL = unlist(RlL[1])
  nlL = length(RlL)

  if(nlL<Tmin | nlR<Tmin){
    ratio = 0
  }else{
    ratio1 = sqrt(sigma2*(sigma2+nl*tau2)/((sigma2+nlL*tau2)*(sigma2+nlR*tau2)))
    ratio2 = tau2/(2*sigma2)*(sum(RlL)^2/(sigma2+nlL*tau2) + sum(RlR)^2/(sigma2+nlR*tau2) - sum(Rl)^2/(sigma2+nl*tau2))
    ratio2 = exp(ratio2)

    ratio = ratio1*ratio2
  }

  ratio

}
