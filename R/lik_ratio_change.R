#'@title Evaluate liklihood ratio of new tree and old tree for CHANGE
#'@export

lik_ratio_change = function(old_t_data,new_t_data,Rj,Tmin,sigma2,tau){
  tau2 = tau^2

  R1New = lapply(new_t_data, function(x) Rj[x])
  R2New = unlist(R1New[2])
  n2New = length(R2New)
  R1New = unlist(R1New[1])
  n1New = length(R1New)

  if(n1New<Tmin | n2New<Tmin){
    ratio = 0
  }else{
    R1Old = lapply(old_t_data, function(x) Rj[x])
    R2Old = unlist(R1Old[2])
    n2Old = length(R2Old)
    R1Old = unlist(R1Old[1])
    n1Old = length(R1Old)

    st2 = sigma2/tau2
    ratio1 = sqrt((st2+n1Old)*(st2+n2Old)/((st2+n1New)*(st2+n2New)))
    ratio2 = 1/(2*sigma2)*(sum(R1New)^2/(n1New+st2)+sum(R2New)^2/(n2New+st2)-sum(R1Old)^2/(n1Old+st2)-sum(R2Old)^2/(n2Old+st2))
    ratio2 = exp(ratio2)
    ratio = ratio1*ratio2
  }

  ratio

}
