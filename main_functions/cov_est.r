##sample covariance estimator
library(huge)
cov_est=function(X_list){
  D=length(X_list)
  p=dim(X_list[[1]])[2]
  S=list()
  for (d in 1:D) {
    S[[d]]=huge.npn(X_list[[d]],npn.func = "skeptic")
  }
  S.hat=matrix(0,p,p)
  for (d in 1:D) {
    S.hat=S.hat+S[[d]]
  }
  S.hat=S.hat/D
  return(list(S.hat=S.hat,
         S.list=S))
  
}