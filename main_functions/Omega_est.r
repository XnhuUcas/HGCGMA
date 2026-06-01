####estimate Omega
Omega_est=function(S.hat,path_sub,model_num){
  ###transform path_sub to a long vector
  p=dim(path_sub)[1]
  path_sub_vec=path_sub[upper.tri(path_sub, diag = FALSE)]
  res_Ome=estimate_Omega(path_sub_vec,S.hat)
  return(list(Omega_hat=res_Ome))
}


