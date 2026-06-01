##copula Gaussian Graphical Model
###Setting 2:Gaussian model,N=150,300,500,800,p=30,50,100
rm(list = ls())
library(MASS)
library(Matrix)
library(igraph)
library(kernlab)
source("candidate.r")
source("data_generate.r")
source("cov_est.r")
source("Omega_est.r")
source("estimate_Omega.r")
##data generate
t=1
Times=100
MSE=rep(0,Times)
KL=rep(0,Times)
while (t<=Times) {
  set.seed(2*t)
  p=30
  N=150
  n1=N*3/10
  n2=N*3/10
  n3=N*4/10
  s=0.125
  max_degree=4
  M=10
  adj_matrix=generate_neighborhood_graph(p,s,max_degree)
  Omega=build_precision_matrix(adj_matrix)
  Omega=ensure_positive_definite(Omega)
  data=generate_nonparanormal_data(N,Omega,mu=rep(0,p),transform_type = "linear")
  X=data$X
  Sigma=data$Sigma
  
  X_list_unscaled=list(X)
  
  ######standardize
  X_list=lapply(X_list_unscaled,function(mat) {
    scale(mat, center = TRUE, scale = TRUE)
  })
  X_matrix=do.call(rbind,X_list)
  ##############################################################
  ###HGCGMA
  res_hete=cv_criteria(X_list = X_list,kfold = 5,model_num = M,method = "glasso")
  Omega_can_hete=res_hete$Omega_can
  S_test_hete=res_hete$S_est_test
  path_hete=res_hete$path
  Omega_oracle_together=res_hete$Omega_hat_MA
  mse=sum((Omega_oracle_together-Omega)^2)
  kl=log(det(Omega)/det(Omega_oracle_together))+
    sum(diag(Omega_oracle_together%*%solve(Omega)))-p
  
  MSE[t]=mse
  KL[t]=kl
  
  print(paste(t,"time finished"))
  t=t+1
}
