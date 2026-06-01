##copula Gaussian Graphical Model
###setting 1: non-paranormal,  N=150,300,500,800;p=30,50,100
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
source("algorithms.r")
t=1
Times=100
MSE=rep(0,Times)
KL=rep(0,Times)
while (t<=Times) {
  set.seed(2*t)
  ##p varies among 30,50,100
  p=30
  ##N varies among 150,300,500,800
  N=150
  n_cdf=N*3/10
  n_power=N*3/10
  n_linear=N*4/10
  s=0.125
  max_degree=4
  mu1=rep(0,p)
  mu2=rep(0,p)
  mu3=rep(0,p)
  M=10
  ##data generate
  adj_matrix=generate_neighborhood_graph(p,s,max_degree)
  Omega=build_precision_matrix(adj_matrix)
  Omega=ensure_positive_definite(Omega)
  ##cdf
  data_cdf=generate_nonparanormal_data(n_cdf,Omega,mu1,transform_type = "cdf")
  X_cdf=data_cdf$X
  Sigma=data_cdf$Sigma
  ##power
  data_power=generate_nonparanormal_data(n_power,Omega,mu2,transform_type = "power")
  X_power=data_power$X
  ##linear
  data_linear=generate_nonparanormal_data(n_linear,Omega,mu3,transform_type = "linear")
  X_linear=data_linear$X
  ####combine
  X_list_unscaled=list(X_cdf=X_cdf,X_power=X_power,X_linear=X_linear)
  
  ######standardize
  X_list=lapply(X_list_unscaled,function(mat) {
    scale(mat, center = TRUE, scale = TRUE)
  })
  X_matrix=do.call(rbind,X_list)
  ##############################################################
  ##HGCMA
  res_hete=cv_criteria(X_list = X_list,kfold = 5,model_num = M,method = "glasso")
  Omega_can_hete=res_hete$Omega_can
  S_test_hete=res_hete$S_est_test
  path_hete=res_hete$path
  Omega_oracle_together=res_hete$Omega_hat_MA
  mse=sum((Omega_oracle_together-Omega)^2)
  kl=log(det(Omega)/det(Omega_oracle_together))+
    sum(diag(Omega_oracle_together%*%solve(Omega)))-p
  MSE[t]=mse
  KL[t]=mse
  
  print(paste(t,"time finished"))
  t=t+1
}


