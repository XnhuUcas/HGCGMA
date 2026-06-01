###candidate models
library(huge)
library(glasso)
library(cluster)
library(factoextra)
library(dbscan)
library(pracma)
library(NlcOptim)
CAN=function(S.hat,model_num,method="glasso"){
  result=huge(x=S.hat,nlambda = model_num,method = method)
  path=result$path
  
  lambda=result$lambda
  Omega_can=result$icov
  
  return(list(path=path,lambda=lambda,Omega_can=Omega_can))
}

##apply CV to calculate the weight
cv_criteria=function(X_list, kfold=5,model_num=10,method="glasso"){
  S.hat.total=cov_est(X_list)$S.hat
  candidate=CAN(S.hat.total,model_num = model_num,method = method)
  path=candidate$path
  lambda=candidate$lambda
  D=length(X_list)
  X_matrix=do.call(rbind,X_list)
  p=dim(X_matrix)[2]
  n.total=dim(X_matrix)[1]
  labels=rep(1:D,times=sapply(X_list,nrow))
  shuffle_index <- sample(nrow(X_matrix))
  X_matrix_shuffled <- X_matrix[shuffle_index, ]
  labels_shuffled <- labels[shuffle_index]
  n.sub=floor(n.total/kfold)
  S_est_test=list()
  S_est_train=list()
  X_sep=list()
  for (k in 1:kfold) {
    X.test=X_matrix[(((k-1)*n.sub+1):(k*n.sub)),]
    X.train=X_matrix[-(((k-1)*n.sub+1):(k*n.sub)),]
    labels_test=labels_shuffled[(((k-1)*n.sub+1):(k*n.sub))]
    labels_train=labels_shuffled[-(((k-1)*n.sub+1):(k*n.sub))]
    ####come back to list
    data_test=vector("list",length=D)
    for (d in 1:D) {
      indices1=which(labels_test == d)
      if(length(indices1) > 0) {
        data_test[[d]] <- X.test[indices1, , drop = FALSE]
      } else {
        data_test[[d]] <- matrix(0, ncol = 2)
      }
    }
    data_train=vector("list",length = D)
    for (d in 1:D) {
      indices2=which(labels_train == d)
      if(length(indices2) > 0) {
        data_train[[d]] <- X.train[indices2, , drop = FALSE]
      } else {
        data_train[[d]] <- matrix(0, ncol = 2)
      }
    }
    data_test=data_test[sapply(data_test, function(mat) { !all(mat == 0)})]
    data_train=data_train[sapply(data_train, function(mat) { !all(mat == 0)})]
    S_est_test[[k]]=cov_est(data_test)$S.hat
    S_est_train[[k]]=cov_est(data_train)$S.hat
  }
  Omega_hat_cv=array(0,dim = c(p,p,kfold,model_num))
  for (m in 1:model_num) {
    path.call=path[[m]]
    for (k in 1:kfold) {
      Omega_hat_cv[,,k,m]=Omega_est(S_est_train[[k]],path.call,model_num = model_num)$Omega_hat
    }
  }
  ###all samples to estimate
  Omega_can_all=array(0,dim = c(p,p,model_num))
  for (m in 1:model_num) {
    path.call=path[[m]]
    #Omega_can_all[,,m]=Omega_est(S.hat.total,path.call,model_num)$Omega_hat
    Omega_can_all[,,m]=candidate$Omega_can[[m]]
  }
  ##weight criteria
  w_0=rep(1/model_num,model_num)
  cv_fun=function(w){
    res=0
    for (kk in c(1: kfold)) {
      Omega_MA <- Reduce('+', lapply(c(1:model_num), function(ii) Omega_hat_cv[ , , kk, ii] * w[ii]))
      res <- res - log(det(Omega_MA)) + sum(diag(Omega_MA %*% S_est_test[[kk]]))  # use covariance of validation samples or population_covariance
    }
    return(res)
  }
  A <- rbind(diag(model_num), -diag(model_num))
  b <- matrix(c(rep(1, model_num), rep(0, model_num)), ncol=1)
  Aeq <- matrix(rep(1, model_num), nrow=1)
  beq <- 1
  w_solution <- pracma::fmincon(x0=w_0, fn=cv_fun, A=A, b=b, Aeq=Aeq, beq=beq)
  w_hat=w_solution$par
  
  Omega_hat_MA=Reduce('+', lapply(c(1:model_num), function(i) Omega_can_all[ , , i] * w_hat[i]))
  
  
  return(list(Omega_can=Omega_can_all,
              S_est_test=S_est_test,
              S_est_train=S_est_train,
              path=path,
              Omega_hat_MA=Omega_hat_MA,
              weight=w_hat,
              lambda=lambda))
}
