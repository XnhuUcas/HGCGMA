#Data generate
generate_neighborhood_graph <- function(p, s = 0.125, max_degree = 4) {
  # 1. 为每个节点生成二维坐标
  Y <- matrix(runif(2 * p), ncol = 2)
  
  # 2. 初始化邻接矩阵
  adj_matrix <- matrix(0, p, p)
  
  # 3. 根据距离概率生成边
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      distance <- sqrt(sum((Y[i,] - Y[j,])^2))
      prob_edge <- (1/sqrt(2*pi)) * exp(-distance^2/(2*s))
      
      if (runif(1) < prob_edge) {
        adj_matrix[i, j] <- 1
        adj_matrix[j, i] <- 1
      }
    }
  }
  
  # 4. 限制最大度数为4
  degrees <- rowSums(adj_matrix)
  while (any(degrees > max_degree)) {
    for (i in which(degrees > max_degree)) {
      # 找到该节点的邻居
      neighbors <- which(adj_matrix[i,] == 1)
      # 随机移除多余的边
      if (length(neighbors) > max_degree) {
        to_remove <- sample(neighbors, length(neighbors) - max_degree)
        adj_matrix[i, to_remove] <- 0
        adj_matrix[to_remove, i] <- 0
      }
    }
    degrees <- rowSums(adj_matrix)
  }
  
  return(adj_matrix)
}

#' 构建逆协方差矩阵
#' @param adj_matrix 邻接矩阵
build_precision_matrix <- function(adj_matrix) {
  p <- nrow(adj_matrix)
  Omega <- adj_matrix * 0.245
  diag(Omega) <- 1
  return(Omega)
}

#' 检查矩阵正定性并调整（如果需要）
ensure_positive_definite <- function(Omega) {
  # 检查特征值
  eigen_vals <- eigen(Omega, only.values = TRUE)$values
  
  if (any(eigen_vals <= 0)) {
    # 如果非正定，添加小的正则化项
    min_eigen <- min(eigen_vals)
    regularization <- ifelse(min_eigen <= 0, abs(min_eigen) + 0.1, 0)
    Omega <- Omega + diag(regularization, nrow(Omega))
    cat("添加正则化项:", regularization, "以确保正定性\n")
  }
  
  return(Omega)
}

#' 生成非paranormal数据
#' @param n 样本量
#' @param Omega 逆协方差矩阵
#' @param mu 均值向量，默认为(0,...,0)
#' @param transform_type 转换类型："cdf" 或 "power"
generate_nonparanormal_data <- function(n, Omega, mu = NULL, transform_type = "cdf") {
  p <- nrow(Omega)
  
  if (is.null(mu)) {
    mu <- rep(0, p)
  }
  
  # 计算协方差矩阵
  Sigma_0 <- solve(Omega)
  Sigma=cov2cor(Sigma_0)
  # 生成多元正态数据 Z ~ N(0, Sigma)
  Z <- mvrnorm(n, mu, Sigma)
  
  # 应用逆变换得到非paranormal数据 X = g(Z)
  X <- matrix(0, n, p)
  
  for (j in 1:p) {
    if (transform_type == "cdf") {
      # 高斯CDF变换（简化版本）
      X[, j] <- g_cdf_transform(Z[, j], mu[j], sqrt(Sigma[j, j]))
    } else if (transform_type == "power") {
      # 对称幂变换（简化版本）
      X[, j] <- g_power_transform(Z[, j], mu[j], sqrt(Sigma[j, j]))
    } else if(transform_type == "linear"){
      X[,j]=Z[,j]
    }
  }
  
  return(list(X=X,Sigma=Sigma))
}

#' 高斯CDF变换（简化实现）
g_cdf_transform <- function(z, mu_j, sigma_j, mu_g0 = 0.05, sigma_g0 = 0.4, n_samples = 10000) {
  # 定义基础变换函数 g0
  g0 <- function(t) {
    pnorm(t, mean = mu_g0, sd = sigma_g0)
  }
  
  # 使用蒙特卡洛方法计算积分
  # 生成服从N(mu_j, sigma_j^2)分布的样本
  t_samples <- rnorm(n_samples, mean = mu_j, sd = sigma_j)
  
  # 计算 E[g0(t)] = ∫ g0(t) φ((t-μ_j)/σ_j) dt
  Eg0 <- mean(g0(t_samples))
  
  # 计算 Var[g0(t)] = ∫ (g0(t) - E[g0(t)])^2 φ((t-μ_j)/σ_j) dt
  g0_samples <- g0(t_samples)
  Varg0 <- mean((g0_samples - Eg0)^2)
  
  # 应用变换
  result <- sigma_j*(g0(z) - Eg0) / sqrt(Varg0) +mu_j
  
  return(result)
}

#' 对称幂变换（简化实现）
g_power_transform <- function(z, mu_j, sigma_j, alpha = 3, n_samples = 10000) {
  # 定义基础变换函数 g0
  g0 <- function(t) {
    sign(t) * abs(t)^alpha
  }
  
  # 使用蒙特卡洛方法计算积分
  # 生成服从N(0, sigma_j^2)分布的样本（因为t-mu_j ~ N(0, sigma_j^2)）
  t_samples <- rnorm(n_samples, mean = 0, sd = sigma_j)
  
  # 计算 E[g0^2(t)] = ∫ g0^2(t) φ(t/σ_j) dt
  g0_sq_samples <- g0(t_samples)^2
  Eg0_sq <- mean(g0_sq_samples)
  
  # 应用变换
  centered_z <- z - mu_j
  result <- sigma_j*g0(centered_z) / sqrt(Eg0_sq)+mu_j
  
  return(result)
}

#' 可视化生成的图
visualize_graph <- function(adj_matrix) {
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  plot(g, 
       vertex.size = 8,
       vertex.label.cex = 0.7,
       vertex.color = "lightblue",
       edge.color = "gray",
       main = "adjacency graph")
}