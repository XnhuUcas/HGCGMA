####estimate Omega
# Constrained MLE
estimate_Omega <- function(graphical_model, covariance_train){
  # Omega is estimated by ELS Algorithm 17.1 in page 634
  p <- dim(covariance_train)[1]               # this p is the dim being chosen
  Omega_hat <- matrix(rep(0, p*p), p)
  # beta_hat is the core parameter to obtain Omega_hat
  beta_hat <- matrix(rep(0, p*p),p) # each col represent a node
  # adjacency_matrix is a 0-1 value matrix with size of p*p
  adjacency_matrix <- matrix(rep(0, p*p), p)
  adjacency_matrix[upper.tri(adjacency_matrix)] <- graphical_model
  adjacency_matrix <- adjacency_matrix + t(adjacency_matrix) + diag(p)
  # set loop
  convergence_status <- rep(FALSE, p)
  current_idx <- 0
  W_current <- W_next <- covariance_train
  iter_num <- 0
  while (!identical(convergence_status, rep(TRUE, p))) {
    current_idx <- current_idx + 1 - (current_idx %/% p)*p
    adjacency_nodes_with_current <- setdiff(which(adjacency_matrix[current_idx,] != 0), current_idx)
    # step (a): Partition W
    W11 <- W_next[-current_idx, -current_idx]
    W11_star <- W_next[adjacency_nodes_with_current, adjacency_nodes_with_current]
    s12_star <- covariance_train[current_idx, adjacency_nodes_with_current]
    if (length(adjacency_nodes_with_current) > 0){
      # it is possible that current node does not connect to any other node,
      # then W11_star is not reversible.
      # step (b)
      beta_star <- solve(as.matrix(W11_star)) %*% s12_star
      beta_hat[adjacency_nodes_with_current, current_idx] <- beta_star
    }
    # step (C), W12 is feasible, because beta_hat is zero vector.
    W12 <- W11 %*% beta_hat[-current_idx, current_idx]
    # Update W_next
    W_next[-current_idx, current_idx] <- W12
    iter_num <- iter_num + 1
    # Termination condition, converge for every idx
    if (norm(W_current-W_next, "F") < 1e-5) {
      # This condition is easy to achieve.
      convergence_status[current_idx] = TRUE
      W_hat <- W_next
    }
    W_current <- W_next
  }
  # Calculate Omega_hat, Omega12 and Omega22
  for (current_idx in c(1:p)) {
    s22 <- covariance_train[current_idx, current_idx]
    W12 <- W_hat[-current_idx, current_idx]
    Omega22 <- 1 / (s22 - W12 %*% beta_hat[-current_idx, current_idx])
    Omega12 <- - beta_hat[-current_idx,current_idx] * as.numeric(Omega22)
    # Omega12 and Omega21 and Omega22 could be updated.
    Omega_hat[-current_idx, current_idx] <- Omega12
    Omega_hat[current_idx, -current_idx] <- Omega12
    Omega_hat[current_idx, current_idx] <- Omega22
  }
  # # debug
  # print(norm(Omega_hat-solve(covariance_train), "F"))
  # # end debug
  gc()    # clean memory
  
  return(Omega_hat)
}

