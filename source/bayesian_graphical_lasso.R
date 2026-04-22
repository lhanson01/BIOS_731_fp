library(brms)
library(glasso) # <-- compare posterior mode to MAP of posterior

eval_class <- function(est_labels, true_labels){
  conf <- table(est_labels, true_labels)
  tn = conf[1,1]
  tp = conf[2,2]
  fp = conf[2,1]
  fn = conf[1,2]
  dt1 = (tp+fp)*(tp+fn)
  dt2 = (tn+fp)*(tn+fn)
  return(list(
      mcc = (tp * tn - fp * fn) / sqrt(exp(log(dt1)+log(dt2))),
      spec = tn / (tn + fp),
      sens = tp / (tp + fn)
    )
  )
}

bayesian_graphical_lasso <- function(data, r = 1, s_par = 0.01, n_iter, n_burn,
                                    simID = 1){
  
  Y <- data$Y
  W_true <- data$W_true
  W_true_vec <- W_true[upper.tri(W_true, diag = TRUE)]
  D <- dim(as.matrix(Y))[1]
  nt <- dim(as.matrix(Y))[2]
  S <- Y%*%t(Y)
  
  comp_time_glasso <- system.time({
    glasso <- do_glasso(Y)
  })

  lambda_init <- glasso$lambda #rgamma(1,
  #                      shape = r,
  #                      scale = s_par)
  W_init <- glasso$Omega
  glasso_labels <- as.numeric(abs(W_init[upper.tri(W_init, diag = TRUE)]) > 0)
  stein_glasso <-
    sum(diag(solve(W_init)%*%W_true)) - log(det(solve(W_init)%*%W_true)) - D
  
  
  comp_time_gibbs <- system.time({
    
  tau_init <- matrix(NA, nrow = D, ncol = D)
  for(j in 1:ncol(tau_init)){
    i = 1
    while(i < j){
      tau_init[i,j] <- tau_init[j,i] <- rexp(n = 1, 
                              rate = lambda_init^2 / 2)
      i = i + 1
    }
  }
  diag(tau_init) <- 0
  
  W <- W_init
  lambda <- lambda_init
  tau <- tau_init
  W_hist <- matrix(NA, nrow = n_iter, ncol = D*(D+1)/2)
  lambda_hist <- rep(NA, n_iter)
  tau_hist <- matrix(NA, nrow = n_iter, ncol = D*(D+1)/2)
  W_hist[1,] <- W[upper.tri(W, diag = TRUE)]
  lambda_hist[1] <- lambda
  tau_hist[1,] <- tau[upper.tri(tau, diag = TRUE)]
  
  for(iter in 2:n_iter){
    print(iter)
    # Sample omega
    for(i in 1:D){
      ### PARTITION ###
      block_w <- W[-i,-i]
      w_col <- W[-i,i]
      w_corner <- W[i,i]
      
      block_s <- S[-i,-i]
      s_col <- S[-i,i]
      s_corner <- S[i,i]
      
      block_tau <- tau[-i,-i]
      tau_col <- tau[-i,i]
      
      C <- solve((s_corner + lambda) * solve(block_w) + solve(diag(tau_col)))
      
      gamma <- rgamma(1, nt/2 + 1, (s_corner+lambda)/2)
      beta <- MASS::mvrnorm(1, -C %*% s_col, C)
      
      W[-i,i] <- W[i,-i] <- beta
      W[i,i] <- gamma + t(beta)%*%solve(block_w)%*%beta
      
    }
    
    lambda <- rgamma(1, shape = r + D*(D+1)/2, 
                     rate = s_par + sum(abs(W)))
    print(lambda)

    ### update tau ###
    W_vec <- as.numeric(W[upper.tri(W)])
    tau_vec <- rep(NA, length(W_vec))
    for(w in 1:length(W_vec)){
      mu_prime <- sqrt(lambda^2/ W_vec[w]^2)
      u <- brms::rinv_gaussian(1, mu_prime, lambda^2)
      tau_vec[w] <- 1/u
    }
    tau[upper.tri(tau)] <- tau[lower.tri(tau)] <- tau_vec
  
    W_hist[iter, ] <- W[upper.tri(W, diag = TRUE)]
    lambda_hist[iter] <- lambda
    tau_hist[iter, ] <- tau[upper.tri(tau, diag = TRUE)]
  }
  
  
  W_post_mean <- colMeans(W_hist[(n_burn+1):n_iter,])
  lambda_post_mean <- mean(lambda_hist[(n_burn+1):n_iter])
  tau_post_mean <- colMeans(tau_hist[(n_burn+1):n_iter,])
  W_cred <- apply(W_hist, MARGIN = 2, FUN = function(x) quantile(x, c(0.025, 0.975)))
  covers <- sapply(1:length(W_post_mean), function(i){
    W_true_vec[i] > W_cred[1,i] & W_true_vec[i] < W_cred[2,i] 
  })
  crude_class <- sapply(1:length(W_post_mean), function(i){
    edge <- !(0 > W_cred[1,i] & 0 < W_cred[2,i]) 
    return(as.numeric(edge))
  })
  
  view_estimates <- rbind(W_cred, W_true_vec)
  
  W_mat <- matrix(NA, D, D)
  W_mat[upper.tri(W_mat, diag = TRUE)] <- W_post_mean
  lower_W <- t(W_mat)
  W_mat[lower.tri(W_mat)] <- lower_W[lower.tri(lower_W)] 
  stein_gibbs <- 
    sum(diag(solve(W_mat)%*%W_true)) - log(det(solve(W_mat)%*%W_true)) - D
  
  #est_adj <- matrix(NA, D, D)
  #est_adj[upper.tri(est_adj, diag = TRUE)] <- evaluate_edge(W_mat, S, nt, D)
  
  #est <- est_adj[upper.tri(est_adj, diag = TRUE)]
  
  gibbs_labels <- crude_class
  true_labels <- data$edge_true
  gibbs_perf <- eval_class(gibbs_labels, true_labels)
  glasso_perf <- eval_class(gibbs_labels, true_labels)
  })
  
  
    
  output <- list(results_tib = tibble(
                   simID = simID,
                   D = D,
                   mcc_gibbs = gibbs_perf$mcc,
                   spec_gibbs = gibbs_perf$spec,
                   sens_gibbs = gibbs_perf$sens,
                   comp_time_gibbs = comp_time_gibbs[3],
                   lambda_gibbs = lambda_post_mean,
                   stein_gibbs = stein_gibbs,
                   omega_coverage = sum(covers)/length(covers),
                   mcc_glasso = glasso_perf$mcc,
                   spec_glasso = glasso_perf$spec,
                   sens_glasso = glasso_perf$sens,
                   comp_time_glasso = comp_time_glasso[3],
                   stein_glasso = stein_glasso,
                   ),
                 W_hist = W_hist,
                 W_true = W_true_vec,
                 g = list(data$g)
              )
                 
  
  return(output)
  
}


##### For testing, delete later ####

