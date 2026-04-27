library(igraph)
library(dplyr)


sim_graph_ts <- function(D, nt, p_mean, pwr = 0.8, rewire_p = 0.004){
  
  g_init <- igraph::sample_pa(D, power = pwr, m = 1, directed = FALSE,
                              out.pref = TRUE)
  #rewire_prob <- ifelse(degree(g_init)==1, 0, rewire_p)
  comp <- 2
  while(comp > 1){
    print(comp)
    g <- g_init %>% rewire(each_edge(p = rewire_p))
    comp <- components(g)$no
  }
  #plot(g)
  adj_g <- matrix(igraph::as_adjacency_matrix(g, type = "both"), nrow = D, ncol = D)
  adj_true <- adj_g[upper.tri(adj_g, diag = TRUE)]
  
  W_true <- matrix(0,D,D)
  p_entries <- adj_g[upper.tri(adj_g)] * 
    rnorm(D*(D-1)/2, mean = p_mean, sd = 2) * sample(c(-1,1), prob = c(0.7,0.3))
  W_true[upper.tri(W_true)] <- p_entries
  lower_tri <- t(W_true)
  W_true[lower.tri(W_true)] <- lower_tri[lower.tri(lower_tri)] 
  diag(W_true) <- rowSums(abs(W_true)) + 1e-02
  print(diag(W_true))

  
    sigma <- solve(W_true) # should be pos_def
    L <- chol(sigma)
    adj_est <- solve(t(L)%*%L)
  
  DxT <- matrix(NA, nrow = D, ncol = nt)
  d2 = 1
  for(t in 1:nt){
    z_d <- rnorm(D)
    DxT[,t] <- L%*%z_d 
  }

  return(list(Y = DxT,
              W_true = W_true,
              edge_true = adj_true,
              g = g)
  )
  
}



