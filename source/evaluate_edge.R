library(dplyr)

evaluate_edge <- function(W, S, nt, D){
  W_frm <- W %>% reshape2::melt() %>%
    filter(Var1 <= Var2) 
  
  w_0 <- 3 
  U_0 <- diag(1, D)
  w_1 <- w_0 + nt
  U_1 <- solve(nt*S + solve(U_0))
  
  E_w <- w_1 * U_1 

  
  both_frame <- W_frm %>% rename(value_post = value) %>%
    mutate(Exp_w = E_w %>% reshape2::melt() %>%
             filter(Var1 <= Var2) %>% select(-Var1, -Var2) %>% 
             rename(value_E = value),
           diag_i_post = mapply(function(i) W[i,i], Var1),
           diag_j_post = mapply(function(j) W[j,j], Var2),
           diag_i_E = mapply(function(i) E_w[i,i], Var1),
           diag_j_E = mapply(function(j) E_w[j,j], Var2),
           rho_post =  -value_post / sqrt(diag_i_post * diag_j_post),
           rho_E = -Exp_w$value_E / sqrt(diag_i_E * diag_j_E)
          ) %>% 
    filter(Var1 < Var2) %>%
    mutate(edge = as.numeric(rho_post / rho_post > 0.5))
    
  return(both_frame$edge)
   
}
