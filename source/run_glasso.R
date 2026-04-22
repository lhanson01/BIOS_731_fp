library(CVglasso)

do_glasso <- function(Y){
  results <- CVglasso(
    X = t(Y),
    K = 10
  )
  
  return(list(Omega = results$Omega,
              lambda = results$Tuning[2]))
}
