W_post <- readRDS(
  file = "output/small_run_W.rds")
tau_post <- readRDS(
  file = "output/small_run_tau.rds")

ymax <- max(W) -10
ymin <- min(W)
plot(W_post[,1,2], type = "l", ylim = c(ymin,ymax))
for(i in 1:100){
  lines(W_post[,sample(1:50, size = 1),sample(1:50, size = 1)], col = i)
}

plot(tau_post[,4,5], type = "l")
for(i in 5:100){
  lines(tau_post[,sample(1:50),sample(1:50)], col = i)
}

W_final <- apply(W_hist, MARGIN = c(2,3), FUN = mean)

plot(W_hist[,3], type = "l")
lines(W_true_vec[3])
