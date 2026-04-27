#### Get graphs, trace plots, and lag plots for lambda, omega off diag, omega diag

these_sims <- sample(1:399,10)
plot_tib <- tibble()

for(sim in these_sims){
  print(sim)
  load(results_paths[sim])
  W_post_means <- colMeans(bgl_test$W_hist[200:400,])
  D <- 100#chain_res[[1]]$results_tib$D[1]
  diag_entries <- seq_len(D) * (seq_len(D)+1)/2
  random_diag <- sample(diag_entries,1)
  non_zero_ind <-  which(bgl_test$W_true != 0)
  non_diag <- non_zero_ind[!(non_zero_ind %in% diag_entries)]
  random_non_diag <- sample(non_diag,1)
  zero_ind <- which(bgl_test$W_true == 0)
  zero_random_non_diag<- sample(zero_ind,1)
  index <- 1:length(W_post_means)
  entry_type <- ifelse(index %in% diag_entries, "red", 
                       ifelse(index %in% non_diag, "blue", 
                              "green"))

  plot( bgl_test$W_true, W_post_means,main = paste("D =", D), pch = 1, 
       col= entry_type, ylab = "W post", xlab = "W true")
  abline(a=0,b=1)
  legend(
    x = 16, y = max(W_post_means)-3, cex = 0.5,
    legend = c("Zero-Off Diag","Non-Zero Off-Diag",  "Diag" ),
    fill = c("green", "blue", "red")
  )

  for(i in c(1,2)){
    tib <- tibble(
      chain = i,
      D = D,
      diag_true = rep(chain_res[[i]]$W_true[random_diag],1000),
      diag_hist = chain_res[[i]]$W_hist[,random_diag],
      non_diag_true = rep(chain_res[[i]]$W_true[random_non_diag],1000),
      non_diag_hist = chain_res[[i]]$W_hist[,random_non_diag],
      zero_non_diag_true = rep(chain_res[[i]]$W_true[zero_random_non_diag],1000),
      zero_non_diag_hist = chain_res[[i]]$W_hist[,zero_random_non_diag],
      g = list(chain_res[[1]]$g[[1]])
    ) %>% mutate(x = 1:1000)
    plot_tib <- bind_rows(plot_tib,tib) 
  }
}

plot_tib_grouped <- plot_tib %>%
  pivot_longer(
    cols = matches("diag|non_diag|zero_non_diag"),
    names_to = c("location",".value"),
    names_pattern = "(diag|non_diag|zero_non_diag)_(.*)"
  )%>%
  group_by(D, chain, location) %>% 
  mutate(chain = factor(chain),
         Dimension = factor(D),
         ESS = round(coda::effectiveSize(hist),1)) 
  



diag_plot <- ggplot(plot_tib_grouped %>% filter(location == "diag"), 
       aes(x = x, y = hist, col = Dimension)) +
  geom_line(linewidth = 0.05) + 
  geom_hline(aes(yintercept = true, col = Dimension, fill = true), linetype = "dashed") +
  facet_wrap(vars(chain)) + 
  theme(plot.title = element_text(size = 10, hjust = 0.5))  +
  labs(x = "Iteration", y = "Value", 
       title = "Trace Plots for Diagonal Element across 4 different
       sims, 2 chains") 

off_diag_plot <- ggplot(plot_tib_grouped %>% filter(location == "non_diag"), 
       aes(x = x, y = hist, col = Dimension)) +
  geom_line(linewidth = 0.05, alpha = 1) + 
  geom_hline(aes(yintercept = true, col = Dimension), linetype = "dashed") +
  facet_wrap(vars(chain)) + 
  theme(plot.title = element_text(size = 10, hjust = 0.5))  +
  labs(x = "Iteration", y = "Value", 
       title = "Trace Plots for Non-Zero Off-Diagonal Element across 4 different
       sims, 2 chains")

zero_off_diag_plot <- ggplot(plot_tib_grouped %>% filter(location == "zero_non_diag"), 
       aes(x = x, y = hist, col = Dimension)) +
  geom_line(linewidth = 0.05,  alpha= 1) + 
  geom_hline(aes(yintercept = true, col = Dimension), linetype = "dashed") +
  facet_wrap(vars(chain)) + 
  theme(plot.title = element_text(size = 10, hjust = 0.5))  +
  labs(x = "Iteration", y = "Value", 
      title = "Trace Plots for Zero Off-Diagonal Element across 4 different
       sims, 2 chains")

ggsave(
  here::here("figures", "diag_trace.png"),
  diag_plot
)

ggsave(
  here::here("figures", "off_diag_trace.png"),
  off_diag_plot
)

ggsave(
  here::here("figures", "zero_off_diag_trace.png"),
  zero_off_diag_plot
)

graph_tib <- plot_tib %>% select(D,g) %>%
  group_by(D,g) %>% distinct(D,g)

for(i in 1:nrow(graph_tib)){
  print(i)
  g <- graph_tib$g[[i]]
  print(g)
  V(g)$size <- degree(g)+1
  V(g)$color <- rainbow(100)[cut(degree(g), breaks = 100)]
  png(file.path("figures",paste0("graph_plot_D_",graph_tib$D[[i]],".png")))
  plot(g,
        vertex.label = NA,
        vertex.col = "blue",
        main = paste("D =",graph_tib$D[[i]]))
  dev.off()
}


  zero_non_diag_lag <- forecast::ggAcf(
  plot_tib_grouped %>% filter(location == "zero_non_diag", chain == 1) %>% 
    select(hist) %>% ungroup(),
  type = "correlation") + 
  labs(title = "Mu 2 Chain 2 lag plot") + ylab("Correlation")

### Omega coverage ###

omega_cov <- res_tib %>%
  group_by(D) %>%
  summarise(Omega_Coverage = mean(omega_coverage)) %>%
  mutate(Dimension = factor(D))

omega_cov_plot <- ggplot(data = omega_cov, aes(y = Omega_Coverage, x = Dimension, 
                             fill = Dimension)) +
  geom_bar(position = "dodge", stat = "identity") +
  coord_cartesian(ylim=c(0.8,1)) + 
  geom_hline(yintercept = 0.95, linetype = "dashed") + 
  labs(title = "Coverage of Omega by Dimension", y = "Coverage")

mcc_plot <- ggplot(data = summarize_res, aes(fill = Dimension, y = MCC, x = Method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = MCC_LL, ymax = MCC_UL), 
                position = position_dodge(0.9), width = 0.5, color = "grey25")+
  labs(title = "MCC by method", y = "Coverage")

spec_plot <- ggplot(data = summarize_res, aes(fill = Dimension, y = Specificity, x = Method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = Spec_LL, ymax = Spec_UL), 
                position = position_dodge(0.9), width = 0.5, color = "grey25") +
  labs(title = "Specificity by method", y = "Specificity")

sens_plot <- ggplot(data = summarize_res, aes(fill = Dimension, y = Sensitivity, x = Method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = Sens_LL, ymax = Sens_UL), 
                position = position_dodge(0.9), width = 0.5, color = "grey25") +
  labs(title = "Sensitivity by method", y = "Sensitivity")

stein_plot <- ggplot(data = summarize_res, aes(fill = Dimension, y = `Stein's Loss`, x = Method)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = Stein_LL, ymax = Stein_UL), 
                position = position_dodge(0.9), width = 0.5, color = "grey25") +
  labs(title = "Stein's loss by method", y = "Stein's loss")

comp_plot <- ggplot(data = summarize_res, 
       aes(y = `Log Comp Time (s)`, x = D, col = Method)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = Comp_LL, ymax = Comp_UL), 
                position = position_dodge(0.9), width = 0.5, color = "grey25") +
  labs(title = "Log Comp Time by method", y = "Log Comp Time", x = "Dimension")

lambda_plot <- ggplot(data = summarize_res, 
       aes(y = Lambda, x = D, col = Method)) +
  geom_line() + geom_point() +
  labs(title = "Lambda by Method", y = "Lambda", x = "Dimension")

ggsave(
  here::here("figures", "omega_cov.png"),
  omega_cov_plot
)

ggsave(
  here::here("figures", "mcc.png"),
  mcc_plot
)

ggsave(
  here::here("figures", "sens.png"),
  sens_plot
)

ggsave(
  here::here("figures", "spec.png"),
  spec_plot
)

ggsave(
  here::here("figures", "stein.png"),
  stein_plot
)

ggsave(
  here::here("figures", "comp.png"),
  comp_plot
)

ggsave(
  here::here("figures", "lambda.png"),
  lambda_plot
)