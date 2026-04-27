extract_true_edge <- function(g){
  adj <- as_adjacency_matrix(g)
  A <- as.matrix(adj)
  labels <- A[upper.tri(A, diag = TRUE)]
  return(labels)
}


results_paths <- file.path(here::here("results"), 
                           list.files("results"))
res_list <- vector(mode="list", length = length(results_paths))
for(i in 1:length(results_paths)){
  print(i)
  load(results_paths[i])
  g <- chain_res[[1]]$g[[1]]
  chain1 <- chain_res[[1]]$results_tib %>% 
    nest(gibbs_labels = labels_gibbs, glasso_labels = labels_glasso) %>%
    mutate(true_edge = list(extract_true_edge(g))) 
    
  chain2 <- chain_res[[2]]$results_tib %>% 
    nest(gibbs_labels = labels_gibbs, glasso_labels = labels_glasso) %>%
    mutate(true_edge = list(extract_true_edge(g))) 
  
  res_list[[i]] <- bind_rows(chain1, chain2)
}



res_tib <- bind_rows(res_list)

## fix glasso metrics
fix_glasso <- res_tib %>%
  mutate(
    mcc_glasso = as.numeric(pmap(
      list(glasso_labels,
           true_edge),
      function(est, true){
        est <- as.numeric(unlist(est))
        true <- as.numeric(unlist(true))
        as.numeric(eval_class(est,true)$mcc)
      }
    )),
    sens_glasso = as.numeric(pmap(
      list(glasso_labels,
           true_edge),
      function(est, true){
        est <- as.numeric(unlist(est))
        true <- as.numeric(unlist(true))
        as.numeric(eval_class(est,true)$sens)
      }
    )),
    spec_glasso = as.numeric(pmap(
      list(glasso_labels,
           true_edge),
      function(est, true){
        est <- as.numeric(unlist(est))
        true <- as.numeric(unlist(true))
        as.numeric(eval_class(est,true)$spec)
      }
  ))
)

##

summarize_res <- fix_glasso %>%
  select(-c(gibbs_labels,glasso_labels))%>%
  pivot_longer(
    cols = matches("gibbs|glasso"),
    names_to = c(".value", "Method"),
    names_pattern = "(.*)_(gibbs|glasso)"
  ) %>%
  group_by(D, Method) %>%
  summarise(
    `Stein's Loss` = mean(stein),
    Stein_SE = sd(stein) / sqrt(n()),
    Stein_LL = `Stein's Loss` - 1.96*Stein_SE,
    Stein_UL = `Stein's Loss` + 1.96*Stein_SE,
    MCC = mean(mcc),
    MCC_SE = sd(mcc)/ sqrt(n()),
    MCC_LL = MCC - 1.96*MCC_SE,
    MCC_UL = MCC + 1.96*MCC_SE,
    Specificity = mean(spec),
    Sensitivity = mean(sens),
    `Log Comp Time (s)` = log(mean(comp_time)),
    Comp_SE = sd(comp_time)/ sqrt(n()),
    Comp_LL = log(mean(comp_time) - 1.96*Comp_SE),
    Comp_UL = log(mean(comp_time) + 1.96*Comp_SE),
    Spec_SE = sd(spec)/ sqrt(n()),
    Sens_SE = sd(sens)/ sqrt(n()),
    Spec_LL = Specificity - 1.96*Spec_SE,
    Spec_UL = Specificity + 1.96*Spec_SE,
    Sens_LL = Sensitivity - 1.96*Sens_SE,
    Sens_UL = Sensitivity + 1.96*Sens_SE,
    Lambda = mean(lambda)
  ) %>%
  mutate(
    Method = factor(Method),
    Dimension = factor(D)
  )






