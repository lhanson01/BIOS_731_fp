suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(CVglasso))


source(here::here("source","bayesian_graphical_lasso.R"))
source(here::here("source","simulate_graph_data.R"))
source(here::here("source","run_glasso.R"))
source(here::here("source","evaluate_edge.R"))

combined_ctrls_txd_full <- readRDS(
  file = here::here("output", "combined_control_TxD.rds")
)

combined_trt_txd_full <- readRDS(
  file = here::here("output", "combined_trt_TxD.rds")
)

combined_ctrls_txd <- as.matrix(combined_ctrls_txd_full)
combined_trt_txd <- as.matrix(combined_trt_txd_full)


n_iter <- 1000
n_burn <- 100
n_chain <- 2
job <- as.numeric(commandArgs(trailingOnly = TRUE))

job_list <- expand.grid(chain = 1:n_chain, group = c("trt","ctrl"))

data_list <- list(trt = list(Y=combined_ctrls_txd),
                  ctrl = list(Y=combined_trt_txd)
  )

group <- job_list$group[job]
chain <- job_list$chain[job]
data <- data_list[[group]]


results <- bayesian_graphical_lasso(
    data = data,
    n_iter = n_iter,
    n_burn = n_burn,
    is_sim = FALSE
  )

file_name <- paste0("real_data_group_", group, "_chain_",chain, ".rds")
file_path <- here::here("output", file_name)

saveRDS(results,
        file = file_path)




