suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(CVglasso))


source(here::here("source","bayesian_graphical_lasso.R"))
source(here::here("source","simulate_graph_data.R"))
source(here::here("source","run_glasso.R"))
source(here::here("source","evaluate_edge.R"))

combined_ctrls_txd <- readRDS(
  file = here::here("output", "combined_control_TxD.rds")
)

combined_trt_txd <- readRDS(
  file = here::here("output", "combined_trt_TxD.rds")
)

n_iter <- 2000
n_burn <- 400
n_chain <- 2

Y_trt_list <- list(Y=as.matrix(combined_trt_txd))
Y_ctrl_list <- list(Y=as.matrix(combined_trt_txd))

ctrl_chain_res <- vector(mode="list", length = n_chain)
trt_chain_res <- vector(mode="list", length = n_chain)

for(i in 1:n_chain){
  ctrl_chain_res[[i]] <- bayesian_graphical_lasso(
    data = Y_trt_list,
    n_iter = n_iter,
    n_burn = n_iter,
    is_sim = FALSE
  )
  
  trt_chain_res[[i]] <- bayesian_graphical_lasso(
    data = Y_trt_list,
    n_iter = n_iter,
    n_burn = n_iter,
    is_sim = FALSE
  )
}

saveRDS(ctrl_chain_res,
        file = "output/ctrl_chain_res.rds")

saveRDS(trt_chain_res,
        file = "output/trt_chain_res.rds")


