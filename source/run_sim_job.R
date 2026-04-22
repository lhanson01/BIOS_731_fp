suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))

source(here::here("source","bayesian_graphical_lasso.R"))
source(here::here("source","simulate_graph_data.R"))
source(here::here("source","run_glasso.R"))
source(here::here("source","evaluate_edge.R"))


D_scen <- c(25,50,100,150)
n_sim <- 100
jobs_table <- expand.grid(simID = 1:n_sim, D = D_scen)

job <- as.numeric(commandArgs(trailingOnly = TRUE))
seed <- floor(runif(n_sim, 1, 10000))
n_iter <- 1000
n_burn <- 200
n_chain <- 2
n_t <- 200
p_mean <- 2

D <- jobs_table$D[job]
simID <- jobs_table$simID[job]

#generate data
data <- sim_graph_ts(D = D,
                     nt = n_t,
                     p_mean = p_mean,
                     )

chain_res <- vector(mode = "list", length = n_chain)
for(chain in 1:n_chain){
  chain_res[[chain]] <- bayesian_graphical_lasso(
    data = data,
    n_iter = n_iter,
    n_burn = n_burn,
    simID = simID
  )
}

date <- gsub("-", "", Sys.Date())
results_path <- file.path(here::here("results"))
if(!file.exists(results_path)){
  dir.create(results_path, Date)
}

filename <- paste0(results_path, "/",
                   "D_", D,"_sim_", simID, ".RDA")

save(chain_res,
     file = filename)






