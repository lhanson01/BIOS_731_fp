  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(here))
  suppressPackageStartupMessages(library(igraph))
  suppressPackageStartupMessages(library(brms))
  suppressPackageStartupMessages(library(CVglasso))
  
  source(here::here("source","bayesian_graphical_lasso.R"))
  source(here::here("source","simulate_graph_data.R"))
  source(here::here("source","run_glasso.R"))
  source(here::here("source","evaluate_edge.R"))
  
  
  D_scen <- c(25,50,100,300)
  n_sim <- 100
  jobs_table <- expand.grid(simID = 1:n_sim, D = D_scen)
  
  job <- as.numeric(commandArgs(trailingOnly = TRUE))
  seed <- floor(runif(nrow(jobs_table), 1, 10000))
  n_iter <- 100
  n_burn <- 20
  n_chain <- 2
  n_t <- 200
  p_mean <- 2
  
  D <- jobs_table$D[job]
  simID <- jobs_table$simID[job]
  
  #generate data
  set.seed(seed[job])
  data <- sim_graph_ts(D = D,
                       nt = n_t,
                       p_mean = p_mean,
                       )
  
  date <- gsub("-", "", Sys.Date())
  results_path <- file.path(here::here("results"))
  if(!file.exists(results_path)){
    dir.create(results_path, Date)
  }
  
  filename <- file.path(results_path, 
                        paste0("D_", D,"_sim_", simID, ".rds")
                        )
  
  use_glasso_init <- c(TRUE, FALSE)
  if(!file.exists(filename)){
  
    chain_res <- vector(mode = "list", length = n_chain)
    i = 1
    for(chain in 1:n_chain){
    gl_init <- use_glasso_init[i]  
      chain_res[[chain]] <- bayesian_graphical_lasso(
        data = data,
        n_iter = n_iter,
        n_burn = n_burn,
        simID = simID,
        gl_init = gl_init
      )
    i = i+1
    }
    
    
    save(chain_res,
         file = filename)
    
    print(paste0("D_", D,"_sim_", simID, "_done"))
  
  } else {
    print("File already exists")
  }
  
  
  
  
