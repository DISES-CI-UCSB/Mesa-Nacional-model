## script: Prioritizr Colombia (matrix)
## Purpose: Set parameters and run prioritization models for Colombia

# ========== SETTING UP =====================================================
## load packages
library(prioritizr)  # modeling package
library(gurobi)      # solver
library(sf)          # vector data
library(terra)       # raster/GIS data
library(tidyverse)   # always
library(Matrix)      # using sparse matrices
library(here)        # easier filepaths
library(purrr)       # run models over list

## Set seed and directories
set.seed(500)

ipt_dir <- here("data/model_inputs")
opt_dir <- here("results")

for (dir in c(ipt_dir, opt_dir)){
  if (!dir.exists(dir)) dir.create(dir)  # Create directories if needed
}; rm(dir)


## Get the list of scenarios
source(here("scripts/utils.R"))
scenarios_df <- filter(scenarios_df, cost == "IHEH2022") # For now, only evaluating one cost

# ========== GET PLANNING UNITS ================================================

## Get list of all non-NA cells (each cell == planning unit)
ids <- cells(template)
n_pus <- length(ids) # number of planning units

## For now, just using IHEH 2022 as sole cost
pus <- readRDS(file.path(ipt_dir, "IHEH_2022.rds"))
pus <- pus[ids, ]
pus[is.na(pus)] <- 0 #shouldn't be any NAs but can use in case


# ========== PRIORITZATION FUNCTION ============================================
## Wrapping all the model building and running inside a function 
## to easily run over a list of scenarios.

prioritizr_model <- function(target, cost, features, includes, model_name, ids, pus) {
  ## Print scenario and time of start
  message("Running scenario: ", model_name, 
          "\nRun start: ", format(Sys.time(), "%H:%M:%S"))
  
  ## Unlist variables
  features <- unlist(features)
  includes <- unlist(includes)
  
  # --------- CONSTRAINTS (LOCK INS/OUTS) -----------------------------
  ## RUNAP is always included
  locked_in <- readRDS(file.path(ipt_dir, "runap.rds"))[ids, ] == 1
  
  ## Add to locked_in matrix depending on scenario
  if ("OMEC" %in% includes) {
    ## Update to make any cell either condition as TRUE
    locked_in <- locked_in | (readRDS(file.path(ipt_dir, "omec.rds"))[ids, ] == 1)
  }
  if ("comunidades" %in% includes) {
    ## Update to make any cell either condition as TRUE
    locked_in <- locked_in | (readRDS(file.path(ipt_dir, "comunidades.rds"))[ids, ] == 1)
  }
  if ("resguardos" %in% includes) {
    ## Update to make any cell either condition as TRUE
    locked_in <- locked_in | (readRDS(file.path(ipt_dir, "resguardos.rds"))[ids, ] == 1)
  }
  

  
  # --------- FEATURES -------------------------------------------------
  # Now read in feature data for the scenario
  ## Start with empty list to add features to
  feature_list <- list()
  
  ## -------- Ecosystems ------------------------------------------
  ## All ecosystems
  if ("ecosystems" %in% features) {
    ## Read in matrix
    ecosys_v <- readRDS(file.path(ipt_dir, "ecosistemas_IAVH_2024.rds"))
    ecosys_v <- ecosys_v[ids, ]
    ## Read in df and filter to only ecosystems not meeting targets
    cons_type <- ifelse("OMEC" %in% includes, "OMEC", "RUNAP") 
    ecosys_df <- read_csv(file.path(ipt_dir, "ecosys_filtered.csv")) %>% 
      filter (targets == target,
              conservation_type == cons_type)
    
    ## Add to list
    col_idx <- which(colnames(ecosys_v) %in% ecosys_df$feature)
    feature_list[["ecosystems"]] <- ecosys_v[ , col_idx]
    ## Remove to keep env memory low
    rm(ecosys_v)
  }
  
  ## Strategic ecosystems
  if ("strategic ecosystems" %in% features) {
    strat_ecos_v <- readRDS(file.path(ipt_dir, "strategic_ecosystems.rds"))
    strat_ecos_v <- strat_ecos_v[ids, ]
    strat_ecos_v[is.na(strat_ecos_v)] <- 0  # fix NAs before converting
    feature_list[["strategic ecosystems"]] <- as(strat_ecos_v, "dgCMatrix")
    rm(strat_ecos_v)
  }
  
  ## Combine if any/both were loaded
  if (length(feature_list) > 0) {
    ecosystems <- do.call(cbind, feature_list)
    ecosystems[is.na(ecosystems)] <- 0  # Remove lingering NAs
    
    ## Transpose matrix and make sparse to match problem format
    ecosys_sparse <- as(t(ecosystems), "sparseMatrix"); rm(ecosystems)
  }
  
  
  ## -------- Species ------------------------------------------
  if ("species" %in% features) {
    mat <- readRDS(file.path(ipt_dir, "biomod_filtered.rds"))
    species_rij <- mat %>% t() %>% as("dgCMatrix"); rm(mat)    # transpose [rows == spp, columns == cell]
    species_rij <- species_rij[, ids]
    
    #If OMEC are in includes, filter by that. Otherwise, just RUNAP
    species_cons_type <- ifelse("OMEC" %in% includes, "OMEC", "RUNAP")  
    
    ## Filter dataframe and then matrix by determined goals
    species_df <- read_csv(file.path(ipt_dir, "biomod_spp_ranges_filtered.csv")) %>%
      filter(targets == target,                        # match target
             conservation_type == species_cons_type,   # match RUNAP/OMEC
             class != "Actinopteri")                   # for now, don't consider fish (Elkin's recommendation) 
    
    row_idx <- setNames(seq_len(nrow(species_rij)), rownames(species_rij))
    idx <- row_idx[species_df$scientific_name]
    idx <- idx[!is.na(idx)]
    species_filtered <- species_rij[idx, ]; rm(species_rij)
  }
  
  
  ## -------- Combine Features ------------------------------------------
  features_mat <- switch(
    paste(c("ecosystems" %in% features | "strategic ecosystems" %in% features,  # if either includes, returns TRUE
            "species" %in% features), collapse = "_"),
    "TRUE_TRUE"  = rbind(ecosys_sparse, species_filtered),  # bind both
    "TRUE_FALSE" = ecosys_sparse,                           # only ecosystems
    "FALSE_TRUE" = species_filtered                         # only species
  )
  
  ## Get number and names of features
  n_features <- as.numeric(features_mat@Dim[1])
  feature_names <- features_mat@Dimnames[[1]]
  
  features_df <- data.frame(
    id = 1:n_features,
    name = feature_names
  )
  

  # --------- SET PROBLEM -------------------------------------------------
  ## Set relative targets for features. For now, these are all the same!
  targets_df <- rep(target/100, n_features)
  
  ## Boundary penalties
  boundaries <- prioritizr::boundary_matrix(template)[ids, ids]
  boundaries <- boundaries/max(boundaries) #scaling issue
  
  
  ## Build problem 
  p <- problem(
    x = pus,
    features = features_df,
    rij_matrix = features_mat) %>% 
    add_min_set_objective() %>% 
    add_relative_targets(targets_df) %>% 
    add_locked_in_constraints(locked_in) %>%
    add_binary_decisions() %>% 
    add_boundary_penalties(penalty = 0.001, data = boundaries) %>% 
    add_gurobi_solver(gap = 0.05, threads = 15, verbose = FALSE)
  
  ## If problem fails presolve check, note it and skip to next
  log_file <- file.path(opt_dir, "failed_scenarios.txt")
  
  s <- tryCatch({
    if (!presolve_check(p))
      stop(paste("Presolve check failed for scenario:", model_name))
    solve(p)
  }, error = function(e) {
    message("Skipping ", model_name, ": ", e$message)
    write(paste(Sys.time(), model_name, e$message, sep = " | "), 
          file = log_file, append = TRUE)
    return(NULL)
  })
  
  ## Exit early if solve failed
  if (is.null(s)) return (NULL)
  
  ## Rasterize solution and save
  s_rast <- rasterize_soln(s, template, locked_in, ids)
  writeRaster(s_rast, 
              file.path(opt_dir, paste0(model_name, ".tif")), 
              overwrite = TRUE) 
  
  # --------- EVALUATE RESULTS ---------------------------------------------
  ## Save target-specific coverage for each scenario
  target_coverage <- eval_target_coverage_summary(p, s) %>%
    mutate(scenario = model_name)
  # filter(absolute_target > 0)
  write_csv(target_coverage, file.path(opt_dir, paste0(model_name, "_summary.csv")))
  
  
  ## Get overview stats and add to running solution df
  freq_tbl <- freq(s_rast)
  cost_summary <- eval_cost_summary(p, s) # Is this even meaningful? 
  eval_summary <- data.frame(
    run = model_name,
    n_total = sum(freq_tbl$count),
    n_new_protection = get_freq(freq_tbl, "Priority area"),
    n_locked_in = get_freq(freq_tbl, "Locked in"),
    cost = cost_summary$cost,
    pct_targets_met = mean(target_coverage$met) * 100
  )
  
  # Append eval_summary row to master CSV 
  csv_path <- file.path(opt_dir, "master_eval_summary.csv")
  
  if (file.exists(csv_path)) {
    summary_df <- read_csv(csv_path, show_col_types = FALSE)
    
    ## Does the run exist? If so, overwrite
    if (any(summary_df$run == eval_summary$run)) {
      summary_df[summary_df$run == eval_summary$run, ] <- eval_summary
    } else {
    ## If not, append to table
      summary_df <- rbind(summary_df, eval_summary)
    }
    
    ## If file hasn't yet been created, then create it
  } else {
    summary_df <- eval_summary
  }
  
  ## Save
  write_csv(summary_df, csv_path)
  
  gc()
}



# ========== RUN PRIORITIZATION ==============================================
## If process stopped part-way, use this code to remove scenarios already completed
# completed <- read_csv(file.path(opt_dir, "master_eval_summary.csv"))
# failed <- file.path(opt_dir, "failed_scenarios.txt")
# failed <- read.table(failed, sep = "|", 
#                          col.names = c("time", "model_name", "error"),
#                          strip.white = TRUE)
# failed_list <- unique(failed$model_name)
# scenarios_df <- scenarios_df %>%
#   filter(!model_name %in% completed$run) %>% 
#   filter(!model_name %in% failed_list)
# rm(completed); rm(failed); rm(failed_list)


## Generate model over list of scenarios
purrr::pmap(scenarios_df, prioritizr_model, ids = ids, pus = pus)

