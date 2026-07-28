## script: SIRAP Orinoquia Model
## Purpose: Set parameters and run prioritization model for Orinoquia region in Colombia

# ========== SETTING UP =======================================================
## Get functions and data
source("scripts/utils.R")

## Load/install packages 
pacman::p_load(  # automatically installs packages if needed
  prioritizr,    # modeling package
  gurobi,        # solver
  Matrix)        # Matrices

## Set seed and directories
set.seed(500)

ipt_dir <- here("data/model_inputs/sirap/orinoquia")
opt_dir <- here("results/sirap/orinoquia")

if (!dir.exists(opt_dir)) dir.create(opt_dir, recursive = TRUE)

## Use specific template for model
template <- template_ori

# ========== PRIORITZATION FUNCTION ============================================
# Wrapping all the model building and running inside a function to easily run over a list of scenarios.

#' @param strat_ecos_target Numeric. Target percentage (0-100) for strategic 
#'   ecosystems; 0 to exclude.
#' @param cong_target Numeric. Target percentage (0-100) for congriales; 
#'   0 to exclude.
#' @param sab_target Numeric. Target percentage (0-100) for savannas within
#'   Orinoquia; 0 to exclude. 
#' @param includes Character vector. Which layers should be "locked-in" to the 
#'   solution (e.g. "RUNAP", "OMEC").
#' @param cost Character. Which cost data to use ("IHEH2022" or 
#'   "IHEH2030"). Also determines available planning units.
#' @param model_name Character. Unique identifier for this scenario, used 
#'   for output file names and logging.
#' @param skip_presolve Logical. If TRUE, skip the presolve check and 
#'   attempt to solve regardless. Default FALSE.
#' @param force_s Logical. If TRUE, force gurobi to return a solution even 
#'   if the presolve/solve process raises non-fatal warnings. Passed to 
#'   solve(p, force = force_s). Default FALSE.
#'
#' @return NULL (invisibly). Writes solution raster and summary CSVs to 
#'   opt_dir for each solution. Also creates and appends log of failed scenarios.

eje_model <- function(strat_ecos_target, cong_target, sab_target, includes, cost, 
                      model_name, skip_presolve = FALSE, force_s = FALSE) { 
  ## Print scenario and time of start
  message("Running scenario: ", model_name, 
          "\nRun start: ", format(Sys.time(), "%H:%M:%S"))
  
  # --------- PLANNING UNITS ------------------------------------------
  if (cost == "IHEH2022") {
    pus <- readRDS(file.path(ipt_dir, "IHEH_orinoquia_2022.rds"))
  } else if (cost == "IHEH2030") {
    pus <- readRDS(file.path(ipt_dir, "IHEH_orinoquia_2030.rds"))
  }
  
  ## Get list of all non-NA cells (each cell == planning unit)
  ids <- cells(template)
  n_pus <- length(ids)
  
  pus <- pus[ids, ]
  pus[is.na(pus)] 
  
  
  # --------- CONSTRAINTS (LOCK INS/OUTS) -----------------------------
  ## Unlist variable
  includes <- unlist(includes)
  
  ## RUNAP is always included
  locked_in <- readRDS(file.path(ipt_dir, "runap_orinoquia.rds"))[ids, ] == 1
  
  ## Add to locked_in matrix depending on scenario
  if ("OMEC" %in% includes) {
    ## Update to make any cell either condition as TRUE
    locked_in <- locked_in | (readRDS(file.path(ipt_dir, "omec_orinoquia.rds"))[ids, ] == 1)
  }
  
  
  # --------- FEATURES & TARGETS ----------------------------------------------
  # Add features (and their targets) to empty lists if they are evaluated
  # in the specific scenario
  features_list <- list()
  targets_list <- list()
  
  ## -------- Strategic Ecosystems --------------------------------
  ## Read in matrix and add to features list if evaluated
  if (strat_ecos_target != 0) {
    ## Read in matrix
    strat_ecos_v <- readRDS(file.path(ipt_dir, "ecosistemas_estrategicos_orinoquia.rds"))
    strat_ecos_v <- t(strat_ecos_v) %>% as("dgCMatrix") # transpose [rows == ecosystem, columns == cell]
    strat_ecos_v <- strat_ecos_v[, ids]
    strat_ecos_v[is.na(strat_ecos_v)] <- 0  # fix NAs before converting
    
    ## Add to features and targets lists
    features_list[["strategic ecosystems"]] <- strat_ecos_v
    targets_list[["strategic ecosystems"]] <- rep(strat_ecos_target/100, nrow(strat_ecos_v))
    rm(strat_ecos_v)
  }
  
  
  ## -------- Congriales ----------------------------------
  if (cong_target != 0) {
    cong_v <- readRDS(file.path(ipt_dir, "congriales.rds"))
    cong_v <- t(cong_v) %>% as("dgCMatrix")
    cong_v <- cong_v[, ids]
    cong_v[is.na(cong_v)] <- 0
    
    ## Add to features and targets list
    features_list[["congriales"]] <- cong_v
    targets_list[["congriales"]] <- (cong_target/100)
    rm(cong_v)
  }
  
  
  ## -------- Savannas ----------------------------------
  if (sab_target != 0) {
    sab_v <- readRDS(file.path(ipt_dir, "sabana_orinoquia.rds"))
    sab_v <- t(sab_v) %>% as("dgCMatrix")
    sab_v <- sab_v[, ids]
    sab_v[is.na(sab_v)] <- 0
    
    ## Add to features and targets list
    features_list[["savannas"]] <- sab_v
    targets_list[["savannas"]] <- (sab_target/100)
    rm(sab_v)
  }
  
  
  ## -------- Combine Features ------------------------------------------
  ## Bind everything into one matrix
  features_mat <- do.call(rbind, features_list)
  targets_df <- round(unlist(targets_list, use.names = FALSE), 4)
  
  ## Get number and names of features
  n_features <- as.numeric(features_mat@Dim[1])
  feature_names <- features_mat@Dimnames[[1]]
  
  features_df <- data.frame(
    id = 1:n_features,
    name = feature_names,
    target = targets_df
  )

  # --------- SET PROBLEM -------------------------------------------------
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
    add_gurobi_solver(gap = 0.01, threads = 10)
  
  ## If problem fails presolve check, note it and skip to next
  log_file <- file.path(opt_dir, "failed_scenarios.txt")
  
  s <- tryCatch({
    if (!skip_presolve && !presolve_check(p))
      stop(paste("Presolve check failed for scenario:", model_name))
    solve(p, force = force_s)
  }, error = function(e) {
    message("Skipping ", model_name, ": ", e$message)
    write(paste(Sys.time(), model_name, e$message, sep = " | "), 
          file = log_file, append = TRUE)
    return(NULL)
  })
  
  ## Exit early if solve failed
  if (is.null(s)) return (NULL)
  
  
  ## ------- Summary statistics -------------------------------------
  ## Get coverage summary & save
  target_coverage <- eval_target_coverage_summary(p, s) %>%
    mutate(scenario = model_name,           # Add the scenario info
           evaluated = "prioritizr_model")  # These features explicitly evaluated in model
  
  write_csv(target_coverage, file.path(opt_dir, paste0(model_name, "_summary.csv")))
  
  ## Rasterize solution and save
  s_rast <- rasterize_soln(s, template, locked_in, ids)
  writeRaster(s_rast,
              file.path(opt_dir, paste0(model_name, ".tif")),
              overwrite = TRUE)
  
  
  ## Get overview stats and add to running list of solutions
  freq_tbl <- freq(s_rast)
  cost_summary <- eval_cost_summary(p, s) # Is this even meaningful?
  eval_summary <- data.frame(
    scenario = model_name,
    n_total = sum(freq_tbl$count),
    n_new_protection = get_freq(freq_tbl, "Priority area"),
    n_locked_in = get_freq(freq_tbl, "Locked in"),
    cost = cost_summary$cost,
    pct_targets_met = mean(target_coverage_full$met, na.rm = TRUE) * 100
  )
  
  # Append eval_summary row to master CSV
  csv_path <- file.path(opt_dir, "master_eval_summary.csv")
  
  if (file.exists(csv_path)) {
    summary_df <- read_csv(csv_path, show_col_types = FALSE)
    
    if (any(summary_df$scenario == eval_summary$scenario)) {
      ## Does the scenario exist? If so, overwrite
      summary_df[summary_df$scenario == eval_summary$scenario, ] <- eval_summary
      
    } else {
      ## If not, append to table
      summary_df <- rbind(summary_df, eval_summary)
    }
    
  } else {
    summary_df <- eval_summary   # If file hasn't yet been created, then create it
  }
  
  ## Save master summary
  write_csv(summary_df, csv_path)
  gc()
  
  
} # END PRIORITIZR FUNCTION


# ========== RUN PRIORITIZATION ==============================================
# ## If process stopped part-way, remove scenarios already completed or permanently failed
# completed_list <- file.path(opt_dir, "master_eval_summary.csv")
# failed_list    <- file.path(opt_dir, "failed_scenarios.txt")
# 
# if (file.exists(completed_list)) {
#   completed <- read_csv(completed_list)
#   scenarios_ec_df <- scenarios_ec_df %>%
#     filter(!model_name %in% completed$scenario)
#   rm(completed)
# }
# 
# if (file.exists(failed_list)) {
#   failed <- read.table(failed_list, sep = "|",
#                        col.names = c("time", "model_name", "error"),
#                        strip.white = TRUE)
#   failed_list <- unique(failed$model_name)
#   scenarios_ec_df <- scenarios_ec_df %>%
#     filter(!model_name %in% failed_list)
#   rm(failed); rm(failed_list)
# }

## Generate model over list of scenarios
purrr::pmap(scenarios_ec_df, eje_model, skip_presolve = TRUE, force_s = TRUE)

