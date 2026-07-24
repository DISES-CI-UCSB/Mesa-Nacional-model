## script: National Marine Model
## Purpose: Set parameters and run marine prioritization national models for Colombia

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

ipt_dir <- here("data/model_inputs/national")
opt_dir <- here("results/national/marine")

for (dir in c(ipt_dir, opt_dir)){
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE) 
}; rm(dir)


# ========== PRIORITZATION FUNCTION ============================================
## Wrapping all the model building and running inside a function 
## to easily run over a list of scenarios.

#' @param target Numeric. Target percentage (0-100) for all the features. 
#' @param includes Character vector. Which layers should be "locked-in" to the 
#'   solution (e.g. "RUNAP", "OMEC").
#' @param features Character vector. Which layers are include as a feature to
#'   be evaluated in the problem (i.e. "ecosystems", "mangroves")
#' @param model_name Character. Unique identifier for this scenario, used 
#'   for output file names and logging.
#' @param skip_presolve Logical. If TRUE, skip the presolve check and 
#'   attempt to solve regardless. Default FALSE.
#' @param force_s Logical. If TRUE, force gurobi to return a solution even 
#'   if the presolve/solve process raises non-fatal warnings. Passed to 
#'   solve(p, force = force_s). Default TRUE
#'
#' @return NULL (invisibly). Writes solution raster and summary CSVs to 
#'   opt_dir for each solution. Also creates and appends log of failed scenarios.
marine_model <- function(target, includes, features, model_name, 
                         skip_presolve = FALSE, force_s = TRUE) {
  ## Print scenario and time of start
  message("Running scenario: ", model_name, 
          "\nRun start: ", format(Sys.time(), "%H:%M:%S"))
  
  # --------- PLANNING UNITS ------------------------------------------
  ## For now, only one cost layer
  pus <- readRDS(file.path(ipt_dir, "huella_humana_marina.rds"))
  
  ## Get list of all non-NA cells (each cell == planning unit)
  ids <- cells(template_mar)
  n_pus <- length(ids) # number of planning units
  
  pus <- pus[ids, ]
  pus[is.na(pus)] <- 0 #shouldn't be any NAs but can use in case
  
  # --------- CONSTRAINTS (LOCK INS/OUTS) -----------------------------
  ## Unlist variable
  includes <- unlist(includes)
  
  ## RUNAP is always included
  locked_in <- readRDS(file.path(ipt_dir, "runap_marinos.rds"))[ids, ] == 1
  
  ## Add to locked_in matrix depending on scenario
  if ("OMEC" %in% includes) {
    ## Update to make any cell either condition as TRUE
    locked_in <- locked_in | (readRDS(file.path(ipt_dir, "omec_marinos.rds"))[ids, ] == 1)
  }
  
  
  # --------- FEATURES ----------------------------------------------
  # Add features to empty lists if they are evaluated
  # NOTE: currently both features (ecosystems and mangroves) always evaluated.
  # This may change in future?
  features_list <- list()
  
  features <- unlist(features)
  
  ## -------- Ecosystems ------------------------------------------
  if ("ecosystems" %in% features) {
    ## Read in matrix
    ecosys_v <- readRDS(file.path(ipt_dir, "ecosistemas_marinos.rds"))
    ecosys_v <- t(ecosys_v) %>% as("dgCMatrix") # transpose [rows == ecosystem, columns == cell]
    ecosys_v <- ecosys_v[, ids] # only keep cells in PUs
    
    ## Read in df and filter to only ecosystems not meeting targets
    cons_type <- ifelse("OMEC" %in% includes, "RUNAP_OMEC", "RUNAP") 
    ecosys_df <- read_csv(file.path(ipt_dir, "marine_ecosys_filtered.csv"), show_col_types = FALSE) %>% 
      filter(targets == target,
             conservation_type == cons_type)
    
    ## Only add ecosystems that need to be evaluated to features list
    row_idx <- which(rownames(ecosys_v) %in% ecosys_df$feature)
    features_list[["ecosystems"]] <- ecosys_v[row_idx, ]
    
    ## Remove matrix to keep env memory low
    rm(ecosys_v)
  }
  
  
  ## -------- Mangroves ------------------------------------------
  # NOTE: this may expand to "strategic ecosystems" if more are added in future
  if ("mangroves" %in% features) {
    manglares_v <- readRDS(file.path(ipt_dir, "manglares.rds"))
    manglares_v <- t(manglares_v) %>% as("dgCMatrix") # transpose [rows == ecosystem, columns == cell]
    manglares_v <- manglares_v[, ids]
    manglares_v[is.na(manglares_v)] <- 0  # fix NAs before converting
    
    ## Add to features list
    features_list[["Manglares"]] <- manglares_v
    rm(manglares_v)
  }
  
  
  ## -------- Combine Features ------------------------------------------
  ## Bind everything into one matrix
  features_mat <- do.call(rbind, features_list)
  
  ## Get number and names of features
  n_features <- as.numeric(features_mat@Dim[1])
  feature_names <- features_mat@Dimnames[[1]]
  
  ## Create list of targets for prioritizr model
  targets_df <- rep((target/100), n_features)
  
  features_df <- data.frame(
    id = 1:n_features,
    name = feature_names,
    target = targets_df
  )
  
  # --------- SET PROBLEM -------------------------------------------------
  ## Boundary penalties
  boundaries <- prioritizr::boundary_matrix(template_mar)[ids, ids]
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
  
  
  ## Get coverage summary & save
  target_coverage <- eval_target_coverage_summary(p, s) %>%
    mutate(scenario = model_name,           # Add the scenario info
           evaluated = "prioritizr_model")  # These features explicitly evaluated in model
  
  write_csv(target_coverage, file.path(opt_dir, paste0(model_name, "_summary.csv")))
  
  ## Rasterize solution and save
  s_rast <- rasterize_soln(s, template_mar, locked_in, ids)
  writeRaster(s_rast,
              file.path(opt_dir, paste0(model_name, ".tif")),
              overwrite = TRUE)
  
  
  # --------- POST-HOC EVALUATION ---------------------------------------------
  # Because some ecosystems (that already met baseline target) 
  # were not included, need to get their coverage stats at a national level as well
  message("Running post-hoc evaluation for scenario: ", model_name)
  
  ## ------- Ecosystems ----------------------------------------
  eco_coverage <- NULL
  
  if ("ecosystems" %in% features) {  # For now, ecosystems always included. But make it flexible in case
    ecosys_mat <- readRDS(file.path(ipt_dir, "ecosistemas_marinos.rds"))[ids, ]
    ecosys_totals <- colSums(ecosys_mat)
    ecosys_names <- colnames(ecosys_mat)
    
    evaluated <- target_coverage$feature
    unevaluated <- setdiff(ecosys_names, evaluated)
    if (length(unevaluated) == 0) return(NULL)
    
    totals <- ecosys_totals[unevaluated]
    in_soln <- colSums(ecosys_mat[s == 1, unevaluated, drop = FALSE])
    abs_target <- totals * (target / 100)
    rel_held <- in_soln / totals
    
    eco_coverage <- data.frame(
      feature = unevaluated,
      met = rel_held >= (target / 100),
      total_amount = totals,
      absolute_target = abs_target,
      absolute_held = in_soln,
      absolute_shortfall = pmax(0, abs_target - in_soln),
      relative_target = target / 100,
      relative_held = rel_held,
      relative_shortfall = pmax(0, (target / 100) - rel_held),
      scenario = model_name,
      evaluated = "post-hoc"
    )
    
    rm(ecosys_mat); gc()
    
  } # END ecosystem post-hoc
  
  
  
  ## ------- Combine and save -------------------------------------
  ## Put all the post-hoc evaluations in one dataframe
  post_hoc_coverage <- eco_coverage # ONLY ONE FOR NOW
  rownames(post_hoc_coverage) <- NULL
  
  ## Combine with solution target coverages, then save
  target_coverage_full <- rbind(target_coverage, post_hoc_coverage)
  
  write_csv(target_coverage_full, file.path(opt_dir, paste0(model_name, "_summary.csv")))
  
  
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
#   scenarios_mar_df <- scenarios_mar_df %>%
#     filter(!model_name %in% completed$scenario)
#   rm(completed)
# }
# 
# if (file.exists(failed_list)) {
#   failed <- read.table(failed_list, sep = "|",
#                        col.names = c("time", "model_name", "error"),
#                        strip.white = TRUE)
#   failed_list <- unique(failed$model_name)
#   scenarios_mar_df <- scenarios_mar_df %>%
#     filter(!model_name %in% failed_list)
#   rm(failed); rm(failed_list)
# }

## Generate model over list of scenarios
purrr::pmap(scenarios_mar_df, marine_model, force_s = TRUE)



