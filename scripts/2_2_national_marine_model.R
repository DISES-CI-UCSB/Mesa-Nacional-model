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

if (!dir.exists(opt_dir)) dir.create(opt_dir, recursive = TRUE)

## Use specific template for model
template <- template_mar


# ========== PRIORITZATION FUNCTION ============================================
## Wrapping all the model building and running inside a function 
## to easily run over a list of scenarios.

#' @param ecos_target Numeric. Target percentage (30 or 50) for ecosystem
#'   representation; some ecosystems will get specific targets added later
#' @param strat_ecos_target Numeric. Target percentage (30 or 50) for strategic 
#'   ecosystems; currently only mangroves evaluated.
#' @param includes Character vector. Which layers should be "locked-in" to the 
#'   solution (e.g. "RUNAP", "OMEC").
#' @param cost Character. Which cost data to use. Currently only one layer used. 
#'   This also determines available planning units.
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
marine_model <- function(ecos_target, strat_ecos_target, includes, cost, model_name, 
                         skip_presolve = FALSE, force_s = TRUE) {
  ## Print scenario and time of start
  message("Running scenario: ", model_name, 
          "\nRun start: ", format(Sys.time(), "%H:%M:%S"))
  
  # --------- PLANNING UNITS ------------------------------------------
  ## Read in correct cost layer
  ## **NOTE: only one cost layer now. Adding this for future flexibility**
  if (cost == "HHM") {
    pus <- readRDS(file.path(ipt_dir, "huella_humana_marina.rds"))
  } else {
    message("Please select an appropriate cost layer.")
  }
  
  ## Get list of all non-NA cells (each cell == planning unit)
  ids <- cells(template)
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
  # Add features (and their targets) to empty lists.
  # If they are evaluated in the specific scenario, then it's in these lists
  features_list <- list()
  targets_list <- list()
  
  ## -------- Ecosystems ------------------------------------------
  if (ecos_target != 0) {
    ## First read in ecosystem IDs
    mar_df <- read_csv(file.path(geo_dir, "national", "ecosistemas_IDs_marinos.csv"))
    
    ## Read in ecosystem matrix
    ecosys_v <- readRDS(file.path(ipt_dir, "ecosistemas_marinos.rds"))
    ecosys_v <- t(ecosys_v) %>% as("dgCMatrix") # transpose [rows == ecosystem, columns == cell]
    ecosys_v <- ecosys_v[, ids] # only keep cells in PUs
    
    ## Read in df and filter to only ecosystems not meeting targets
    cons_type <- ifelse("OMEC" %in% includes, "RUNAP_OMEC", "RUNAP") 
    ecosys_df <- read_csv(file.path(ipt_dir, "marine_ecosys_filtered.csv"), show_col_types = FALSE) %>% 
      filter(targets == ecos_target,
             conservation_type == cons_type)
    
    ## Only add ecosystems that need to be evaluated to features list
    row_idx <- which(rownames(ecosys_v) %in% ecosys_df$feature)
    features_list[["ecosystems"]] <- ecosys_v[row_idx, ]
    
    ## Add unique ecosystem targets to list
    targets_list[["ecosystems"]] <- 
      ## Pull the targets from the correct variable (e.g. 30 or 50)
      ecosys_df[[sprintf("umbral_del_%s_percent", ecos_target)]] / 100 #Convert to percentage
  
    
    ## Remove objects to keep env memory low
    rm(ecosys_v, ecosys_df)
  }
  
  
  ## -------- Strategic Ecosystems --------------------------------------
  # NOTE: this only contains mangroves now. Will need to be updated (along with `1_data_prep.R`) once more are added
  if (strat_ecos_target != 0) {
    manglares_v <- readRDS(file.path(ipt_dir, "manglares.rds"))
    manglares_v <- t(manglares_v) %>% as("dgCMatrix") # transpose [rows == ecosystem, columns == cell]
    manglares_v <- manglares_v[, ids]
    manglares_v[is.na(manglares_v)] <- 0  # fix NAs before converting
    
    ## Add to features and targets list
    features_list[["manglares"]] <- manglares_v
    targets_list[["manglares"]] <- (strat_ecos_target/100)
    rm(manglares_v)
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
  
  ## Assign appropriate feature types. Easier for interpreting summary CSV outputs
  feature_lookup <- bind_rows(
    if (!is.null(features_list[["ecosystems"]]))
      tibble(feature = rownames(features_list[["ecosystems"]]),
             feature_type = "ecosystem", class = NA_character_),
    
    if (!is.null(features_list[["manglares"]]))
      tibble(feature = rownames(features_list[["manglares"]]),
             feature_type = "strategic ecosystem", class = NA_character_)
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
    ## **NOTE: change threads and node_file_start depending on computer config**
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
    left_join(feature_lookup, by = "feature") %>%  # Add feature types
    mutate(scenario = model_name,           # Add the scenario info
           evaluated = "prioritizr_model")  # These features explicitly evaluated in model
  
  write_csv(target_coverage, file.path(opt_dir, paste0(model_name, "_summary.csv")))
  
  ## Rasterize solution and save
  s_rast <- rasterize_soln(s, template, locked_in, ids)
  writeRaster(s_rast,
              file.path(opt_dir, paste0(model_name, ".tif")),
              overwrite = TRUE)
  
  
  # --------- POST-HOC EVALUATION ---------------------------------------------
  # Because some ecosystems (that already met baseline target) 
  # were not included, need to get their coverage stats at a national level as well
  message("Running post-hoc evaluation for scenario: ", model_name)
  
  ## ------- Ecosystems ----------------------------------------
  eco_coverage <- NULL  #Create empty dataframe to store info in below
  
  if (ecos_target != 0) {  # For now, ecosystems always included. But make it flexible in case
    
    ecosys_mat <- readRDS(file.path(ipt_dir, "ecosistemas_marinos.rds"))[ids, ]
    ecosys_totals <- colSums(ecosys_mat)
    ecosys_names <- colnames(ecosys_mat)
    
    evaluated <- target_coverage$feature
    unevaluated <- setdiff(ecosys_names, evaluated)
    
    if (length(unevaluated) > 0) {
      ## Read in spreadsheet provided by INVEMAR
      ecos_target_df <- 
        read_excel(file.path(ipt_dir, "Umbrales Profundo somero_V2.xlsx"),
                   .name_repair = "unique_quiet") %>% 
        janitor::clean_names() %>% 
        ## Clean biome names to match more easily
        mutate(feature = clean(definicion)) %>% 
        filter(feature %in% unevaluated) %>% 
        ## Make sure no duplicates
        distinct(feature, .keep_all = TRUE)
      
      ## Pull the correct target column (e.g. 30 or 50)
      target_col <- sprintf("umbral_del_%s_percent", ecos_target)
      stopifnot(target_col %in% names(ecos_target_df))
      
      ## Align targets to `unevaluated` order, then convert to proportion
      rel_target <- ecos_target_df[[target_col]][match(unevaluated, ecos_target_df$feature)] / 100
      
      ## Warn about (and handle) ecosystems with no match in the spreadsheet
      if (anyNA(rel_target)) {
        warning(model_name, ": no unique target found for ", sum(is.na(rel_target)),
                " ecosystem(s); using general target (", ecos_target, "%).")
        rel_target[is.na(rel_target)] <- ecos_target / 100
      }
      
      totals     <- ecosys_totals[unevaluated]
      in_soln    <- colSums(ecosys_mat[s == 1, unevaluated, drop = FALSE])
      abs_target <- totals * rel_target
      rel_held   <- in_soln / totals
      
      eco_coverage <- data.frame(
        feature = unevaluated,
        met = rel_held >= rel_target,
        total_amount = totals,
        absolute_target = abs_target,
        absolute_held = in_soln,
        absolute_shortfall = pmax(0, abs_target - in_soln),
        relative_target = rel_target,
        relative_held = rel_held,
        relative_shortfall = pmax(0, rel_target - rel_held),
        scenario = model_name,
        evaluated = "post-hoc",
        feature_type = "ecosystem",
        class = NA_character_,
        row.names = NULL
      )
    }
    
    rm(ecosys_mat); gc()
    
  } # END ecosystem post-hoc
  
  
  
  ## ------- Combine and save -------------------------------------
  # ## Put all the post-hoc evaluations in one dataframe
  # post_hoc_coverage <- eco_coverage # ONLY ONE FOR NOW
  # rownames(post_hoc_coverage) <- NULL
  # 
  # ## Combine with solution target coverages, then save
  # target_coverage_full <- rbind(target_coverage, post_hoc_coverage)
  target_coverage_full <- dplyr::bind_rows(target_coverage, eco_coverage)
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
## -------- First Pass --------------------------------------
## First time running all scenarios.
## If process interrupted part-way, remove scenarios already completed or failed.
completed_list <- file.path(opt_dir, "master_eval_summary.csv")
failed_list    <- file.path(opt_dir, "failed_scenarios.txt")

if (file.exists(completed_list)) {
  completed <- read_csv(completed_list)
  scenarios_mar_df <- scenarios_mar_df %>%
    filter(!model_name %in% completed$scenario)
  rm(completed)
}

if (file.exists(failed_list)) {
  failed <- read.table(failed_list, sep = "|",
                       col.names = c("time", "model_name", "error"),
                       strip.white = TRUE)
  failed_list <- unique(failed$model_name)
  scenarios_mar_df <- scenarios_mar_df %>%
    filter(!model_name %in% failed_list)
  rm(failed); rm(failed_list)
}

## Iterate model over list of scenarios
purrr::pmap(scenarios_mar_df, marine_model, force_s = TRUE)


## -------- Second Pass --------------------------------------
## Once complete, evaluate the list of failed scenarios and errors. 
## Here, we'll force any scenarios to run that failed a presolve check
completed_list <- file.path(opt_dir, "master_eval_summary.csv")
failed_list    <- file.path(opt_dir, "failed_scenarios.txt")

if (file.exists(failed_list)) {
  ## Get full list
  failed <- read.table(failed_list, sep = "|",
                       col.names = c("time", "model_name", "error"),
                       strip.white = TRUE) %>% 
    mutate(time = as.POSIXct(time))
  
  ## Only the latest failure entry per scenario
  ## (e.g. maybe first time failed presolve, but later failed from memory issue)
  latest_failed <- failed %>%
    arrange(model_name, time) %>%
    group_by(model_name) %>%
    slice_tail(n = 1) %>%
    ungroup()
  
  ## Only rerun those that failed from presolve check
  presolve_fails <- latest_failed %>% 
    filter(grepl("^Presolve check failed", error)) %>%
    pull(model_name) %>%
    unique()
  
  rerun_df <- scenarios_mar_df %>%
    filter(model_name %in% presolve_fails) 
  
  ## Make sure no previously completed models are included
  if (file.exists(completed_list)){
    completed <- read_csv(completed_list)
    rerun_df <- rerun_df %>%
      filter(!model_name %in% completed$scenario)
  }
}


## Rerun and force solutions
purrr::pmap(rerun_df, marine_model, skip_presolve = TRUE, force_s = TRUE)


# NOTE: If there are still unrun scenarios after the second attemp, 
# they are likely due to another error (such as Memory related).
# Investigate error messages carefully, but it is likely is due to higher RAM requirements. 
# A separate virtual machine may be needed. 
