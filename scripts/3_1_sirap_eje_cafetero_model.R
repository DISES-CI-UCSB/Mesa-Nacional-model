## script: SIRAP Eje Cafetero Model
## Purpose: Set parameters and run prioritization model for Eje Cafetero region in Colombia

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

ipt_dir <- here("data/model_inputs/sirap/eje_cafetero")
opt_dir <- here("results/sirap/eje_cafetero")

if (!dir.exists(opt_dir)) dir.create(opt_dir, recursive = TRUE)

## Use specific template for model
template <- template_ec

## Get cell area in km2
## **NOTE: although EPSG9377 is not equal-area, checked distortion to be -0.08%. Basically negligible, so using simple constant.**
cell_res_m <- terra::res(template)
cell_area_km2 <- (cell_res_m[1] * cell_res_m[2]) / 1e6

# ========== PRIORITZATION FUNCTION ============================================
# Wrapping all the model building and running inside a function to easily run over a list of scenarios.

#' @param strat_ecos_target Numeric. Target percentage (0-100) for strategic 
#'   ecosystems; 0 to exclude.
#' @param bs_target Numeric. Target percentage (0-100) for bosque seco specifically. 
#'   If NA, evaluated at strat_ecos_target. If not NA, evaluated separately.
#' @param hum_target Numeric. Target percentage (0-100) for the humedales data
#'   unique to Eje Cafetero. 
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

eje_model <- function(strat_ecos_target, bs_target, hum_target, includes, cost, 
                      model_name, skip_presolve = FALSE, force_s = FALSE) { 
  ## Print scenario and time of start
  message("Running scenario: ", model_name, 
          "\nRun start: ", format(Sys.time(), "%H:%M:%S"))
  
  # --------- PLANNING UNITS ------------------------------------------
  if (cost == "IHEH2022") {
    pus <- readRDS(file.path(ipt_dir, "IHEH_EC_2022.rds"))
  } else if (cost == "IHEH2030") {
    pus <- readRDS(file.path(ipt_dir, "IHEH_EC_2030.rds"))
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
  locked_in <- readRDS(file.path(ipt_dir, "runap_EC.rds"))[ids, ] == 1
  
  ## Add to locked_in matrix depending on scenario
  if ("OMEC" %in% includes) {
    ## Update to make any cell either condition as TRUE
    locked_in <- locked_in | (readRDS(file.path(ipt_dir, "omec_EC.rds"))[ids, ] == 1)
  }
  
  
  # --------- FEATURES & TARGETS ----------------------------------------------
  # Add features (and their targets) to empty lists if they are evaluated
  # in the specific scenario
  features_list <- list()
  targets_list <- list()
  
  ## -------- Strategic Ecosystems --------------------------------
  ## Add if either all strategic ecosystems or bosque seco are evaluated
  if (strat_ecos_target != 0 | !is.na(bs_target)) {
    ## Read in matrix
    strat_ecos_v <- readRDS(file.path(ipt_dir, "ecosistemas_estrategicos_EC.rds"))
    strat_ecos_v <- t(strat_ecos_v) %>% as("dgCMatrix") # transpose [rows == ecosystem, columns == cell]
    strat_ecos_v <- strat_ecos_v[, ids]
    strat_ecos_v[is.na(strat_ecos_v)] <- 0  # fix NAs before converting
    
    ## If bosque seco has distinct target, separate feature
    if (!is.na(bs_target)) {
      ## Pull bosque seco into own matrix, and remove from strategic ecosystems
      bs_v <- strat_ecos_v[rownames(strat_ecos_v) == "bosque_seco", ]
      strat_ecos_v <- strat_ecos_v[rownames(strat_ecos_v) != "bosque_seco", ]
      
      ## Add both to features and targets list
      features_list[["strategic ecosystems"]] <- strat_ecos_v
      features_list[["bosque seco"]] <- bs_v
      
      targets_list[["strategic ecosystems"]] <- rep(strat_ecos_target/100, nrow(strat_ecos_v))
      targets_list[["bosque seco"]] <- (bs_target/100)
      
      rm(strat_ecos_v, bs_v)
      
    ## Otherwise, keep all three together and add to list
    } else {
      features_list[["strategic ecosystems"]] <- strat_ecos_v
      targets_list[["strategic ecosystems"]] <- rep(strat_ecos_target/100, nrow(strat_ecos_v))
      rm(strat_ecos_v)
    }
  }
  
  
  ## -------- Eje Cafetero Wetlands ----------------------------------
  # Read in matrix and add to features list if evaluated
  if (hum_target != 0) {
    hum_v <- readRDS(file.path(ipt_dir, "humedales_EC.rds"))
    hum_v <- t(hum_v) %>% as("dgCMatrix")
    hum_v <- hum_v[, ids]
    hum_v[is.na(hum_v)] <- 0
    
    ## Add to features and targets list
    features_list[["EC wetlands"]] <- hum_v
    targets_list[["EC wetlands"]] <- (hum_target/100)
    rm(hum_v)
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
    mutate(scenario = model_name,  # Add the scenario info
           ## Was explicitly included in the model
           evaluated = "prioritizr_model",
           ## Translate number of PUs into area
           total_amount_km2 = total_amount * cell_area_km2,
           absolute_held_km2 = absolute_held * cell_area_km2,
           feature_type = case_when(
             feature %in% rownames(features_list[["strategic ecosystems"]]) ~ "strategic ecosystem",
             feature %in% rownames(features_list[["bosque seco"]]) ~ "bosque seco",
             feature %in% rownames(features_list[["EC wetlands"]]) ~ "EC wetlands",
             TRUE ~ NA_character_
           ),
           class = NA_character_)
  
  write_csv(target_coverage, file.path(opt_dir, paste0(model_name, "_summary.csv")))
  
  ## Rasterize solution and save
  s_rast <- rasterize_soln(s, template, locked_in, ids)
  writeRaster(s_rast,
              file.path(opt_dir, paste0(model_name, ".tif")),
              overwrite = TRUE)
  
  
  # --------- POST-HOC EVALUATION (SPECIES & ECOSYSTEMS) -----------------------
  # Species and ecosystems aren't explicit features in this model, so evaluate
  # their coverage post-hoc against BOTH the 17% and 30% thresholds, using only
  # cells within the Orinoquia region (ids). These are for providing statistics
  # in the webtool.
  message("Running post-hoc evaluation for scenario: ", model_name)
  
  ## ------- Ecosystems ----------------------------------------
  ## NOTE: adjust filename to match your Orinoquia ecosystems matrix
  ecosys_mat <- readRDS(file.path(ipt_dir, "ecosistemas_IAVH_2024_EC.rds"))[ids, ]
  ecosys_totals  <- colSums(ecosys_mat)
  ecosys_in_soln <- colSums(ecosys_mat[s == 1, , drop = FALSE])
  rm(ecosys_mat); gc()
  
  eco_coverage <- tibble(
    feature = names(ecosys_totals),
    total_amount = ecosys_totals,
    absolute_held = ecosys_in_soln
  ) %>% 
    # filter(total_amount > 0) %>%
    mutate(
      relative_held     = absolute_held / total_amount,
      total_amount_km2  = total_amount * cell_area_km2,
      absolute_held_km2 = absolute_held * cell_area_km2,
      met               = NA,
      relative_target   = NA_real_,
      scenario     = model_name,
      evaluated    = "post-hoc",
      feature_type = "ecosystem",
      class        = NA_character_
    )
  
  ## ------- Species --------------------------------------------
  taxon_names <- c("Aves", "Amphibia", "Mammalia", "Crocodylia",
                   "Squamata", "Magnoliopsida_1", "Magnoliopsida_2")
  
  taxon_files <- list.files(ipt_dir, pattern = "_EC\\.rds$", full.names = TRUE) %>%
    keep(~ tools::file_path_sans_ext(basename(.x)) %>%
           str_remove("_EC$") %in% taxon_names)
  
  spp_coverage <- list()
  
  for (f in taxon_files) {
    taxon_name <- tools::file_path_sans_ext(basename(f)) %>% str_remove("_EC$")
    message("Processing species group: ", taxon_name)
    
    mat <- readRDS(f)[ids, ]
    spp_totals  <- colSums(mat)
    spp_in_soln <- colSums(mat[s == 1, , drop = FALSE])
    rm(mat); gc()
    
    spp_coverage[[taxon_name]] <- tibble(
      feature = names(spp_totals),
      total_amount = spp_totals,
      absolute_held = spp_in_soln
    ) %>% 
      # filter(total_amount > 0) %>%  # drop species not present in the region
      mutate(
        relative_held     = absolute_held / total_amount,
        total_amount_km2  = total_amount * cell_area_km2,
        absolute_held_km2 = absolute_held * cell_area_km2,
        met               = NA,
        relative_target   = NA_real_,
        scenario     = model_name,
        evaluated    = "post-hoc",
        feature_type = "species",
        class        = taxon_name
      )
  }
  
  post_hoc_coverage <- bind_rows(c(list(ecosystems = eco_coverage), spp_coverage))
  
  
  ## ------- Combine with explicit target coverage and save ---------------
  ## Using bind_rows (not rbind) since post-hoc rows carry met_17/met_30
  ## instead of a single "met" column — mismatched columns fill as NA.
  target_coverage_full <- bind_rows(target_coverage, post_hoc_coverage) %>%
    mutate(class = case_when(
      class %in% c("Magnoliopsida_1", "Magnoliopsida_2") ~ "Magnoliopsida",
      .default = class
    ))
  
  write_csv(target_coverage_full, file.path(opt_dir, paste0(model_name, "_summary.csv")))
  
  
  ## Get overview stats and add to running list of solutions
  freq_tbl <- freq(s_rast)
  cost_summary <- eval_cost_summary(p, s) # Is this even meaningful?
  eval_summary <- data.frame(
    scenario = model_name,
    n_total = sum(freq_tbl$count),
    n_new_protection = get_freq(freq_tbl, "Priority area"),
    n_locked_in = get_freq(freq_tbl, "Locked in"),
    area_total_km2 = sum(freq_tbl$count) * cell_area_km2,
    area_new_protection_km2 = get_freq(freq_tbl, "Priority area") * cell_area_km2,
    area_locked_in_km2 = get_freq(freq_tbl, "Locked in") * cell_area_km2,
    cost = cost_summary$cost,
    pct_targets_met = mean(target_coverage$met, na.rm = TRUE) * 100
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
## If process stopped part-way, remove scenarios already completed or permanently failed
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



## quick code for combining all results
cols_to_round <- c("total_amount_km2", "absolute_held_km2")

csv_files <- list.files(opt_dir, pattern = "\\.csv$", full.names = T)
csv_files <- csv_files[-41]
master_df <- csv_files %>%
  map_dfr(~ read_csv(.x, show_col_types = FALSE) %>%
            mutate(across(all_of(cols_to_round), ~ round(.x, 4)))) %>% 
  relocate(scenario, .before = everything()) %>% 
  relocate(total_amount_km2, .before = absolute_target) %>% 
  relocate(absolute_held_km2, .before = absolute_shortfall) %>% 
  mutate(absolute_target_km2 = absolute_target * cell_area_km2, .before = absolute_held)

# ---- Write out the master CSV ----
write_csv(master_df, file.path(opt_dir, "resultados_todos.csv"))

df <- read_csv(file.path(opt_dir, "master_eval_summary.csv"))

