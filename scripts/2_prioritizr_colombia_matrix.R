## script: Prioritizr Colombia (matrix)
## Purpose: Set parameters and run prioritization models for Colombia

# ========== SETTING UP =======================================================
## load packages
library(prioritizr)  # modeling package
library(gurobi)      # solver
library(sf)          # vector data
library(terra)       # raster/GIS data
library(tidyverse)   # always
library(Matrix)      # using sparse matrices
library(readxl)      # read .xls format
library(here)        # easier filepaths
library(purrr)       # run models over list
source(here("scripts/utils.R"))

## Set seed and directories
set.seed(500)

ipt_dir <- here("data/model_inputs")
opt_dir <- here("results_new")

for (dir in c(ipt_dir, opt_dir)){
  if (!dir.exists(dir)) dir.create(dir)  # Create directories if needed
}; rm(dir)

## Testing subset for now
scenarios_df <- scenarios_df[1:5, ]


# ========== PRIORITZATION FUNCTION ============================================
## Wrapping all the model building and running inside a function 
## to easily run over a list of scenarios.

# prioritizr_model <- function(target, cost, features, includes, model_name, ids, pus) {
prioritizr_model <- function(ecos_target, strat_ecos_target, sp_rep_target, 
                             sp_rn_target, ecos_serv_target, includes, cost, model_name) {
  ## Print scenario and time of start
  message("Running scenario: ", model_name, 
          "\nRun start: ", format(Sys.time(), "%H:%M:%S"))
  
  # --------- PLANNING UNITS ------------------------------------------
  if (cost == "IHEH2022") {
    pus <- readRDS(file.path(ipt_dir, "IHEH_2022.rds"))
  } else if (cost == "IHEH2030") {
    pus <- readRDS(file.path(ipt_dir, "IHEH_2030.rds"))
  }
  
  ## Get list of all non-NA cells (each cell == planning unit)
  ids <- cells(template_terra)
  n_pus <- length(ids) # number of planning units
  
  pus <- pus[ids, ]
  pus[is.na(pus)] <- 0 #shouldn't be any NAs but can use in case
  
  # --------- CONSTRAINTS (LOCK INS/OUTS) -----------------------------
  ## Unlist variable
  includes <- unlist(includes)
  
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
  

  
  # --------- FEATURES & TARGETS ----------------------------------------------
  # Add features (and their targets) to empty lists if they are evaluated
  # in the specific scenario
  features_list <- list()
  targets_list <- list()
  
  ## -------- Ecosystems ------------------------------------------
  ## All ecosystems
  if (ecos_target != 0) {
    ## Read in matrix
    ecosys_v <- readRDS(file.path(ipt_dir, "ecosistemas_IAVH_2024.rds"))
    ecosys_v <- t(ecosys_v) %>% as("dgCMatrix") # transpose [rows == ecosystem, columns == cell]
    ecosys_v <- ecosys_v[, ids] # only keep cells in PUs
    
    ## Read in df and filter to only ecosystems not meeting targets
    cons_type <- ifelse("OMEC" %in% includes, "OMEC", "RUNAP") 
    ecosys_df <- read_csv(file.path(ipt_dir, "ecosys_filtered.csv")) %>% 
      filter(targets == ecos_target,
             conservation_type == cons_type)
    
    ## Only add ecosystems that need to be evaluated to features list
    row_idx <- which(rownames(ecosys_v) %in% ecosys_df$feature)
    features_list[["ecosystems"]] <- ecosys_v[row_idx, ]
    
    ## Add to targets list
    targets_list[["ecosystems"]] <- rep(ecos_target/100, 
                                       nrow(features_list[["ecosystems"]]))
    ## Remove matrix to keep env memory low
    rm(ecosys_v)
  }
  
  ## Strategic ecosystems
  if (strat_ecos_target != 0) {
    strat_ecos_v <- readRDS(file.path(ipt_dir, "ecosistemas_estrategicos_terrestres.rds"))
    strat_ecos_v <- t(strat_ecos_v) %>% as("dgCMatrix") # transpose [rows == ecosystem, columns == cell]
    strat_ecos_v <- strat_ecos_v[, ids]
    strat_ecos_v[is.na(strat_ecos_v)] <- 0  # fix NAs before converting
    
    ## Add to features and targets list
    features_list[["strategic ecosystems"]] <- strat_ecos_v
    targets_list[["strategic ecosystems"]] <- rep(strat_ecos_target/100, nrow(strat_ecos_v))
    rm(strat_ecos_v)
  }
  
  
  ## -------- Ecosystem Services ----------------------------------
  # Read in matrix and add to features list if evaluated
  if (ecos_serv_target != 0) {
    ecos_serv_v <- readRDS(file.path(ipt_dir, "servicios_ecosistemicos.rds"))
    ecos_serv_v <- t(ecos_serv_v) %>% as("dgCMatrix")
    ecos_serv_v <- ecos_serv_v[, ids]
    ecos_serv_v[is.na(ecos_serv_v)] <- 0
    
    ## Add to features and targets list
    features_list[["ecosystem services"]] <- ecos_serv_v
    targets_list[["ecosystem services"]] <- rep(ecos_serv_target/100, nrow(ecos_serv_v))
    rm(ecos_serv_v)
  }
  
  
  ## -------- Species ------------------------------------------
  # Representativeness and national responsbility are mutually exclusive
  # Only run following code if either are evaluated in the scenario
  if (sp_rep_target != 0 | sp_rn_target == TRUE) {
    ## Are we "including" OMEC+RUNAP, or just RUNAP?
    species_cons_type <- ifelse("OMEC" %in% includes, "OMEC", "RUNAP")  
    
    ## Species representativeness
    if (sp_rep_target != 0) {
      ## First read in matrix
      mat <- readRDS(file.path(ipt_dir, "biomod_filtered_representatividad.rds"))
      species_rij <- mat %>% t() %>% as("dgCMatrix"); rm(mat)    # transpose [rows == spp, columns == cell]
      species_rij <- species_rij[, ids]
      
      ## Filter dataframe and then matrix by targets
      species_df <- read_csv(file.path(ipt_dir, "biomod_spp_filtered_representatividad.csv")) %>%
        filter(targets == sp_rep_target,                # match target
               conservation_type == species_cons_type,  # match RUNAP/RUNAP+OMEC
               class != "Actinopteri")                  # for now, don't consider fish (Elkin's recommendation) 
      
      row_idx <- setNames(seq_len(nrow(species_rij)), rownames(species_rij))
      idx <- row_idx[species_df$scientific_name]
      idx <- idx[!is.na(idx)]
      species_filtered <- species_rij[idx, ]; rm(species_rij)
      
      ## Add to features and targets list
      features_list[["species representativeness"]] <- species_filtered
      targets_list[["species representativeness"]] <- rep(sp_rep_target/100, nrow(species_filtered))
      rm(species_filtered)
      
    ## Species national responsibility
    } else if (sp_rn_target == TRUE) {
      mat <- readRDS(file.path(ipt_dir, "biomod_filtered_responsibilidad_nacional.rds"))
      species_rij <- mat %>% t() %>% as("dgCMatrix"); rm(mat)    # transpose [rows == spp, columns == cell]
      species_rij <- species_rij[, ids]
      
      ## Filter dataframe
      species_df <- read_csv(file.path(ipt_dir, "biomod_spp_responsibilidad_nacional.csv")) %>% 
        filter(target_met == FALSE,                    # hasn't met target yet
               conservation_type == species_cons_type) # match RUNAP/RUNAP+OMEC
      
      row_idx <- setNames(seq_len(nrow(species_rij)), rownames(species_rij))
      idx <- row_idx[species_df$scientific_name]
      
      ## Keep species_df synched
      matched <- !is.na(idx)
      idx <- idx[matched]
      species_df <- species_df[matched, ]
      
      species_filtered <- species_rij[idx, ]; rm(species_rij)
      
      ## Add to features and targets list
      features_list[["species national responsbility"]] <- species_filtered
      targets_list[["species national responsibility"]] <- species_df$responsibility # Using the percentage for targets rather than raw 'target_km2' since all other targets are "relative" 
      rm(species_filtered)
    }
  }
  
  
  ## -------- Combine Features ------------------------------------------
  features_mat <- do.call(rbind, features_list)
  targets_df <- round(unlist(targets_list, use.names = FALSE),4)
  
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
  boundaries <- prioritizr::boundary_matrix(template_terra)[ids, ids]
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
    add_gurobi_solver(gap = 0.05, threads = 15, verbose = TRUE)
  
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
  s_rast <- rasterize_soln(s, template_terra, locked_in, ids)
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
    
    if (any(summary_df$run == eval_summary$run)) {
      ## Does the run exist? If so, overwrite
      summary_df[summary_df$run == eval_summary$run, ] <- eval_summary
      
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
# # If process stopped part-way, use this code to remove scenarios already completed
# completed <- read_csv(file.path(opt_dir, "master_eval_summary.csv"))
# failed <- file.path(opt_dir, "failed_scenarios.txt")
# failed <- read.table(failed, sep = "|",
#                          col.names = c("time", "model_name", "error"),
#                          strip.white = TRUE) %>%
#   ## remove memory errors for now
#   filter(error == "Error 10001: Out of memory")
# failed_list <- unique(failed$model_name)
# scenarios_df <- scenarios_df %>%
#   filter(!model_name %in% completed$run)
#   # filter(!model_name %in% failed_list)
# rm(completed); rm(failed); rm(failed_list)

## Generate model over list of scenarios
purrr::pmap(scenarios_df, prioritizr_model)




# ========== FULL SUMMARY STATS ==============================================
# Because some species and ecosystems (that already met baseline target) 
# were not included, need to get their coverage stats at a national level as well
ids <- cells(template)

## --------- 1. Get which cells were part of each solution -------------------
## Get file paths
solution_files <- list.files(opt_dir, pattern = "\\.tif$", full.names = TRUE)

## Table of completed solutions
soln_lookup <- tibble(path = solution_files,
                      model_name = tools::file_path_sans_ext(basename(solution_files))) %>%
  filter(model_name %in% scenarios_df$model_name)

## For each scenario, get the selected cell indices and evaluated feature names
scenario_meta <- purrr::pmap(soln_lookup, function(path, model_name) {
  s_rast <- rast(path)
  soln_vec <- values(s_rast, dataframe = FALSE)[ids]
  selected_ids <- ids[!is.na(soln_vec) & soln_vec > 0]
  
  evaluated <- read_csv(
    file.path(opt_dir, paste0(model_name, "_summary.csv")),
    show_col_types = FALSE
  )$feature
  
  ## Pull target value for this scenario from scenarios_df
  target <- scenarios_df$target[scenarios_df$model_name == model_name]
  
  list(model_name = model_name, selected_ids = selected_ids, 
       evaluated = evaluated, target = target)
})


## --------- 2. Loop through species taxons -----------------------------------
## Load each species taxon class matrix, compute coverage across ALL scenarios, then discard it
class_names <- c("Aves", "Mammalia", "Crocodylia", 
                 "Squamata", "Magnoliospida_1", "Magnoliospida_2")

class_files <- list.files(ipt_dir, pattern = "\\.rds$", full.names = TRUE) %>% 
  keep(~ tools::file_path_sans_ext(basename(.x)) %in% class_names)

## Will collect one data frame per class, rbind at the end
spp_coverage <- list()

for (f in class_files) {
  class_name <- tools::file_path_sans_ext(basename(f))
  message("Processing class: ", class_name)
  
  mat <- readRDS(f)  # [cells x species], sparse binary
  
  ## Total range per species (denominator). Compute once per matrix
  spp_totals <- colSums(mat)
  spp_names  <- colnames(mat)
  
  ## For each scenario, subset to selected cells and compute coverage
  class_coverage <- purrr::map_dfr(scenario_meta, function(scen) {
    ## Which features were not explicitly evaluated?
    unevaluated <- setdiff(spp_names, scen$evaluated)
    if (length(unevaluated) == 0) return(NULL)
    
    ## Get stats
    totals <- spp_totals[unevaluated]
    in_soln <- colSums(mat[scen$selected_ids, unevaluated, drop = FALSE])
    abs_target <- totals * (scen$target / 100)
    rel_held <- in_soln / totals
    
    ## Create table matching format of prioritizr summary outputs
    tibble(
      feature = unevaluated,
      met = rel_held >= (scen$target / 100),
      total_amount = totals,
      absolute_target = abs_target,
      absolute_held = in_soln,
      absolute_shortfall = pmax(0, abs_target - in_soln),
      relative_target = scen$target / 100,
      relative_held = rel_held,
      relative_shortfall = pmax(0, (scen$target / 100) - rel_held),
      scenario = scen$model_name,
      type = "species",
      class = class_name
    )
  })
  
  spp_coverage[[class_name]] <- class_coverage
  rm(mat); gc()
}


## --------- 3. Loop of ecosystems -------------------------------------
ecosys_mat <- readRDS(file.path(ipt_dir, "ecosistemas_IAVH_2024.rds"))
ecosys_totals <- colSums(ecosys_mat)
ecosys_names <- colnames(ecosys_mat)

ecosys_coverage <- purrr::map_dfr(scenario_meta, function(scen) {
  unevaluated <- setdiff(ecosys_names, scen$evaluated)
  if (length(unevaluated) == 0) return(NULL)
  
  totals <- ecosys_totals[unevaluated]
  in_soln <- colSums(ecosys_mat[scen$selected_ids, unevaluated, drop = FALSE])
  abs_target <- totals * (scen$target / 100)
  rel_held <- in_soln / totals
  
  tibble(
    feature = unevaluated,
    met = rel_held >= (scen$target / 100),
    total_amount = totals,
    absolute_target = abs_target,
    absolute_held = in_soln,
    absolute_shortfall = pmax(0, abs_target - in_soln),
    relative_target = scen$target / 100,
    relative_held = rel_held,
    relative_shortfall = pmax(0, (scen$target / 100) - rel_held),
    scenario = scen$model_name,
    type = "ecosystem",
    class = NA_character_  #can replace this with specific ecosys class later
  )
})

rm(ecosys_mat); gc()

## --------- 4. Combine and save -------------------------------------
all_coverage <- bind_rows(c(spp_coverage, list(ecosys_coverage))) %>% 
  ## Merge all plants back 
  mutate(class = case_when(
    class %in% c("Magnoliospida_1", "Magnoliospida_2") ~ "Magnoliospida",
    .default = class
  ))
write_csv(all_coverage, file.path(opt_dir, "unevaluated_feature_coverage.csv"))



## NOTE: Ask Will if one master CSV vs individual ones are best

## Get the "type" of feature and "class" of already evaluated features
species_lookup <- purrr::map_dfr(class_files, ~ {
  mat <- readRDS(.x)
  tibble(
    feature = colnames(mat),
    type = "species",
    class = tools::file_path_sans_ext(basename(.x))
  )
})

ecosys_lookup <- tibble(
  feature = ecosys_names,  # already in memory from step 3
  type = "ecosystem",
  class = NA_character_)

feature_lookup <- bind_rows(species_lookup, ecosys_lookup) %>% 
  ## Merge all plants back 
  mutate(class = case_when(
    class %in% c("Magnoliospida_1", "Magnoliospida_2") ~ "Magnoliospida",
    .default = class
  ))


### If individual: 
all_coverage %>% 
  group_by(scenario) %>% 
  group_walk(~ {
    ## Make sure scenario summary CSV exists
    csv_path <- file.path(opt_dir, paste0(.y$scenario, "_summary.csv"))
    if (!file.exists(csv_path)) {
      message("Skipping (no CSV found): ", .y$scenario)
      return()
    }
    
    existing <- read_csv(csv_path, show_col_types = FALSE) %>%
      left_join(feature_lookup, by = "feature")
    bind_rows(existing, .x) %>%
      write_csv(csv_path)
  })

# 
# ### If one master one:
# all_summaries <-
#   list.files(opt_dir, 
#              pattern = "_summary\\.csv$", 
#              full.names = TRUE) %>%
#   map_dfr(read_csv, show_col_types = FALSE)
# 
# write_csv(all_summaries, file.path(opt_dir, "all_feature_coverage.csv"))
# 


## Exploring how many filtered out spp didn't meet targets
library(tidyverse)

# ============================================================
# Configuration
# ============================================================

# Path to the folder containing your CSVs
csv_dir <- here("results")

# Column names — adjust if yours differ
col_scenario    <- "scenario"       # NA for post-hoc (filtered-out) species
col_met_target  <- "met"            # logical/binary: did prioritizr meet the target?
col_coverage    <- "relative_held"       # numeric: proportion of range covered (post-hoc)
col_feature     <- "feature"        # species/feature name
col_relative_target <- "relative_target"

# ============================================================
# 1. Read CSVs from "Esp" scenarios only
# ============================================================

all_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
esp_files <- all_files[grepl("Esp", basename(all_files))]

message("Found ", length(esp_files), " 'Esp' scenario CSVs out of ",
        length(all_files), " total CSVs.")

results_raw <- map_dfr(esp_files, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    mutate(source_file = basename(f))
}) %>% 
  filter(type == "species")

# ============================================================
# 2. Classify species into two groups
# ============================================================

results <- results_raw %>%
  mutate(
    species_type = if_else(is.na(.data[[col_scenario]]),
                           "post_hoc",     # LC/NT filtered out, range coverage added post-hoc
                           "in_solution")  # explicitly included in prioritizr run
  )
# ============================================================
# 3. Determine whether each species met its target
#    - In-solution species: use the `met` column from prioritizr
#    - Post-hoc species:    compare `coverage` column against threshold
# ============================================================

# ============================================================
# 4. Summary: per scenario × species type
# ============================================================

summary_by_scenario <- results %>%
  group_by(source_file, species_type, .data[[col_relative_target]]) %>%
  summarise(
    n_species   = n(),
    n_met       = sum(met, na.rm = TRUE),
    n_not_met   = sum(!met, na.rm = TRUE),
    pct_not_met = round(100 * n_not_met / n_species, 1),
    .groups = "drop"
  ) %>%
  rename(relative_target = .data[[col_relative_target]])
# ============================================================
# 5. Summary: overall across all Esp scenarios
# ============================================================
summary_overall <- results %>%
  group_by(species_type, .data[[col_relative_target]]) %>%
  summarise(
    n_records        = n(),
    n_unique_species = n_distinct(.data[[col_feature]]),
    n_met            = sum(met, na.rm = TRUE),
    n_not_met        = sum(!met, na.rm = TRUE),
    pct_not_met      = round(100 * n_not_met / n_records, 1),
    .groups = "drop"
  ) %>%
  rename(relative_target = .data[[col_relative_target]])
# ============================================================
# 6. Which post-hoc species fail most often?
# ============================================================

posthoc_failures <- results %>%
  filter(species_type == "post_hoc", !met) %>%
  count(.data[[col_feature]], name = "n_scenarios_failed") %>%
  arrange(desc(n_scenarios_failed))

species_df <- read_csv(file.path(temp_dir, "biomod_spp_ranges_updatedIUCN.csv"))

posthoc_failures2 <- left_join(posthoc_failures, species_df, join_by("feature" == "scientific_name"))

by_class <- posthoc_failures2 %>% 
  group_by(class) %>% 
  summarize(count = n())
# ============================================================
# 7. Print results
# ============================================================

cat("\n===== OVERALL SUMMARY =====\n")
print(summary_overall)

cat("\n===== PER-SCENARIO SUMMARY (first 20 rows) =====\n")
print(head(summary_by_scenario, 20))

cat("\n===== TOP POST-HOC SPECIES FAILING ACROSS SCENARIOS =====\n")
print(head(posthoc_failures, 20))

# ============================================================
# 8. Quick visualisation: % of species not meeting target
# ============================================================


p <- summary_by_scenario %>%
  filter(species_type == "post_hoc") %>%
  ggplot(aes(x = reorder(source_file, pct_not_met),
             y = pct_not_met,
             fill = factor(relative_target))) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Especies LC/NT que no alcanzan los objetivos según el escenario",
    x     = "Escenario (CSV file)",
    y     = "% de especies que no alcanzan el objetivo",
    fill  = "objetivo relativo"
  ) +
  theme_minimal(base_size = 11)
p


p2 <- summary_by_scenario %>%
  filter(species_type == "post_hoc") %>%
  ggplot(aes(x = reorder(source_file, n_not_met),
             y = n_not_met,
             fill = factor(relative_target))) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Post-hoc species not meeting coverage target per scenario",
    x     = "Scenario (CSV file)",
    y     = "% of species not meeting target",
    fill  = "Relative target"
  ) +
  theme_minimal(base_size = 11)
p2


ggsave("coverage_failures_by_scenario.png", p,
       width = 10, height = max(6, 0.2 * length(esp_files)),
       dpi = 150)

message("Plot saved to coverage_failures_by_scenario.png")

# ============================================================
# 9. Export summaries
# ============================================================

write_csv(summary_by_scenario, "summary_by_scenario.csv")
write_csv(summary_overall,     "summary_overall.csv")
write_csv(posthoc_failures,    "posthoc_species_failures.csv")

message("CSV summaries written.")

