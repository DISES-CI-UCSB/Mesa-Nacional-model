## script: National Terrestrial Model
## Purpose: Set the parameters and run the terrestrial national prioritization models for Colombia

# ========== SETTING UP =======================================================
## Load packages
library(prioritizr)  # modeling package
library(gurobi)      # solver
library(Matrix)      # using sparse matrices
library(parallel)    # parLapply / cluster management

## Load required functions and objects
source("scripts/utils.R")

## Set seed and directories
set.seed(500)
ipt_dir <- here("data/model_inputs/national")
opt_dir <- here("results/national/terrestrial/NEW")

for (dir in c(ipt_dir, opt_dir)){
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)  # Create directories if needed
}; rm(dir)


# ========== PRIORITZATION FUNCTION ============================================
# Wrapping all the model building and running inside a function to easily run over a list of scenarios.

#' @param ecos_target Numeric. Target percentage (0-100) for ecosystem
#'   representation; 0 to exclude ecosystems from this scenario. 
#' @param strat_ecos_target Numeric. Target percentage (0-100) for strategic 
#'   ecosystems; 0 to exclude.
#' @param sp_rep_target Numeric. Target percentage (0-100) for species 
#'   representativeness; 0 to exclude. Mutually exclusive with sp_rn_target.
#' @param sp_rn_target Logical. If TRUE, evaluates species national 
#'   responsibility targets (species-specific) instead of representativeness (17% or 30%).
#' @param ecos_serv_target Numeric. Target percentage (0-100) for ecosystem 
#'   services; 0 to exclude.
#' @param includes Character vector. Which layers should be "locked-in" to the 
#'   solution (e.g. "RUNAP", "OMEC", "comunidades", "resguardos").
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

terrestrial_model <- function(ecos_target, strat_ecos_target, sp_rep_target, 
                              sp_rn_target, ecos_serv_target, includes, cost, model_name,
                              skip_presolve = FALSE, force_s = FALSE) { 
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
  locked_in <- readRDS(file.path(ipt_dir, "runap_terrestres.rds"))[ids, ] == 1
  
  ## Add to locked_in matrix depending on scenario
  if ("OMEC" %in% includes) {
    ## Update to make any cell either condition as TRUE
    locked_in <- locked_in | (readRDS(file.path(ipt_dir, "omec_terrestres.rds"))[ids, ] == 1)
  }
  if ("comunidades" %in% includes) {
    ## Update to make any cell either condition as TRUE
    locked_in <- locked_in | (readRDS(file.path(ipt_dir, "comunidades_national.rds"))[ids, ] == 1)
  }
  if ("resguardos" %in% includes) {
    ## Update to make any cell either condition as TRUE
    locked_in <- locked_in | (readRDS(file.path(ipt_dir, "resguardos_national.rds"))[ids, ] == 1)
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
    ecosys_v <- readRDS(file.path(ipt_dir, "ecosistemas_IAVH_2024_terrestres.rds"))
    ecosys_v <- t(ecosys_v) %>% as("dgCMatrix") # transpose [rows == ecosystem, columns == cell]
    ecosys_v <- ecosys_v[, ids] # only keep cells in PUs
    
    ## Read in df and filter to only ecosystems not meeting targets
    cons_type <- ifelse("OMEC" %in% includes, "RUNAP_OMEC", "RUNAP") 
    ecosys_df <- read_csv(file.path(ipt_dir, "terrestrial_ecosys_filtered.csv"), show_col_types = FALSE) %>% 
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
    ecos_serv_v <- readRDS(file.path(ipt_dir, "servicios_ecosistemicos_terrestres.rds"))
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
      species_df <- read_csv(file.path(ipt_dir, "biomod_spp_filtered_representatividad.csv"), show_col_types = FALSE) %>% 
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
      mat <- readRDS(file.path(ipt_dir, "biomod_filtered_responsibilidad_national.rds"))
      species_rij <- mat %>% t() %>% as("dgCMatrix"); rm(mat)    # transpose [rows == spp, columns == cell]
      species_rij <- species_rij[, ids]
      
      ## Filter dataframe
      species_df <- read_csv(file.path(ipt_dir, "biomod_spp_responsibilidad_national.csv"), show_col_types = FALSE) %>% 
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
  
  ## Assign appropriate feature types to solution target summaries later
  feature_lookup <- bind_rows(
    if (!is.null(features_list[["ecosystems"]]))
      tibble(feature = rownames(features_list[["ecosystems"]]),
             feature_type = "ecosystem", class = NA_character_),
    
    if (!is.null(features_list[["strategic ecosystems"]]))
      tibble(feature = rownames(features_list[["strategic ecosystems"]]),
             feature_type = "strategic ecosystem", class = NA_character_),
    
    if (!is.null(features_list[["ecosystem services"]]))
      tibble(feature = rownames(features_list[["ecosystem services"]]),
             feature_type = "ecosystem service", class = NA_character_),
    
    if (!is.null(features_list[["species representativeness"]]))
      tibble(feature = rownames(features_list[["species representativeness"]])) %>%
      left_join(species_df %>% select(scientific_name, class),
                by = c("feature" = "scientific_name")) %>%
      mutate(feature_type = "species"),
    
    if (!is.null(features_list[["species national responsbility"]]))
      tibble(feature = rownames(features_list[["species national responsbility"]])) %>% 
      left_join(species_df %>% select(scientific_name, class),
                by = c("feature" = "scientific_name")) %>% 
      mutate(feature_type = "species")
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
    add_gurobi_solver(gap = 0.05, threads = 10, verbose = TRUE)
  
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
    left_join(feature_lookup, by = "feature") %>%  # Add feature type (and species class)
    mutate(scenario = model_name,           # Add the scenario info
           evaluated = "prioritizr_model")  # These features explicitly evaluated in model
  
  write_csv(target_coverage, file.path(opt_dir, paste0(model_name, "_summary.csv")))

  ## Rasterize solution and save
  s_rast <- rasterize_soln(s, template_terra, locked_in, ids)
  writeRaster(s_rast,
              file.path(opt_dir, paste0(model_name, ".tif")),
              overwrite = TRUE)
  
  
  # --------- EVALUATE SPECIES (& RERUN) ---------------------------------------------
  # If a scenario evaluated species representativeness, double-check that 
  # all filtered-out species met the target (17% or 30%). If any species did NOT,
  # "lock in" solution and re-run for just those species.
  
  if (sp_rep_target != 0) {
    
    ## Load each taxon class names (as matching matrices)
    taxon_names <- c("Aves", "Amphibia", "Mammalia", "Crocodylia", 
                     "Squamata", "Magnoliopsida_1", "Magnoliopsida_2")
    
    taxon_files <- list.files(ipt_dir, pattern = "\\.rds$", full.names = TRUE) %>% 
      keep(~ tools::file_path_sans_ext(basename(.x)) %>% 
             str_remove("_national$") %in% taxon_names)
    
    ## Will collect one matrix per class, rbind at the end
    unmet_spp_list <- list()
    
    ## Loop through each taxonomic class:
    ## Determine if any species haven't met targets and then add to list.
    for (f in taxon_files) {
      taxon_name <- tools::file_path_sans_ext(basename(f)) %>% str_remove("_national$")
      message("Double-checking: ", taxon_name)
      
      ## Read in taxonomic class matrix
      mat <- readRDS(f) %>% t() 
      mat <- mat[, ids]
      
      ## Total range per species
      spp_totals <- rowSums(mat)
      spp_names  <- rownames(mat)
      
      ## Which species were not explicitly evaluated?
      evaluated <- target_coverage$feature
      unevaluated <- setdiff(spp_names, evaluated)
      if (length(unevaluated) == 0) next
      
      ## Get coverage stats
      totals <- spp_totals[unevaluated]
      in_soln <- rowSums(mat[unevaluated, s == 1, drop = FALSE]) # only look at cells in solution (s==1)
      abs_target <- totals * (sp_rep_target / 100)
      rel_held <- in_soln / totals
      
      ## Create table of species that don't meet target
      unmet_df <-
        tibble(species = unevaluated, 
               met = rel_held >= (sp_rep_target / 100)) %>%
        filter(met == FALSE)
      
      ## If any unmet species, filter the matrix and add to list
      if (nrow(unmet_df) > 0) {
        ## Filter matrix by species
        row_idx <- setNames(seq_len(nrow(mat)), rownames(mat))
        idx <- row_idx[unmet_df$species]
        idx <- idx[!is.na(idx)]
        mat_filtered <- mat[idx, ]
        
        ## Add to list
        unmet_spp_list[[taxon_name]] <- mat_filtered
        
        ## Track type/class alongside, keyed by feature name
        feature_lookup <- bind_rows(
          feature_lookup,
          tibble(feature = rownames(mat_filtered), feature_type = "species", class = taxon_name)
        )
        
        rm(mat_filtered)
      }
      
      rm(mat); gc() # clear memory
      
    } # END TAXON COVERAGE LOOP
    
    
    ## If any species didn't meet targets, re-run prioritization!
    unmet_spp_mat <- do.call(rbind, unmet_spp_list); rm(unmet_spp_list)
    
    if (!is.null(unmet_spp_mat)) {
      ## Keep all PUs previously selected, so add to locked_in
      locked_in_p2 <- locked_in | (s == 1)
      
      ## Update targets and features dataframes
      targets_p2 <- rep((sp_rep_target/100), unmet_spp_mat@Dim[1])
      
      features_p2 <- data.frame(
        id = 1:as.numeric(unmet_spp_mat@Dim[1]),
        name = unmet_spp_mat@Dimnames[[1]],
        target = targets_p2
      )
      
      ## Build problem, using same PUs and boundary data
      p2 <- problem(
        x = pus,
        features = features_p2,
        rij_matrix = unmet_spp_mat) %>% 
        add_min_set_objective() %>% 
        add_relative_targets(targets_p2) %>% 
        add_locked_in_constraints(locked_in_p2) %>%
        add_binary_decisions() %>% 
        add_boundary_penalties(penalty = 0.001, data = boundaries) %>% 
        add_gurobi_solver(gap = 0.05, threads = 10, verbose = FALSE)
      
      ## Always force this problem
      s2 <- solve(p2, force = TRUE)
      
      ## Get updated summary coverage & overwrite
      summary_1 <- eval_target_coverage_summary(p, s2)   # original features in new solution
      summary_2 <- eval_target_coverage_summary(p2, s2)  # new features in new solution
      
      target_coverage <- rbind(summary_1, summary_2) %>% 
        left_join(feature_lookup, by = "feature") %>% 
        mutate(scenario = model_name,           # Add the scenario info
               evaluated = "prioritizr_model")  # These features explicitly evaluated in model
      
      write_csv(target_coverage, file.path(opt_dir, paste0(model_name, "_summary.csv")))
      
      ## Replace original solution; Rasterize & save!
      s <- s2; rm(p2, s2, summary_1, summary_2)
      
      s_rast <- rasterize_soln(s, template_terra, locked_in, ids)
      writeRaster(s_rast,
                  file.path(opt_dir, paste0(model_name, ".tif")),
                  overwrite = TRUE)
      
    } # END CONDITIONAL RERUN 
    
  } # END CONDITIONAL SPECIES EVALUATION

  
  # --------- POST-HOC EVALUATION ---------------------------------------------
  # Because some species and ecosystems (that already met baseline target) 
  # were not included, need to get their coverage stats at a national level as well
  message("Running post-hoc evaluation for scenario: ", model_name)
  
  ## ------- Species ----------------------------------------
  spp_coverage <- list()
  
  ## Only run if species was evaluate
  if (sp_rep_target != 0 | sp_rn_target == TRUE) {
    taxon_names <- c("Aves", "Amphibia", "Mammalia", "Crocodylia", 
                     "Squamata", "Magnoliopsida_1", "Magnoliopsida_2")
    
    taxon_files <- list.files(ipt_dir, pattern = "\\.rds$", full.names = TRUE) %>% 
      keep(~ tools::file_path_sans_ext(basename(.x)) %>% 
             str_remove("_national$") %in% taxon_names)
    
    ## Loop through each group 
    for (f in taxon_files) {
      taxon_name <- tools::file_path_sans_ext(basename(f)) %>% str_remove("_national$")
      message("Processing species group: ", taxon_name)
      
      ## Total range per species (denominator). Compute once per matrix
      mat <- readRDS(f)[ids, ]  
      spp_totals <- colSums(mat)
      spp_names  <- colnames(mat)
      
      ## Which species were not explicitly evaluated?
      evaluated <- target_coverage$feature
      unevaluated <- setdiff(spp_names, evaluated)
      if (length(unevaluated) == 0) next
      
      ## Get coverage stats
      totals <- spp_totals[unevaluated]
      in_soln <- colSums(mat[s == 1, unevaluated, drop = FALSE]) # only look at cells in solution (s==1)
      rel_held <- in_soln / totals
      
      if (sp_rep_target != 0) {
        abs_target <- totals * (sp_rep_target / 100)
        
        class_coverage <- data.frame(
          feature = unevaluated,
          met = rel_held >= (sp_rep_target / 100),
          total_amount = totals,
          absolute_target = abs_target,
          absolute_held = in_soln,
          absolute_shortfall = pmax(0, abs_target - in_soln),
          relative_target = sp_rep_target / 100,
          relative_held = rel_held,
          relative_shortfall = pmax(0, (sp_rep_target / 100) - rel_held),
          scenario = model_name,
          evaluated = "post-hoc",
          feature_type = "species",
          class = taxon_name
        )
      } else {
        stats_df <- tibble(
          feature = unevaluated,
          total_amount = totals,
          absolute_held = in_soln,
          relative_held = rel_held
        )
        
        ## Get species-specific targets for national responsibility
        species_df <- 
          read_csv(file.path(ipt_dir, "biomod_spp_responsibilidad_national.csv"), show_col_types = FALSE) %>% 
          filter(conservation_type == species_cons_type,
                 scientific_name %in% unevaluated) %>% 
          select(scientific_name, responsibility) %>% 
          rename(feature = scientific_name,
                 relative_target = responsibility)
        
        class_coverage <- species_df %>% 
          inner_join(stats_df, by = "feature") %>% 
          mutate(
            met = rel_held >= relative_target,
            absolute_target = (totals * relative_target),
            absolute_shortfall = pmax(0, (totals*relative_target) - in_soln),
            relative_shortfall = pmax(0, relative_target - rel_held),
            scenario = model_name,
            evaluated = "post-hoc",
            feature_type = "species",
            class = taxon_name
          )
        
        rm(stats_df, species_df)
      }
      
      spp_coverage[[taxon_name]] <- class_coverage
      rm(mat); gc()
    } # END taxon loop
  } # END species post-hoc
  
  
  ## ------- Ecosystems ----------------------------------------
  eco_coverage <- NULL
  
  if (ecos_target != 0) {  # For now, ecosystems always included. But make it flexible in case
    ecosys_mat <- readRDS(file.path(ipt_dir, "ecosistemas_IAVH_2024_terrestres.rds"))[ids, ]
    ecosys_totals <- colSums(ecosys_mat)
    ecosys_names <- colnames(ecosys_mat)
    
    evaluated <- target_coverage$feature
    unevaluated <- setdiff(ecosys_names, evaluated)
    if (length(unevaluated) == 0) return(NULL)
    
    totals <- ecosys_totals[unevaluated]
    in_soln <- colSums(ecosys_mat[s == 1, unevaluated, drop = FALSE])
    abs_target <- totals * (ecos_target / 100)
    rel_held <- in_soln / totals
    
    eco_coverage <- data.frame(
      feature = unevaluated,
      met = rel_held >= (ecos_target / 100),
      total_amount = totals,
      absolute_target = abs_target,
      absolute_held = in_soln,
      absolute_shortfall = pmax(0, abs_target - in_soln),
      relative_target = ecos_target / 100,
      relative_held = rel_held,
      relative_shortfall = pmax(0, (ecos_target / 100) - rel_held),
      scenario = model_name,
      evaluated = "post-hoc",
      feature_type = "ecosystem",
      class = NA_character_  #can replace this with specific ecosys class later
    )
    
    rm(ecosys_mat); gc()
    
  } # END ecosystem post-hoc

  
  
  ## ------- Combine and save -------------------------------------
  ## Put all the post-hoc evaluations in one dataframe
  post_hoc_coverage <- bind_rows(c(spp_coverage, list(eco_coverage)))
  rownames(post_hoc_coverage) <- NULL
  
  ## Combine with solution target coverages, then save
  target_coverage_full <- rbind(target_coverage, post_hoc_coverage) %>% 
    ## Merge all plants back 
    mutate(class = case_when(
      class %in% c("Magnoliopsida_1", "Magnoliopsida_2") ~ "Magnoliopsida",
      .default = class
    ))
  
  write_csv(target_coverage_full, file.path(opt_dir, paste0(model_name, "_summary.csv")))

  
  ## Get overview stats to hand back to the caller.
  ## NOTE: this used to read-modify-write master_eval_summary.csv right here.
  ## That's fine when scenarios run one at a time, but once this function is
  ## called from multiple parallel workers, two workers can read the file at
  ## the same time and each write back a version that's missing the other's
  ## row -- a classic race condition. So instead we just return the row, and
  ## a helper (update_master_summary(), below) appends it to the CSV once,
  ## serially, in the main session after every worker has finished.
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

  gc()
  return(eval_summary)
  
} # END PRIORITIZR FUNCTION


#' Append (or update) scenario summary rows in the shared master_eval_summary.csv.
#' Call this ONCE, serially, from the main session after a batch of parallel
#' terrestrial_model() calls finishes -- never from inside a parallel worker
#' (see note above for why).
#'
#' @param results List returned by parLapply/parLapplyLB -- one element per
#'   scenario, each either a one-row data.frame (from a successful run) or
#'   NULL (scenario failed and was logged to failed_scenarios.txt).
#' @param opt_dir Character. Output directory containing master_eval_summary.csv.
update_master_summary <- function(results, opt_dir) {
  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) == 0) {
    message("No successful scenarios to add to master_eval_summary.csv")
    return(invisible(NULL))
  }
  new_df <- do.call(rbind, results)

  csv_path <- file.path(opt_dir, "master_eval_summary.csv")

  if (file.exists(csv_path)) {
    summary_df <- read_csv(csv_path, show_col_types = FALSE)

    for (i in seq_len(nrow(new_df))) {
      match_idx <- which(summary_df$scenario == new_df$scenario[i])
      if (length(match_idx) > 0) {
        ## Scenario already in master file -- overwrite its row
        summary_df[match_idx, ] <- new_df[i, ]
      } else {
        summary_df <- rbind(summary_df, new_df[i, ])
      }
    }

  } else {
    summary_df <- new_df  # If file hasn't yet been created, then create it
  }

  write_csv(summary_df, csv_path)
  invisible(summary_df)
}


# ========== RUN PRIORITIZATION ==============================================
## **NOTE: RERUNNING FOR AMPHIBIAN MODELS ONLY** -- mirrors the update made in
## 2_1_national_terrestrial_model.R. We're no longer re-running everything
## logged in failed_scenarios.txt; instead we rerun exactly the scenarios
## listed in rerun_scenarios.csv (the ones where amphibians didn't reach
## their target), now that "Amphibia" has been added to the species checks
## above.
rerun_scenarios <- read_csv("rerun_scenarios.csv")

if (nrow(rerun_scenarios) == 0) {
  stop("rerun_scenarios.csv was found but contained no rows.")
}

## ---- Parallel setup ----------------------------------------------------
## Each gurobi solve uses 15 threads. With 36 physical/logical cores (no
## hyperthreading), that leaves room for 2 concurrent solves (2 x 15 = 30
## cores used, 6 held back for the OS / R's own data prep work) without
## oversubscribing. RAM doesn't push this number up: 2 workers x 50GB peak
## = 100GB, far under the 1TB available, so CPU threads -- not memory -- are
## the binding constraint here.
##
## If you want to test a 3-wide run instead, drop threads = 15 to
## threads = 12 in BOTH add_gurobi_solver() calls inside terrestrial_model()
## (3 x 12 = 36) and compare total wall-clock time against the 2-wide/15-thread
## version -- gurobi's speedup per thread is sublinear, so more concurrent
## solves with fewer threads each sometimes wins on throughput, but it has to
## be benchmarked on your own problems/scenarios to know which is faster.
n_workers <- 3

cl <- makeCluster(n_workers, type = "FORK")
## FORK (Linux/macOS only) forks the current R session, so every worker
## already has prioritizr/gurobi/Matrix loaded and already has template_terra,
## rerun_scenarios, ipt_dir, opt_dir, and terrestrial_model() in scope --
## no clusterExport needed. On Windows, FORK isn't available; use
## type = "PSOCK" instead and add, before the parLapplyLB() call below:
##   clusterEvalQ(cl, {library(prioritizr); library(gurobi); library(Matrix)
##                      library(dplyr); library(readr); library(here)})
##   clusterExport(cl, c("template_terra", "rerun_scenarios", "ipt_dir",
##                        "opt_dir", "terrestrial_model", "presolve_check",
##                        "rasterize_soln", "get_freq"))
##   # plus any other object/function from utils.R that terrestrial_model() uses

## NOTE: unlike the old failed_scenarios.txt rerun, these scenarios are not
## being forced past a presolve failure -- they're being re-evaluated because
## Amphibia was added to the species checks inside terrestrial_model(), so
## skip_presolve/force_s are left at their defaults (FALSE), same as the
## purrr::pmap(rerun_scenarios, terrestrial_model) call in
## 2_1_national_terrestrial_model.R.
##
## Each row of rerun_scenarios is one task. parLapplyLB (load-balanced) hands
## a worker its next scenario as soon as it finishes its current one, which
## matters here since scenarios with sp_rep_target != 0 can trigger a second
## solve internally and run noticeably longer than others -- a plain
## parLapply would pre-split tasks evenly up front and could leave a worker
## idle while another is still stuck on a slow scenario.
scenario_args <- split(rerun_scenarios, seq_len(nrow(rerun_scenarios)))

results <- parLapplyLB(cl, scenario_args, function(row) {
  tryCatch(
    do.call(terrestrial_model, as.list(row)),
    error = function(e) {
      ## terrestrial_model() already catches errors from presolve_check()/
      ## solve() internally and logs them to failed_scenarios.txt. This is
      ## just a safety net for anything upstream of that (e.g. a bad
      ## readRDS path), so one bad scenario can't take down the whole batch.
      message("Uncaught error in ", row$model_name, ": ", e$message)
      NULL
    }
  )
})

stopCluster(cl)

## ---- Force any of THESE scenarios that failed the presolve check --------
## The first pass above uses the default skip_presolve/force_s = FALSE, so
## anything that fails presolve_check() gets logged to failed_scenarios.txt
## and skipped rather than solved. This second pass rereads that log,
## filters to only scenarios that are actually part of rerun_scenarios (so
## it doesn't pick up stale failures logged by some earlier, unrelated run),
## keeps just the ones that failed presolve specifically (not some other
## error), and reruns exactly those with skip_presolve = TRUE, force_s = TRUE.
##
## Reminder: forcing bypasses the presolve gate, it doesn't fix whatever made
## the problem infeasible in the first place -- Gurobi may return an
## infeasible or low-quality solution for a forced scenario. Worth spot-
## checking the *_summary.csv for anything that gets force-solved here.
failed_file <- file.path(opt_dir, "failed_scenarios.txt")

if (file.exists(failed_file)) {
  failed <- read.table(failed_file, sep = "|",
                       col.names = c("time", "model_name", "error"),
                       strip.white = TRUE) %>%
    mutate(time = as.POSIXct(time))

  ## Only the latest failure entry per scenario (e.g. maybe it failed
  ## presolve once, then a different error the next time it was tried).
  presolve_fails <- failed %>%
    filter(model_name %in% rerun_scenarios$model_name) %>%
    arrange(model_name, time) %>%
    group_by(model_name) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    filter(grepl("^Presolve check failed", error)) %>%
    pull(model_name) %>%
    unique()

  force_df <- rerun_scenarios %>%
    filter(model_name %in% presolve_fails)

  if (nrow(force_df) > 0) {
    message("Forcing ", nrow(force_df), " scenario(s) that failed presolve: ",
            paste(force_df$model_name, collapse = ", "))

    cl2 <- makeCluster(n_workers, type = "FORK")

    force_args <- split(force_df, seq_len(nrow(force_df)))

    force_results <- parLapplyLB(cl2, force_args, function(row) {
      tryCatch(
        do.call(terrestrial_model,
                c(as.list(row), list(skip_presolve = TRUE, force_s = TRUE))),
        error = function(e) {
          message("Uncaught error in ", row$model_name, ": ", e$message)
          NULL
        }
      )
    })

    stopCluster(cl2)

    ## Combine with the first pass's results so both end up in the master file
    results <- c(results, force_results)
  }
}

## Combine every returned summary row into master_eval_summary.csv, once,
## serially, here in the main session (see the note on update_master_summary()
## above for why this can't safely happen inside the parallel workers).
update_master_summary(results, opt_dir)
