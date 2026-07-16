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
  locked_in <- readRDS(file.path(ipt_dir, "runap_terrestres.rds"))[ids, ] == 1
  
  ## Add to locked_in matrix depending on scenario
  if ("OMEC" %in% includes) {
    ## Update to make any cell either condition as TRUE
    locked_in <- locked_in | (readRDS(file.path(ipt_dir, "omec_terrestres.rds"))[ids, ] == 1)
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
      mat <- readRDS(file.path(ipt_dir, "biomod_filtered_responsibilidad_nacional.rds"))
      species_rij <- mat %>% t() %>% as("dgCMatrix"); rm(mat)    # transpose [rows == spp, columns == cell]
      species_rij <- species_rij[, ids]
      
      ## Filter dataframe
      species_df <- read_csv(file.path(ipt_dir, "biomod_spp_responsibilidad_nacional.csv"), show_col_types = FALSE) %>% 
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
    taxon_names <- c("Aves", "Mammalia", "Crocodylia", 
                     "Squamata", "Magnoliopsida_1", "Magnoliopsida_2")
    
    taxon_files <- list.files(ipt_dir, pattern = "\\.rds$", full.names = TRUE) %>% 
      keep(~ tools::file_path_sans_ext(basename(.x)) %in% taxon_names)
    
    ## Will collect one matrix per class, rbind at the end
    unmet_spp_list <- list()
    
    ## Loop through each taxonomic class:
    ## Determine if any species haven't met targets and then add to list.
    for (f in taxon_files) {
      taxon_name <- tools::file_path_sans_ext(basename(f))
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
        add_gurobi_solver(gap = 0.05, threads = 15, verbose = TRUE)
      
      ## Solve
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
    taxon_names <- c("Aves", "Mammalia", "Crocodylia", 
                     "Squamata", "Magnoliopsida_1", "Magnoliopsida_2")
    
    taxon_files <- list.files(ipt_dir, pattern = "\\.rds$", full.names = TRUE) %>% 
      keep(~ tools::file_path_sans_ext(basename(.x)) %in% taxon_names)
    
    ## Loop through each group 
    for (f in taxon_files) {
      taxon_name <- tools::file_path_sans_ext(basename(f))
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
          read_csv(file.path(ipt_dir, "biomod_spp_responsibilidad_nacional.csv"), show_col_types = FALSE) %>% 
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
    ecosys_mat <- readRDS(file.path(ipt_dir, "ecosistemas_IAVH_2024.rds"))[ids, ]
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
      class %in% c("Magnoliospida_1", "Magnoliospida_2") ~ "Magnoliospida",
      .default = class
    ))
  
  write_csv(target_coverage_full, file.path(opt_dir, paste0(model_name, "_summary.csv")))

  
  ## Get overview stats and add to running list of solutions
  freq_tbl <- freq(s_rast)
  cost_summary <- eval_cost_summary(p, s) # Is this even meaningful?
  eval_summary <- data.frame(
    run = model_name,
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
## Testing subset for now
# scenarios_df <- scenarios_df[1:5, ]


## If process stopped part-way, use this code to remove scenarios already completed
completed <- read_csv(file.path(opt_dir, "master_eval_summary.csv"))
# failed <- file.path(opt_dir, "failed_scenarios.txt")
# failed <- read.table(failed, sep = "|",
#                          col.names = c("time", "model_name", "error"),
#                          strip.white = TRUE) %>%
#   ## remove memory errors for now
#   filter(error == "Error 10001: Out of memory")
# failed_list <- unique(failed$model_name)
scenarios_df <- scenarios_df %>%
  filter(!model_name %in% completed$run)
  # filter(!model_name %in% failed_list)
# rm(completed); rm(failed); rm(failed_list)

## Generate model over list of scenarios
purrr::pmap(scenarios_df, prioritizr_model)



