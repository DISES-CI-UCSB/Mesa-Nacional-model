

# ========== SETTING UP =====================================================
## load packages
library(prioritizr)
library(sf)
library(terra)
library(vegan)
library(cluster)
library(tidyverse)
library(doParallel)
library(raster)
library(ggplot2)
library(gurobi)
library(Matrix)
library(here)

## Set seed and directories
set.seed(500)

ipt_dir <- here("data/model_inputs")
opt_dir <- here("results")

for (dir in c(ipt_dir, opt_dir)){
  if (!dir.exists(dir)) dir.create(dir)  # Create directories if needed
}

## Load cost raster as template
template <- rast(here("data/costs/human_footprint_2022.tif"))

# ========== TARGETS ==========================================================

## Determine which targets to use
target <-
  17
  # 30
  
includes <-
  "RUNAP"
  # "OMEC"


# ========== GET PLANNING UNITS ================================================

## Get list of all non-NA cells (each cell == planning unit)
ids <- cells(template)
n_pus <- length(ids) # number of planning units

## For now, just using IHEH 2022 as sole cost
pus <- readRDS(file.path(ipt_dir, "IHEH_2022.rds"))
pus <- pus[ids, ]
pus[is.na(pus)] <- 0 #shouldn't be any NAs but can use in case


# ========== CONSTRAINTS (LOCK INS/OUTS) =====================================

## Create lock_in matrix based on Targets selection
if (includes == "RUNAP") {
  runap <- readRDS(file.path(ipt_dir, "runap.rds"))
  runap <- runap[ids, ]
  
  locked_in <- (runap == 1)
  rm(runap)
  
} else if (includes == "OMEC") {
  runap <- readRDS(file.path(ipt_dir, "runap.rds"))
  runap <- runap[ids, ]
  omec <- readRDS(file.path(ipt_dir, "omec.rds"))
  omec <- omec[ids, ]
  
  locked_in <- (runap == 1) | (omec == 1)
  rm(runap); rm(omec)
}



# ========== FEATURES ========================================================

## Can change...
# features_list <- c("ecosistemas", "paramos", "manglares", "humedales", "bosque_seco")

## -------- Ecosystems -----------------------------------------------------
## All ecosystems
## NOTE: still need to fix this approach!!!!
## Issueing with scaling, so just omit for now
# ecosys_v <- readRDS(file.path(ipt_dir, "ecosistemas.rds"))

## Strategic ecosystems
strat_ecos_v <- readRDS(file.path(ipt_dir, "strategic_ecosystems.rds"))


## Combine all
# ecosystems <- cbind(ecosys_v, strat_ecos_v); rm(ecosys_v); rm(strat_ecos_v)
ecosystems <- strat_ecos_v; rm(strat_ecos_v)
ecosystems <- ecosystems[ids,]      # Only keep PUs
ecosystems[is.na(ecosystems)] <- 0  # Remove lingering NAs

## Transpose matrix and make sparse to match problem format
ecosystems <- t(ecosystems)
ecosys_sparse <- as(ecosystems, "sparseMatrix"); rm(ecosystems)


## -------- Species -----------------------------------------------------
mat <- readRDS(file.path(ipt_dir, "biomod_filtered.rds"))
## transpose (rows == spp, columns == cell)
species_rij <- mat %>% t() %>% as("dgCMatrix"); rm(mat)
species_rij <- species_rij[, ids]

## Filter matrix by determined goals
species_df <- read_csv(file.path(ipt_dir, "biomod_spp_ranges_filtered.csv")) %>% 
  filter(targets == target,
         conservation_type == includes)

row_idx <- setNames(seq_len(nrow(species_rij)), rownames(species_rij))
idx <- row_idx[species_df$scientific_name]
idx <- idx[!is.na(idx)]

species_filtered <- species_rij[idx, ]; rm(species_rij)


## Combine all features into one mat
features_mat <- rbind(ecosys_sparse, species_filtered)

# ## try to fix scaling issue?
# features_scalar <- max(features_mat)
# features_mat <- features_mat / features_scalar

n_features <- as.numeric(features_mat@Dim[1])
feature_names <- features_mat@Dimnames[[1]]

features_df <- data.frame(
  id = 1:n_features,
  name = feature_names
)

# ========== SET PROBLEM ======================================================

# Set relative targets for features
## For now, these are all the same!
targets_df <- rep(target/100, n_features)
boundaries <- prioritizr::boundary_matrix(template)[ids, ids]
boundaries <- boundaries/max(boundaries) #scaling issue


# Build problem with locked_in constraints
## NOTE: eventually change this to a function to run over multiple targets and input scenarios...
p1 <- problem(
  x = pus,
  features = features_df,
  rij_matrix = features_mat) %>% 
  add_min_set_objective() %>% 
  add_relative_targets(targets_df) %>% 
  add_locked_in_constraints(locked_in) %>% 
  # add_locked_in_constraints(comunidades) |>
  # add_locked_in_constraints(resguardo) |>
  add_binary_decisions() %>% 
  add_boundary_penalties(penalty = 0.001, data = boundaries) %>% 
  # add_highs_solver(gap = 0.1, threads = 4)
  add_gurobi_solver(gap = 0.05, threads = 6)

print(p1)
presolve_check(p1)


# ========== SOLVE PROBLEM ====================================================
## Get solution
s1 <- solve(p1)

## Create fxn to rasterize solution (outputs as matrix)
rasterize_soln <- function(s, template) {
  # Create output raster from template
  rast <- template
  rast[] <- NA
  
  # Assign solution values to planning unit cells
  rast[ids] <- s
  
  # Mark existing PAs (locked-in units that were selected)
  # 1 = new cells selected; 2 = existing PA; NA = not selected
  rast[ids[which(locked_in == 1)]] <- 2
  
  # Set 0s to NA (not selected)
  rast[rast == 0] <- NA
  
  # Add category labels
  levels(rast) <- data.frame(
    value = 1:2,
    layer = c("Selected", 
              "Locked in") #NOTE!! : change this to just locked-in bc sometimes includes communities..
  )
  
  return(rast)
}

s_rast <- rasterize_soln(s1, template)
writeRaster(s_rast, 
            file.path(opt_dir, sprintf("solution_%s_%s.tif", target, includes)), overwrite = TRUE) ## ALSO CHANGE THIS TO SPRINTF to bring in targets


# ========== EVALUATE RESULTS =================================================
run_name <- "Ecos30+ESTR30+RUNAP+Comunidades_HF"  # NOTE: build this into workflow!!
target_coverage <- eval_target_coverage_summary(p1, s1) %>% 
  mutate(scenario = run_name)
# filter(absolute_target > 0)

## Save summary df
write_csv(target_coverage,
          file.path(opt_dir, sprintf("solution_summary_%s_%s.csv", target, includes)))



## Get overview stats and add to running solution df
freq <- freq(s_rast)
cost_summary <- eval_cost_summary(p1, s1) # Is this even meaninful? 
eval_summary <- data.frame(
  run = run_name,
  n_total = sum(freq$count),
  n_new_protection = pluck(subset(freq, value == "Selected"), "count"),
  n_locked_in = pluck(subset(freq, value == "Locked in"), "count"),
  cost = cost_summary$cost,
  pct_targets_met = mean(p1_target_coverage$met) * 100
)

# Append eval_summary row to master CSV 
if (file.exists(file.path(opt_dir, "master_eval_summary.csv"))) {   # Does this file exist?
  ## Read in
  summary_df <- read_csv(file.path(opt_dir, "master_eval_summary.csv"))
  
  ## Does the run exist? If so, overwrite
  if (summary_df$run == eval_summary$run) {
    summary_df[summary_df$run == eval_summary$run, ] <- eval_summary
  } else {
    ## If not, append to table
    summary_df <- rbind(summary_df, eval_summary)
  }
  
  ## Save
  write_csv(summary_df, file.path(opt_dir, "master_eval_summary.csv"))
  
  ## If file hasn't yet been created, then create it
} else {
  write.csv(eval_summary, file.path(opt_dir, "master_eval_summary.csv"), row.names = FALSE)
}


