# load packages
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

#set seed and working directory
set.seed(500)
setwd("C:/Users/kevinramos/Year1/Colombia_Prioritization/mesa_prioridades/mesa_prioridades_kevin")
# =====================================================
# PROTECTED AREA AND INCLUDES RASTER STACK
# =====================================================

#path to raster stack
PA_RASTER_STACK_PATH <- "./costs_and_constraints/cost_constraints_stack_1km.tif"

# =====================================================
# LOAD FEATURES
# =====================================================
ecosistemas <- rast("./features/1/ecosistemas.tif")
paramos <- rast("./features/4/paramos.tif")
manglar_inv <- rast("./features/24/Manglares INVEMAR.tif")
humedales <- rast("./features/6/humedales.tif")
bosque_seco <- rast("./features/7/bosque_seco.tif")

#change 'features' based on scenario
features <- c(paramos, manglar_inv, humedales, bosque_seco, ecosistemas)

# =====================================================
# LOAD PROTECTED AREAS (1km RASTER STACK)
# =====================================================
raster_stack <- rast(PA_RASTER_STACK_PATH)
hf <- raster_stack[[1]]  # Huella 2030 conservacionista
runap <- raster_stack[[2]]  # RUNAP_23_mode
comunidades <- raster_stack[[3]]  # comunidades_mode
renta_agropecuaria <- raster_stack[[4]]  # renta_agropecuaria_mode
resguardo <- raster_stack[[5]]  # resguardo_mode
omecs <- raster_stack[[6]]  # OMECs_mode
IHEH_2022 <- raster_stack[[7]]  # IHEH_2022
coca_muertes <- raster_stack[[8]]  # coca_muertes
climate_refugia <- raster_stack[[9]]  # climate_refugia

# =====================================================
# BUILD PRIORITIZATION PROBLEM (WITH locked_in)
# =====================================================
#set relative targets for features
# Check feature order
names(features)

# Start with 0.30 for all features
targets <- rep(0.30, nlyr(features))

# Set relative target for any specific feature
targets[names(features) == "ecosistemas"] <- 0.30

# Build problem with locked_in constraints
p1 <- problem(IHEH_2022, features = features) |>
  add_min_set_objective() |>
  add_relative_targets(targets) |>
  add_locked_in_constraints(runap) |>
  #add_locked_in_constraints(omecs) |>
  #add_locked_in_constraints(renta_agropecuaria) |>
  add_locked_in_constraints(comunidades) |>
  add_locked_in_constraints(resguardo) |>
  add_binary_decisions() |>
  add_boundary_penalties(penalty = 0.001, edge_factor = 0.5) |>
  add_highs_solver(gap = 0.1, threads = 4)
print(p1)
# =====================================================
# SOLVE PROBLEM
# =====================================================
s1 <- solve(p1, force = TRUE)

# Plot map of prioritization
dev.new()
plot(
  s1,
  main = "Prioritization",
  col = c("grey90", "darkgreen"),
  axes = FALSE
)

# =====================================================
# EVALUATE RESULTS
# =====================================================

# Store evaluation summaries in a table
run_name <- "Ecos30+ESTR30+RUNAP+Comunidades_HF"  # label for this run
n_summary <- eval_n_summary(p1, s1)
cost_summary <- eval_cost_summary(p1, s1)
p1_target_coverage <- eval_target_coverage_summary(p1, s1)

# Build a one-row summary table for this run
eval_summary <- data.frame(
  run = run_name,
  n_selected = n_summary$n,
  cost = cost_summary$cost,
  pct_targets_met = mean(p1_target_coverage$met) * 100
)
print(eval_summary)

# Target coverage detail per feature
p1_target_coverage$run <- run_name
print(p1_target_coverage)

# =====================================================
# APPEND EVALUATION RESULTS TO MASTER CSVS
# =====================================================
master_eval_path <- "./Nacional_1km_solutions/master_eval_summary.csv"
master_target_path <- "./Nacional_1km_solutions/master_target_coverage.csv"

# Append eval_summary row to master CSV (create with header if first run)
if (file.exists(master_eval_path)) {
  write.table(eval_summary, file = master_eval_path,
              sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
} else {
  write.csv(eval_summary, file = master_eval_path, row.names = FALSE)
}

# Append target_coverage rows to master CSV
if (file.exists(master_target_path)) {
  write.table(p1_target_coverage, file = master_target_path,
              sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
} else {
  write.csv(p1_target_coverage, file = master_target_path, row.names = FALSE)
}

cat("Evaluation tables appended to master CSVs.\n")

# save solution raster
writeRaster(s1,
  filename = "./Nacional_1km_solutions/Ecos30+ESTR30+RUNAP+Comunidades_HF.tif",
  overwrite = FALSE
)



# =====================================================
# IGNORE BELOW: TESTING CODE FOR LOADING RASTERS, PREPARING STACKS, AND OTHER DATA EXPLORATION
# =====================================================
pozos <- rast("C:/Users/kevinramos/Year1/Colombia_Prioritization/mesa_prioridades/mesa_prioridades_kevin/Taller_Expertos/original_data/rasters/pozos.tif")

minerales <- rast("C:/Users/kevinramos/Year1/Colombia_Prioritization/mesa_prioridades/mesa_prioridades_kevin/Taller_Expertos/original_data/rasters/minerales_transicion.tif")

other_stack <- rast("C:/Users/kevinramos/Year1/courses/pstats_231_ML/Final/cost_constraints_stack_1km.tif")

coberturas <- rast("C:/Users/kevinramos/Year1/Colombia_Prioritization/mesa_prioridades/mesa_prioridades_kevin/Taller_Expertos/original_data/rasters/coberturas.tif")

climate_refugia <- rast("C:/Users/kevinramos/Year1/Colombia_Prioritization/mesa_prioridades/mesa_prioridades_kevin/costs_and_constraints/cambio_climatico/Refugios_2030/riqueza_total_2021-2040_ssp585_10.tif")

# =====================================================
# MORE TEST CODE
# =====================================================


#import total_rent
renta_agropecuaria <- rast("C:/Users/kevinramos/Year1/Colombia_Prioritization/mesa_prioridades/mesa_prioridades_kevin/Taller_Expertos/original_data/rasters/total_rent_km2_norm.tif")
plot(renta_agropecuaria, main = "Total Rent", axes = FALSE)

names(raster_stack)

raster_stack_will <- rast("C:/Users/kevinramos/Year1/Colombia_Prioritization/mesa_prioridades/data/cost_constraints_stack_1km.tif")
names(raster_stack_will)
plot(raster_stack_will[[4]], main = "agricultural rent", axes = FALSE)

#resample normalized renta to match raster_stack_will grid before replacing
renta_agropecuaria <- resample(renta_agropecuaria, raster_stack_will, method = "bilinear")

#replace agricultural rent layer in-place (keeps layer order intact)
raster_stack_will[[4]] <- renta_agropecuaria
names(raster_stack_will)[4] <- "Renta_agropecuaria"
names(raster_stack_will)
plot(raster_stack_will[[1]], main = "HF_2030", axes = FALSE)
plot(raster_stack_will[[2]], main = "RUNAP", axes = FALSE)
plot(raster_stack_will[[3]], main = "comunidades", axes = FALSE)
plot(raster_stack_will[[4]], main = "normalized agricultural rent", axes = FALSE)
plot(raster_stack_will[[5]], main = "resguardo", axes = FALSE)
plot(raster_stack_will[[6]], main = "omecs", axes = FALSE)
plot(raster_stack_will[[7]], main = "IHEH_2022", axes = FALSE)
plot(raster_stack_will[[8]], main = "coca_muertes", axes = FALSE)
plot(raster_stack_will[[9]], main = "climate_refugia", axes = FALSE)

#write updated stack
writeRaster(raster_stack_will, filename = "./costs_and_constraints/cost_constraints_stack_1km.tif", overwrite = TRUE)

raster_stack <- rast(PA_RASTER_STACK_PATH)
hf <- raster_stack[[1]]  # Huella 2030 conservacionista
runap <- raster_stack[[2]]  # RUNAP_23_mode
comunidades <- raster_stack[[3]]  # comunidades_mode
renta_agropecuaria <- raster_stack[[4]]  # renta_agropecuaria_mode
resguardo <- raster_stack[[5]]  # resguardo_mode
omecs <- raster_stack[[6]]  # OMECs_mode
IHEH_2022 <- raster_stack[[7]]  # IHEH_2022
coca_muertes <- raster_stack[[8]]  # coca_muertes
climate_refugia <- raster_stack[[9]]  # climate_refugia

conserve_colombia <- rast("C:/Users/kevinramos/Year1/courses/pstats_231_ML/Final/Ramos_Final/conserve_colombia.tif")
names(conserve_colombia)
coberturas <- conserve_colombia[[11]]
plot(coberturas)

writeRaster(coberturas, filename = "./costs_and_constraints2/coberturas.tif", overwrite = TRUE)

#====================================================
# TEST CODE FOR MASTERS WORK
#====================================================
# Import raster stack using for masters work
masters_stack <- rast("C:/Users/kevinramos/Year1/Projects/OMECs/masters-thesis/data/conserve_colombia.tif")
names(masters_stack)

#replace current renta_agropecuaria layer in stack with new one for masters work
 masters_stack[[4]] <- renta_agropecuaria

#make pozos layer continous by cell density
density_window <- 9

# Local neighborhood count of wells.
w <- matrix(1, nrow = density_window, ncol = density_window)
pozos_neighbors <- focal(
  pozos,
  w = w,
  fun = "sum",
  na.rm = TRUE,
  na.policy = "omit",
  fillvalue = 0
)

# Scale well-present cells to 1-2.
# A well with no neighboring wells has count = 1, so it gets value 1.
max_neighbors <- global(pozos_neighbors, "max", na.rm = TRUE)[1, 1]

if (!is.na(max_neighbors) && max_neighbors > 1) {
  pozos_score_1_2 <- 1 + (pozos_neighbors - 1) / (max_neighbors - 1)
  pozos_continuous <- ifel(is.na(pozos), NA, ifel(pozos == 1, pozos_score_1_2, 0))
} else {
  pozos_continuous <- ifel(is.na(pozos), NA, ifel(pozos == 1, 1, 0))
}

names(pozos_continuous) <- "pozos"
plot(pozos_continuous, main = "Pozos Density (1-2)", axes = FALSE)

#replace pozos layer in stack with continuous version
masters_stack[[8]] <- pozos_continuous
names(masters_stack)
plot(masters_stack[[8]], main = "pozos density", axes = FALSE)

#make sure band name is called "pozos" for masters work
names(masters_stack)[8] <- "pozos"
#write updated stack for masters work
writeRaster(masters_stack, filename = "C:/Users/kevinramos/Year1/Projects/OMECs/masters-thesis/data/conserve_colombia2.tif", overwrite = TRUE)
