## script: Data Prep
## Purpose: Prep the input data into matrix formats to easily read into priotizr script

## Load/install required libraries
if (!require("pacman")) install.packages("pacman")

## Load required packages
pacman::p_load(       # automatically installs packages if needed
  tidyverse,          # always
  terra,              # GIS functions
  sf,                 # vector functions that play nicer w/tidyverse
  here,               # easier file paths
  svMisc,             # progress bar
  purrr,              # faster lapply
  Matrix)             # Matrices

## Local directories
temp_dir <- here("data/temp_outputs")
ipt_dir <- here("data/model_inputs")

for (dir in c(temp_dir, ipt_dir)){
  if (!dir.exists(dir)) dir.create(dir)
} ; rm(dir)

## Get functions and data
source(here("scripts/utils.R"))


#---------------------------- Costs & Constraints ------------------------------
## Human footprint
## IHEH 2022 was used as template, so just vectorize and save
iheh_v <- as.matrix(template)
saveRDS(iheh_v, file.path(ipt_dir, "IHEH_2022.rds"))


## IHEH 2030
iheh_2030_r <- rast(here("data/costs/human_footprint_2030.tif")) %>% 
  resample(template, method = "bilinear")
iheh_2030_v <- as.matrix(iheh_2030_r)
saveRDS(iheh_2030_v, file.path(ipt_dir, "IHEH_2030.rds"))


## Agricultural rent
## NOTE: this isn't the correct, normalized layer??
# renta_ag_r <- rast(here("data/costs/net_benefit.tif")) %>% 
#   resample(template, method = "bilinear")


## RUNAP
runap_r <- rast("data/includes/runap_protected_areas.tif") %>% 
  resample(template, method = "near") %>% 
  classify(matrix(ncol = 2, c(NA, 0))) %>% 
  mask(template)
## For now, just make it binary. Later can decide whether categories matter.
runap_r[runap_r > 1] <- 1
writeRaster(runap_r, file.path(temp_dir, "runap.tif"))

## Save as matrix
runap_v <- as.matrix(runap_r)
runap_v[is.na(runap_v)] <- 0
saveRDS(runap_v, file.path(ipt_dir, "runap.rds"))



## OMECs
omec_r <- rast(here("data/includes/omecs.tif")) %>% 
  resample(template, method = "near") %>% 
  classify(matrix(ncol = 2, c(NA, 0))) %>% 
  mask(template)
## For now, make binary. Later can decide whether categories matter.
omec_r[omec_r > 1] <- 1
writeRaster(omec_r, file.path(temp_dir, "omec.tif"))

## Save as matrix
omec_v <- as.matrix(omec_r)
omec_v[is.na(omec_v)] <- 0
saveRDS(omec_v, file.path(ipt_dir, "omec.rds"))



## Afro-Colombian communities
comunidades_r <- rast(here("data/includes/comunidades.tif")) %>% 
  resample(template, method = "near")
writeRaster(comunidades_r, file.path(temp_dir, "comunidades.tif"))

comunidades_v <- as.matrix(comunidades_r)
comunidades_v[is.na(comunidades_v)] <- 0
saveRDS(comunidades_v, file.path(ipt_dir, "comunidades.rds"))



## Indigenous reserves
resguardos_r <- rast(here("data/includes/resguardos.tif")) %>% 
  resample(template, method = "near")
writeRaster(resguardos_r, file.path(temp_dir, "resguardos.tif"))

resguardos_v <- as.matrix(resguardos_r)
resguardos_v[is.na(resguardos_v)] <- 0
saveRDS(resguardos_v, file.path(ipt_dir, "resguardos.rds"))



#-------------------------------- Features -------------------------------------
##------------------------- Strategic Ecosystems -------------------------------
## Paramos 
paramos_r <- rast(here("data/features/paramos.tif")) %>%
  ## match
  resample(template, method = "near") %>%
  ## reclassify
  classify(matrix(ncol = 2, c(2, 0)))
writeRaster(paramos_r, file.path(temp_dir, "paramos.tif"), overwrite = TRUE)

paramos_v <- as.matrix(paramos_r)


## Bosque seco
bosque_seco_r <- rast(here("data/features/bosque_seco.tif")) %>% 
  resample(template, method = "near") %>% 
  ## reclassify.. starts with only values of 1 and NA for some reason
  classify(matrix(ncol = 2, c(NA, 0))) %>% 
  mask(template)

writeRaster(bosque_seco_r, file.path(temp_dir, "bosque_seco.tif"), overwrite = TRUE)

bosque_seco_v <- as.matrix(bosque_seco_r)


## Mangroves
manglares_r <- rast(here("data/features/mangroves.tif")) %>% 
  ## match
  resample(template, method = "near")

writeRaster(manglares_r, file.path(temp_dir, "mangroves.tif"), overwrite = TRUE)

manglares_v <- as.matrix(manglares_r)


## Wetlands
humedales_r <- rast(here("data/features/humedales.tif")) %>% 
  ## match
  resample(template, method = "near") %>%
  ## reclassify
  classify(matrix(ncol = 2, c(2, 0)))

writeRaster(humedales_r, file.path(temp_dir, "humedales.tif"), overwrite = TRUE)

humedales_v <- as.matrix(humedales_r)

strat_ecos_v <- cbind(paramos_v, bosque_seco_v, manglares_v, humedales_v)
strat_ecos_v <- as(strat_ecos_v, "dgCMatrix")
saveRDS(strat_ecos_v, file.path(ipt_dir, "strategic_ecosystems.rds"))

##------------------------------- Ecosystems -----------------------------------
## First, convert shapefile to raster. 
## Second, turn into gridded matrix matching other inputs.

## Read in shapefile
ecosys_sf <- read_sf(
  dsn = here("data/features/Mapa_Ecosistemas_Continentales_Costeros_Marinos_100K_2024/SHAPE/e_eccmc_100K_2024.shp")) %>%
  st_transform(crs(template)) 

## Get "code list" of biomes
## Could choose several different attributes, but for now using IAVH Biomes
ecosys_df <- data.frame(
  biome = unique(ecosys_sf$bioma_IAvH)) %>% 
  filter(biome != "N.A.") %>%           # weird manual NA in dataset
  mutate(biome_id = seq_len(nrow(.)))

## Add biome codes to shapefile, then rasterize
ecosys_r <- ecosys_sf %>% 
  left_join(ecosys_df, join_by("bioma_IAvH" == "biome")) %>% 
  vect() %>%  # make terra obj
  rasterize(template, field = "biome_id") %>% 
  mask(template) # remove coastal areas

## Save the intermediate raster and biome code df
writeRaster(ecosys_r, file.path(temp_dir, "ecosistemas_IAVH_2024.tif"), overwrite = TRUE)
write_csv(ecosys_df, file.path(temp_dir, "ecosistemas_IDs_IAVH_2024.csv"))


## Some ecosystems lost after masking, so only keep those with cells left
valid_ids <- freq(ecosys_r)$value

ecosys_df_valid <- ecosys_df %>% 
  filter(biome_id %in% valid_ids)

## Conver raster in to df of values
vals <- as.data.frame(ecosys_r, cells = TRUE, na.rm = TRUE) %>% 
  filter(biome_id %in% valid_ids) %>% 
  ## sparseMatrix requires sequential numbers, so remap w/temporary variable
  mutate(biome_id_j = match(biome_id, ecosys_df_valid$biome_id))

## Create sparse matrix
ecosys_mat <- sparseMatrix(
  i = vals$cell,         # each row = raster cell
  j = vals$biome_id_j,   # each column = ecosystem
  x = 1,
  dims = c(ncell(ecosys_r), nrow(ecosys_df_valid)),
  dimnames = list(
    NULL,
    ecosys_df_valid$biome
  )
)

## Save
saveRDS(ecosys_mat, file.path(ipt_dir, "ecosistemas_IAVH_2024.rds"))


## How many ecosystems already are meeting targets under RUNAP and OMEC?
ecosys_mat <- readRDS(file.path(ipt_dir, "ecosistemas_IAVH_2024.rds"))
ids <- cells(template)

ecosys_mat <- ecosys_mat[ids, ]
ecosys_mat[is.na(ecosys_mat)] <- 0

runap <- readRDS(file.path(ipt_dir, "runap.rds"))[ids, ] == 1
runap[is.na(runap)] <- FALSE

omec <- readRDS(file.path(ipt_dir, "omec.rds"))[ids, ] == 1
omec[is.na(omec)] <- FALSE

locked_runap <- runap
locked_runap_omec <- runap | omec


# Build summary dataframe
ecosys_summary <- data.frame(
  feature        = colnames(ecosys_mat),
  total_cells    = colSums(ecosys_mat),
  pct_runap      = colSums(ecosys_mat[locked_runap, ])      / colSums(ecosys_mat) * 100,
  pct_runap_omec = colSums(ecosys_mat[locked_runap_omec, ]) / colSums(ecosys_mat) * 100
) %>%
  mutate(
    meets_17_runap      = pct_runap      >= 17,
    meets_30_runap      = pct_runap      >= 30,
    meets_17_runap_omec = pct_runap_omec >= 17,
    meets_30_runap_omec = pct_runap_omec >= 30
  )

rownames(ecosys_summary) <- NULL
## Save for reference
write_csv(ecosys_summary, file.path(temp_dir, "ecosystem_coverage.csv"))

## Now make one in tidy format so we can filter matrix in main script
targets <- c(17, 30)
conservation_types <- c("RUNAP", "OMEC")
combos <- expand_grid(target = targets, conservation_type = conservation_types)

filter_ecos <- function(target, conservation_type) {
  pct_col <- if (conservation_type == "RUNAP") "pct_runap" else "pct_runap_omec"
  
  ecosys_summary %>%
    # Remove ecosystems already meeting the target
    filter(.data[[pct_col]] < target) %>%
    mutate(
      targets = target,
      conservation_type = conservation_type
    )
}

ecosys_filtered_df <- map2_dfr(combos$target, combos$conservation_type, filter_ecos) %>% 
  select(!c(meets_17_runap:meets_30_runap_omec))

write_csv(ecosys_filtered_df, file.path(ipt_dir, "ecosys_filtered.csv"))


##------------------------------ BioModelos -----------------------------------
## Will generate one matrix that includes all potential spp
## Then can further filter (based on name) for specific 17 or 30 target runs in prioritizr script

## Path to locally stored BioModelos data
biomod_fp <- "C:/Users/nmcmanus/OneDrive - Conservation International Foundation/Documents/Projects/DISES/biomodelos/BioModelos_Data/NatGeo_NGS-86896T-21"

## Read in dataframe of spp list with updated IUCN ranges (generated in 0_biomodelos_exploration.R)
spp_ranges_updated_df <- read_csv(file.path(temp_dir, "biomod_spp_ranges_updatedIUCN.csv"))

## Set parameters
targets <- c(17, 30)
conservation_types <- c("RUNAP", "OMEC")
combos <- expand_grid(target = targets, conservation_type = conservation_types)

## Function that applies all filters for a given target + conservation type
filter_spp <- function(target, conservation_type) {
  ## Select the right coverage % column based on conservation type
  pct_col <- if (conservation_type == "RUNAP") {
    "range_pct_runap"
  } else {
    "range_pct_omec"
  }
  
  spp_ranges_updated_df %>%
    ## 1 Remove species already meeting the area target
    filter(.data[[pct_col]] < target) %>%
    ## 2: Remove species with range > 50% of country
    filter(range_pct_country < 50) %>%
    ## 3: Remove species with range under 1km
    filter(range_km2 > 1) %>%
    ## 4: remove LC and NT species
    filter(!iucn_status %in% c("LC", "NT")) %>%
    ## Tag rows with the scenario info
    mutate(
      targets = target,
      conservation_type = conservation_type,
      n_species = n()
    )
}

## Apply across all combos and stack into one tidy df
spp_filtered_df <- map2_dfr(combos$target, combos$conservation_type, filter_spp) %>%
  ## get filename of each spp
  mutate(file_name = file.path(
    biomod_fp, 
    "presente",
    paste0(sub(" ", "_", scientific_name), "_10_MAXENT.tif")
  ))

## Save this df for filtering in prioritizr run script
write_csv(spp_filtered_df, file.path(ipt_dir, "biomod_spp_ranges_filtered.csv"))


## Next, create the sparse matrix for filtered potential spp
## So filter for lowest threshold, returning most spp
spp_list <- spp_filtered_df %>% 
  filter(targets == 30,
         conservation_type == "RUNAP") %>% 
  select(scientific_name, class, file_name)


## Loop through remaining spp
for (i in 1:nrow(spp_list)){
  if(i == 1) {
    ## Read in and resample raster to match res and ext exactly (slight diff); already same CRS
    r <- rast(spp_list$file_name[1]) %>% 
      resample(., template, method = "near")
    
    ## change name to just spp
    names(r) <- spp_list$scientific_name[1]
    
    ## Get binary values, then turn into sparse matrix
    v <- values(r)
    v[is.na(v)]<-0
    vmat <- as(v,'sparseMatrix')
    
    ## Loop through the rest and add to the matrix each time
  } else {
    r <- rast(spp_list$file_name[i]) %>% 
      resample(., template, method = "near")
    names(r) <- spp_list$scientific_name[i]
    
    v2 <- values(r)
    v2[is.na(v2)]<-0
    vmat2 <- as(v2,'sparseMatrix')
    vmat <- cbind(vmat, vmat2)
    
    progress(i-1, max.value=nrow(spp_list))
  }
}

## Export
saveRDS(vmat, file.path(ipt_dir, "biomod_filtered.rds"))



## Finally, create a sparse matrix for ALL of the species data, which will be used for 
## generating metrics after running the prioritization models
spp_df <- read_csv(file.path(temp_dir, "biomod_spp_ranges_updatedIUCN.csv")) %>% 
  ## Remove fish for now
  # filter(class != "Actinopteri") %>% 
  filter(class == "Magnoliopsida") %>% #run later
  ## 3: Remove species with range under 1km
  filter(range_km2 > 1) %>%
  mutate(file_name = file.path(
    biomod_fp, 
    "presente",
    paste0(sub(" ", "_", scientific_name), "_10_MAXENT.tif")
  ))

classes <- unique(spp_df$class)

## Loop through spp and make matrices for each taxonomic class
for (current_class in classes) {
  message("Working on: ", current_class)
  class_list <- spp_df %>% filter(class == current_class)
  
  for (i in 1:nrow(class_list)){
    if(i == 1) {
      ## Read in and resample raster to match res and ext exactly (slight diff); already same CRS
      r <- rast(class_list$file_name[1]) %>% 
        resample(., template, method = "near")
      
      ## change name to just spp
      names(r) <- class_list$scientific_name[1]
      
      ## Get binary values, then turn into sparse matrix
      v <- values(r)
      v[is.na(v)]<-0
      vmat <- as(v,'sparseMatrix')
      
      ## Loop through the rest and add to the matrix each time
    } else {
      r <- rast(class_list$file_name[i]) %>% 
        resample(., template, method = "near")
      names(r) <- class_list$scientific_name[i]
      
      v2 <- values(r)
      v2[is.na(v2)]<-0
      vmat2 <- as(v2,'sparseMatrix')
      vmat <- cbind(vmat, vmat2)
      
      progress(i-1, max.value=nrow(class_list))
    }
  }
  ## Export
  saveRDS(vmat, file.path(ipt_dir, sprintf("%s.rds", current_class)))
  rm(vmat)
}


# for (i in 1:nrow(spp_df)){
#   if(i == 1) {
#     ## Read in and resample raster to match res and ext exactly (slight diff); already same CRS
#     r <- rast(spp_df$file_name[1]) %>% 
#       resample(., template, method = "near")
#     
#     ## change name to just spp
#     names(r) <- spp_df$scientific_name[1]
#     
#     ## Get binary values, then turn into sparse matrix
#     v <- values(r)
#     v[is.na(v)]<-0
#     vmat <- as(v,'sparseMatrix')
#     
#     ## Loop through the rest and add to the matrix each time
#   } else {
#     r <- rast(spp_df$file_name[i]) %>% 
#       resample(., template, method = "near")
#     names(r) <- spp_df$scientific_name[i]
#     
#     v2 <- values(r)
#     v2[is.na(v2)]<-0
#     vmat2 <- as(v2,'sparseMatrix')
#     vmat <- cbind(vmat, vmat2)
#     
#     progress(i-1, max.value=nrow(spp_df))
#   }
# }
# 
# ## Export
# saveRDS(vmat, file.path(ipt_dir, "biomod_full.rds"))

