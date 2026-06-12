## script: Data Prep
## Purpose: Prepare all the input data into the same CRS and resolution for analysis. 
## Outputs include prepped rasters/shapefiles, as well as matrices used for prioritizr inputs. 

#-------------------------------- Set up ---------------------------------------
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
temp_dir <- here("data/temp_outputs")     # Store intermediate/temporary outputs
geo_dir <- here("data/model_input_lyrs")  # GeoTIFs of input layers used
ipt_dir <- here("data/model_inputs")      # Inputs directly used in prioritizr model

for (dir in c(temp_dir, geo_dir, ipt_dir)){
  if (!dir.exists(dir)) dir.create(dir)
} ; rm(dir)

## Get functions and data
source(here("scripts/utils.R"))


#---------------------------- Costs & Constraints ------------------------------
## File paths for these datasets
costs <- here("data/costs")
includes <- here("data/includes")

## Human footprint ---------------------------------------
## IHEH 2022 was used as template, so just vectorize and save
iheh_v <- as.matrix(template_terra)
iheh_v[is.na(iheh_v)] <- 0
saveRDS(iheh_v, file.path(ipt_dir, "IHEH_2022.rds"))


## IHEH 2030
## NOTE: need source - just pulled from Mesa Google Drive
iheh_2030_r <- rast(file.path(costs, "IHEH_2030/Huella_2030.tif")) %>% 
  ## First put into CRS of interest, mostly preserving native resolution
  project(., crs(template_terra), method = "bilinear") %>% 
  ## Second, aggregate to get cells closer to desired resolution (1km)
  aggregate(
    fact = floor(1000 / res(.)[1]), #factor must be integer, so round down
    fun = "mean", 
    na.rm = TRUE
  ) %>% 
  ## Finally, make sure it's exactly matching template
  resample(template_terra, method = "bilinear")

## Save raster
writeRaster(iheh_2030_r, file.path(geo_dir, "IHEH_2030.tif"), overwrite = TRUE)

## Turn into matrix and save
iheh_2030_v <- as.matrix(iheh_2030_r)
iheh_2030_v[is.na(iheh_2030_v)] <- 0
saveRDS(iheh_2030_v, file.path(ipt_dir, "IHEH_2030.rds"))


## Agricultural rent
## NOTE: this isn't the correct, normalized layer??
# renta_ag_r <- rast(here("data/costs/net_benefit.tif")) %>% 
#   resample(template, method = "bilinear")


## RUNAP ---------------------------------------
runap_vect <- vect(file.path(includes, "RUNAP/runap.shp")) %>% 
  project(., crs(template_combined))
writeVector(runap_vect, file.path(geo_dir, "runap.shp"))

### For now, make all categories the same (can change later)
runap_r <- rasterize(runap_vect, template_combined)
names(runap_r) <- "RUNAP"
writeRaster(runap_r, file.path(geo_dir, "runap.tif"), overwrite = TRUE)

## Save as matrix
runap_v <- as.matrix(runap_r)
runap_v[is.na(runap_v)] <- 0
saveRDS(runap_v, file.path(ipt_dir, "runap.rds"))


## OMECs ---------------------------------
omec_fp <- file.path(includes, "WDOECM_Jun2026_Public_shp")

### Read in latest polygons, filter for Colombia, and join them together
omec0_sf <- read_sf(file.path(omec_fp, "WDOECM_Jun2026_Public_shp_0/WDOECM_Jun2026_Public_shp-polygons.shp")) %>% 
  filter(ISO3 == "COL")
omec1_sf <- read_sf(file.path(omec_fp, "WDOECM_Jun2026_Public_shp_1/WDOECM_Jun2026_Public_shp-polygons.shp")) %>% 
  filter(ISO3 == "COL") #NOTE: returns 0
omec2_sf <- read_sf(file.path(omec_fp, "WDOECM_Jun2026_Public_shp_2/WDOECM_Jun2026_Public_shp-polygons.shp")) %>% 
  filter(ISO3 == "COL")

omec_vect <- rbind(omec0_sf, omec1_sf, omec2_sf) %>% 
  st_transform(., crs = crs(template_combined)) %>% 
  vect()

## Save filtered shapefile
writeVector(omec_vect, file.path(geo_dir, "omecs_col.shp"), overwrite = TRUE)

## Now rasterize
## NOTE: for now, treating all categories the same. Can change later
omec_r <- rasterize(omec_vect, template_combined)
names(omec_r) <- "OMEC"
writeRaster(omec_r, file.path(temp_dir, "omec.tif"), overwrite = TRUE)

## Save as matrix
omec_v <- as.matrix(omec_r)
omec_v[is.na(omec_v)] <- 0
saveRDS(omec_v, file.path(ipt_dir, "omec.rds"))



## Afro-Colombian communities ----------------------------------
## Read in raw shapefile, transform and save
comunidades_vect <- vect(file.path(includes, "Consejo_Comunitario_Titulado/Consejo_Comunitario_Titulado.shp")) %>% 
  project(., crs(template_terra))
writeVector(comunidades_vect, file.path(geo_dir, "comunidades.shp"), overwrite = TRUE)

## Now rasterize and save
comunidades_r <- rasterize(comunidades_vect, template_terra)
names(comunidades_r) <- "comunidades"
writeRaster(comunidades_r, file.path(geo_dir, "comunidades.tif"), overwrite = TRUE)

## Finally, turn to matrix and save
comunidades_v <- as.matrix(comunidades_r)
comunidades_v[is.na(comunidades_v)] <- 0
saveRDS(comunidades_v, file.path(ipt_dir, "comunidades.rds"))



## Indigenous reserves ---------------------------------------
resguardos_vect <- vect(file.path(includes, "Resguardo_Indigena_Formalizado/Resguardo_Indígena_Formalizado.shp")) %>% 
  project(., crs(template_terra))
writeVector(resguardos_vect, file.path(geo_dir, "resguardos.shp"), overwrite = TRUE)

resguardos_r <- rasterize(resguardos_vect, template_terra)
names(resguardos_r) <- "resguardos"
writeRaster(resguardos_r, file.path(geo_dir, "resguardos.tif"), overwrite = TRUE)

resguardos_v <- as.matrix(resguardos_r)
resguardos_v[is.na(resguardos_v)] <- 0
saveRDS(resguardos_v, file.path(ipt_dir, "resguardos.rds"))



#-------------------------------- Features -------------------------------------
## File path for all feature datasets
features <- here("data/features")

##------------------------- Strategic Ecosystems -------------------------------
###---------------------------- Terrestrial ------------------------------------
#### Paramos ---------------------------------------
## Rasterize shapefile
paramos_r <- read_sf(
  dsn = file.path(features, "Paramos/Paramo_delimitado.shp")) %>%
  st_transform(crs(template_terra)) %>%
  vect() %>%
  rasterize(., template_terra) 
names(paramos_r) <- "paramos"

## Save raster
writeRaster(paramos_r, file.path(geo_dir, "paramos.tif"), overwrite = TRUE)

## Convert to matrix
paramos_v <- as.matrix(paramos_r)
paramos_v[is.na(paramos_v)] <- 0


#### Bosque seco ----------------------------------
bosque_seco_r <- read_sf(
  file.path(features, "Bosque_Seco_Tropical/Bosque_Seco_Tropical.shp")) %>% 
  st_transform(crs(template_terra)) %>% 
  vect() %>% 
  rasterize(., template_terra)
names(bosque_seco_r) <- "bosque_seco"

writeRaster(bosque_seco_r, file.path(geo_dir, "bosque_seco.tif"), overwrite = TRUE)

## Convert to matrix
bosque_seco_v <- as.matrix(bosque_seco_r)
bosque_seco_v[is.na(bosque_seco_v)] <- 0

#### Wetlands ------------------------------------
humedales_r <- read_sf(
  file.path(features, "Humedales/Humedales.shp")) %>% 
  st_transform(crs(template_terra)) %>% 
  vect() %>% 
  rasterize(., template_terra)
names(humedales_r) <- "humedales"

writeRaster(humedales_r, file.path(geo_dir, "humedales.tif"), overwrite = TRUE)

## Convert to matrix
humedales_v <- as.matrix(humedales_r)
humedales_v[is.na(humedales_v)] <- 0


## Now combine terrestrial strategic ecosystems to one matrix
strat_ecos_ter_v <- cbind(paramos_v, bosque_seco_v, humedales_v)
strat_ecos_ter_v <- as(strat_ecos_ter_v, "dgCMatrix")
saveRDS(strat_ecos_ter_v, file.path(ipt_dir, "ecosistemas_estrategicos_terrestres.rds"))


###----------------------------- Marine ---------------------------------------
# NOTE: adding coral and seagrass to this once we get data

# #### Mangroves
# manglares_r <- read_sf(
#   file.path(features, "MANGLARES_COLOMBIA/MANGLARES_COLOMBIA.shp")) %>% 
#   st_transform(crs(template_combined)) %>% 
#   vect() %>% 
#   rasterize(., template_combined)
# 
# writeRaster(manglares_r, file.path(geo_dir, "mangroves.tif"), overwrite = TRUE)
# 
# manglares_v <- as.matrix(manglares_r)
# manglares_v[is.na(manglares_v)] <- 0
# 
# #### Corales
# 
# 
# #### Sea grass
# 
# 
# 
# ## Now combine marine strategic ecosystems to one matrix
# strat_ecos_mar_v <- cbind(manglares_v, corales_v, aq_pastos_v)
# strat_ecos_mar_v <- as(strat_ecos_mar_v, "dgCMatrix")
# saveRDS(strat_ecos_mar_v, file.path(ipt_dir, "ecosistemas_estrategicos_marinos.rds"))



##------------------------------- Ecosystems -----------------------------------
###----------------------------- Terrestrial -----------------------------------
# First, convert shapefile to raster of matching CRS & resolution. 
# Second, turn into gridded matrix matching other inputs.

## Read in shapefile
ecosys_sf <- read_sf(
  dsn = file.path(features, "Mapa_Ecosistemas_Continentales_Costeros_Marinos_100K_2024/SHAPE/e_eccmc_100K_2024.shp")) %>%
  st_transform(crs(template_terra)) 

## Get "code list" of biomes
## Many different classification levels, but for now using the "IAVH Biomes"
ecosys_df <- data.frame(
  biome = unique(ecosys_sf$bioma_IAvH)) %>% 
  filter(biome != "N.A.") %>%           # weird manual NA in dataset
  mutate(biome_id = seq_len(nrow(.)))

## Add biome codes to shapefile, then rasterize
ecosys_r <- ecosys_sf %>% 
  left_join(ecosys_df, join_by("bioma_IAvH" == "biome")) %>% 
  vect() %>%     # make terra obj
  rasterize(template_terra, field = "biome_id") %>% 
  mask(template_terra) # remove coastal areas; won't be included in PUs anyways

## Save the intermediate raster and biome code df
writeRaster(ecosys_r, file.path(geo_dir, "ecosistemas_IAVH_2024.tif"), overwrite = TRUE)
write_csv(ecosys_df, file.path(temp_dir, "ecosistemas_IDs_IAVH_2024.csv"))


## Some ecosystems won't be evaluated (outside PUs -- defined by IHEH2022)
## Only keep ecosystems that remain after masking.
## First, get list of which ecosystems remain
valid_ids <- freq(ecosys_r)$value

## Then, update df with remaining ecosystems
ecosys_df_valid <- ecosys_df %>% 
  filter(biome_id %in% valid_ids)

## Convert raster in to df of values
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
ids <- cells(template_terra)

ecosys_mat <- ecosys_mat[ids, ]
ecosys_mat[is.na(ecosys_mat)] <- 0

runap <- readRDS(file.path(ipt_dir, "runap.rds"))[ids, ] == 1
runap[is.na(runap)] <- FALSE

omec <- readRDS(file.path(ipt_dir, "omec.rds"))[ids, ] == 1
omec[is.na(omec)] <- FALSE

locked_runap <- runap
locked_runap_omec <- runap | omec


## Build summary dataframe
ecosys_summary <- data.frame(
  feature = colnames(ecosys_mat),
  total_cells = colSums(ecosys_mat),
  pct_runap = colSums(ecosys_mat[locked_runap, ]) / colSums(ecosys_mat) * 100,
  pct_runap_omec = colSums(ecosys_mat[locked_runap_omec, ]) / colSums(ecosys_mat) * 100
) %>%
  mutate(
    meets_17_runap = pct_runap >= 17,
    meets_30_runap = pct_runap >= 30,
    meets_17_runap_omec = pct_runap_omec >= 17,
    meets_30_runap_omec = pct_runap_omec >= 30
  )
rownames(ecosys_summary) <- NULL

## Save for reference
write_csv(ecosys_summary, file.path(temp_dir, "ecosystem_coverage.csv"))

# ## Visualize size of ecosystems
# ggplot(ecosys_summary, aes (x = total_cells)) +
#   geom_histogram(bins = 30, color = "black", fill = "steelblue") +
#   theme_minimal() +
#   stat_bin(
#     bins = 30,
#     geom = "text",
#     aes(label = after_stat(count)),
#     vjust = -0.5
#   )+
#   labs(x = "Ecosystem Area (km2)",
#        y = "Number of Ecosystems")


## Now make one in tidy format so we can filter matrix in main script
targets <- c(17, 30)
conservation_types <- c("RUNAP", "OMEC")
combos <- expand_grid(target = targets, 
                      conservation_type = conservation_types)

## Fxn to filter ecosystems
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

## Return tidy df using fxn
ecosys_filtered_df <- map2_dfr(combos$target, combos$conservation_type, filter_ecos) %>%  
  select(!c(meets_17_runap:meets_30_runap_omec))

## Save dataframe for use in main prioritizr script
write_csv(ecosys_filtered_df, file.path(ipt_dir, "ecosys_filtered.csv"))


###----------------------------- Marine ---------------------------------------
# Follow same approach as terrestrial ecosystems.

## Read in shapefile
ecosys_mar_sf <- read_sf(
  dsn = file.path(features, "Union_Profundo_Somero/Union_Profundo_Somero.shp")) %>% 
  st_transform(crs(template_mar)) 

## Get "code list" of marine biomes
## As directed, using the "consolidated" attribute
ecosys_mar_df <- data.frame(
  biome = unique(ecosys_mar_sf$Consolidad)) %>%  
  mutate(biome_id = seq_len(nrow(.)))

## Add biome codes to shapefile, then rasterize
ecosys_mar_r <- ecosys_mar_sf %>% 
  left_join(ecosys_mar_df, join_by("Consolidad" == "biome")) %>% 
  vect() %>%         # make terra obj
  rasterize(template_mar, field = "biome_id") %>% 
  mask(template_mar) 

## Save the intermediate raster and biome code df
writeRaster(ecosys_mar_r, 
            file.path(geo_dir, "ecosistemas_marinos.tif"), 
            overwrite = TRUE)
write_csv(ecosys_mar_df, file.path(temp_dir, "ecosistemas_IDs_marinos.csv"))


##------------------------------ BioModelos -----------------------------------
# This section contains several steps, some of which take a long time to run!
# For more explanation and visualizations behind decisions here, see the
# biomodelos_exploration.qmd

## Path to locally stored BioModelos data
biomod_fp <- "C:/Users/nmcmanus/OneDrive - Conservation International Foundation/Documents/Projects/DISES/biomodelos/BioModelos_Data/NatGeo_NGS-86896T-21"

## Get list of species
spp_df <- read_csv(file.path(biomod_fp, "listas_spp_natgeo_sib_2023.csv")) %>% 
  janitor::clean_names()


## Read in RUNAP and OMEC data (if needed) prepped above. 
## Mask to only terrestrial areas and erase any overlapping polygons to avoid double-counting
runap_sf <- vect(file.path(geo_dir, "runap.shp")) %>% 
  mask(., outline) %>%           # Read in as terra vect bc easier for masking
  st_as_sf() %>%                 # Convert to sf obj for exactextract
  select(c(objectid, geometry))  # only need polygons for this analysis; easier to join with other datasets

omec_sf <- vect(file.path(geo_dir, "omecs_col.shp")) %>% 
  mask(., outline) %>% 
  st_as_sf() %>% 
  rename(objectid = SITE_ID) %>% 
  select(c(objectid, geometry)) %>% 
  rbind(., runap_sf) %>%   # Combine with RUNAP
  st_make_valid() %>% 
  st_union()               # Join to remove overlaps

## After combining with other boundaries, now remove any overlaps from RUNAP too
runap_sf <- st_union(runap_sf)


### -------------- Species range size and conservation coverage ----------------
# To filter out species from prioritizr run, we first need to calculate
# how much of each species' range is covered by either RUNAP or RUNAP+OMEC

## Evaluating only present ranges for now
present_ranges <- file.path(biomod_fp, "presente")

## Update list with info we need
spp_ranges_df <- spp_df %>%
  select(scientific_name, class, endemic, threat_status_uicn, threat_status_mads) %>%
  ## file names to read in rasters easier
  mutate(file_name = paste0(sub(" ", "_", scientific_name), "_10_MAXENT.tif"))

## Example raster to get area of country
r <- rast(file.path(present_ranges, "Nothocercus_julius_10_MAXENT.tif")) %>% 
  project(., my_crs)
r[!is.na(r)] <- 1    # make everything value 1

country_area_km2 <- global(cellSize(r, unit = "km")*r, "sum", na.rm = TRUE)[[1]]
rm(r)


## How much of the country's area does existing conservation cover? 
runap_km2 <- as.numeric(st_area(runap_sf)) / 1e6 
omec_km2 <- as.numeric(st_area(omec_sf)) / 1e6 

runap_pct_country <- (runap_km2/country_area_km2) * 100
omec_pct_country <- (omec_km2/country_area_km2) * 100


## Create function for getting range size and RUNAP coverage
spp_ranges_fxn <- function(file_name, ...) {
  ## Read in raster
  r <- rast(file.path(present_ranges, file_name)) %>% 
    project(., my_crs)
  
  ## Get total range area from the binary raster.
  ## Not equal area proj, so get corrected sum of range habitat
  range_km2 <- global(cellSize(r, unit = "km")*r, "sum", na.rm = TRUE)[[1]]
  
  ## Get range coverage within RUNAP
  range_runap_km2 <- sum(exact_extract(r, runap_sf, "sum", coverage_area = TRUE, progress = FALSE), na.rm = TRUE) / 1e6 
  
  ## Get coverage within OMECs
  range_omec_km2 <- sum(exact_extract(r, omec_sf, "sum", coverage_area = TRUE, progress = FALSE), na.rm = TRUE) / 1e6 
  
  ## Return list of values
  list(
    range_km2 = range_km2,
    range_pct_country = round((range_km2/country_area_km2)*100,2),
    range_runap_km2 = range_runap_km2,
    range_pct_runap = round((range_runap_km2/range_km2)*100, 2),
    range_omec_runap_km2 = range_omec_km2,
    range_pct_omec_runap = round((range_omec_km2/range_km2)*100, 2)
  )
}

## Update df with function values
## NOTE: THIS WILL TAKE A WHILE! (~8hrs) Run overnight or on a virtual machine
spp_ranges_df <- spp_ranges_df %>%
  mutate(results = pmap(pick(everything()), spp_ranges_fxn, .progress = TRUE)) %>% 
  unnest_wider(results) %>% 
  select(!file_name)

## Save intermediate output
write_csv(spp_ranges_df, file.path(temp_dir, "biomod_spp_ranges.csv"))

### -------------------- Update threatened status -----------------------------
# Many species missing IUCN threatened status, which we need for proiritization methods.
# Download newest version of RedList and manually update below. 

## Read in species range df if needed
spp_ranges_df <- read_csv(file.path(temp_dir, "biomod_spp_ranges.csv"))

## Read in RedList data and match BioModelos df
iucn_df <- read_csv(here("data/redlist_species_data_20260424/assessments.csv")) %>% 
  janitor::clean_names() %>% 
  select(scientific_name, redlist_category) %>% 
  ## only keep spp matching biomodelos
  filter(scientific_name %in% spp_ranges_df$scientific_name) %>% 
  ## match notation of other df
  mutate(iucn_status = case_when(
    redlist_category == "Least Concern" ~ "LC",
    redlist_category == "Near Threatened" ~ "NT",
    redlist_category == "Vulnerable" ~ "VU",
    redlist_category == "Endangered" ~ "EN",
    redlist_category == "Data Deficient" ~ "DD",
    redlist_category == "Lower Risk/near threatened" ~ "LR/nt",
    redlist_category == "Lower risk/least concern" ~ "LR/lc",
    redlist_category == "Extinct" ~ "EX",
    redlist_category == "Critically Endangered" ~ "CR",
    redlist_category == "Extinct in the wild" ~ "EW",
  )) %>% 
  select(!redlist_category)


## Update the main df with these new statuses
spp_ranges_updated_df <-
  left_join(spp_ranges_df, iucn_df, by = "scientific_name") %>%
  ## Only keep old status where new list is NA
  mutate(
    updated_status = case_when(
      !is.na(iucn_status) ~ iucn_status,
      is.na(iucn_status) ~ threat_status_uicn
    ),
    .before = threat_status_uicn
    ) %>% 
  ## Manually change some statuses due to taxonomic mismatches (updated names):
  mutate(updated_status = case_when(
    scientific_name == "Mazama sanctaemartae" ~ "LC", #Mazama gouazoubira
    scientific_name == "Mazama zetta" ~ "DD", #Mazama americana
    scientific_name == "Odocoileus cariacou" ~ "LC", #Odocoileus virginianus
    scientific_name == "Odocoileus goudotii" ~ "LC", #Odocoileus virginianus
    scientific_name == "Dicotyles tajacu" ~ "LC", #Pecari tajacu
    scientific_name == "Puma yagouaroundi" ~ "LC", #Herpailurus yagouaroundi
    scientific_name == "Eptesicus andinus" ~ "LC", #Neoeptesicus andinus
    scientific_name == "Eptesicus brasiliensis" ~ "LC", #Neoeptesicus brasiliensis
    scientific_name == "Eptesicus furinalis" ~ "LC", #Neoeptesicus furinalis
    scientific_name == "Dasypus pastasae" ~ "LC", #match to D. kappleri
    scientific_name == "Marmosa isthmica" ~ "LC", #match to M. robinsoni
    scientific_name == "Marmosops chucha" ~ "LC", #match to M. parvidens
    scientific_name == "Metachirus myosuros" ~ "LC", #match to M. nudicaudatus
    scientific_name == "Philander melanurus" ~ "LC", #match to P. opossum
    scientific_name == "Sylvilagus salentus" ~ "DD", #S. andinus
    scientific_name == "Saguinus geoffroyi" ~ "NT", #Oedipomidas geoffroyi
    scientific_name == "Saguinus inustus" ~ "LC", #Tamarinus inustus
    scientific_name == "Saguinus leucopus" ~ "VU", #Oedipomidas leucopus
    scientific_name == "Saguinus oedipus" ~ "CR", #Oedipomidas oedipus
    scientific_name == "Cavia porcellus" ~ "LC", #domesticated spp
    scientific_name == "Nephelomys childi" ~ "LC", #match to N. meridensis
    scientific_name == "Olallamys albicaudus" ~ "DD", #O. albicauda
    scientific_name == "Hadrosciurus igniventris" ~ "LC", #Sciurus igniventris
    scientific_name == "Syntheosciurus granatensis" ~ "LC", #Sciurus granatensis
    .default = updated_status
  )) %>% 
  select(!c(iucn_status, threat_status_uicn, threat_status_mads)) %>% 
  rename(iucn_status = updated_status) %>% 
  ## Make endemic status binary
  mutate(endemic = case_when(
    is.na(endemic) ~ 0,  # 0 if not endemic
    !is.na(endemic) ~ 1  # 1 if endemic
  ))

## save intermediate output!
write_csv(spp_ranges_updated_df, file.path(temp_dir, "biomod_spp_ranges_updatedIUCN.csv"))


### ------------------------- Filter species list ------------------------------
# Finally, we can determine which species meet certain criteria and determine
# whether they need to be explicitly evaluated in the prioritizr problem. 
# Those filtered out will be evaluated post-hoc.

## Read in updated df if needed
spp_ranges_updated_df <- read_csv(file.path(temp_dir, "biomod_spp_ranges_updatedIUCN.csv"))

## Set parameters
targets <- c(17, 30)                     # 17% and 30% representivity targets
conservation_types <- c("RUNAP", "OMEC") # Existing coverage by RUNAP and OMEC+RUNAP
combos <- expand_grid(                   # Create tibble of all the combinations
  target = targets,  
  conservation_type = conservation_types)

## Function that applies all filters for a given target + conservation type
filter_spp <- function(target, conservation_type) {
  ## Select the right coverage % column based on conservation type
  pct_col <- if (conservation_type == "RUNAP") {
    "range_pct_runap"
  } else {
    "range_pct_omec_runap"
  }
  
  spp_ranges_updated_df %>%
    ## 1 Remove species already meeting the area target (17% or 30%)
    filter(.data[[pct_col]] < target) %>%
    ## 2: Remove species with range >= 50% of country
    filter(range_pct_country < 50) %>%
    ## 3: Remove species with range under 1km
    filter(range_km2 > 1) %>%
    ## 4: remove LC and NT species
    filter(!iucn_status %in% c("LC", "NT")) %>%
    ## At Elkin's sugggestion, remove fish for now
    ## NOTE: CHANGE THIS?? 
    filter(class != "Actinopteri") %>% 
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



### -------------------------- Sparse Matrices --------------------------------
# Here, we create one sparse matrix of "filtered" species for use in prioritizaiton,
# and a series of matrices for all species that can be used in other post-hoc analysis.

#### ---------------------------- Filtered ---------------------------------
# Use the lowest filter threshold, which returns the highest number of species

## List of species for "filtered" matrix
spp_m_list <- spp_filtered_df %>% 
  filter(targets == 30,
         conservation_type == "RUNAP") %>% 
  ## Only need these variables
  select(scientific_name, class, file_name)


## Create matrix by looping through species list
for (i in 1:nrow(spp_m_list)){
  ## For the first species, create the matrix
  if(i == 1) {
    ## Read in SDM and match template
    r <- rast(spp_m_list$file_name[1]) %>% 
      project(., my_crs, method = "near") %>%
      ## Already close to template resolution, so no need to aggregate
      resample(., template_terra, method = "near")
    
    ## Change name to just species
    names(r) <- spp_m_list$scientific_name[1]
    
    ## Get binary values, then turn into sparse matrix
    v <- values(r)
    v[is.na(v)]<-0
    vmat <- as(v,'sparseMatrix')
    
  ## For the rest, add to the matrix each time
  } else {
    r <- rast(spp_m_list$file_name[i]) %>%
      project(., my_crs, method = "near") %>% 
      resample(., template_terra, method = "near")
    names(r) <- spp_m_list$scientific_name[i]
    
    v2 <- values(r)
    v2[is.na(v2)]<-0
    vmat2 <- as(v2,'sparseMatrix')
    vmat <- cbind(vmat, vmat2)
    
    progress(i-1, max.value=nrow(spp_m_list))
  }
}

## Export
saveRDS(vmat, file.path(ipt_dir, "biomod_filtered.rds"))



#### ------------------------- All BioModelos -------------------------------
# Create a sparse matrix for ALL of the species data, which will be used for 
# generating metrics after running the prioritization models

## Read in list of all species
spp_df <- read_csv(file.path(temp_dir, "biomod_spp_ranges_updatedIUCN.csv"))

## Plants too large to fit into one matrix, so must split into two
plants_df <- spp_df %>% 
  filter(class == "Magnoliopsida") %>% 
  mutate(
    id = seq_len(nrow(.)),
    class = case_when(
      id <= 3000 ~ "Magnoliopsida_1",
      id > 3000 ~ "Magnoliopsida_2"
  )) %>% 
  select(!id)

## Merge back with rest
spp_df <-
  filter(spp_df, class != "Magnoliopsida") %>%
  rbind(., plants_df) %>%
  filter(range_km2 > 1) %>%
  mutate(file_name = file.path(biomod_fp, "presente", paste0(
    sub(" ", "_", scientific_name), "_10_MAXENT.tif"
  )))

classes <- unique(spp_df$class)

## Loop through species and make matrices for each taxonomic class
for (current_class in classes) {
  message("Working on: ", current_class)
  class_list <- spp_df %>% filter(class == current_class)
  
  for (i in 1:nrow(class_list)){
    if(i == 1) {
      ## Read in and resample raster to match res and ext exactly (slight diff); already same CRS
      r <- rast(class_list$file_name[1]) %>% 
        project(., my_crs, method = "near") %>% 
        resample(., template_terra, method = "near")
      
      ## change name to just spp
      names(r) <- class_list$scientific_name[1]
      
      ## Get binary values, then turn into sparse matrix
      v <- values(r)
      v[is.na(v)]<-0
      vmat <- as(v,'sparseMatrix')
      
      ## Loop through the rest and add to the matrix each time
    } else {
      r <- rast(class_list$file_name[i]) %>% 
        project(., my_crs, method = "near") %>% 
        resample(., template_terra, method = "near")
      
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

