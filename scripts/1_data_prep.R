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
  svMisc,             # progress bar for looping
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

## Human footprint ---------------------------------------------------
### Terrestrial -------------------
## IHEH 2022 already rasterized, so just vectorize and save
iheh_v <- as.matrix(rast(file.path(geo_dir, "IHEH_2022.tif")))
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
names(iheh_2030_r) <- "IHEH_2030"

## Save raster
writeRaster(iheh_2030_r, file.path(geo_dir, "IHEH_2030.tif"), overwrite = TRUE)

## Turn into matrix and save
iheh_2030_v <- as.matrix(iheh_2030_r)
iheh_2030_v[is.na(iheh_2030_v)] <- 0
saveRDS(iheh_2030_v, file.path(ipt_dir, "IHEH_2030.rds"))


### Marine --------------------------------
# Any cells that overlap with terrestrial IHEH, replace with those values.
# NOTE: Requested to do this, but may skew results for mangroves...

## Read in marine human footprint data and match CRS and resolution of template
hm_r <- rast(file.path(costs, "total 2.tif")) %>% 
  project(., my_crs, "bilinear") %>% 
  resample(template_mar, "bilinear")

## Read in terrestrial IHEH and match to marine extent
iheh_r <- rast(file.path(geo_dir, "IHEH_2022.tif")) %>% 
  resample(., template_mar)

## Return values of IHEH where they overlap
hm_cover_r <- mosaic(hm_r, iheh_r, fun = "last") %>% 
  mask(hm_r)  # Remove excess terrestrial portion

## Save raster output
writeRaster(hm_cover_r, 
            file.path(geo_dir, "huella_humana_marina.tif"),
            overwrite = TRUE)

## Convert to matrix and save
hm_v <- as.matrix(hm_cover_r)
hm_v[is.na(hm_v)] <- 0
saveRDS(hm_v, file.path(ipt_dir, "huella_humana_marina.rds"))

rm(iheh_r, hm_cover_r, hm_r)

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
writeRaster(omec_r, file.path(geo_dir, "omec.tif"), overwrite = TRUE)

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

##-------------------------- Ecosystem Services --------------------------------
### Carbon -----------------------------------
# Data from Spawn et. al 2020. 
# NOTE: Using prepped raster from Jaime. Need to figure out how this was created.

carbono_r <- rast(file.path(features, "agb_plus_bgb_spawn_2020_fixed_1km.tif")) %>% 
  project(., my_crs, method = "bilinear") %>%
  resample(., template_terra, method = "bilinear")

names(carbono_r) <- "carbono"

## Save raster
writeRaster(carbono_r, file.path(geo_dir, "carbono.tif"), overwrite = TRUE)

## Turn into matrix
carbono_v <- as.matrix(carbono_r)
carbono_v[is.na(carbono_v)] <- 0

### Freshwater -------------------------------
# Data from IDEAM 2018 ENA. Already provided as categorical GeoTiff. 
# Only want to evaluate "moderate" and "high" value areas in model. 
# NOTE: check this is correct??

## Read in raster
agua_r <- rast(
  file.path(
    features,
    "Zonas Potenciales de Recarga de Agua Subterraneas Ena2018",
    "GeoTiff/ass_h_ena2018_rcg.tif"
  ))

## See full attribute table. Terra automatically loads "cantidad"
# cats(agua_r)

## We only want the categorical values ("vals")
levels(agua_r) <- NULL

## Reproject in our CRS
agua_r <- project(agua_r, my_crs, method = "near")

## Aggregate and resample to match 1km resolution
fact <- floor(res(template_terra)[[1]]/res(agua_r)[[1]]) # Original is ~90m resolution

agua_r <- aggregate(agua_r, fact, fun = "modal") %>% 
  resample(., template_terra, method = "near")

## Only keep moderate and high levels (vals 3 and 4)
agua_binary_r <- classify(
  agua_r, 
  matrix(c(1, 0,  # muy bajo
           2, 0,  # bajo
           3, 1,  # moderado
           4, 1), # alto
         ncol = 2, byrow = TRUE))

names(agua_binary_r) <- "agua_dulce"

## Save rasters
writeRaster(agua_r, 
            file.path(temp_dir, "recarga_agua_subterranea.tif"), 
            overwrite = TRUE)
writeRaster(agua_binary_r,
            file.path(geo_dir, "recarga_agua_subterranea_moderado_alto.tif"),
            overwrite = TRUE)

## Turn into matrix
agua_v <- as.matrix(agua_binary_r)
agua_v[is.na(agua_v)] <- 0



## Compile all ecosystem services into one matrix and save
ecosys_serv_v <- cbind(carbono_v, agua_v)
ecosys_serv_v <- as(ecosys_serv_v, "dgCMatrix")
saveRDS(ecosys_serv_v, file.path(ipt_dir, "servicios_ecosistemicos.rds"))


##------------------------- Strategic Ecosystems -------------------------------
###---------------------------- Terrestrial ------------------------------------
#### Paramos -----------------------------------
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
# NOTE: In future, may have coral and sea grass also. 
# If so, aggregate into one matrix using same approach as terrestrial above

#### Mangroves ---------------------------------------
# Already rasterized and saved in `utils.R`, so just make a matrix
manglares_r <- rast(file.path(geo_dir, "manglares.tif"))

manglares_v <- as.matrix(manglares_r)
manglares_v[is.na(manglares_v)] <- 0

saveRDS(manglares_v, file.path(ipt_dir, "manglares.rds"))


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
write_csv(ecosys_df, file.path(geo_dir, "ecosistemas_IDs_IAVH_2024.csv"))


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

## Run function to determine which terrestrial ecosystems have already met targets.
eco_terra_filtered <- ecosys_coverage(ecosys_mat,
                                      targets = c(17, 30),
                                      ids, "terrestrial")

## Read in intermediate fxn product to visualize ecosystems size if useful
# ecosys_summary <- read_csv(file.path(temp_dir, "terrestrial_ecosystem_coverage.csv"))
# 
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


###----------------------------- Marine ---------------------------------------
# The shapefile is already rasterized and saved in `utils.R`
# So just read that in and create sparse matrix. 

## Get raster
ecosys_mar_r <- rast(file.path(geo_dir, "ecosistemas_marinos.tif"))

## Read in dataframe of biome ids
ecosys_mar_df <- read_csv(file.path(geo_dir, "ecosistemas_IDs_marinos.csv"))

## Convert raster in to df of values
vals <- as.data.frame(ecosys_mar_r, cells = TRUE, na.rm = TRUE)

## Create sparse matrix
ecosys_mar_mat <- sparseMatrix(
  i = vals$cell,                # each row = raster cell
  j = vals$ecosistemas_marino,  # each column = ecosystem
  x = 1,
  dims = c(ncell(ecosys_mar_r), nrow(ecosys_mar_df)),
  dimnames = list(
    NULL,
    ecosys_mar_df$biome
  )
)

## Save
saveRDS(ecosys_mar_mat, file.path(ipt_dir, "ecosistemas_marinos.rds"))


## How many ecosystems already are meeting targets under RUNAP and OMEC?
ecosys_mar_mat <- readRDS(file.path(ipt_dir, "ecosistemas_marinos.rds"))
ids <- cells(template_mar)  # which cells are PUs for marine model?

## Run function to determine which marine ecosystems have already met goals.
## Automatically saves tidy df in input_dir, but can store as object here to visualize
eco_mar_filtered <- ecosys_coverage(ecosys_mar_mat, 
                                    targets = c(30, 50), 
                                    ids, "marine")


## Read in intermediate fxn product to visualize ecosystems size if useful
# ecosys_summary <- read_csv(file.path(temp_dir, "marine_ecosystem_coverage.csv"))
# 
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
  ## Remove old IUCN status
  select(!c(iucn_status, threat_status_uicn)) %>% 
  ## Clean up some names
  rename(iucn_status = updated_status,
         mads_status = threat_status_mads) %>% 
  mutate(
    ## Make endemic status binary
    endemic = case_when(
      is.na(endemic) ~ 0,  # 0 if not endemic
      !is.na(endemic) ~ 1  # 1 if endemic
    ),
    ## Combine threat statuses
    threat_status = case_when(
      !is.na(mads_status) ~ mads_status,  # Use national (MADS) status when available
      .default = iucn_status),            # Where NA, use updated IUCN
    .before = iucn_status
    )

## save intermediate output
write_csv(spp_ranges_updated_df, file.path(temp_dir, "biomod_spp_ranges_updatedIUCN.csv"))


### ------------------------- Filter species list ------------------------------
# Finally, we can determine which species meet certain criteria and determine
# whether they need to be explicitly evaluated in the prioritizr problem. 
# Those filtered out will be evaluated post-hoc.

## Read in updated df if needed
spp_ranges_updated_df <- read_csv(file.path(temp_dir, "biomod_spp_ranges_updatedIUCN.csv"))


####---- Representativeness --------------------------------
# One method for setting species targets is general representativeness,
# or ensuring at least 17% or 30% of their range is within a conserved area.

## Set parameters
targets <- c(17, 30)                     # 17% and 30% representivity targets
conservation_types <- c("RUNAP", "OMEC") # Existing coverage by RUNAP and OMEC+RUNAP
combos <- expand_grid(                   # Create tibble of all the combinations
  target = targets,  
  conservation_type = conservation_types)

## Function that applies all filters for a given target + conservation type
filter_spp_rep <- function(target, conservation_type) {
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
    filter(!threat_status %in% c("LC", "NT")) %>%
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
spp_filtered_rep_df <- map2_dfr(combos$target, combos$conservation_type, filter_spp_rep) %>%
  ## get filename of each spp
  mutate(file_name = file.path(
    biomod_fp, 
    "presente",
    paste0(sub(" ", "_", scientific_name), "_10_MAXENT.tif")
  ))

## Save this df for filtering in prioritizr run script
write_csv(spp_filtered_rep_df, file.path(ipt_dir, "biomod_spp_filtered_representatividad.csv"))


####---- National Responsibility --------------------------------
# The other method for targets is national responsibility, which sets individual
# species targets based on their range size, threatened status, and endemism.

spp_nr_df <- spp_ranges_updated_df %>% 
  ## Filter out any species with ranges of 0km2
  filter(range_km2 > 1) %>% 
  ## At Elkin's sugggestion, remove fish for now
  ## NOTE: CHANGE THIS?? 
  filter(class != "Actinopteri") %>% 
  
  ## Set min and max values for country percentage covered by range 
  mutate(
    range_pct_country_adjusted = pmin(pmax(range_pct_country, 0.1), 90),
    .after = range_pct_country
  ) %>% 
  
  ## Also set "floor" and "ceiling" for range values to match (NOTE: maybe change in future...)
  mutate(
    range_adjusted_km2 = pmin(pmax(range_km2, country_area_km2*0.001), country_area_km2*0.9),
    .after = range_km2
  ) %>% 
  
  mutate(
    ## Create weights based on threatened status (NOTE: for now using IUCN bc not enough national data)
    w_amenaza = case_when(
      threat_status == "CR" ~ 1.0,
      threat_status == "EN" ~ 0.78,
      threat_status == "VU" ~ 0.33,
      threat_status == "NT" ~ 0.14,
      .default = 0.05 # LC, LR, and DD 
    ),
    
    ## Get the max and min ranges across all species
    min_range = min(range_adjusted_km2, na.rm = TRUE),
    max_range = max(range_adjusted_km2, na.rm = TRUE),
    
    ## Determine "range effect" for each spp (NOTE: methods for this formula??)
    range_effect = 1 - (log10(range_adjusted_km2) - log10(min_range)) / (log10(max_range) - log10(min_range)), 
    range_effect = 0.1 + range_effect * 0.9,
    
    ## Finally, get the overall national responsibility 
    responsibility = case_when(
      endemic == 1 ~ (w_amenaza * range_effect),
      endemic == 0 ~ (w_amenaza * range_effect * ((100 - range_pct_country_adjusted)/100))
    ),
    
    ## What is the area based target for each species?
    ## [NOTE: Here we need to use REAL range (not adjusted), since
    ## we can't have a target that is higher than the actually range]
    target_km2 = range_km2 * responsibility,
    
    ## Also cannot have target less than 1km2 (resolution of solution)
    target_km2 = pmax(target_km2, 1),
    
    ## So, which species already meet their target within existing conservation areas?
    target_met_runap = case_when(
      range_runap_km2 >= target_km2 ~ TRUE,
      range_runap_km2 < target_km2 ~ FALSE
    ),
    
    target_met_omec_runap = case_when(
      range_omec_runap_km2 >= target_km2 ~ TRUE,
      range_omec_runap_km2 < target_km2 ~ FALSE
    )
  ) %>% 
  
  ## Put into longer format to more easily filter
  pivot_longer(
    cols = c(target_met_runap, target_met_omec_runap),
    names_to = c(".value", "conservation_type"),
    names_pattern = "(target_met)_(runap|omec_runap)"
  ) %>% 
  mutate(
    conservation_type = recode(conservation_type,
                               "runap" = "RUNAP",
                               "omec_runap" = "OMEC")
  )

## Save
write_csv(spp_nr_df, file.path(ipt_dir, "biomod_spp_responsibilidad_nacional.csv"))

# ## Visualize responsibility spread
# thres <- spp_nr_df %>% filter(responsibility <= 0.3) #Cut off long tail to visualize easier
# 
# ggplot(thres, aes(x = responsibility)) +
#   geom_histogram(bins = 30, color = "black", fill = "steelblue4", alpha = 0.6) +
#   theme_minimal()+
#   stat_bin(
#     bins = 30,
#     geom = "text",
#     size = 3.5,
#     aes(label = after_stat(count)),
#     vjust = -0.5
#   ) +
#   labs(
#     x = "Responsibilidad",
#     y = "n especies"
#   )
# 
# ## How many spp met target?
# targets_met <- spp_nr_df %>%
#   pivot_longer(cols = target_met_runap:target_met_omec_runap, names_to = "conservation_type") %>%
#   group_by(conservation_type, value, class) %>%
#   summarize(n = n()) %>%
#   mutate(conservation_type = recode(conservation_type,
#                                       "target_met_runap" = "RUNAP",
#                                       "target_met_omec_runap"= "OMEC"),
#          value = case_when(
#            value == TRUE ~ "NR target met",
#            value == FALSE ~ "NR target not met"))
# 
# ggplot(targets_met, aes(x = class, y = n, fill = conservation_type)) +
#   geom_col(position = "dodge",
#            alpha = 0.55,
#            linewidth = 1,
#            aes(color = conservation_type)) +
#   facet_wrap(~value)+
#   geom_text(
#     aes(label = n),
#     position = position_dodge(width = 0.9),
#     vjust = -.8,
#     color = "black",
#     fontface = "bold",
#     size = 3
#   ) +
#   theme_bw() +
#   labs(y = "n especies",
#        fill = "Conservation Type",
#        color = "Conservation Type")+
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(angle = 25, vjust = 0.6),
#   )
# 
# targets_met %>%
#   filter(value == "NR target not met") %>%
# ggplot(aes(x = class, y = n, fill = conservation_type)) +
#   geom_col(position = "dodge",
#            alpha = 0.55,
#            linewidth = 1,
#            aes(color = conservation_type)) +
#   # facet_wrap(~value)+
#   geom_text(
#     aes(label = n),
#     position = position_dodge(width = 0.9),
#     vjust = -.8,
#     color = "black",
#     fontface = "bold",
#     size = 4
#   ) +
#   theme_bw() +
#   labs(y = "n especies",
#        title = "National Responsibility Not Met",
#        fill = "Conservation Type",
#        color = "Conservation Type")+
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(angle = 25, vjust = 0.6),
#   )
# 

### -------------------------- Sparse Matrices --------------------------------
# Here, we create two sparse matrices for "filtered" species for use in prioritizaiton,
# and a series of matrices for all species that can be used in other post-hoc analysis.

#### ------------------------ Representativeness -------------------------------
# One matrix created for "representativeness", using the lowest filter threshold
# (30% and just RUNAP) to return the highest number of species. 

## List of species for "filtered" matrix
spp_filtered_rep_df <- read_csv(file.path(ipt_dir, "biomod_spp_filtered_representatividad.csv"))

spp_m_list <- spp_filtered_rep_df %>% 
  filter(targets == 30,
         conservation_type == "RUNAP") %>% 
  ## Only need these variables
  select(scientific_name, class, file_name)

## Function for created matrix from filtered species list
filtered_matrix <- function(spp_m_list) {
  ## Looping through species list
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
  return(vmat)
}

## Run for representativeness
vmat_rep <- filtered_matrix(spp_m_list)

## Export
saveRDS(vmat_rep, file.path(ipt_dir, "biomod_filtered_representatividad.rds"))
rm(spp_m_list, vmat_rep)

#### ------------------------ National Responsibility --------------------------
# A second matrix filtering species based on "national responsibility" thresholds.
spp_nr_df <- read_csv(file.path(ipt_dir, "biomod_spp_responsibilidad_nacional.csv"))

spp_m_list <- spp_nr_df %>% 
  filter(target_met == FALSE,
         conservation_type == "RUNAP") %>% 
  ## Only need these variables
  select(scientific_name, class) %>% 
  mutate(file_name = file.path(biomod_fp, "presente", paste0(
    sub(" ", "_", scientific_name), "_10_MAXENT.tif")))

vmat_nr <- filtered_matrix(spp_m_list)

## Export
saveRDS(vmat_nr, file.path(ipt_dir, "biomod_filtered_responsibilidad_nacional.rds"))
rm(spp_m_list, vmat_nr)


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

