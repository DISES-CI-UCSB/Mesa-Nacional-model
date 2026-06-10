## script: utils
## Purpose: Functions and dataframes created for use in other scripts that can be loade

## Load/install required libraries
if (!require("pacman")) install.packages("pacman")

## Load required packages
pacman::p_load(       # automatically installs packages if needed
  tidyverse,          # always
  here,               # easier file paths
  terra,              # GIS 
  sf,                 # vector functions
  purrr)              # faster lapply

## Local directories
temp_dir <- here("data/temp_outputs")     # Store intermediate/temporary outputs
geo_dir <- here("data/model_input_lyrs")  # GeoTIFs of input layers used
ipt_dir <- here("data/model_inputs")      # Inputs directly used in prioritizr model

for (dir in c(temp_dir, geo_dir, ipt_dir)){
  if (!dir.exists(dir)) dir.create(dir)
} ; rm(dir)

# ========== TEMPLATES ==============================================
## Use the MAGNA-SIRGAS/CTM-12 as it's official projection for Colombia
my_crs <- "EPSG:9377"

## ------ Terrestrial -------------------------------------------
## Use Humboldt-produced raster as base for terrestrial Colombian extent. 
## In a geodatabase, so find the correct layer (IHEH 2022)
# info <- describe("data/costs/HEH_2022.gdb")
# print(info)


## If template hasn't yet been created, run this code. 
## Otherwise, save time and read in existing layer
if (!file.exists(file.path(geo_dir, "template_terrestre.tif"))) {
  iheh_r <- rast('OpenFileGDB:"data/costs/HEH_2022.gdb":IHEH') %>% 
    ## First put into CRS of interest, mostly preserving native resolution
    project(., my_crs, method = "bilinear") %>%
    ## Second, aggregate to get cells closer to desired resolution (1km)
    aggregate(
      fact = floor(1000 / res(.)[1]), #factor must be integer, so round down
      fun = "mean", 
      na.rm = TRUE
    ) %>% 
    ## Finally, make sure it's exactly 1km resolution
    project(., my_crs, method = "bilinear", res = 1000)
  
  ## Save template raster as updated IHEH2022
  writeRaster(iheh_r, file.path(geo_dir, "IHEH_2022.tif"), overwrite = TRUE)
  
  ## Also save as terrestrial template of one value
  template_terra <- iheh_r
  template_terra[!is.na(template_terra)] <- 1
  names(template_terra) <- "template_terrestre"
  writeRaster(template_terra, 
              file.path(geo_dir, "template_terrestre.tif"), 
              overwrite = TRUE)
  
} else {
  template_terra <- rast(file.path(geo_dir, "template_terrestre.tif"))
}


## ------ Marine -----------------------------------
## Using marine ecosystems to generate template raster
if (!file.exists(file.path(geo_dir, "template_marino.tif"))) {
  mar_v <- vect(here("data/features/Union_Profundo_Somero/Union_Profundo_Somero.shp")) %>% 
    project(., my_crs)
  
  mar_r <- rast(
    ext(mar_v),
    resolution = 1000,
    crs = my_crs)
  
  ## Just create template for now, marine ecosys layer fully processed later in 1_data_prep.R
  template_mar <- rasterize(mar_v, mar_r)
  names(template_mar) <- "template_marino"
  writeRaster(template_mar, 
              file.path(geo_dir, "template_marino.tif"), 
              overwrite = TRUE)
} else {
  template_mar <- rast(file.path(geo_dir, "template_marino.tif"))
  
}

## ------ Combined template -----------------------------------
## Combined template doesn't have values, since it will not be used
## to generate PUs for models. Just used for PAs and OMECs
combined_ext <- ext(
  min(xmin(template_terra), xmin(template_mar)), # xmin
  max(xmax(template_terra), xmax(template_mar)), # xmax
  min(ymin(template_terra), ymin(template_mar)), # ymin
  max(ymax(template_terra), ymax(template_mar))  # ymax
)

template_combined <- rast(
  ext = combined_ext,
  resolution = 1000,
  crs = my_crs
)


# ========== FUNCTIONS ==============================================
## Create fxn to rasterize solution (outputs as matrix)
rasterize_soln <- function(s, template, locked_in, ids) {
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
    layer = c("Priority area", 
              "Locked in")
  )
  
  return(rast)
}


## Get cell counts from rasterized outputs
get_freq <- function(freq_df, val) {
  x <- subset(freq_df, as.character(value) == val)$count
  if (length(x) == 0) 0 else x
}



# ========== MODEL SCENARIOS =================================================
# Create a dataframe with all model scenario permutations
## Fxn to create combos
get_combos <- function(x, min_size = 1) {
  unlist(
    lapply(min_size:length(x), function(k) combn(x, k, simplify = FALSE)),
    recursive = FALSE
  )
}

## All valid combinations of features (at least 1)
feature_combos <- get_combos(c("ecosystems", "strategic ecosystems", "species"))

## All valid combinations of includes (must include RUNAP + any combo of others)
include_combos <- get_combos(c("OMEC", "comunidades", "resguardos")) %>%
  lapply(function(x) c("RUNAP", x)) %>%   # prepend RUNAP to each
  c(list("RUNAP"))                        # add RUNAP-only option

## Lookup tables for abbreviations (for run name)
feature_abbr <- c(
  "ecosystems" = "Ecos",
  "strategic ecosystems" = "ESTR",
  "species" = "Esp")

include_abbr <- c(
  "RUNAP" = "RUNAP",
  "OMEC" = "OMEC",
  "comunidades" = "Com",
  "resguardos" = "Res")

cost_abbr <- c(
  "IHEH2022" = "IHEH",
  "net benefit" = "Agr")



## Create all permutations
scenarios_df <- expand.grid(
  target = c(17, 30),
  cost   = c("IHEH2022", 
             "net benefit"),   # not using net benefit/ag rent currently, but including for future
  KEEP.OUT.ATTRS = FALSE
) %>%
  merge(tibble(features = feature_combos)) %>%
  merge(tibble(includes = include_combos)) %>%
  mutate(
    model_name = paste0(
      map2_chr(features, target, ~ paste(paste0(feature_abbr[.x], .y), collapse = "+")),
      "+",
      map_chr(includes, ~ paste(include_abbr[.x], collapse = "+")),
      "_",
      cost_abbr[cost]
    )
  )

## Only want dataframe
rm(feature_combos); rm(include_combos)
rm(cost_abbr); rm(feature_abbr); rm(include_abbr)
rm(get_combos)

