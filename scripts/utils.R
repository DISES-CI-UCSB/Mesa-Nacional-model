## script: utils
## Purpose: Generates necessary directories, functions, and data used for other scripts. 
## Also loads and installs most of the required packages.


# ========== PACKAGES & DIRECTORIES ==========================================
## Load/install required libraries
if (!require("pacman")) install.packages("pacman") # Installs pacman package if needed

## Load required packages
pacman::p_load(       # automatically installs packages if needed
  tidyverse,          # always
  here,               # easier file paths
  janitor,            # cleans dataframe variables
  readxl,             # read .xls format
  terra,              # GIS functions
  sf,                 # vector functions (plays nicer with tidyverse)
  purrr)              # faster lapply-type functions


## Create local directories
temp_dir <- here("data/temp_outputs")     # Store intermediate/temporary outputs from data preparation
geo_dir <- here("data/model_input_lyrs")  # GeoTIFs and shapefiles of input layers used in models
ipt_dir <- here("data/model_inputs")      # Inputs directly used in prioritizr model (typically matrices)

base_dirs <- c(temp_dir, geo_dir, ipt_dir)

## Sub-directories based on model level and region
sub_dirs  <- c("national", 
               "sirap/eje_cafetero",
               "sirap/orinoquia") 

## Combine base and sub-directories
all_dirs <- c(
  base_dirs,
  as.vector(outer(base_dirs, sub_dirs, file.path))
)

## Create all directories locally (if needed)
for (dir in all_dirs){
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

## Remove unneeded variables from environment
rm(dir, base_dirs, all_dirs)


# ========== TEMPLATES ==============================================
# This section creates the template rasters for each region, which are used to prepare all other model inputs.
# If you'd like to change a particular model's resolution or extent, this is where to do so!

## Use MAGNA-SIRGAS/CTM-12 for all regions as it's official projection for Colombia
my_crs <- "EPSG:9377"


## ------ Terrestrial -------------------------------------------
# Use Humboldt-produced IHEH raster as base for terrestrial Colombian extent. 
# This is how planning units will be defined for the terrestrial model runs.

## Raw data in a geodatabase, so first find the correct layer (IHEH 2022)
# info <- describe("data/costs/HEH_2022.gdb")
# print(info)

## Where to save/read the template
template_path <- file.path(geo_dir, "national", "template_terrestre.tif")

## If the template hasn't yet been created, then the following section will make it. 
## Otherwise, save time and by skipping and reading in the existing template. 
if (!file.exists(template_path)) {
  ## Open raster from geodatabase
  iheh_r <- rast('OpenFileGDB:"data/costs/HEH_2022.gdb":IHEH') %>% 
    ## First put into CRS of interest, mostly preserving native resolution
    project(., my_crs, method = "bilinear") %>%
    ## Second, aggregate to get cells closer to desired resolution (1km)
    aggregate(
      fact = floor(1000 / res(.)[1]), # factor must be integer, so round down
      fun = "mean", 
      na.rm = TRUE
    ) %>% 
    ## Finally, make sure it's exactly 1km resolution
    project(., my_crs, method = "bilinear", res = 1000) %>% 
    setNames("IHEH_2022") # specify raster name in file metadata
  
  ## Save raster as one of the model costs (IHEH2022)
  writeRaster(iheh_r, 
              file.path(geo_dir, "national", "IHEH_2022.tif"), 
              overwrite = TRUE)
  
  ## Save the terrestrial template as binary version of cost layer
  template_terra <- iheh_r
  template_terra[!is.na(template_terra)] <- 1
  names(template_terra) <- "template_terrestre"
  
  writeRaster(template_terra, template_path, overwrite = TRUE)
  
} else {
  # Reads in file if already created before
  template_terra <- rast(template_path)
}

## Country outline for mapping
mask <- as.numeric(template_terra)
mask[mask > -1] <- 1
outline <- as.polygons(mask); rm(mask)


## ------ Marine -----------------------------------
# Use marine ecosystems, mangroves, and marine human footprint to generate the 
# template raster; this makes sure all planning units have a cost.
# Same as terrestrial, this template will define the PUs for marine model.

template_path <- file.path(geo_dir, "national", "template_marino.tif")

## Only run all this code if needed (the first time)
if (!file.exists(template_path)) {
  
  ### 1. Start with ecosystems -----------------------------------
  ## Read in shapefile
  mar_sf <- read_sf(
    dsn = file.path("data/features", 
                    "Union_Profundo_Somero/Union_Profundo_Somero.shp")) %>% 
    st_transform(my_crs) 
  
  ## Get "code list" of marine biomes, then update shapefile with attribute.
  ## As directed, using the "consolidad" attribute
  mar_df <- data.frame(
    biome = unique(mar_sf$Consolidad)) %>%  
    mutate(biome_id = seq_len(nrow(.)))
  
  mar_sf <- left_join(mar_sf, mar_df, join_by("Consolidad" == "biome"))
  
  ## Create initial (empty) marine template
  mar_r <- rast(
    ext(mar_sf),
    resolution = 1000,
    crs = my_crs
  )
  
  template_mar <- rasterize(vect(mar_sf), mar_r)
  
  ## Use template to rasterize marine ecosystems
  ecosys_mar_r <- mar_sf %>% 
    vect() %>%    
    rasterize(template_mar, field = "biome_id") %>% 
    mask(template_mar)
  
  ## Some small ecosystems are lost during rasterization process. 
  ## If we want to force these to be included (even if one pixel), can identify and include below. 
  present <- unique(na.omit(values(ecosys_mar_r)))
  missing_ids <- setdiff(mar_df$biome_id, present) # 4 biomes get lost
  
  ## Fix raster to include additional pixels of missing biomes (if needed)
  if (length(missing_ids) > 0) {
    message("Patching ", length(missing_ids)," missing biome(s): ",
            paste(missing_ids, collapse = ", "))
    
    ## For each missing biome, find its largest polygon and use its centroid
    missing_patches <- mar_sf %>%
      filter(biome_id %in% missing_ids) %>%
      group_by(biome_id) %>%
      slice_max(st_area(geometry), n = 1, with_ties = FALSE) %>%  # largest polygon per biome
      ungroup() %>%
      st_centroid() %>%  # centroid of that polygon
      vect()
    
    ## Stamp each centroid into the raster (overwrites whatever was there)
    ecosys_mar_r <- rasterize(missing_patches,
                              ecosys_mar_r,
                              field = "biome_id",
                              update = TRUE) %>% # preserves existing values
      setNames("ecosistemas_marino")
    
    ## Verify all biomes now present
    final_ids <- unique(na.omit(values(ecosys_mar_r)))
    still_missing <- setdiff(mar_df$biome_id, final_ids)
    
    message("Now missing ", length(still_missing), " biomes(s).")
  }
  
  ## Update marine template to include these additional pixels
  template_mar <- ecosys_mar_r
  
  ## Save the marine ecosystems raster, and marine IDs CSV
  writeRaster(ecosys_mar_r, 
              file.path(geo_dir, "national", "ecosistemas_marinos.tif"), 
              overwrite = TRUE)
  
  write_csv(mar_df, file.path(geo_dir, "national", "ecosistemas_IDs_marinos.csv"))
  
  
  ## 2. Add in Mangroves -----------------------------------
  ## Read in shapefile and rasterize using existing template
  manglares_r <- read_sf(
    file.path("data/features/MANGLARES_COLOMBIA/MANGLARES_COLOMBIA.shp")) %>%
    st_transform(crs(template_mar)) %>%
    vect() %>%
    rasterize(., template_mar) %>% 
    setNames("manglares")
  
  ## Add pixels with values to marine template, then make all values of 1
  template_mar <- cover(template_mar, manglares_r)
  template_mar[!is.na(template_mar)] <- 1  # all one value
  names(template_mar) <- "template_marino"
  
  ## Save mangroves raster
  writeRaster(manglares_r, 
              file.path(geo_dir, "national", "manglares.tif"), 
              overwrite = TRUE)
  
  
  ## 3. Mask by cost (marine footprint) -----------------------------------
  # The model cannot evaluate cells with no cost. Therefore, need to make sure
  # the final marine template doesn't extend beyond marine footprint raster. 
  
  ## Read in raster and match CRS and resolution of template
  hm_r <- rast(here("data/costs/total 2.tif")) %>% 
    project(., my_crs, method = "bilinear") %>% 
    resample(., template_mar, method = "bilinear")
  
  ## Remove template cells with no cost value
  template_mar <- mask(template_mar, hm_r)
  
  ## Save final marine template
  writeRaster(template_mar, template_path, overwrite = TRUE)
  
  
  ## Remove all the extra variables from environment
  rm(present, missing_ids, final_ids, still_missing, missing_patches, 
     mar_df, mar_r, mar_sf, ecosys_mar_r, manglares_r, hm_r)
  
} else {
  # If the template has already been made, then simply read it in
  template_mar <- rast(template_path)
}

rm(template_path)



# # NOTE: This section below is not currently used, but could be adapted in the future
# # if a combined marine + terrestrial model is desired. 
# 
# ## Get the combined extent from marine and terrestrial templates
# combined_ext <- ext(
#   min(xmin(template_terra), xmin(template_mar)), # xmin
#   max(xmax(template_terra), xmax(template_mar)), # xmax
#   min(ymin(template_terra), ymin(template_mar)), # ymin
#   max(ymax(template_terra), ymax(template_mar))  # ymax
# )
# 
# ## Generate the empty raster and correct resolution and CRS
# template_combined <- rast(
#   ext = combined_ext,
#   resolution = 1000,
#   crs = my_crs
# )
# 
# rm(combined_ext)


## -------- SIRAPs ---------------------------------------------------
# Make templates (serving as planning units) for regional SIRAP models.
# Each may differ in resolution. 

### ------- Eje Cafetero ---------------
# Using 300m resolution, roughly matching native res of IHEH2022 cost layer

template_path <- file.path(geo_dir, "sirap/eje_cafetero", "template_eje_cafetero.tif")

## If template hasn't yet been created, run following code
if (!file.exists(template_path)) {
  ## First get the boundary of region from the provided SIRAP shapefile data
  ec_v <- read_sf(file.path("data/sirap_actualizado/Territoriales_y_SIRAPs.shp")) %>%
    filter(Tematico == "Eje Cafetero") %>%  # Only keep EC
    st_union() %>%                          # Dissolve municipal boundaries
    vect() %>%                              # transform to a terra object to match the rasters
    project(., my_crs)                      # match projection
    
  ## Save this if needed
  writeVector(ec_v, 
              file.path(geo_dir, "sirap/eje_cafetero", "eje_cafetero.shp"),
              overwrite = TRUE)
  
  
  ## Now read in raw cost raster 
  iheh_r <- rast('OpenFileGDB:"data/costs/HEH_2022.gdb":IHEH') %>% 
    ## Match CRS, and slightly change resolution (~309m -> 300m)
    project(., my_crs, method = "bilinear", res = 300) %>% 
    ## Crop and mask to EC
    crop(., ec_v, mask = TRUE) %>% 
    setNames("IHEH_EC_2022")
  
  ## Save raster as one of the model costs (IHEH2022)
  writeRaster(iheh_r, 
              file.path(geo_dir, "sirap/eje_cafetero", "IHEH_EC_2022.tif"), 
              overwrite = TRUE)
  
  ## Save template as the binary version of the cost raster
  template_ec <- iheh_r
  template_ec[!is.na(template_ec)] <- 1
  names(template_ec) <- "template_eje_cafetero"
  
  writeRaster(template_ec, template_path, overwrite = TRUE)
  
} else {
  template_ec <- rast(template_path)
}


### ------- Orinoquia ------------------
# Using 500m resolution to balance computational load with better regional resolution

template_path <- file.path(geo_dir, "sirap/orinoquia", "template_orinoquia.tif")

## Only create template if needed
if (!file.exists(template_path)) {
  ## First get boundary of region
  ori_v <- read_sf(file.path("data/sirap_actualizado/Territoriales_y_SIRAPs.shp")) %>% 
    filter(Territoria == "DTOR") %>%  # Only keep Orinoquia
    st_union() %>%                    # Dissolve municipal boundaries
    vect() %>%                        # terra obj to match rasters
    project(., my_crs)                # match projection
  
  ## Save this if needed
  writeVector(ori_v, 
              file.path(geo_dir, "sirap/orinoquia", "orinoquia.shp"),
              overwrite = TRUE)

  
  ## Now read in raw cost raster 
  iheh_r <- rast('OpenFileGDB:"data/costs/HEH_2022.gdb":IHEH') %>% 
    ## Match CRS and change resolution to 500m
    project(., my_crs, method = "bilinear", res = 500) %>% 
    ## Crop and mask to orinoquia
    crop(., ori_v, mask = TRUE) %>% 
    setNames("IHEH_orinoquia_2022")
  
  ## Save raster as one of the model costs (IHEH2022)
  writeRaster(iheh_r, 
              file.path(geo_dir, "sirap/orinoquia", "IHEH_orinoquia_2022.tif"), 
              overwrite = TRUE)
  
  ## Save template as binary version of cost raster
  template_ori <- iheh_r
  template_ori[!is.na(template_ori)] <- 1
  names(template_ori) <- "template_orinoquia"
  
  writeRaster(template_ori, template_path, overwrite = TRUE)
  
} else {
  template_ori <- rast(template_path)
}

rm(template_path)



# ========== FUNCTIONS ==============================================
# This section creates some functions that are utilized multiple times 
# in other scripts. Feel free to update to continue streamlining scripts. 


## Function to rasterize the prioritizr solution
#' @param s This is the prioritizr solution output (as a matrix) 
#' @param template The region's template raster will be used to convert  
#'   matrix outputs into a GeoTIFF.
#' @param locked_in Matrix of all the areas "locked in" to the solution (e.g. RUNAP y/o OMECs).
#'   They are given a different value in raster to easily differentiate.
#' @param ids Numeric list of which cells in the raster are part of the solution.
#'
#' @return rast. Creates a raster of the solution.

rasterize_soln <- function(s, template, locked_in, ids) {
  ## Create empty raster from region's template (terrestrial, marine, or SIRAP)
  rast <- template
  rast[] <- NA
  
  ## Assign solution values to planning unit cells
  rast[ids] <- s
  
  ## Mark existing conservation (locked-in planning units)
  ## 1 = new cells selected; 2 = existing conservation; NA = not selected
  rast[ids[which(locked_in == 1)]] <- 2
  rast[rast == 0] <- NA   # Set 0s to NA (not selected)
  
  ## Add category labels
  levels(rast) <- data.frame(
    value = 1:2,
    layer = c("Priority area",  # NOTE: can change these later if needed
              "Locked in")
  )
  
  return(rast)
}



## Quick function that gets cell counts from rasterized solution outputs
get_freq <- function(freq_df, val) {
  x <- subset(freq_df, as.character(value) == val)$count
  if (length(x) == 0) 0 else x
}



## Get summary coverage stats for ecosystems
#' @param ecosys_m binary matrix of ecosystem values, where each column is a
#'   unique ecosystem category and each row is a cell. Values of 1 indicate a presence.
#' @param targets Numeric list of model target coverages. Defaults to 17% and 30%.
#' @param ecosys_type String. Either "terrestrial" or "marine", used to 
#'   differentiate process and naming conventions.
#' @param ids Numeric list of which cells in the raster are part of the solution.
#'
#' @return A tidy data frame with one row per feature/target/conservation-type
#'   combination that has not yet met its target. Columns: `feature` (ecosystem
#'   name), `total_cells` (cell count for that feature), `pct_runap` and
#'   `pct_runap_omec` (percent coverage under each conservation scenario),
#'   `targets` (target coverage being evaluated), and `conservation_type`
#'   ("RUNAP" or "RUNAP_OMEC"). Also writes summary and filtered CSVs to disk
#'   as a side effect.

ecosys_coverage <- function(ecosys_m, targets = c(17, 30), ids, ecosys_type) {
  ## Only eval cells in PUs
  ecosys_m <- ecosys_m[ids, ]
  ecosys_m[is.na(ecosys_m)] <- 0
  
  ## Read in RUNAP and OMEC matrices
  if (ecosys_type == "terrestrial") {
    runap <- readRDS(file.path(ipt_dir, "national", "runap_terrestres.rds"))[ids, ] == 1
    runap[is.na(runap)] <- FALSE
    
    omec <- readRDS(file.path(ipt_dir, "national", "omec_terrestres.rds"))[ids, ] == 1
    omec[is.na(omec)] <- FALSE
    
  } else if (ecosys_type == "marine") {
    runap <- readRDS(file.path(ipt_dir, "national", "runap_marinos.rds"))[ids, ] == 1
    runap[is.na(runap)] <- FALSE
    
    omec <- readRDS(file.path(ipt_dir, "national", "omec_marinos.rds"))[ids, ] == 1
    omec[is.na(omec)] <- FALSE
  }
  
  ## Define locked in areas for both conservation scenarios
  locked_runap <- runap
  locked_runap_omec <- runap | omec
  
  ## Build summary dataframe
  ecosys_summary <- data.frame(
    feature = colnames(ecosys_m),
    total_cells = colSums(ecosys_m),
    pct_runap = colSums(ecosys_m[locked_runap, ]) / colSums(ecosys_m) * 100,
    pct_runap_omec = colSums(ecosys_m[locked_runap_omec, ]) / colSums(ecosys_m) * 100
  ) 
  
  ## Dynamically add targets information to df
  for (t in targets) {
    ecosys_summary[[paste0("meets_", t, "_runap")]] <- ecosys_summary$pct_runap >= t
    ecosys_summary[[paste0("meets_", t, "_runap_omec")]] <- ecosys_summary$pct_runap_omec >= t
  }
  
  rownames(ecosys_summary) <- NULL # fix row names
  
  ## Save for reference
  write_excel_csv(ecosys_summary, 
                  file.path(temp_dir, "national", paste0(ecosys_type, "_ecosystem_coverage.csv")))
  
  
  ## Now make one in tidy format so we can filter matrix in main script
  conservation_types <- c("RUNAP", "RUNAP_OMEC")
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
  
  ## Drop all dynamically-created meets_* columns before returning
  meets_cols <- grep("^meets_", names(ecosys_summary), value = TRUE)
  
  ## Return tidy dataframe
  ecosys_filtered_df <- map2_dfr(combos$target, combos$conservation_type, filter_ecos) %>%
    select(-any_of(meets_cols))
  
  ## Save dataframe for use in main prioritizr script
  write_excel_csv(ecosys_filtered_df,
                  file.path(ipt_dir, "national", paste0(ecosys_type, "_ecosys_filtered.csv")))
  
  ## Return df if useful to examine or visualize
  return(ecosys_filtered_df)
}


# ========== MODEL SCENARIOS =================================================
# This section creates a dataframe of scenarios for each model, 
# where each column is a modeling parameter and each row is a different scenario.
# These are used to iteratively run the prioritizr models!


## ------ Terrestrial -----------------------------------
# NOTE: For now, just using excel sheet provided by Mesa Nacional de Prioridades.

scenarios_terra_df <- 
  ## Read in excel file from Mesa Nacional
  read_excel(file.path(ipt_dir, "national", "corridas_05062026.xlsx"), 
             sheet = "Hoja1", skip = 1, .name_repair = "unique_quiet") %>% 
  janitor::clean_names() %>% 
  ## Specify parameter names  to call in modeling function
  rename(
    ecos_target = umbral_3,
    strat_ecos_target = umbral_5,
    sp_rep_target = umbral_7,
    sp_rn_target = umbral_9,
    ecos_serv_target = umbral_11,
    includes = figuras_de_manejo, 
    cost = costo
  ) %>% 
  ## Combine locked-in parameter as a list within one cell
  mutate(
    includes = map(includes, ~ if (.x == "RUNAP y OMEC") c("RUNAP", "OMEC") else c("RUNAP"))
    ) %>% 
  ## Don't need these columns
  select(-starts_with("atributo"), -id) %>% 
  ## Species national responsibility targets are variable, 
  ## so correct to make these targets binary (evaluate or not)
  mutate(sp_rn_target = sp_rn_target > 0) %>% 
  ## Remove unnecessary duplicated scenarios
  distinct()

## Abbreviated naming convention
feature_abbr <- c(
  ecos_target = "Eco",         # Ecosistemas
  strat_ecos_target = "Estr",  # Ecosistemas estratégicos
  ecos_serv_target = "Serv",   # Ecosistemas servicios
  sp_rep_target = "EspRep"     # Especies (representividad)
)

## Function to build out model name (and file name) to reflect scenario parameters
build_model_name <- function(ecos_target, strat_ecos_target, sp_rep_target,
                             sp_rn_target, ecos_serv_target, includes, cost) {
  ## Only include a component in the name if it was considered in the specific scenario
  parts <- c(
    if (ecos_target != 0) paste0(feature_abbr["ecos_target"], ecos_target),
    if (strat_ecos_target != 0) paste0(feature_abbr["strat_ecos_target"], strat_ecos_target),
    if (ecos_serv_target != 0) paste0(feature_abbr["ecos_serv_target"], ecos_serv_target),
    if (sp_rep_target != 0) paste0(feature_abbr["sp_rep_target"], sp_rep_target),
    if (isTRUE(sp_rn_target)) "EspRN" # Especies (responsibilidad nacional)
  )
  feature_str <- paste(parts, collapse = "+")
  
  includes_names <- unlist(includes)
  includes_str <- if ("OMEC" %in% includes_names) {
    "RUNAP+OMEC"
  } else {
    "RUNAP"
  }
  paste0(feature_str, "+", includes_str, "_", cost)
}

## Update dataframe with the model names
scenarios_terra_df <- scenarios_terra_df %>%
  mutate(
    model_name = pmap_chr(
      list(ecos_target, strat_ecos_target, sp_rep_target,
           sp_rn_target, ecos_serv_target, includes, cost),
      build_model_name
    )
  )


## ------ Marine -----------------------------------
# No spreadsheet provided, so manually build out dataframe for modeling scenarios
# matching format of other geographic regions. 

## Abbreviations for modeling parameters
feature_abbr <- c(
  "ecosystems" = "Ecos",  # Ecosistemas
  "mangroves"  = "Mang")  # Manglares

include_abbr <- c(
  "RUNAP" = "RUNAP",
  "OMEC"  = "OMEC")

# NOTE: only one cost for now, so not currently needed
# cost_abbr <- c(
#   "huella_marina" = "HHM")

## Build out the dataframe
scenarios_mar_df <- expand.grid(
  ## Targets currently either 30% or 50% coverage
  target = c(30, 50),
  includes = list(c("RUNAP"), c("RUNAP", "OMEC"))
  ) %>%
  ## Ecosystems and mangroves both always considered, with matching targets
  mutate(
    features = list(c("ecosystems", "mangroves")),
    model_name = paste0(
      map_chr(target, ~ paste(paste0(feature_abbr, .x), collapse = "+")),
      "+",
      map_chr(includes, ~ paste(include_abbr[.x], collapse = "+")),
      "_",
      "HHM" #cost layer
    )
  )

rm(include_abbr)


## ------ Eje Cafetero -----------------------------------
# Read in spreadsheet shared by Mesa national
scenarios_ec_df <- 
  read_excel(file.path(ipt_dir, "sirap/eje_cafetero", "corridas_SIRAP_EC_16072026.xlsx"),
             .name_repair = "unique_quiet") %>% 
  janitor::clean_names() %>% 
  ## Rename variables to reflect parameters
  rename(
    strat_ecos_target = umbral_3,
    bs_target = umbral_5,
    hum_target = umbral_7,
    includes = figuras_de_manejo, 
    cost = costo
  ) %>% 
  mutate(
    ## Make locked-in areas a list
    includes = map(includes, ~ if (.x == "RUNAP y OMEC") c("RUNAP", "OMEC") else c("RUNAP")),
  ) %>% 
  ## Don't need these columns
  select(-starts_with("atributo"), -id) %>% 
  ## Remove duplicate scenarios
  distinct()

## Abbreviated naming convention
feature_abbr <- c(
  strat_ecos_target = "Estr", # Ecosistemas estratégicos
  bs_target = "Bs",           # Bosque seco
  hum_target = "HuEC"         # Humedales (de Eje Cafetero)
)

## Function to build out the model name for each scenario
build_model_name <- function(strat_ecos_target, bs_target, hum_target, includes, cost) {
  ## Only include components in name if considered
  parts <- c(
    if (strat_ecos_target != 0) paste0(feature_abbr["strat_ecos_target"], strat_ecos_target),
    if (!is.na(bs_target)) paste0(feature_abbr["bs_target"], bs_target),
    if (hum_target != 0) paste0(feature_abbr["hum_target"], hum_target)
  )
  feature_str <- paste(parts, collapse = "+")
  
  includes_names <- unlist(includes)
  includes_str <- if ("OMEC" %in% includes_names) {
    "RUNAP+OMEC"
  } else {
    "RUNAP"
  }
  
  paste0(feature_str, "+", includes_str, "_", cost)
}

## Update dataframe with model names
scenarios_ec_df <- scenarios_ec_df %>%
  mutate(
    model_name = pmap_chr(
      list(strat_ecos_target, bs_target, hum_target, includes, cost),
      build_model_name
    )
  )


## ------ Orinoquia -----------------------------------
# Read in spreadsheet shared by Mesa national
scenarios_ori_df <- 
  read_excel(file.path(ipt_dir, "sirap/orinoquia", "corridas_SIRAP_ORI_16072026.xlsx"),
             .name_repair = "unique_quiet") %>% 
  janitor::clean_names() %>% 
  rename(
    strat_ecos_target = umbral_3,
    cong_target = umbral_5,
    includes = figuras_de_manejo, 
    cost = costo
  ) %>% 
  mutate(
    ## make list
    includes = map(includes, ~ if (.x == "RUNAP y OMEC") c("RUNAP", "OMEC") else c("RUNAP")),
    ) %>% 
  ## Don't need these columns
  select(-starts_with("atributo"), -id) %>% 
  ## Cross each scenario with both savanna targets
  tidyr::expand_grid(sab_target = c(17, 30)) %>% 
  relocate(sab_target, .before = includes)


## Abbreviated naming convention
feature_abbr <- c(
  strat_ecos_target = "Estr",  # Ecosistemas estratégicos
  cong_target = "Cong",        # Congriales
  sab_target = "Sab"           # Sabanas
)

## Create model (scenario) names
build_model_name <- function(strat_ecos_target, cong_target, sab_target, includes, cost) {
  parts <- c(
    if (strat_ecos_target != 0) paste0(feature_abbr["strat_ecos_target"], strat_ecos_target),
    if (cong_target != 0) paste0(feature_abbr["cong_target"], cong_target),
    if (sab_target != 0) paste0(feature_abbr["sab_target"], sab_target)
  )
  feature_str <- paste(parts, collapse = "+")
  
  includes_names <- unlist(includes)
  includes_str <- if ("OMEC" %in% includes_names) {
    "RUNAP+OMEC"
  } else {
    "RUNAP"
  }
  
  paste0(feature_str, "+", includes_str, "_", cost)
}

## Update dataframe with names
scenarios_ori_df <- scenarios_ori_df %>%
  mutate(
    model_name = pmap_chr(
      list(strat_ecos_target, cong_target, sab_target, includes, cost),
      build_model_name
    )
  )

## Remove unneeded objects/fxns from environment
rm(feature_abbr, build_model_name)
