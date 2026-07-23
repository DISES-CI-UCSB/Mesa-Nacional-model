## script: utils
## Purpose: Functions and dataframes created for use in other scripts that can be loade

## Load/install required libraries
if (!require("pacman")) install.packages("pacman")

## Load required packages
pacman::p_load(       # automatically installs packages if needed
  tidyverse,          # always
  here,               # easier file paths
  readxl,             # read .xls format
  terra,              # GIS 
  sf,                 # vector functions
  purrr)              # faster lapply

## Create local directories
temp_dir <- here("data/temp_outputs")     # Store intermediate/temporary outputs
geo_dir <- here("data/model_input_lyrs")  # GeoTIFs of input layers used
ipt_dir <- here("data/model_inputs")      # Inputs directly used in prioritizr model

base_dirs <- c(temp_dir, geo_dir, ipt_dir)

## Sub-directories based on model level/region
sub_dirs  <- c("nacional", 
               "sirap/eje_cafetero",
               "sirap/orinoquia") 

all_dirs <- c(
  base_dirs,
  as.vector(outer(base_dirs, sub_dirs, file.path))
)

## Creates all directories (if needed)
for (dir in all_dirs){
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

## Remove unneeded variables
rm(dir, base_dirs, all_dirs)


# ========== TEMPLATES ==============================================
## Use the MAGNA-SIRGAS/CTM-12 as it's official projection for Colombia
my_crs <- "EPSG:9377"

## ------ Terrestrial -------------------------------------------
# Use Humboldt-produced raster as base for terrestrial Colombian extent. 
# This is how planning units will be defined for the terrestrial model runs.

## Data in a geodatabase, so first find the correct layer (IHEH 2022)
# info <- describe("data/costs/HEH_2022.gdb")
# print(info)


## If the template hasn't yet been created, then it will be. 
## Otherwise, save time and read in the existing template
template_path <- file.path(geo_dir, "nacional", "template_terrestre.tif")

if (!file.exists(template_path)) {
  ## Open raster from geodatabase
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
  
  names(iheh_r) <- "IHEH_2022"
  
  ## Save raster as one of the model costs (IHEH2022)
  writeRaster(iheh_r, 
              file.path(geo_dir, "nacional", "IHEH_2022.tif"), 
              overwrite = TRUE)
  
  ## Save terrestrial template as binary version
  template_terra <- iheh_r
  template_terra[!is.na(template_terra)] <- 1
  names(template_terra) <- "template_terrestre"
  writeRaster(template_terra, 
              template_path, 
              overwrite = TRUE)
  
} else {
  template_terra <- rast(template_path)
}

## Country outline for mapping
mask <- as.numeric(template_terra)
mask[mask > -1] <- 1
outline <- as.polygons(mask); rm(mask)


## ------ Marine -----------------------------------
# Use marine ecosystems, mangroves, and marine human footprint
# to generate a template raster (that makes sure all PUs have a cost).
# Same as terrestrial, this template will define the PUs for marine model.

template_path <- file.path(geo_dir, "nacional", "template_marino.tif")

## Only run all this code if needed (the first time)
if (!file.exists(template_path)) {
  
  ### 1. Start with ecosystems -----------------------------------
  ## Read in shapefile
  mar_sf <- read_sf(
    dsn = file.path("data/features", 
                    "Union_Profundo_Somero/Union_Profundo_Somero.shp")) %>% 
    st_transform(my_crs) 
  
  ## Get "code list" of marine biomes, then update shapefile with attribute
  ## As directed, using the "consolidated" attribute
  mar_df <- data.frame(
    biome = unique(mar_sf$Consolidad)) %>%  
    mutate(biome_id = seq_len(nrow(.)))
  
  mar_sf <- left_join(mar_sf, mar_df, join_by("Consolidad" == "biome"))
  
  ## Create initial marine template
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
                              update = TRUE)  # preserves existing values
    names(ecosys_mar_r) <- "ecosistemas_marino"
    
    ## Verify all biomes now present
    final_ids <- unique(na.omit(values(ecosys_mar_r)))
    still_missing <- setdiff(mar_df$biome_id, final_ids)
    
    message("Now missing ", length(still_missing), " biomes(s).")
  }
  
  ## Update marine template to include these additional pixels
  template_mar <- ecosys_mar_r
  
  ## Save the marine ecosystems raster, and marine IDs CSV
  writeRaster(ecosys_mar_r, 
              file.path(geo_dir, "nacional", "ecosistemas_marinos.tif"), 
              overwrite = TRUE)
  
  write_csv(mar_df, file.path(geo_dir, "nacional", "ecosistemas_IDs_marinos.csv"))
  
  
  ## 2. Add in Mangroves -----------------------------------
  ## Read in shapefile and rasterize using existing template
  manglares_r <- read_sf(
    file.path("data/features/MANGLARES_COLOMBIA/MANGLARES_COLOMBIA.shp")) %>%
    st_transform(crs(template_mar)) %>%
    vect() %>%
    rasterize(., template_mar)
  names(manglares_r) <- "manglares"
  
  ## Add pixels with values to marine template, then make all values of 1
  template_mar <- cover(template_mar, manglares_r)
  template_mar[!is.na(template_mar)] <- 1  # all one value
  names(template_mar) <- "template_marino"
  
  ## Save mangroves raster
  writeRaster(manglares_r, 
              file.path(geo_dir, "nacional", "manglares.tif"), 
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
  writeRaster(template_mar,
              template_path,
              overwrite = TRUE)
  
  
  ## Remove all the extra variables from environment
  rm(present, missing_ids, final_ids, still_missing, missing_patches, 
     mar_df, mar_r, mar_sf, ecosys_mar_r, manglares_r, hm_r)
  
} else {
  template_mar <- rast(template_path)
}

rm(template_path)

## ------ Combined template -----------------------------------
# # The combined template is empty (no values), since it is not used to define
# # planning units for the model. It is currently only used for data spanning
# # both marine and terrestrial (e.g. RUNAP and OMECs)
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
# Make templates (serivng as planning units) for regional SIRAP models.
# Each may differ in resolution. 

### ------- Eje Cafetero ---------------
# Using 300m resolution, roughly matching native res of IHEH2022 cost layer

template_path <- file.path(geo_dir, "sirap/eje_cafetero", "template_eje_cafetero.tif")

## If template hasn't yet been created, run
if (!file.exists(template_path)) {
  ## First get boundary of region
  ec_v <- read_sf(file.path("data/sirap_actualizado/Territoriales_y_SIRAPs.shp")) %>%
    filter(Tematico == "Eje Cafetero") %>%  # Only keep EC
    st_union() %>%                          # Dissolve municipal boundaries
    vect() %>%                              # terra obj to match rasters
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
    crop(., ec_v, mask = TRUE)
  
  names(iheh_r) <- "IHEH_EC_2022"
  
  ## Save raster as one of the model costs (IHEH2022)
  writeRaster(iheh_r, 
              file.path(geo_dir, "sirap/eje_cafetero", "IHEH_EC_2022.tif"), 
              overwrite = TRUE)
  
  ## Save template as binary version
  template_ec <- iheh_r
  template_ec[!is.na(template_ec)] <- 1
  names(template_ec) <- "template_eje_cafetero"
  writeRaster(template_ec, 
              template_path, 
              overwrite = TRUE)
  
} else {
  template_ec <- rast(template_path)
}


### ------- Orinoquia ------------------
# Using 500m resolution to balance computational load with better regional resolution

template_path <- file.path(geo_dir, "sirap/orinoquia", "template_orinoquia.tif")

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
    crop(., ori_v, mask = TRUE)
  
  names(iheh_r) <- "IHEH_orinoquia_2022"
  
  ## Save raster as one of the model costs (IHEH2022)
  writeRaster(iheh_r, 
              file.path(geo_dir, "sirap/orinoquia", "IHEH_orinoquia_2022.tif"), 
              overwrite = TRUE)
  
  ## Save template as binary version
  template_ori <- iheh_r
  template_ori[!is.na(template_ori)] <- 1
  names(template_ori) <- "template_orinoquia"
  writeRaster(template_ori, 
              template_path, 
              overwrite = TRUE)
  
} else {
  template_ori <- rast(template_path)
}

rm(template_path)



# ========== FUNCTIONS ==============================================
## Create fxn to rasterize solution (outputs as matrix)
rasterize_soln <- function(s, template, locked_in, ids) {
  ## Create output raster from template (either marine or terrestrial)
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


## Get cell counts from rasterized outputs
get_freq <- function(freq_df, val) {
  x <- subset(freq_df, as.character(value) == val)$count
  if (length(x) == 0) 0 else x
}


## Get summary coverage stats for ecosystems
ecosys_coverage <- function(ecosys_m,            # ecosystem matrix
                            targets = c(17, 30), # list of targets (defaults to 17% and 30%)
                            ids,                 # list of cells that are PUs
                            ecosys_type          # terrestrial or marine (string)
                            ) { 
  ## Only eval cells in PUs
  ecosys_m <- ecosys_m[ids, ]
  ecosys_m[is.na(ecosys_m)] <- 0
  
  ## Read in RUNAP and OMEC matrices
  if (ecosys_type == "terrestrial") {
    runap <- readRDS(file.path(ipt_dir, "nacional", "runap_terrestres.rds"))[ids, ] == 1
    runap[is.na(runap)] <- FALSE
    
    omec <- readRDS(file.path(ipt_dir, "nacional", "omec_terrestres.rds"))[ids, ] == 1
    omec[is.na(omec)] <- FALSE
    
  } else if (ecosys_type == "marine") {
    runap <- readRDS(file.path(ipt_dir, "nacional", "runap_marinos.rds"))[ids, ] == 1
    runap[is.na(runap)] <- FALSE
    
    omec <- readRDS(file.path(ipt_dir, "nacional", "omec_marinos.rds"))[ids, ] == 1
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
                  file.path(temp_dir, "nacional", paste0(ecosys_type, "_ecosystem_coverage.csv")))
  
  
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
                  file.path(ipt_dir, "nacional", paste0(ecosys_type, "_ecosys_filtered.csv")))
  
  ## Return df if useful to examine or visualize
  return(ecosys_filtered_df)
}


# ========== MODEL SCENARIOS =================================================
# Create a dataframe with parameters for each model scenario


## ------ Terrestrial -----------------------------------
# NOTE: For now, just using excel sheet provided by Mesa Nacional.
# Once confirmed, may re-introduce code to more flexibly generate dataframe.

scenarios_terra_df <- 
  ## Read in excel file from Mesa Nacional
  read_excel(file.path(ipt_dir, "nacional", "corridas_05062026.xlsx"), 
             sheet = "Hoja1", skip = 1) %>% 
  janitor::clean_names() %>% 
  rename(
    ecos_target = umbral_3,
    strat_ecos_target = umbral_5,
    sp_rep_target = umbral_7,
    sp_rn_target = umbral_9,
    ecos_serv_target = umbral_11,
    includes = figuras_de_manejo, 
    cost = costo
  ) %>% 
  mutate(
    includes = map(includes, ~ if (.x == "RUNAP y OMEC") c("RUNAP", "OMEC") else c("RUNAP"))
    ) %>% 
  ## Don't need these columns
  select(-starts_with("atributo"), -id) %>% 
  ## Species national responsibility targets are variable, 
  ## so correct to make these targets binary (evaluate or not)
  mutate(sp_rn_target = sp_rn_target > 0) %>% 
  ## Remove unnecessary scenarios
  distinct()

feature_abbr <- c(
  ecos_target = "Eco",
  strat_ecos_target = "Estr",
  ecos_serv_target = "Serv",
  sp_rep_target = "EspRep"
)

build_model_name <- function(ecos_target, strat_ecos_target, sp_rep_target,
                             sp_rn_target, ecos_serv_target, includes, cost) {
  parts <- c(
    if (ecos_target != 0) paste0(feature_abbr["ecos_target"], ecos_target),
    if (strat_ecos_target != 0) paste0(feature_abbr["strat_ecos_target"], strat_ecos_target),
    if (ecos_serv_target != 0) paste0(feature_abbr["ecos_serv_target"], ecos_serv_target),
    if (sp_rep_target != 0) paste0(feature_abbr["sp_rep_target"], sp_rep_target),
    if (isTRUE(sp_rn_target)) "EspRN"
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

scenarios_terra_df <- scenarios_terra_df %>%
  mutate(
    model_name = pmap_chr(
      list(ecos_target, strat_ecos_target, sp_rep_target,
           sp_rn_target, ecos_serv_target, includes, cost),
      build_model_name
    )
  )

rm(feature_abbr)


## ------ Marine -----------------------------------
feature_abbr <- c(
  "ecosystems" = "Ecos",
  "mangroves"  = "Mang")

include_abbr <- c(
  "RUNAP" = "RUNAP",
  "OMEC"  = "OMEC")

# cost_abbr <- c(
#   "huella_marina" = "HHM") #only one cost for now

scenarios_mar_df <- expand.grid(
  target = c(30, 50),
  includes = list(c("RUNAP"), c("RUNAP", "OMEC"))
  ) %>%
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

## Only want dataframe
rm(feature_abbr); rm(include_abbr)


## ------ Eje Cafetero -----------------------------------
# Read in spreadsheet shared by Mesa Nacional
scenarios_ec_df <- 
  read_excel(file.path(ipt_dir, "sirap/eje_cafetero", "corridas_SIRAP_EC_16072026.xlsx")) %>% 
  janitor::clean_names() %>% 
  rename(
    strat_ecos_target = umbral_3,
    bs_target = umbral_5,
    hum_target = umbral_7,
    includes = figuras_de_manejo, 
    cost = costo
  ) %>% 
  mutate(
    ## make list
    includes = map(includes, ~ if (.x == "RUNAP y OMEC") c("RUNAP", "OMEC") else c("RUNAP")),
    ## remove NAs
    bs_target = case_when(
      is.na(bs_target) ~ 0,
      .default = bs_target
    )
  ) %>% 
  ## Don't need these columns
  select(-starts_with("atributo"), -id) %>% 
  ## Remove duplicate scenarios
  distinct()


## ------ Orinoquia -----------------------------------
# Read in spreadsheet shared by Mesa Nacional
scenarios_ori_df <- 
  read_excel(file.path(ipt_dir, "sirap/orinoquia", "corridas_SIRAP_ORI_16072026.xlsx")) %>% 
  janitor::clean_names() %>% 
  rename(
    strat_ecos_target = umbral_3,
    cong_target = umbral_5,
    includes = figuras_de_manejo, 
    cost = costo
  ) %>% 
  mutate(
    ## make list
    includes = map(includes, ~ if (.x == "RUNAP y OMEC") c("RUNAP", "OMEC") else c("RUNAP"))
    ) %>% 
  ## Don't need these columns
  select(-starts_with("atributo"), -id)
