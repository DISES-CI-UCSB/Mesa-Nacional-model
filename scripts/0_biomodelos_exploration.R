## script: BioModelos Exploration
## Purpose: Filter BioModelos species to those required for prioritization, and visualize results.

## Load/install required libraries
if (!require("pacman")) install.packages("pacman")

## Load required packages
pacman::p_load(       # automatically installs packages if needed
  tidyverse,          # always
  terra,              # GIS functions
  sf,                 # vector functions that play nicer w/tidyverse
  here,               # easier file paths
  purrr,              # iterative fxns
  exactextractr,      # better coverage data
  kableExtra,         # tables
  svMisc,             # progress bar
  Matrix)             # Matrices

## Path to locally stored BioModelos data
biomod_fp <- "C:/Users/nmcmanus/OneDrive - Conservation International Foundation/Documents/Projects/DISES/biomodelos/BioModelos_Data/NatGeo_NGS-86896T-21"

## Store outputs
temp_dir <- here("data/temp_outputs")
ipt_dir <- here("data/model_inputs")

for (dir in c(temp_dir, ipt_dir)){
  if (!dir.exists(dir)) dir.create(dir)
}

#-------------------------- BioModelos Exploration ----------------------------
## In this section, investigating how many spp already hit area-based conservation
## targets within current PAs and OMECs. Then filtering spp to see how many remain to evaluate.

## Get list of spp
spp_df <- read_csv(file.path(biomod_fp, "listas_spp_natgeo_sib_2023.csv")) %>% 
  janitor::clean_names()

## Read in RUNAP data
runap_sf <- read_sf(
  dsn = here("data/RUNAP/runap.shp"))


## Read in and join OMEC data for Colombia
omec_fpath <- here("data/WDPA_WDOECM_Nov2025_Public_COL_shp")

omec0 <- read_sf(
  file.path(
    omec_fpath,
    "WDPA_WDOECM_Nov2025_Public_COL_shp_0",
    "WDPA_WDOECM_Nov2025_Public_COL_shp-polygons.shp"))
omec1 <- read_sf(
  file.path(
    omec_fpath,
    "WDPA_WDOECM_Nov2025_Public_COL_shp_1",
    "WDPA_WDOECM_Nov2025_Public_COL_shp-polygons.shp"))
omec2 <- read_sf(
  file.path(
    omec_fpath,
    "WDPA_WDOECM_Nov2025_Public_COL_shp_2",
    "WDPA_WDOECM_Nov2025_Public_COL_shp-polygons.shp"))

omec_sf <- rbind(omec0, omec1, omec2)
rm(omec0, omec1, omec2)


## ---------------- Get species range size and current coverage ----------------
## Evaluating only present ranges for now
present_ranges <- file.path(biomod_fp, "presente")

## Update list with info we need
spp_ranges_df <- spp_df %>%
  select(scientific_name, class, endemic, threat_status_uicn, threat_status_mads) %>%
  ## file names to read in rasters easier
  mutate(file_name = paste0(sub(" ", "_", scientific_name), "_10_MAXENT.tif"))

## example raster for reprojecting
r <- rast(file.path(present_ranges, "Nothocercus_julius_10_MAXENT.tif"))

## Change sf CRS to match rasters
runap_sf <- st_transform(runap_sf, crs(r))
omec_sf <- st_transform(omec_sf, crs(r))

## Crop the vectors by the country
## Converting to terra to mask, since st_crop isn't working with weird omec geometries
outline <- as.numeric(r) 
outline[outline > -1] <- 1
outline <- as.polygons(outline) 
outline <- project(outline, crs(omec_sf))

omec_sf <- vect(omec_sf) %>% 
  mask(., outline) %>%  
  st_as_sf

runap_sf <- vect(runap_sf) %>% 
  mask(., outline) %>% 
  st_as_sf


## Get total area of country (from raster)
r[!is.na(r)] <- 1 #make everything 1
country_area_km2 <- global(cellSize(r, unit = "km")*r, "sum", na.rm = TRUE)[[1]]
rm(r)


## How much of the country's area does existing conservation cover? 
runap_km2 <- as.numeric(st_area(st_union(runap_sf))) / 1e6 #st_union bc slight overlap in some polygons
omec_km2 <- as.numeric(st_area(st_make_valid(st_union(st_make_valid(omec_sf))))) / 1e6 #problems with polygons so need to use make_valid twice to correct

runap_pct_country <- (runap_km2/country_area_km2) * 100
omec_pct_country <- (omec_km2/country_area_km2) * 100


## Create function for getting range size and RUNAP coverage
spp_ranges_fxn <- function(file_name, ...) {
  ## Read in raster
  r <- rast(file.path(present_ranges, file_name))
  
  ## Get total range area from binary raster
  ## not equal area proj, so get corrected sum of range habitat
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
    range_omec_km2 = range_omec_km2,
    range_pct_omec = round((range_omec_km2/range_km2)*100, 2)
  )
}

## Update df with function values
## NOTE: THIS WILL TAKE A WHILE! Run overnight or on separate machine
spp_ranges_df <- spp_ranges_df %>%
  mutate(results = pmap(pick(everything()), spp_ranges_fxn, .progress = TRUE)) %>% 
  unnest_wider(results) %>% 
  select(!file_name)

## save intermediate output!
write_csv(spp_ranges_df, file.path(temp_dir, "biomod_spp_ranges.csv"))


## -------------------------- Explore spp stats --------------------------------
## Read in file if needed
spp_ranges_df <- read_csv(file.path(temp_dir, "biomod_spp_ranges.csv"))

### ------1. How many spp meet area targets within existing conservation?-------
targets <- c(17, 30, 34)

## Overall summary
target_summary <- map_dfr(targets, function(t) {
  spp_ranges_df %>%
    summarise(
      target_pct = t,
      runap_n = sum(range_pct_runap >= t, na.rm = TRUE),
      runap_pct_spp = (runap_n / n()) * 100,
      omec_n = sum(range_pct_omec  >= t, na.rm = TRUE),
      omec_pct_spp = (omec_n  / n()) * 100
    )
})

## By class
target_summary_by_class <- map_dfr(targets, function(t) {
  spp_ranges_df %>%
    group_by(class) %>%
    summarise(
      target_pct = t,
      runap_n = sum(range_pct_runap >= t, na.rm = TRUE),
      runap_pct_spp = (runap_n / n()) * 100,
      omec_n = sum(range_pct_omec  >= t, na.rm = TRUE),
      omec_pct_spp = (omec_n  / n()) * 100,
      .groups = "drop"
    )
}) %>%
  arrange(class, target_pct)

## Table outputs
target_summary %>%
  mutate(
    target_pct = paste0(target_pct, "%"),
    runap_pct_spp = paste0(round(runap_pct_spp, 1), "%"),
    omec_pct_spp = paste0(round(omec_pct_spp,  1), "%")
  ) %>%
  kable(
    col.names = c("Area Target", "N Species (RUNAP)", "% Species (RUNAP)", "N Species (OMEC)", "% Species (OMEC)"),
    align = "c",
    # caption = "Species meeting area targets"
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#1B6CA8")  # header styling


target_summary_by_class %>%
  mutate(
    target_pct = paste0(target_pct, "%"),
    runap_pct_spp = paste0(round(runap_pct_spp, 1), "%"),
    omec_pct_spp = paste0(round(omec_pct_spp,  1), "%")
  ) %>%
  kable(
    col.names = c("Class", "Area Target", "N Species (RUNAP)", "% Species (RUNAP)", "N Species (OMEC)", "% Species (OMEC)"),
    align = "c",
    # caption = "Species meeting area targets"
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE
  ) %>%
  collapse_rows(
    columns = 1,
    valign = "middle"
  ) %>% 
  row_spec(0, bold = TRUE, color = "white", background = "#1B6CA8")  # header styling


## Plot outputs
target_summary_by_class %>% 
  mutate(target_pct = paste0(target_pct, "%") %>% fct_inorder()) %>%  # keep target order
  pivot_longer(
    cols = c(runap_n:omec_pct_spp),
    names_to = c("conservation_system", ".value"),
    names_sep = "_(?=pct|n)"
  ) %>% 
  mutate(conservation_system = recode(conservation_system,
                                      "runap" = "RUNAP",
                                      "omec"= "OMEC")) %>%
  ggplot(aes(x = target_pct, y = pct_spp, fill = conservation_system)) +
  geom_col(position = "dodge", linewidth = 1, alpha = 0.8, aes(color = conservation_system)) +
  facet_wrap(~ class, ncol = 2) +                          
  scale_fill_manual(values = c("RUNAP" = "#2CA25F",
                               "OMEC" = "dodgerblue3")) +
  scale_color_manual(values = c("RUNAP" = "forestgreen",
                                "OMEC" = "dodgerblue4"),
                     guide = "none") +
  labs(x = "Area Target",
       y = "% of Species",
       fill = NULL,
       title = "Percent of species meeting conservation area targets by taxonomic class") +
  geom_text(
    aes(label = n, group = conservation_system),
    position = position_dodge(0.9),
    vjust = -0.5,
    size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = c(0.75, 0.1),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold", size = 10)
  )


###-------------2. What are the range sizes of each spp? ----------------------
## How many have a range of 0? Under 1,000?
no_range <- spp_ranges_df %>% 
  group_by(class) %>% 
  mutate(class_total = n()) %>% 
  ungroup %>% 
  filter(range_km2 == 0) %>% 
  group_by(class, class_total) %>% 
  summarize(count = n(),
            pct_class = (count/class_total)*100) %>% 
  distinct()

## How many above/under 30% of country? 
n_spp_under30 <- nrow(filter(spp_ranges_df, range_pct_country < 30))
n_spp_over30 <- nrow(filter(spp_ranges_df, range_pct_country >= 30))

ggplot(data = spp_ranges_df, aes(x = range_pct_country)) +
  geom_histogram(binwidth = 5, color = "steelblue4", fill = "steelblue1", alpha = 0.8) +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  geom_vline(xintercept = 30, color = "darkslateblue", linewidth = 1, linetype = "dashed") +
  annotate("text",
           x = 15, 
           y = 1300,
           label = n_spp_under30,
           size = 4, fontface = "bold") +
  annotate("text",
           x = 70, 
           y = 1300,
           label = n_spp_over30,
           size = 4, fontface = "bold") +
  theme_minimal()+
  labs(x = "Range percent of country",
       y = "Count",
       title = "Distribution of BioModelos total spp ranges") +
  theme(
    axis.title = element_text(face = "bold")
  )


###------------------3. Break down by IUCN status? ----------------------------
## How many spp have a threatened status in MADS?
nrow(filter(spp_ranges_df, !is.na(threat_status_mads)))

## How many spp have a threatened status in IUCN?
nrow(filter(spp_ranges_df, !is.na(threat_status_uicn)))

## Where do IUCN and MADS disagree?
spp_ranges_df %>% 
  filter(!is.na(threat_status_mads),
         !is.na(threat_status_uicn),
         threat_status_mads != threat_status_uicn) %>% 
  count(threat_status_mads, threat_status_uicn, sort = TRUE) %>% 
  kable(
    align = "c"
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#1B6CA8")



## Plot current spp breakdown by class
spp_ranges_df %>%
  ## remove NAs
  # filter(!is.na(threat_status_uicn)) %>%
  group_by(class) %>%
  mutate(class_total = n()) %>%
  ungroup %>%
  group_by(class, threat_status_uicn) %>%
  summarise(count = n(), per = (count / class_total) * 100) %>%
  distinct() %>%
  ungroup() %>%
  ggplot(aes(x = threat_status_uicn, y = per)) +
  geom_col(fill = "steelblue2", color = "steelblue4") +
  facet_wrap(~ class) +
  labs(
    y = "Percent species (by class)",
    x = "IUCN Threat Status"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold", size = 10)
  )


## Fill in gaps with updated IUCN data? 
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

## How many disagree?
left_join(spp_ranges_df, iucn_df, by = "scientific_name") %>% 
  filter(!is.na(iucn_status),
         !is.na(threat_status_uicn),
         iucn_status != threat_status_uicn) %>% 
  count(iucn_status, threat_status_uicn, sort = TRUE) %>% 
  rename(new_redlist = iucn_status) %>% 
  kable(
    align = "c"
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#1B6CA8")

## How many old NAs does new list fill in?
left_join(spp_ranges_df, iucn_df, by = "scientific_name") %>% 
  filter(!is.na(iucn_status),
         is.na(threat_status_uicn)) %>% 
  count(iucn_status, sort = TRUE) %>% 
  rename(new_redlist = iucn_status) %>% 
  kable(
    align = "c"
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE
  ) %>%
  row_spec(0, bold = TRUE, color = "white", background = "#1B6CA8")


## Let's update the main spreadsheet with the new redlist statuses
spp_ranges_updated_df <- 
  left_join(spp_ranges_df, iucn_df, by = "scientific_name") %>% 
  ## Only keep old status where new list is NA
  mutate(updated_status = case_when(
    !is.na(iucn_status) ~ iucn_status,
    is.na(iucn_status) ~ threat_status_uicn
  ), .before = threat_status_uicn) %>% 
  ## Manually change some status due to taxonomic mismatch:
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
  select(!c(iucn_status, threat_status_uicn, threat_status_mads, endemic)) %>% 
  rename(iucn_status = updated_status)


## save intermediate output!
write_csv(spp_ranges_updated_df, file.path(temp_dir, "biomod_spp_ranges_updatedIUCN.csv"))


## Now visualize breakdown with updated list
spp_ranges_updated_df %>% 
  group_by(class) %>%
  mutate(class_total = n()) %>%
  ungroup %>%
  group_by(class, iucn_status) %>%
  summarise(count = n(), 
            per = (count / class_total) * 100) %>%
  distinct() %>%
  ungroup() %>%
  ggplot(aes(x = iucn_status, y = per)) +
  geom_col(fill = "steelblue2", color = "steelblue4") +
  facet_wrap(~ class) +
  labs(
    y = "Percent species (by class)",
    x = "Updated IUCN Threat Status"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold", size = 10)
  )


###---------------------- 4. Apply filters to spp list -------------------------
## read in df if neeed
spp_ranges_updated_df <- read_csv(file.path(dises, "temp_outputs/biomod_spp_ranges_updatedIUCN.csv"))

## Parameters
targets <- c(17, 30, 34)
conservation_types <- c("RUNAP", "OMEC")

## Create all combinations of target x conservation type
combos <- expand_grid(
  target = targets,
  conservation_type = conservation_types
)

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
    ## 2: Remove species with range > 30% of country
    filter(range_pct_country < 30) %>%
    ## 3: Remove species with range of 0
    filter(range_km2 > 0) %>%
    ## 4: remove LC and NT species
    filter(!iucn_status %in% c("LC", "NT")) %>%
    ## Tag rows with the scenario info
    mutate(
      targets = target,
      conservation_type = conservation_type,
      n_species = n())
}

## Apply across all combos and stack into one tidy df
spp_filtered_df <- map2_dfr(combos$target, combos$conservation_type, filter_spp)

## Save intermediate
write_csv(spp_filtered_df, file.path(temp_dir, "biomod_spp_ranges_filtered.csv"))



## Get summary
spp_filtered_df %>% 
  select(target:n_species) %>% 
  distinct() %>% 
  pivot_wider(
    names_from = conservation_type,
    values_from = n_species
  )


## Plot of results
spp_filtered_df %>%
  group_by(target, conservation_type, class) %>%
  summarize(n_species = n(), .groups = "drop") %>% 
  mutate(target = fct_inorder(as.character(target))) %>% 
  ggplot(aes(x = target, y = n_species, fill = class)) +
  geom_col() +
  geom_text(
    aes(label = n_species),
    position = position_stack(vjust = 0.5),
    color = "white",
    fontface = "bold",
    size = 3
  ) +
  # facet_wrap(~conservation_type) +
  facet_grid(class ~ conservation_type, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    x = "Area Target",
    y = "N Species",
    fill = "Class") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

## What are the ranges of left over spp?
ggplot(spp_filtered_df, aes(x = range_pct_country)) +
  geom_histogram(bins = 15, color = "steelblue4", fill = "steelblue1", alpha = 0.8) +
  facet_grid(class ~ conservation_type, scales = "free_y") +
  theme_minimal()+
  labs(x = "Range percent of country",
       y = "Count",
       title = "Distribution of remaining spp ranges (after filtered)") +
  theme(
    axis.title = element_text(face = "bold")
  )


##---------------------- Generate spp richness layers --------------------------
## We want to compare how the species richness changes between all spp
## and after filter are applied. Do areas of high biodiversity change?

###----------------- 1. Spp richness for all biomodelos ----------------------
## Function for creating richness rasters
spp_richness <- function(df, dir_out) {
  classes <- unique(df$class)
  ## First create richness for each class
  for (taxon in classes) {
    message("Processing: ", taxon)
    ## Filter for class
    class_df <- filter(df, class == taxon)
    
    ## Reach in all the rasters as a stack
    r_stack <- rast(class_df$file_name)
    
    richness <- terra::app(r_stack, "sum", na.rm = TRUE)
    
    ## Save output
    writeRaster(richness, file.path(dir_out, paste0(taxon, "_richness.tif")))
  }
  
  files <- list.files(dir_out,
                      pattern = "^.*_richness\\.tif$",
                      full.names = TRUE)
  
  richness_stack <- rast(files)
  total_richness <- terra::app(richness_stack, "sum", na.rm = TRUE)
  
  writeRaster(total_richness, file.path(dir_out,  "total_richness.tif"))
}


## First get results for all BioModelos spp
biomod_spp <- read_csv(file.path(temp_dir, "biomod_spp_ranges_updatedIUCN.csv")) %>% 
  ## get filepaths directly in df
  mutate(file_name = file.path(
    present_ranges, 
    paste0(sub(" ", "_", scientific_name), "_10_MAXENT.tif")
    ))

dir_out <- file.path(temp_dir, "total_biomod_richness")
if (!dir.exists(dir_out)) dir.create(dir_out)

spp_richness(biomod_spp, dir_out)

###---------------------- 2. Spp richness filtered ----------------------------
## Now look at filtered richness layers
biomod_filtered_spp <- read_csv(file.path(temp_dir, "biomod_spp_ranges_filtered.csv")) %>% 
  mutate(file_name = file.path(
    present_ranges, 
    paste0(sub(" ", "_", scientific_name), "_10_MAXENT.tif")
    )) %>% 
  ## MANY different ways to slice this, but for most conservative outcome,
  ## looking at 17% and 30% targets within RUNAP only
  filter(conservation_type == "RUNAP")

dir_out <- file.path(temp_dir, "filtered_biomod_richness")
if (!dir.exists(dir_out)) dir.create(dir_out)


spp_richness_filtered <- function(df, dir_out, targets) {
  classes <- unique(df$class)
  target_df <- df %>% filter(target == targets)
  
  ## First create richness for each class
  for (taxon in classes) {
    message("Processing: ", taxon)
    ## Filter for class
    class_df <- filter(target_df, class == taxon)
    
    ## Reach in all the rasters as a stack then sum
    r_stack <- rast(class_df$file_name)
    richness <- terra::app(r_stack, "sum", na.rm = TRUE)
    
    ## Save output
    writeRaster(richness, file.path(dir_out, paste0(taxon,"_",targets,"_richness.tif")))
  }
  
  files <- list.files(dir_out,
                      pattern = paste0(targets, "_richness\\.tif$"),
                      full.names = TRUE)
  
  richness_stack <- rast(files)
  total_richness <- terra::app(richness_stack, "sum", na.rm = TRUE)
  
  writeRaster(total_richness, file.path(dir_out, paste0("total_", targets, "_richness.tif")))
}

## Run for two different targets
spp_richness_filtered(biomod_filtered_spp, dir_out, target = 17)
spp_richness_filtered(biomod_filtered_spp, dir_out, target = 30)


