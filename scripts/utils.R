## script: utils
## Purpose: Functions and dataframes created for use in other scripts that can be loade

## Load/install required libraries
if (!require("pacman")) install.packages("pacman")

## Load required packages
pacman::p_load(       # automatically installs packages if needed
  tidyverse,          # always
  here,               # easier file paths
  purrr)              # faster lapply

## Set template for creating datasets 
template <- rast(here("data/costs/human_footprint_2022.tif"))
## NOTE: Can change CRS in this script later if neededd

## Create outline polygon if neede for masking too
mask <- as.numeric(template); mask[mask > -1] <- 1
outline <- as.polygons(mask); rm(mask)

# ========== RASTERIZE SOLUTION ==============================================
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

