## script: Scenario builder
## Purpose: Build a dataframe with all the possible scenarios to run

## Load/install required libraries
if (!require("pacman")) install.packages("pacman")

## Load required packages
pacman::p_load(       # automatically installs packages if needed
  tidyverse,          # always
  here,               # easier file paths
  purrr)              # faster lapply



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
  cost   = c("IHEH2022", "net benefit"),
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
rm(cost_abbr); rm(feature_abbr); rm(includes_abbr)
rm(get_combos)

