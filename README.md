**English** \| [Español](README.es.md)

# Conservation Prioritization Models

## Overview

This repository contains scripts for preparing spatial data and running systematic conservation prioritization analyses across four geographic regions in Colombia: two national-scale analyses (terrestrial and marine) and two regional analyses within priority conservation landscapes (SIRAP Eje Cafetero and SIRAP Orinoquía) — covering the full workflow from raw data preparation through to solved prioritization outputs. All models are built using the [`prioritizr`](https://prioritizr.net/) R package, formulating spatial conservation planning as mixed integer linear programs (MILPs). Problems are solved with Gurobi by default, though open-source solvers are also supported (see [Requirements](#requirements)).

**Regions currently evaluated:**

| Region | Realm | Scale | Resolution | Script |
|---------------|---------------|---------------|---------------|---------------|
| Colombia (national) | Terrestrial | National | 1 km | `2_1_national_terrestrial_model.R` |
| Colombia (national) | Marine | National | 1 km | `2_2_national_marine_model.R` |
| SIRAP Eje Cafetero | Terrestrial | Regional | 300 m | `3_1_sirap_eje_cafetero_model.R` |
| SIRAP Orinoquía | Terrestrial | Regional | 500 m | `3_2_sirap_orinoquia_model.R` |

Discrete models were built for each of the four regions, utilizing different input data, and generating different resolution solutions. The goal is to create models for each of the territorial and thematic SIRAPs that utilize locally relevant data and modeling parameters; currently, only two have been created.

<!-- TODO: add maps of each region here to visualize areas -->

## Repository structure

```         
.
└── scripts/
    ├── utils.R                              # Shared helper functions, planning unit templates, and scenario tables sourced by every script
    ├── 1_data_prep.R                        # Data cleaning, cropping, and layer preparation (run first, required by all models)
    ├── 2_1_national_terrestrial_model.R     # National terrestrial prioritization model
    ├── 2_2_national_marine_model.R          # National marine prioritization model
    ├── 3_1_sirap_eje_cafetero_model.R       # Regional model: SIRAP Eje Cafetero
    └── 3_2_sirap_orinoquia_model.R          # Regional model: SIRAP Orinoquía
```

All scripts assume the working directory is the **repository root** (calling `source("scripts/utils.R")`), and `utils.R` itself builds several other paths off that. Run scripts from an R project/`.Rproj` at the root, or `setwd()` there first.

<!-- TODO: if there end up being 1-2 more scripts, add them here (e.g. a post-processing / visualization script) -->

## Data

Input data are **not included in this repository** due to file size. All datasets must be downloaded/requested separately and placed in the expected local directory structure (below) before running `1_data_prep.R`.

| Dataset | Used for | Source / notes |
|------------------------|------------------------|------------------------|
| IHEH (Índice de Huella Humana) 2022 & projected 2030 | Cost layer (terrestrial and SIRAPs) | [Produced by Humboldt Institute](https://geonetwork.humboldt.org.co/geonetwork/srv/spa/catalog.search#/metadata/7d8f0aeb-8136-45a7-a469-f0016f618250) (IAVH, 2024. contrato de prestación de servicios No. 23-22/187-23/0017-197PS. Instituto Humboldt. Bogotá, D. C., Colombia.) |
| IHEH 2030 | Cost layer (terrestrial and SIRAPs) | [Produced by Humboldt Institute](https://reporte.humboldt.org.co/biodiversidad/2019/cap2/203/#seccion10); data provided directly by report authors. |
| Marine human footprint | Cost layer (marine) | Prepared by INVEMAR and reviewed by the Mesa Nacional. This layer combines pressures from factors contributing to biodiversity loss across the Colombian Caribbean and Pacific marine and coastal zones, using 2019–2020 baseline data at 300 m resolution. It considers nine variables: 1) coastal habitat transformation, 2) human access to the coastal zone, 3) human occupation, 4) nighttime light, 5) marine pollution, 6) industrial fishing, 7) artisanal fishing, 8) marine traffic, and 9) port activity. |
| RUNAP (Registro Único Nacional de Áreas Protegidas) | Locked-in constraint (all) | [Registro Unico Nacional AP](https://www.datos.gov.co/dataset/runap-Registro-Unico-Nacional-AP/n9kx-xwgg/about_data) (last updated April 9, 2026) |
| WDOECM — Other Effective Conservation Measures (OMECs) | Locked-in constraint (all) | [Protected Planet / WDOECM](https://www.protectedplanet.net/en/thematic-areas/oecms), filtered to Colombia (downloaded June 2026) |
| Marine ecosystems | Feature layer (marine) | Consolidated by INVEMAR, combining its ecosystem unit analyses for deep and shallow zones into a single layer; not publicly hosted. |
| Mangroves | Feature layer (marine) | [Produced by INVEMAR](https://acceso-datos-ambientales-invemar.hub.arcgis.com/datasets/INVEMAR::manglares-colombia/about) (last updated Mar 22, 2023) |
| Terrestrial and coastal ecosystems, (IAvH biomes) 2024 | Feature layer (terrestrial) | [Produced by IDEAM](https://www.ideam.gov.co/ecosistemas) — *Mapa de Ecosistemas Continentales, Costeros y Marinos (MEC)*, 100K, 2024 |
| Páramos | Feature layer (terrestrial and SIRAPs) | Produced by MADS and [sourced from SIAC](https://siac-datosabiertos-mads.hub.arcgis.com/maps/a85c9d818c9640649318a33b9fa115fd/about) (last updated Aug 27, 2026) |
| Bosque seco tropical | Feature layer (terrestrial and SIRAPs) | Produced by MADS and [sourced from SIAC](https://siac-datosabiertos-mads.hub.arcgis.com/maps/a85c9d818c9640649318a33b9fa115fd/about) (last updated Aug 27, 2026) |
| Humedales (national) | Feature layer (terrestrial and SIRAPs) | Produced by MADS and [sourced from SIAC](https://siac-datosabiertos-mads.hub.arcgis.com/maps/a85c9d818c9640649318a33b9fa115fd/about) (last updated Aug 27, 2026) |
| Humedales (Eje Cafetero specific) | Feature layer (Eje Cafetero) | Supplementary data for wetlands within Eje Cafetero were provided directly, not publicly hosted |
| Sabanas (SIRAP Orinoquía) | Feature layer (Orinoquía) | Supplementary data for savannas in Orinoquía provided directly, not publicly hosted |
| Congriales (SIRAP Orinoquía) | Feature layer (Orinoquía) | Data provided directly, not publicly hosted |
| Carbon (AGB + BGB) | Feature layer (terrestrial) | Spawn et al. 2020 — [dataset](https://doi.org/10.3334/ORNLDAAC/1763); repo currently uses a version pre-processed by a collaborator ("Jaime") — **method of that pre-processing needs documenting** |
| Freshwater recharge zones | Feature layer (terrestrial) | [Sourced by IDEAM](https://bart.ideam.gov.co/cneideam/Capasgeo/#:~:text=Zonas%5FPotenciales%5Fde%5FRecarga%5Fde%5FAgua%5FSubterraneas%5FEna2018), *Zonas Potenciales de Recarga de Aguas Subterráneas*, ENA 2018 |
| Land cover (IDEAM CLC 2022) | Reference/visualization | IDEAM, Cobertura de la Tierra 100K, 2022 — [ArcGIS Open Data portal](https://experience.arcgis.com/experience/568ddab184334f6b81a04d2fe9aac262/page/Datos-Abiertos-Geográficos-/) |
| Species distributions models (BioModelos) | Feature layers (terrestrial) | [BioModelos dataset](https://geonetwork.humboldt.org.co/geonetwork/srv/spa/catalog.search#/metadata/0a1a6bdf-3231-4a77-8031-0dc3fa40f21b) hosted by IAvH |
| IUCN Red List assessments | Species feature layer (terrestrial) and conservation target setting | [IUCN Red List](https://www.iucnredlist.org/) bulk download. Nationally determined species threat status were used where available, and data gaps were filled with IUCN Red List. |
| Afro-Colombian community lands | Reference/visualization | [Produced by ANT](https://data-agenciadetierras.opendata.arcgis.com/datasets/agenciadetierras::consejo-comunitario-titulado/about) (Agencia Nacional de Tierras) |
| Indigenous reserves | Reference/visualization | [Produced by ANT](https://data-agenciadetierras.opendata.arcgis.com/datasets/agenciadetierras::resguardo-indigena-formalizado/about) (Agencia Nacional de Tierras) |
| RAMSAR sites | Reference/visualization | [Colombian RAMSAR sites](https://siac-datosabiertos-mads.hub.arcgis.com/maps/3b63a187479543eb858028ecaa9d068b/about) produced by MADS |
| UNESCO Biosphere Reserves | Reference/visualization | [Colombian Biosphere Reserves](https://siac-datosabiertos-mads.hub.arcgis.com/maps/3b63a187479543eb858028ecaa9d068b/about) produced by MADS |
| Reservas Forestales (Ley 2ª de 1959) | Reference/visualization | [Data produced by MADS](https://siac-datosabiertos-mads.hub.arcgis.com/maps/3b63a187479543eb858028ecaa9d068b/about) |
| Zonas de Reserva Campesina Constituida | Reference/visualization | [Produced by ANT](https://data-agenciadetierras.opendata.arcgis.com/datasets/zonas-de-reserva-campesina-constituida/explore?location=3.091663%2C-71.907458%2C6) |
| ECC Eje Cafetero (*Estructura Ecológica Complementaria*) | Reference/visualization (Eje Cafetero) | Data provided directly, not publicly hosted |
| Model scenario and target definitions (`corridas_*.xlsx`) | Defines each model run's conservation features, targets, locked-in layers, and cost choice | Spreadsheets provided by Mesa Nacional — one per region: `corridas_05062026.xlsx` (national terrestrial), `corridas_SIRAP_EC_16072026.xlsx` (Eje Cafetero), `corridas_SIRAP_ORI_16072026.xlsx` (Orinoquía). The national marine scenarios are built manually in `utils.R`. These are **required inputs**, not raw GIS data — without them the model scripts have nothing to iterate over. <!-- TODO: pending partner permission, these will be pushed directly to the repo (e.g. under data/model_inputs/) so they come with `git clone` — update this row and the Setup section once that happens --> |
| AICAs / KBAs | Planned feature, not yet implemented | Data access pending as of the current script version |

### Expected local data directory structure {#expected-local-data-directory-structure}

While not part of the repo, the scripts assume this layout relative to the project root (paths come directly from `utils.R` / `1_data_prep.R`):

```         
data/
├── costs/                        # IHEH rasters (2022 & 2030), marine human footprint/
├── includes/                     # RUNAP, OMEC, community & indigenous reserve shapefiles
├── features/                     # Ecosystems, páramos, bosque seco, humedales, mangroves, carbon, freshwater recharge, Orinoquía-specific layers
├── sirap_actualizado/            # SIRAP/territorial boundary shapefiles (used to build SIRAP planning-unit templates)
├── visualization/                # Layers used only for display: ECC, ZRC, RAMSAR, biosphere reserves, reservas forestales, IDEAM CLC
├── redlist_species_data_.../     # IUCN Red List bulk download (assessments.csv)
├── temp_outputs/                 # Intermediate outputs written by scripts (auto-created)
│   ├── national/
│   └── sirap/{eje_cafetero,orinoquia}/
├── model_input_lyrs/             # Prepared GeoTIFFs/shapefiles per region (auto-created). Used for visualizing prepared model inputs.
│   ├── national/
│   └── sirap/{eje_cafetero,orinoquia}/
└── model_inputs/                 # Matrices (.rds) + scenario Excel files fed directly into prioritizr models (auto-created, except the Excel files)
    ├── national/                 #   -> also holds corridas_05062026.xlsx
    └── sirap/
        ├── eje_cafetero/         #   -> also holds corridas_SIRAP_EC_16072026.xlsx
        └── orinoquia/            #   -> also holds corridas_SIRAP_ORI_16072026.xlsx

results/                          # Model outputs (auto-created by each model script)
├── national/{terrestrial,marine}/
└── sirap/{eje_cafetero,orinoquia}/
```

`temp_outputs/`, `model_input_lyrs/`, `model_inputs/`, and `results/...` are created automatically by the scripts if missing — you only need to seed `costs/`, `includes/`, `features/`, `sirap_actualizado/`, `visualization/`, `redlist_species_data_.../`, and the three `corridas_*.xlsx` scenario files.

## Requirements

-   **R** with the [`pacman`](https://cran.r-project.org/package=pacman) package (installed automatically by the scripts if missing; used to auto-install everything else).
-   **Key R packages** (auto-installed via `pacman::p_load()` where possible): `tidyverse`, `readxl`, `terra`, `sf`, `purrr`, `Matrix`, `prioritizr`, `exactextractr`.
-   **Gurobi** (with a valid license) — the model scripts currently call `add_gurobi_solver()`. This is a commercial solver: you'll need a license (free academic licenses are available from [Gurobi](https://www.gurobi.com/academia/academic-program-and-licenses/)) and the `gurobi` R package installed from the Gurobi installation itself (it is **not** on CRAN, so `pacman`/`p_load` can't fetch it automatically — install it manually first per Gurobi's R setup instructions). `prioritizr` also supports open-source solvers (e.g. `Rsymphony`, `HiGHS` via `highs`, or `cbc` via `rcbc`) — swap out the `add_gurobi_solver()` calls in the four model scripts for one of these if you'd rather avoid the license requirement.

## Setup

1.  Clone this repository and open it as the R working directory (e.g. via an `.Rproj` at the repo root).
2.  Create the `data/` subfolders listed under [Expected local data directory structure](#expected-local-data-directory-structure) and populate them with the datasets in [Data](#data), including the three `corridas_*.xlsx` scenario files.
3.  Install Gurobi and its R package manually (if using, see [Requirements](#requirements)) — this is the one dependency the scripts can't install for you.
4.  Run `scripts/1_data_prep.R`. This is required before running **any** model script — it builds the planning-unit templates and all prepared feature/cost matrices used across all four regions. Note: this script, particularly the BioModelos species-range section, can take **many (8+) hours** to run; consider running it overnight or on a separate machine.

## Usage / run order

```         
1_data_prep.R                       # Required first — prepares data for all regions
   │
   ├── 2_1_national_terrestrial_model.R
   ├── 2_2_national_marine_model.R
   ├── 3_1_sirap_eje_cafetero_model.R
   └── 3_2_sirap_orinoquia_model.R
```

After `1_data_prep.R` has been run once, the four model scripts are independent of one another and can be run in any order, individually or all together, depending on which region(s) you need. `utils.R` is sourced automatically by every other script and doesn't need to be run directly.

Each model script:

1.  Sources `utils.R` and loads the region's prepared inputs and its scenario table (from the corresponding `corridas_*.xlsx`, parsed into a data frame in `utils.R`).
2.  Iterates over every scenario (row) in that table with `purrr::pmap()`, building and solving a `prioritizr` problem for each — locking in the layers specified (e.g. RUNAP, RUNAP+OMEC), applying the specified feature targets, and using the specified cost layer.
3.  Solves with Gurobi (`gap = 0.05`) and writes a rasterized solution and summary stats for each scenario.
4.  **Runs are resumable**: at the top of each script's run section, it checks the region's `master_eval_summary.csv` and `failed_scenarios.txt` and skips any scenario already completed or already logged as a non-presolve failure. A second pass then retries scenarios that failed only their presolve check, forcing a solve.

## Outputs

Written to `results/<region>/` (auto-created), one folder per region:

| File | Contents |
|------------------------------------|------------------------------------|
| `<model_name>.tif` | Rasterized solution for that scenario: selected new priority cells vs. existing locked-in (already-protected) cells. Cells not part of the solution are NA. |
| `<model_name>_summary.csv` | Per-scenario target coverage detail, listing each feature as a row (including individual species and ecosystems, where applicable) with details on how much is held in the solution. |
| `master_eval_summary.csv` | One row per scenario across the whole region: total/new/locked-in cell counts, cost, % of targets met. Appended to (not overwritten) as scenarios complete, which is what makes reruns resumable. |
| `failed_scenarios.txt` | Log of scenario name + timestamp + error message for any scenario that failed to solve. Used to re-run model if interrupted partway. |
| `resultados_todos.csv` *(SIRAP scripts only)* | All per-scenario summary CSVs in the region concatenated into one file |

`model_name` encodes the scenario's targets and settings (e.g. feature targets, which locked-in layers were included, and the cost layer used) — see the `build_model_name()` functions in `utils.R` for the exact naming logic per region.

## Notes / caveats

-   Gurobi run time scales with the number of scenarios and planning units; national-scale runs will take considerably longer than the SIRAP regions. Some national terrestrial scenarios evaluate all \~8,000 species targets and can require **over 64 GB of RAM** — running the national scripts on a virtual machine with higher memory/processing power is recommended rather than a typical laptop. `failed_scenarios.txt` will document which scenarios failed due to memory limitations.
-   Several data layers (Orinoquía sabanas and congriales, Eje Cafetero wetlands, marine ecosystems and human footprint) were shared directly by partners and are not from a public portal — flagged in [Data](#data) where that's the case.
-   The BioModelos species section currently only runs for the **national terrestrial** model; SIRAP models do not yet include species as a feature, but a post-hoc evaluation for these scripts still evaluates species-level coverage in the solutions and appends to their `<model_name>_summary.csv`.

<!-- TODO: add anything else worth flagging — e.g. known data gaps, assumptions in target-setting, who to contact for access to the shared/internal datasets -->

## Citation / contact

<!-- TODO: write this out crediting folks. -->
