# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 3: Download and process AFAI Sargassum data
#
# Dataset: noaa_aoml_atlantic_oceanwatch_AFAI_7D
# Source:  NOAA AOML / USF Optical Oceanography Lab
# ERDDAP:  https://cwcgom.aoml.noaa.gov/erddap/
# Years:   2017–2025 (9 complete years)
#
# AFAI (Alternative Floating Algae Index):
#   Positive values = floating Sargassum detected at surface
#   Negative / near-zero = open ocean, no Sargassum
#   NA = cloud cover
#
# What we build:
#   1. Treatment: Sargassum coverage per DR municipality per month
#      - afai_mean:     mean AFAI across coastal pixels
#      - afai_coverage: fraction of pixels with positive AFAI (0–1)
#   2. Instrument: mean open-ocean AFAI per month
#      (western Atlantic, exogenous to DR local economic conditions)
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(rerddap)
library(sf)
library(geodata)

path_processed <- "data/processed/"
path_raw_sat   <- "data/raw/satellite/"
dir.create(path_raw_sat,   showWarnings = FALSE, recursive = TRUE)
dir.create(path_processed, showWarnings = FALSE, recursive = TRUE)

ERDDAP_URL  <- "https://cwcgom.aoml.noaa.gov/erddap/"
DATASET_ID  <- "noaa_aoml_atlantic_oceanwatch_AFAI_7D"
STUDY_YEARS <- 2017:2025

# Dominican Republic bounding box
DR_LAT_MIN <- 17.0; DR_LAT_MAX <- 21.0
DR_LON_MIN <- -73.0; DR_LON_MAX <- -68.0

# Open-ocean instrument box (western Atlantic, upstream of DR)
OC_LAT_MIN <- 10.0; OC_LAT_MAX <- 25.0
OC_LON_MIN <- -65.0; OC_LON_MAX <- -45.0

# ── Helper: download one year of AFAI for a bounding box ──────
download_afai_year <- function(year, lat_min, lat_max, lon_min, lon_max, label) {
  cat("  Downloading", label, year, "...\n")
  tryCatch({
    raw <- griddap(DATASET_ID,
                   url       = ERDDAP_URL,
                   time      = c(paste0(year, "-01-01"), paste0(year, "-12-31")),
                   latitude  = c(lat_min, lat_max),
                   longitude = c(lon_min, lon_max),
                   fields    = "AFAI")
    raw$data %>%
      mutate(
        year  = as.integer(year),
        month = lubridate::month(time)
      ) %>%
      filter(!is.na(AFAI)) %>%
      select(longitude, latitude, year, month, AFAI)
  }, error = function(e) {
    cat("    WARNING: failed for", label, year, "–", conditionMessage(e), "\n")
    NULL
  })
}

# ── Step 1: Load DR municipality boundaries ───────────────────
cat("=== Loading DR municipality boundaries ===\n")
dr_gadm     <- gadm(country = "DOM", level = 2, path = path_raw_sat)
dr_sf       <- st_as_sf(dr_gadm) %>%
  select(municipio = NAME_2, provincia = NAME_1) %>%
  st_transform(4326)
dr_buffered <- st_buffer(dr_sf, dist = 20000)   # 20 km seaward buffer
cat("Municipalities loaded:", nrow(dr_sf), "\n\n")

# ── Step 2: Download + spatially aggregate DR coastal AFAI ────
# Process year by year to keep memory usage manageable.
# Raw pixel data (~1 GB/year) is never saved to disk — only the
# municipality × month aggregates are retained.

out_coastal <- file.path(path_processed, "satellite_coastal.rds")

if (file.exists(out_coastal)) {
  cat("=== Loading cached coastal AFAI ===\n")
  muni_afai <- readRDS(out_coastal)
} else {
  cat("=== Downloading coastal AFAI (2017–2025) ===\n")
  muni_afai <- bind_rows(lapply(STUDY_YEARS, function(yr) {

    pixels <- download_afai_year(yr,
                                 DR_LAT_MIN, DR_LAT_MAX,
                                 DR_LON_MIN, DR_LON_MAX,
                                 "DR coast")
    if (is.null(pixels) || nrow(pixels) == 0) return(NULL)

    # Spatially join pixels to buffered municipalities
    pts    <- st_as_sf(pixels, coords = c("longitude", "latitude"), crs = 4326)
    joined <- st_join(pts, dr_buffered, join = st_within, left = FALSE) %>%
      st_drop_geometry()

    # Aggregate to municipality × month
    joined %>%
      group_by(municipio, provincia, year, month) %>%
      summarise(
        afai_mean     = mean(AFAI, na.rm = TRUE),
        afai_max      = max(AFAI,  na.rm = TRUE),
        afai_coverage = mean(AFAI > 0, na.rm = TRUE),  # fraction of pixels with Sargassum
        n_pixels      = n(),
        .groups = "drop"
      )
  }))

  saveRDS(muni_afai, out_coastal)
  write_csv(muni_afai, file.path(path_processed, "satellite_coastal.csv"))
  cat("Coastal AFAI saved.\n\n")
}

# ── Step 3: Download + aggregate open-ocean AFAI (instrument) ─
out_instr <- file.path(path_processed, "satellite_instrument.rds")

if (file.exists(out_instr)) {
  cat("=== Loading cached open-ocean AFAI (instrument) ===\n")
  instrument <- readRDS(out_instr)
} else {
  cat("=== Downloading open-ocean AFAI instrument (2017–2025) ===\n")
  instrument <- bind_rows(lapply(STUDY_YEARS, function(yr) {

    pixels <- download_afai_year(yr,
                                 OC_LAT_MIN, OC_LAT_MAX,
                                 OC_LON_MIN, OC_LON_MAX,
                                 "Open ocean")
    if (is.null(pixels) || nrow(pixels) == 0) return(NULL)

    # Aggregate all ocean pixels to a single monthly mean
    pixels %>%
      group_by(year, month) %>%
      summarise(
        afai_ocean_mean     = mean(AFAI, na.rm = TRUE),
        afai_ocean_coverage = mean(AFAI > 0, na.rm = TRUE),
        .groups = "drop"
      )
  }))

  saveRDS(instrument, out_instr)
  write_csv(instrument, file.path(path_processed, "satellite_instrument.csv"))
  cat("Instrument AFAI saved.\n\n")
}

# ── Summary ───────────────────────────────────────────────────
cat("====== AFAI SATELLITE SUMMARY ======\n")
cat("Municipality-year-month rows:", nrow(muni_afai), "\n")
cat("Municipalities with data:    ", n_distinct(muni_afai$municipio), "\n")
cat("Years:                       ", paste(sort(unique(muni_afai$year)), collapse = ", "), "\n")
cat("Mean AFAI (coastal):         ", round(mean(muni_afai$afai_mean,     na.rm=TRUE), 5), "\n")
cat("Mean Sargassum coverage:     ", round(mean(muni_afai$afai_coverage, na.rm=TRUE), 3), "(fraction of pixels)\n")
cat("\nInstrument (open-ocean, first 24 rows):\n")
print(head(instrument[, c("year", "month", "afai_ocean_mean", "afai_ocean_coverage")], 24))
cat("====================================\n")
