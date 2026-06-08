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
dr_gadm <- gadm(country = "DOM", level = 2, path = path_raw_sat)
dr_sf   <- st_as_sf(dr_gadm) %>%
  select(municipio = NAME_2, provincia = NAME_1) %>%
  st_transform(4326)

# Build spatial filters in projected CRS (UTM Zone 19N) so distances are in metres
dr_sf_proj  <- st_transform(dr_sf, 32619)
dr_buffered <- st_buffer(dr_sf_proj, dist = 20000) %>% st_transform(4326)  # 20km seaward

# Nearshore exclusion zone: strip pixels within 3km of the coastline.
# The 0–3km zone is shallow water dominated by seagrass and coral bottom
# reflectance that permanently registers positive AFAI regardless of Sargassum.
coastline_proj  <- st_union(dr_sf_proj) %>% st_boundary()
nearshore_excl  <- st_buffer(coastline_proj, dist = 3000) %>% st_transform(4326)

cat("Municipalities loaded:", nrow(dr_sf), "\n")
cat("Nearshore exclusion zone: 3km from coastline\n\n")

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
    joined <- st_join(pts, dr_buffered, join = st_within, left = FALSE)

    # Remove nearshore pixels (0–3km from coastline) — shallow water
    # contamination from seagrass and coral bottom reflectance
    nearshore_flag <- lengths(st_intersects(joined, nearshore_excl)) > 0
    joined <- joined[!nearshore_flag, ] %>% st_drop_geometry()

    # Aggregate to municipality × month
    # Three complementary Sargassum measures:
    #   afai_coverage  — fraction of pixels above the detection threshold (0–1)
    #                    captures spatial extent: "how much of the coast has Sargassum?"
    #   afai_sargassum — mean excess AFAI above threshold, pmax(AFAI - 0.001, 0)
    #                    captures both extent AND density: "how much Sargassum overall?"
    #                    zero when absent, higher = more/denser mats
    #   afai_max       — peak AFAI value (useful for detecting episodic heavy events)
    # Threshold 0.001 from Wang & Hu (2016) AFAI methodology.
    joined %>%
      group_by(municipio, provincia, year, month) %>%
      summarise(
        afai_mean      = mean(AFAI, na.rm = TRUE),
        afai_max       = max(AFAI,  na.rm = TRUE),
        afai_coverage  = mean(AFAI > 0.001,           na.rm = TRUE),
        afai_sargassum = mean(pmax(AFAI - 0.001, 0),  na.rm = TRUE),
        n_pixels       = n(),
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
  # Ocean box is large (~300 sq degrees), so download one month at a time
  # and aggregate to a single mean immediately — never holding a full year in memory.
  instrument <- bind_rows(lapply(STUDY_YEARS, function(yr) {
    cat("  Downloading Open ocean", yr, "...\n")
    bind_rows(lapply(1:12, function(mo) {
      start <- sprintf("%d-%02d-01", yr, mo)
      end   <- format(as.Date(start) + months(1) - lubridate::days(1), "%Y-%m-%d")
      tryCatch({
        raw <- griddap(DATASET_ID, url = ERDDAP_URL,
                       time      = c(start, end),
                       latitude  = c(OC_LAT_MIN, OC_LAT_MAX),
                       longitude = c(OC_LON_MIN, OC_LON_MAX),
                       fields    = "AFAI")
        raw$data %>%
          filter(!is.na(AFAI)) %>%
          summarise(
            year                = yr,
            month               = mo,
            afai_ocean_mean     = mean(AFAI, na.rm = TRUE),
            afai_ocean_coverage = mean(AFAI > 0, na.rm = TRUE)
          )
      }, error = function(e) {
        cat("    WARNING: Ocean", yr, mo, "–", conditionMessage(e), "\n")
        NULL
      })
    }))
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
