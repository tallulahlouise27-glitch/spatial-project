# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 3: Download and process satellite chlorophyll data
#
# Data source: MODIS Aqua Monthly Chlorophyll-a (erdMBchlamday)
# from NOAA CoastWatch ERDDAP — no login required
#
# What we build:
#   1. Treatment: mean coastal chlorophyll per DR municipality per year
#      (proxy for Sargassum exposure in nearshore waters)
#   2. Instrument: mean open-ocean chlorophyll in the central Atlantic
#      per year (proxy for offshore Sargassum — exogenous to DR economy)
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(rerddap)
library(sf)
library(terra)
library(geodata)

path_processed <- "data/processed/"
path_raw_sat   <- "data/raw/satellite/"
dir.create(path_raw_sat, showWarnings = FALSE, recursive = TRUE)

ERDDAP_URL <- "https://coastwatch.pfeg.noaa.gov/erddap/"
DATASET_ID <- "erdMBchlamday"

# Dominican Republic bounding box (0-360 longitude notation)
# DR spans ~17.5-20.5N, -72 to -68W => 288-292 in 0-360
DR_LAT_MIN <- 17.0; DR_LAT_MAX <- 21.0
DR_LON_MIN <- 287;  DR_LON_MAX <- 293

# Open-ocean instrument box: western Atlantic upstream of Caribbean
# (~10-22N, 45-65W) = 295-315 in 0-360 — far from DR coast
OC_LAT_MIN <- 10; OC_LAT_MAX <- 22
OC_LON_MIN <- 295; OC_LON_MAX <- 315

years <- 2016:2026

# Query dataset to find its actual last available date — avoids hard-coding
ds_info     <- info(DATASET_ID, url = ERDDAP_URL)
dataset_end <- tryCatch({
  nc <- ds_info$alldata$NC_GLOBAL
  raw <- nc$value[nc$attribute_name == "time_coverage_end"]
  as.Date(substr(raw, 1, 10))
}, error = function(e) Sys.Date())
cat("Dataset available through:", format(dataset_end), "\n\n")

# ── Helper: download one year from ERDDAP ─────────────────────
download_chla <- function(year, lat_min, lat_max, lon_min, lon_max, label) {
  cat("  Downloading", label, year, "...\n")
  start_date <- paste0(year, "-01-01")
  end_date   <- format(min(as.Date(paste0(year, "-12-31")), dataset_end), "%Y-%m-%d")
  tryCatch({
    dat <- griddap(
      DATASET_ID,
      url       = ERDDAP_URL,
      time      = c(start_date, end_date),
      latitude  = c(lat_min, lat_max),
      longitude = c(lon_min, lon_max),
      fields    = "chlorophyll"
    )
    dat$data %>%
      mutate(year = year) %>%
      filter(!is.na(chlorophyll))
  }, error = function(e) {
    cat("    WARNING: failed for", label, year, "-", conditionMessage(e), "\n")
    NULL
  })
}

# ── Step 1: Download DR coastal chlorophyll ────────────────────
dr_raw_file <- file.path(path_raw_sat, "dr_chla_raw.rds")
if (file.exists(dr_raw_file)) {
  cat("=== Loading cached DR coastal chlorophyll ===\n")
  dr_chla <- readRDS(dr_raw_file)
} else {
  cat("=== Downloading DR coastal chlorophyll (2016-2026) ===\n")
  dr_list <- lapply(years, download_chla,
                    lat_min = DR_LAT_MIN, lat_max = DR_LAT_MAX,
                    lon_min = DR_LON_MIN, lon_max = DR_LON_MAX,
                    label   = "DR coast")
  dr_chla <- bind_rows(dr_list)
  saveRDS(dr_chla, dr_raw_file)
}
cat("DR pixels loaded:", nrow(dr_chla), "\n\n")

# ── Step 2: Download open-ocean chlorophyll (instrument) ───────
oc_raw_file <- file.path(path_raw_sat, "oc_chla_raw.rds")
if (file.exists(oc_raw_file)) {
  cat("=== Loading cached open-ocean chlorophyll ===\n")
  oc_chla <- readRDS(oc_raw_file)
} else {
  cat("=== Downloading open-ocean chlorophyll (instrument) ===\n")
  oc_list <- lapply(years, download_chla,
                    lat_min = OC_LAT_MIN, lat_max = OC_LAT_MAX,
                    lon_min = OC_LON_MIN, lon_max = OC_LON_MAX,
                    label   = "Open ocean")
  oc_chla <- bind_rows(oc_list)
  saveRDS(oc_chla, oc_raw_file)
}
cat("Open-ocean pixels loaded:", nrow(oc_chla), "\n\n")

# ── Step 3: Load DR municipality boundaries ───────────────────
cat("=== Loading DR municipality boundaries ===\n")
dr_gadm <- gadm(country = "DOM", level = 2, path = path_raw_sat)
dr_sf   <- st_as_sf(dr_gadm) %>%
  select(municipio = NAME_2, provincia = NAME_1) %>%
  st_transform(4326)

# ── Step 4: Spatially join coastal pixels to municipalities ────
cat("=== Joining satellite pixels to municipalities ===\n")

# Convert DR chlorophyll to sf points (standard -180/180 longitude)
dr_pts <- dr_chla %>%
  mutate(lon_std = longitude - 360) %>%
  st_as_sf(coords = c("lon_std", "latitude"), crs = 4326)

# Buffer municipalities 20km seaward to capture nearshore ocean pixels
dr_buffered <- st_buffer(dr_sf, dist = 20000)

# Spatial join: assign each ocean pixel to the nearest municipality
joined <- st_join(dr_pts, dr_buffered, join = st_within, left = FALSE)

# Extract month from the time column (format: "2021-02-14T12:00:00Z")
joined <- joined %>%
  mutate(month = lubridate::month(time))

# Monthly mean per municipality
muni_chla <- joined %>%
  st_drop_geometry() %>%
  group_by(municipio, provincia, year, month) %>%
  summarise(
    chla_coastal_mean = mean(chlorophyll, na.rm = TRUE),
    chla_coastal_max  = max(chlorophyll, na.rm = TRUE),
    n_pixels          = n(),
    .groups = "drop"
  )

# ── Step 5: Summarise open-ocean instrument ────────────────────
instrument <- oc_chla %>%
  mutate(month = lubridate::month(time)) %>%
  group_by(year, month) %>%
  summarise(
    chla_ocean_mean = mean(chlorophyll, na.rm = TRUE),
    chla_ocean_max  = max(chlorophyll, na.rm = TRUE),
    .groups = "drop"
  )

# ── Step 6: Save ──────────────────────────────────────────────
saveRDS(muni_chla,  file.path(path_processed, "satellite_coastal.rds"))
saveRDS(instrument, file.path(path_processed, "satellite_instrument.rds"))
write_csv(muni_chla,  file.path(path_processed, "satellite_coastal.csv"))
write_csv(instrument, file.path(path_processed, "satellite_instrument.csv"))

# ── Summary ──────────────────────────────────────────────────
cat("\n====== SATELLITE SUMMARY ======\n")
cat("Municipality-year-month rows: ", nrow(muni_chla), "\n")
cat("Municipalities with data:    ", n_distinct(muni_chla$municipio), "\n")
cat("Years:                       ", paste(sort(unique(muni_chla$year)), collapse = ", "), "\n")
cat("Mean coastal chlorophyll:    ", round(mean(muni_chla$chla_coastal_mean, na.rm=TRUE), 3), "mg/m3\n")
cat("\nInstrument (open-ocean monthly mean, first 12 rows):\n")
print(head(instrument[, c("year", "month", "chla_ocean_mean")], 12))
cat("================================\n")
