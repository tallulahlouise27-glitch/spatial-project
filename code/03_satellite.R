# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 3: Download and process satellite chlorophyll data
#
# Two-source strategy to cover 2016-2025 without gaps:
#   - MODIS Aqua (erdMH1chlamday, var: chlorophyll): 2016 – May 2022
#   - VIIRS/NPP  (nesdisVHNSQchlaMonthly, var: chlor_a): Jun 2022 – 2025
#
# What we build:
#   1. Treatment: mean coastal chlorophyll per DR municipality per month
#   2. Instrument: mean open-ocean chlorophyll per month
#      (western Atlantic, exogenous to DR local economic conditions)
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(rerddap)
library(sf)
library(geodata)

path_processed <- "data/processed/"
path_raw_sat   <- "data/raw/satellite/"
dir.create(path_raw_sat, showWarnings = FALSE, recursive = TRUE)

ERDDAP_URL <- "https://coastwatch.pfeg.noaa.gov/erddap/"

# Dominican Republic bounding box
DR_LAT_MIN <- 17.0; DR_LAT_MAX <- 21.0
DR_LON_MIN <- -73;  DR_LON_MAX <- -68

# Open-ocean instrument box (western Atlantic, 10-22N, 65-45W)
OC_LAT_MIN <- 10; OC_LAT_MAX <- 22
OC_LON_MIN <- -65; OC_LON_MAX <- -45

# ── Helper: download one year from a given ERDDAP dataset ─────
download_chla <- function(year, lat_min, lat_max, lon_min, lon_max,
                          label, dataset_id, field, end_cap = NULL) {
  cat("  Downloading", label, year, paste0("(", dataset_id, ")"), "...\n")
  start_date <- paste0(year, "-01-01")
  end_date   <- paste0(year, "-12-31")
  if (!is.null(end_cap)) end_date <- format(min(as.Date(end_date), as.Date(end_cap)), "%Y-%m-%d")
  if (start_date > end_date) return(NULL)
  tryCatch({
    dat <- griddap(dataset_id, url = ERDDAP_URL,
                   time      = c(start_date, end_date),
                   latitude  = c(lat_min, lat_max),
                   longitude = c(lon_min, lon_max),
                   fields    = field)
    dat$data %>%
      rename(chlorophyll = all_of(field)) %>%
      mutate(year  = as.integer(year),
             month = lubridate::month(time)) %>%
      filter(!is.na(chlorophyll)) %>%
      select(longitude, latitude, time, year, month, chlorophyll)
  }, error = function(e) {
    cat("    WARNING: failed for", label, year, "-", conditionMessage(e), "\n")
    NULL
  })
}

# ── Download wrapper: combine MODIS (2016-2022) + VIIRS (2022-2025) ──
download_combined <- function(lat_min, lat_max, lon_min, lon_max, label) {
  # MODIS Aqua monthly: 2016 through May 2022
  modis <- bind_rows(lapply(2016:2022, download_chla,
    lat_min=lat_min, lat_max=lat_max, lon_min=lon_min, lon_max=lon_max,
    label=label, dataset_id="erdMH1chlamday", field="chlorophyll",
    end_cap="2022-05-31"))

  # VIIRS monthly: June 2022 through 2025 (avoids the 2018-2021 gap)
  viirs_years <- 2022:2025
  viirs <- bind_rows(lapply(viirs_years, function(yr) {
    start <- if (yr == 2022) "2022-06-01" else paste0(yr, "-01-01")
    tryCatch({
      dat <- griddap("nesdisVHNSQchlaMonthly", url=ERDDAP_URL,
                     time=c(start, paste0(yr, "-12-31")),
                     latitude=c(lat_min, lat_max),
                     longitude=c(lon_min, lon_max),
                     fields="chlor_a")
      dat$data %>%
        rename(chlorophyll = chlor_a) %>%
        mutate(year=as.integer(yr), month=lubridate::month(time)) %>%
        filter(!is.na(chlorophyll)) %>%
        select(longitude, latitude, time, year, month, chlorophyll)
    }, error=function(e) {
      cat("    WARNING: VIIRS failed for", label, yr, "-", conditionMessage(e), "\n")
      NULL
    })
  }))

  bind_rows(modis, viirs)
}

# ── Step 1: Download DR coastal chlorophyll ────────────────────
dr_raw_file <- file.path(path_raw_sat, "dr_chla_raw.rds")
if (file.exists(dr_raw_file)) {
  cat("=== Loading cached DR coastal chlorophyll ===\n")
  dr_chla <- readRDS(dr_raw_file)
} else {
  cat("=== Downloading DR coastal chlorophyll (2016-2025) ===\n")
  dr_chla <- download_combined(DR_LAT_MIN, DR_LAT_MAX, DR_LON_MIN, DR_LON_MAX, "DR coast")
  saveRDS(dr_chla, dr_raw_file)
}
cat("DR pixels loaded:", nrow(dr_chla), "| Years:", paste(sort(unique(dr_chla$year)), collapse=", "), "\n\n")

# ── Step 2: Download open-ocean chlorophyll (instrument) ───────
oc_raw_file <- file.path(path_raw_sat, "oc_chla_raw.rds")
if (file.exists(oc_raw_file)) {
  cat("=== Loading cached open-ocean chlorophyll ===\n")
  oc_chla <- readRDS(oc_raw_file)
} else {
  cat("=== Downloading open-ocean chlorophyll (instrument, 2016-2025) ===\n")
  oc_chla <- download_combined(OC_LAT_MIN, OC_LAT_MAX, OC_LON_MIN, OC_LON_MAX, "Open ocean")
  saveRDS(oc_chla, oc_raw_file)
}
cat("Open-ocean pixels loaded:", nrow(oc_chla), "| Years:", paste(sort(unique(oc_chla$year)), collapse=", "), "\n\n")

# ── Step 3: Load DR municipality boundaries ───────────────────
cat("=== Loading DR municipality boundaries ===\n")
dr_gadm <- gadm(country = "DOM", level = 2, path = path_raw_sat)
dr_sf   <- st_as_sf(dr_gadm) %>%
  select(municipio = NAME_2, provincia = NAME_1) %>%
  st_transform(4326)

# ── Step 4: Spatially join coastal pixels to municipalities ────
cat("=== Joining satellite pixels to municipalities ===\n")

dr_pts      <- st_as_sf(dr_chla, coords = c("longitude", "latitude"), crs = 4326)
dr_buffered <- st_buffer(dr_sf, dist = 20000)
joined      <- st_join(dr_pts, dr_buffered, join = st_within, left = FALSE)

muni_chla <- joined %>%
  st_drop_geometry() %>%
  group_by(municipio, provincia, year, month) %>%
  summarise(
    chla_coastal_mean = mean(chlorophyll, na.rm = TRUE),
    chla_coastal_max  = max(chlorophyll,  na.rm = TRUE),
    n_pixels          = n(),
    .groups = "drop"
  )

# ── Step 5: Summarise open-ocean instrument ────────────────────
instrument <- oc_chla %>%
  group_by(year, month) %>%
  summarise(
    chla_ocean_mean = mean(chlorophyll, na.rm = TRUE),
    chla_ocean_max  = max(chlorophyll,  na.rm = TRUE),
    .groups = "drop"
  )

# ── Step 6: Save ──────────────────────────────────────────────
saveRDS(muni_chla,  file.path(path_processed, "satellite_coastal.rds"))
saveRDS(instrument, file.path(path_processed, "satellite_instrument.rds"))
write_csv(muni_chla,  file.path(path_processed, "satellite_coastal.csv"))
write_csv(instrument, file.path(path_processed, "satellite_instrument.csv"))

cat("\n====== SATELLITE SUMMARY ======\n")
cat("Municipality-year-month rows: ", nrow(muni_chla), "\n")
cat("Municipalities with data:    ", n_distinct(muni_chla$municipio), "\n")
cat("Years:                       ", paste(sort(unique(muni_chla$year)), collapse = ", "), "\n")
cat("Mean coastal chlorophyll:    ", round(mean(muni_chla$chla_coastal_mean, na.rm=TRUE), 3), "mg/m3\n")
cat("\nInstrument (open-ocean, first 24 rows):\n")
print(head(instrument[, c("year", "month", "chla_ocean_mean")], 24))
cat("================================\n")
