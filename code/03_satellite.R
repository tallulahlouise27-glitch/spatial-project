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
  cat("  Downloading", label, year, "(month by month)...\n")
  months_data <- bind_rows(lapply(1:12, function(mo) {
    start <- sprintf("%d-%02d-01", year, mo)
    end   <- format(as.Date(start) + months(1) - lubridate::days(1), "%Y-%m-%d")
    tryCatch({
      raw <- griddap(DATASET_ID,
                     url       = ERDDAP_URL,
                     time      = c(start, end),
                     latitude  = c(lat_min, lat_max),
                     longitude = c(lon_min, lon_max),
                     fields    = "AFAI")
      cat("    month", mo, ":", nrow(raw$data), "pixels\n")
      raw$data %>%
        filter(!is.na(AFAI)) %>%
        mutate(year = as.integer(year), month = as.integer(mo)) %>%
        select(longitude, latitude, year, month, AFAI)
    }, error = function(e) {
      cat("    WARNING: failed for", label, year, "month", mo, "–", conditionMessage(e), "\n")
      NULL
    })
  }))
  if (is.null(months_data) || nrow(months_data) == 0) return(NULL)
  months_data
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

  process_year <- function(yr) {
    cache_file <- file.path(path_raw_sat, sprintf("coastal_%d.rds", yr))

    # Use per-year cache so crashes don't lose completed years
    if (file.exists(cache_file)) {
      cat("  Loading cached year", yr, "\n")
      return(readRDS(cache_file))
    }

    pixels <- download_afai_year(yr,
                                 DR_LAT_MIN, DR_LAT_MAX,
                                 DR_LON_MIN, DR_LON_MAX,
                                 "DR coast")
    if (is.null(pixels) || nrow(pixels) == 0) return(NULL)

    pts    <- st_as_sf(pixels, coords = c("longitude", "latitude"), crs = 4326)
    joined <- st_join(pts, dr_buffered, join = st_within, left = FALSE)

    nearshore_flag <- lengths(st_intersects(joined, nearshore_excl)) > 0
    joined <- joined[!nearshore_flag, ] %>% st_drop_geometry()

    result <- joined %>%
      group_by(municipio, provincia, year, month) %>%
      summarise(
        afai_mean      = mean(AFAI, na.rm = TRUE),
        afai_max       = max(AFAI,  na.rm = TRUE),
        afai_coverage  = mean(AFAI > 0.001,           na.rm = TRUE),
        afai_sargassum = mean(pmax(AFAI - 0.001, 0),  na.rm = TRUE),
        n_pixels       = n(),
        .groups = "drop"
      )

    saveRDS(result, cache_file)
    cat("  Year", yr, "saved to cache.\n")
    result
  }

  muni_afai <- bind_rows(lapply(STUDY_YEARS, process_year))

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
cat("Mean AFAI (coastal):         ", round(mean(muni_afai$afai_mean,      na.rm=TRUE), 5), "\n")
cat("Mean Sargassum intensity:    ", round(mean(muni_afai$afai_sargassum, na.rm=TRUE), 6), "\n")
cat("Mean Sargassum coverage:     ", round(mean(muni_afai$afai_coverage,  na.rm=TRUE), 3), "(fraction of pixels)\n")
cat("====================================\n\n")

# ── Diagnostic 1: Seasonal pattern ────────────────────────────
# Sargassum in the Caribbean peaks May–September.
# afai_sargassum should be clearly higher in those months.
# A flat seasonal profile would indicate the measure is not
# capturing Sargassum but some other factor (e.g. algal noise).

cat("====== DIAGNOSTIC 1: SEASONAL PATTERN ======\n")
cat("Mean afai_sargassum by month (averaged across all years and municipalities)\n")
cat("Sargassum should peak May-Sep (months 5-9) and be near zero Jan-Feb\n\n")

seasonal <- muni_afai %>%
  group_by(month) %>%
  summarise(
    mean_sargassum = mean(afai_sargassum, na.rm = TRUE),
    mean_coverage  = mean(afai_coverage,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    month_name = month.abb[month],
    peak_flag  = ifelse(month %in% 5:9, "<<< PEAK SEASON", "")
  )

print(seasonal[, c("month_name", "mean_sargassum", "mean_coverage", "peak_flag")],
      row.names = FALSE)

peak_mean    <- mean(seasonal$mean_sargassum[seasonal$month %in% 5:9])
offpeak_mean <- mean(seasonal$mean_sargassum[!seasonal$month %in% 5:9])
cat(sprintf("\nPeak season mean: %.6f  |  Off-peak mean: %.6f  |  Ratio: %.1fx\n",
            peak_mean, offpeak_mean, peak_mean / offpeak_mean))
if (peak_mean > offpeak_mean * 1.5) {
  cat("PASS: Peak season is at least 1.5x higher than off-peak — seasonal signal present.\n")
} else {
  cat("WARNING: Weak seasonal signal. Peak is less than 1.5x off-peak — check for contamination.\n")
}

# ── Diagnostic 2: Known spike years ───────────────────────────
# Documented large Sargassum events from satellite literature and
# NOAA AOML ocean data (visible in our instrument data):
#   2018: Record year (March 2018 ocean coverage 0.065 — all-time high at the time)
#   2021: Sustained bloom Apr–Aug (ocean coverage ~0.031–0.041)
#   2022: Strong summer bloom (June 2022 ocean coverage 0.050)
#   2023: Early-onset event (March 2023 ocean coverage 0.063)
#   2025: All-time record (July 2025 ocean coverage 0.092, 41% above 2018 peak)
# Low year: 2017 (ocean coverage never above 0.012), 2024 (max 0.025)
#
# Our coastal afai_sargassum should reflect these patterns:
# 2018, 2021, 2022, 2023, 2025 >> 2017, 2024

cat("\n====== DIAGNOSTIC 2: KNOWN SPIKE YEARS ======\n")
cat("Comparing annual mean afai_sargassum against documented Sargassum events\n\n")

annual_intensity <- muni_afai %>%
  group_by(year) %>%
  summarise(
    mean_sargassum = mean(afai_sargassum, na.rm = TRUE),
    peak_sargassum = mean(afai_sargassum[month %in% 5:9], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    expected = case_when(
      year %in% c(2018, 2021, 2022, 2023, 2025) ~ "HIGH (documented spike)",
      year %in% c(2017, 2024)                    ~ "LOW  (documented quiet)",
      TRUE                                        ~ "moderate"
    )
  )

print(annual_intensity, row.names = FALSE)

high_mean <- mean(annual_intensity$mean_sargassum[annual_intensity$year %in% c(2018,2021,2022,2023,2025)])
low_mean  <- mean(annual_intensity$mean_sargassum[annual_intensity$year %in% c(2017,2024)])
cat(sprintf("\nMean for spike years (2018/21/22/23/25): %.6f\n", high_mean))
cat(sprintf("Mean for quiet years (2017/24):          %.6f\n", low_mean))
if (high_mean > low_mean * 1.3) {
  cat("PASS: Spike years are at least 1.3x higher than quiet years.\n")
} else {
  cat("WARNING: Spike years not clearly higher than quiet years — measure may not be tracking real events.\n")
}
cat("=============================================\n")
