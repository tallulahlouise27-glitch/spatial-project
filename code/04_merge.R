# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 4: Merge ENCFT survey data with satellite data
#
# Panel structure: municipality × year × quarter
#
# ENCFT is designed to be representative at the quarterly level.
# Using quarterly income preserves this structure and avoids repeating
# an annual figure across months where income doesn't actually vary.
#
# Treatment (afai_sargassum): quarterly mean coastal Sargassum intensity
#   mean of monthly afai_sargassum across the 3 months in each quarter
#
# Instrument: Bartik shift-share (constructed in script 05)
#   z_bartik = afai_ocean_coverage_qt × baseline_sargassum_i
#   Quarterly ocean AFAI × municipality's full-sample mean exposure
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(sf)
library(geodata)

path_processed <- "data/processed/"

# ── Load data ─────────────────────────────────────────────────
cat("Loading ENCFT quarterly municipality panel...\n")
encft <- readRDS(file.path(path_processed, "encft_municipio_quarterly.rds"))

cat("Loading ENCFT annual municipality panel (for robustness)...\n")
encft_annual <- readRDS(file.path(path_processed, "encft_municipio.rds"))

cat("Loading satellite coastal data (monthly)...\n")
sat_coastal <- readRDS(file.path(path_processed, "satellite_coastal.rds"))

cat("Loading satellite instrument data (monthly)...\n")
sat_instr <- readRDS(file.path(path_processed, "satellite_instrument.rds"))

cat("ENCFT quarterly rows:", nrow(encft), "\n")
cat("Satellite coastal rows:", nrow(sat_coastal), "\n")
cat("Instrument rows:", nrow(sat_instr), "\n\n")

# ── Step 1: Standardise municipality names for matching ───────
clean_name <- function(x) {
  x %>%
    tolower() %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    gsub("[^a-z0-9 ]", "", .) %>%
    trimws() %>%
    gsub("\\s+", " ", .)
}

name_crosswalk <- tribble(
  ~encft_name,                 ~gadm_name,
  "azua",                      "azua de compostela",
  "barahona",                  "santa cruz de barahona",
  "bisono",                    "villa bisono",
  "cambita garabitos",         "cambita garabito",
  "castanuelas",               "castanuela",
  "el seibo",                  "santa cruz del seybo",
  "hato mayor",                "hato mayor del rey",
  "higuey",                    "salvaleon de higuey",
  "la vega",                   "concepcion de la vega",
  "monte cristi",              "san fernando de monte cristi",
  "neiba",                     "neyba",
  "peralvillo",                "esperalvillo",
  "puerto plata",              "san felipe de puerto plata",
  "quisqueya",                 "quisquella",
  "samana",                    "santa barbara de samana",
  "san antonio de guerra",     "guerra",
  "san gregorio de nigua",     "nigua",
  "san juan",                  "san juan de la maguana",
  "santiago",                  "santiago de los caballeros",
  "santo domingo de guzman",   "distrito nacional",
  "tabara arriba",             "villa tabara arriba",
  "villa isabela",             "la isabela",
  "villa montellano",          "montellano",
  "villa riva",                "villa rivas",
  "villa vasquez",             "villa vazquez",
  "yaguate",                   "san gregorio de yaguate"
)

# Aggregate satellite to quarterly (Q1=Jan-Mar, Q2=Apr-Jun, Q3=Jul-Sep, Q4=Oct-Dec)
sat_quarterly <- sat_coastal %>%
  mutate(
    muni_key = clean_name(municipio),
    quarter  = ceiling(month / 3)
  ) %>%
  group_by(muni_key, municipio, provincia, year, quarter) %>%
  summarise(
    afai_sargassum = mean(afai_sargassum, na.rm = TRUE),
    afai_coverage  = mean(afai_coverage,  na.rm = TRUE),
    afai_mean      = mean(afai_mean,      na.rm = TRUE),
    n_months       = n(),
    .groups = "drop"
  )

instr_quarterly <- sat_instr %>%
  mutate(quarter = ceiling(month / 3)) %>%
  group_by(year, quarter) %>%
  summarise(
    afai_ocean_coverage = mean(afai_ocean_coverage, na.rm = TRUE),
    afai_ocean_mean     = mean(afai_ocean_mean,     na.rm = TRUE),
    n_months            = n(),
    .groups = "drop"
  )

# ── Step 2: Classify municipalities by coastal proximity ──────
cat("Classifying municipalities by coastal proximity...\n")

dr_gadm <- gadm(country = "DOM", level = 2, path = "data/raw/satellite/")
dr_sf   <- st_as_sf(dr_gadm) %>%
  select(municipio = NAME_2) %>%
  st_transform(32619)

coastline  <- st_union(dr_sf) %>% st_boundary()
coast_zone <- st_buffer(coastline, dist = 500)
touches_coast <- lengths(st_intersects(dr_sf, coast_zone)) > 0
dist_m <- as.numeric(st_distance(st_centroid(dr_sf), coastline))

coastal_lookup <- dr_sf %>%
  st_drop_geometry() %>%
  mutate(
    muni_key = clean_name(municipio),
    coastal_type = case_when(
      touches_coast  ~ "coastal",
      dist_m < 40000 ~ "near_coastal",
      TRUE           ~ "inland"
    )
  ) %>%
  select(muni_key, coastal_type)

cat("Coastal type breakdown:\n")
print(table(coastal_lookup$coastal_type))
cat("\n")

# ── Step 3: Build quarterly panel ─────────────────────────────
# ENCFT quarterly income × quarterly satellite × quarterly ocean instrument.
# The ENCFT is designed to be representative at the quarterly level —
# using quarterly income avoids repeating an annual figure and correctly
# reflects the survey design.

encft_match <- encft %>%
  mutate(
    ANO          = as.integer(ANO),
    TRIMESTRE    = as.integer(TRIMESTRE),
    muni_key_raw = clean_name(DES_MUNICIPIO)
  ) %>%
  left_join(name_crosswalk, by = c("muni_key_raw" = "encft_name")) %>%
  mutate(muni_key = coalesce(gadm_name, muni_key_raw)) %>%
  select(-muni_key_raw, -gadm_name)

panel <- encft_match %>%
  left_join(sat_quarterly,  by = c("muni_key", "ANO" = "year", "TRIMESTRE" = "quarter")) %>%
  left_join(instr_quarterly, by = c("ANO" = "year", "TRIMESTRE" = "quarter")) %>%
  left_join(coastal_lookup,  by = "muni_key")

cat("Quarterly merge result:", nrow(panel), "rows\n")
cat("Muni-year-quarters with satellite data:",
    sum(!is.na(panel$afai_sargassum)), "of", nrow(panel), "\n")

truly_unmatched <- panel %>%
  group_by(DES_MUNICIPIO, DES_PROVINCIA) %>%
  summarise(any_sat = any(!is.na(afai_sargassum)), .groups = "drop") %>%
  filter(!any_sat)
cat("Municipalities with no satellite data at all:", nrow(truly_unmatched), "\n\n")

# ── Step 4: Create analysis variables ─────────────────────────
panel <- panel %>%
  mutate(
    log_income    = log(ingreso_medio + 1),
    log_income_t1 = log(ingreso_T1 + 1),
    log_income_t2 = log(ingreso_T2 + 1),
    log_income_t3 = log(ingreso_T3 + 1),
    muni_fe          = factor(ID_MUNICIPIO),
    year_fe          = factor(ANO),
    quarter_fe       = factor(TRIMESTRE),
    # Year×quarter FE: absorbs all common national shocks in each quarter
    # (Sargassum season, macroeconomic fluctuations, survey-wave effects)
    year_quarter_fe  = factor(paste0(ANO, "_q", TRIMESTRE))
  ) %>%
  filter(!is.na(log_income), !is.na(afai_sargassum), !is.na(afai_ocean_coverage))

cat("Final quarterly panel:", nrow(panel), "observations\n")
cat("Municipalities:       ", n_distinct(panel$ID_MUNICIPIO), "\n")
cat("Years:                ", paste(sort(unique(panel$ANO)), collapse = ", "), "\n")
cat("Quarters per muni-year (median):",
    median(panel %>% count(ID_MUNICIPIO, ANO) %>% pull(n)), "\n")
cat("Log income mean:      ", round(mean(panel$log_income, na.rm = TRUE), 3), "\n")
cat("\nMunicipalities by coastal type:\n")
print(table(panel$coastal_type, useNA = "ifany"))

# ── Step 5: Annual panel for robustness ───────────────────────
encft_annual_match <- encft_annual %>%
  mutate(
    ANO          = as.integer(ANO),
    muni_key_raw = clean_name(DES_MUNICIPIO)
  ) %>%
  left_join(name_crosswalk, by = c("muni_key_raw" = "encft_name")) %>%
  mutate(muni_key = coalesce(gadm_name, muni_key_raw)) %>%
  select(-muni_key_raw, -gadm_name)

sat_annual <- sat_coastal %>%
  mutate(muni_key = clean_name(municipio)) %>%
  group_by(muni_key, year) %>%
  summarise(
    afai_sargassum_annual = mean(afai_sargassum, na.rm = TRUE),
    afai_cov_annual       = mean(afai_coverage,  na.rm = TRUE),
    .groups = "drop"
  )

instr_annual <- sat_instr %>%
  group_by(year) %>%
  summarise(
    afai_ocean_cov      = mean(afai_ocean_coverage, na.rm = TRUE),
    afai_ocean_cov_peak = mean(afai_ocean_coverage[month %in% 5:9], na.rm = TRUE),
    .groups = "drop"
  )

panel_annual <- encft_annual_match %>%
  left_join(sat_annual,      by = c("muni_key", "ANO" = "year")) %>%
  left_join(instr_annual,    by = c("ANO" = "year")) %>%
  left_join(coastal_lookup,  by = "muni_key") %>%
  mutate(
    log_income    = log(ingreso_medio + 1),
    log_income_t1 = log(ingreso_T1 + 1),
    log_income_t2 = log(ingreso_T2 + 1),
    log_income_t3 = log(ingreso_T3 + 1),
    muni_fe = factor(ID_MUNICIPIO),
    year_fe = factor(ANO)
  ) %>%
  filter(!is.na(log_income), !is.na(afai_sargassum_annual), !is.na(afai_ocean_cov_peak))

cat("\nAnnual robustness panel:", nrow(panel_annual), "observations\n")

# ── Step 6: Save ──────────────────────────────────────────────
saveRDS(panel,        file.path(path_processed, "panel_quarterly.rds"))
saveRDS(panel_annual, file.path(path_processed, "panel_analysis.rds"))
write_csv(panel,        file.path(path_processed, "panel_quarterly.csv"))
write_csv(panel_annual, file.path(path_processed, "panel_analysis.csv"))
cat("\nSaved panel_quarterly.rds (main) and panel_analysis.rds (annual robustness)\n")
