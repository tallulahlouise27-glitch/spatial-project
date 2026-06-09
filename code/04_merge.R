# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 4: Merge ENCFT survey data with satellite data
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(sf)
library(geodata)

path_processed <- "data/processed/"

# ── Load data ─────────────────────────────────────────────────
cat("Loading ENCFT municipality panel...\n")
encft <- readRDS(file.path(path_processed, "encft_municipio.rds"))

cat("Loading satellite coastal data...\n")
sat_coastal <- readRDS(file.path(path_processed, "satellite_coastal.rds"))

cat("Loading satellite instrument data...\n")
sat_instr <- readRDS(file.path(path_processed, "satellite_instrument.rds"))

cat("ENCFT rows:", nrow(encft), "\n")
cat("Satellite coastal rows:", nrow(sat_coastal), "\n")
cat("Instrument rows:", nrow(sat_instr), "\n\n")

# ── Step 1: Average satellite data to municipality × year ─────
# (ENCFT is annual at municipality level)
sat_annual <- sat_coastal %>%
  group_by(municipio, provincia, year) %>%
  summarise(
    afai_mean_annual     = mean(afai_mean,       na.rm = TRUE),
    afai_cov_annual      = mean(afai_coverage,   na.rm = TRUE),
    # afai_sargassum: mean excess AFAI above 0.001 threshold — captures both
    # spatial coverage and mat density. Zero = no Sargassum, higher = more/denser.
    afai_sargassum_annual = mean(afai_sargassum, na.rm = TRUE),
    # Peak season (May–Sep) when Sargassum is most intense in Caribbean
    afai_peak             = mean(afai_coverage[month %in% 5:9],   na.rm = TRUE),
    afai_sargassum_peak   = mean(afai_sargassum[month %in% 5:9],  na.rm = TRUE),
    .groups = "drop"
  )

instr_annual <- sat_instr %>%
  group_by(year) %>%
  summarise(
    afai_ocean_annual    = mean(afai_ocean_mean,                                  na.rm = TRUE),
    afai_ocean_cov       = mean(afai_ocean_coverage,                              na.rm = TRUE),
    # Peak-season ocean coverage: restricts to May–Sep when the Atlantic
    # Sargassum belt is active. Annual mean dilutes the signal with winter months.
    afai_ocean_cov_peak  = mean(afai_ocean_coverage[month %in% 5:9],              na.rm = TRUE),
    .groups = "drop"
  )

# ── Step 2: Standardise municipality names for matching ───────
# ENCFT uses DES_MUNICIPIO (uppercase, accented Spanish)
# Satellite uses NAME_2 from GADM (mixed case, may differ)
# Strategy: convert both to lowercase, strip accents, remove punctuation

clean_name <- function(x) {
  x %>%
    tolower() %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    gsub("[^a-z0-9 ]", "", .) %>%
    trimws() %>%
    gsub("\\s+", " ", .)
}

# Crosswalk: ENCFT short/common names → GADM official names
# ENCFT uses popular names; GADM uses full official municipality names
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

encft_match <- encft %>%
  mutate(
    ANO          = as.integer(ANO),
    muni_key_raw = clean_name(DES_MUNICIPIO)
  ) %>%
  left_join(name_crosswalk, by = c("muni_key_raw" = "encft_name")) %>%
  mutate(muni_key = coalesce(gadm_name, muni_key_raw)) %>%
  select(-muni_key_raw, -gadm_name)

sat_match <- sat_annual %>%
  mutate(
    muni_key = clean_name(municipio)
  )

# ── Step 3: Classify municipalities by coastal proximity ──────
# Three categories:
#   "coastal"      — boundary directly touches the sea
#   "near_coastal" — within 40km of the coast but no direct shoreline
#   "inland"       — more than 40km from the coast
#
# Method: dissolve all municipality polygons → DR outline → extract coastline →
# test each municipality against it.

cat("Classifying municipalities by coastal proximity...\n")

dr_gadm <- gadm(country = "DOM", level = 2, path = "data/raw/satellite/")
dr_sf   <- st_as_sf(dr_gadm) %>%
  select(municipio = NAME_2) %>%
  st_transform(32619)   # UTM Zone 19N — metres, covers Dominican Republic

# Build coastline as the outer boundary of the dissolved island polygon
coastline <- st_union(dr_sf) %>% st_boundary()

# 500m buffer around coastline to handle floating-point edge cases
coast_zone <- st_buffer(coastline, dist = 500)

# coastal = municipality polygon overlaps the 500m coastal strip
touches_coast <- lengths(st_intersects(dr_sf, coast_zone)) > 0

# distance (metres) from each municipality centroid to coastline
dist_m <- as.numeric(st_distance(st_centroid(dr_sf), coastline))

coastal_lookup <- dr_sf %>%
  st_drop_geometry() %>%
  mutate(
    muni_key = clean_name(municipio),
    coastal_type = case_when(
      touches_coast      ~ "coastal",
      dist_m < 40000     ~ "near_coastal",
      TRUE               ~ "inland"
    )
  ) %>%
  select(muni_key, coastal_type)

cat("Coastal type breakdown:\n")
print(table(coastal_lookup$coastal_type))
cat("\n")

# ── Step 4: Merge ─────────────────────────────────────────────
panel <- encft_match %>%
  left_join(sat_match,      by = c("muni_key", "ANO" = "year")) %>%
  left_join(instr_annual,   by = c("ANO" = "year")) %>%
  left_join(coastal_lookup, by = "muni_key")

cat("Merge result:", nrow(panel), "rows\n")
cat("Municipalities matched with satellite data:",
    sum(!is.na(panel$afai_cov_annual)), "of", nrow(panel), "\n")

# Municipalities with NO satellite data in any year (truly unmatched)
truly_unmatched <- panel %>%
  group_by(DES_MUNICIPIO, DES_PROVINCIA) %>%
  summarise(any_sat = any(!is.na(afai_cov_annual)), .groups = "drop") %>%
  filter(!any_sat)

cat("\nMunicipalities with no satellite data at all:", nrow(truly_unmatched),
    "(likely landlocked — will be dropped from regression)\n")

# ── Step 4: Create analysis variables ─────────────────────────
panel <- panel %>%
  mutate(
    # Log income (add 1 to handle zeros)
    log_income    = log(ingreso_medio + 1),
    # Log income by tertile — main heterogeneity outcomes
    log_income_t1 = log(ingreso_T1 + 1),   # bottom third
    log_income_t2 = log(ingreso_T2 + 1),   # middle third
    log_income_t3 = log(ingreso_T3 + 1),   # top third
    # Year and municipality as factors for fixed effects
    year_fe  = factor(ANO),
    muni_fe  = factor(ID_MUNICIPIO)
  ) %>%
  # Keep only rows with all key variables present
  filter(!is.na(log_income), !is.na(afai_sargassum_annual), !is.na(afai_ocean_cov_peak))

cat("\nFinal analysis panel:", nrow(panel), "observations\n")
cat("Municipalities:       ", n_distinct(panel$ID_MUNICIPIO), "\n")
cat("Years:                ", paste(sort(unique(panel$ANO)), collapse = ", "), "\n")
cat("Employment rate mean: ", round(mean(panel$tasa_empleo, na.rm=TRUE), 3), "\n")
cat("Log income mean:      ", round(mean(panel$log_income, na.rm=TRUE), 3), "\n")
cat("\nMunicipalities by coastal type:\n")
print(table(panel$coastal_type, useNA = "ifany"))

# ── Step 5: Save ──────────────────────────────────────────────
saveRDS(panel, file.path(path_processed, "panel_analysis.rds"))
write_csv(panel, file.path(path_processed, "panel_analysis.csv"))
cat("\nSaved to data/processed/panel_analysis.rds\n")
