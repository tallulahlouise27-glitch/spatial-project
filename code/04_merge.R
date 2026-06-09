# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 4: Merge ENCFT survey data with satellite data
#
# Panel structure: municipality × year × month
#
# Income: quarterly ENCFT, matched to each month within the quarter.
#   Months 1–3 → Q1 income, months 4–6 → Q2, etc.
#   ENCFT is designed to be representative at quarterly level, so income
#   does not vary within quarters — only between quarters and years.
#   Year×quarter FEs in the regression absorb this mechanical within-quarter
#   constancy and focus identification on month-to-month Sargassum variation.
#
# Treatment: monthly afai_sargassum per municipality
# Instrument: Bartik shift-share (built in script 05)
#   z_bartik = monthly_ocean_coverage_mt × baseline_sargassum_i
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(sf)
library(geodata)

path_processed <- "data/processed/"

# ── Load data ─────────────────────────────────────────────────
cat("Loading ENCFT quarterly municipality panel...\n")
encft_q <- readRDS(file.path(path_processed, "encft_municipio_quarterly.rds"))

cat("Loading ENCFT annual panel (for robustness)...\n")
encft_a <- readRDS(file.path(path_processed, "encft_municipio.rds"))

cat("Loading satellite coastal data (monthly)...\n")
sat_coastal <- readRDS(file.path(path_processed, "satellite_coastal.rds"))

cat("Loading satellite instrument data (monthly)...\n")
sat_instr <- readRDS(file.path(path_processed, "satellite_instrument.rds"))

cat("ENCFT quarterly rows:", nrow(encft_q), "\n")
cat("Satellite coastal rows:", nrow(sat_coastal), "\n")
cat("Instrument rows:", nrow(sat_instr), "\n\n")

# ── Step 1: Standardise municipality names ────────────────────
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
  "villa los almacigos",       "los almacigos",
  "villa vasquez",             "villa vazquez",
  "yaguate",                   "san gregorio de yaguate"
)

apply_crosswalk <- function(df, year_col = "ANO") {
  df %>%
    mutate(
      !!year_col := as.integer(.data[[year_col]]),
      muni_key_raw = clean_name(DES_MUNICIPIO)
    ) %>%
    left_join(name_crosswalk, by = c("muni_key_raw" = "encft_name")) %>%
    mutate(muni_key = coalesce(gadm_name, muni_key_raw)) %>%
    select(-muni_key_raw, -gadm_name)
}

encft_q_match <- apply_crosswalk(encft_q)
encft_a_match <- apply_crosswalk(encft_a)

# ── Step 2: Classify municipalities by coastal proximity ──────
cat("Classifying municipalities by coastal proximity...\n")

dr_gadm <- gadm(country = "DOM", level = 2, path = "data/raw/satellite/")
dr_sf   <- st_as_sf(dr_gadm) %>%
  select(municipio = NAME_2) %>%
  st_transform(32619)

# Ocean coastline only: subtract the Haiti land border from DR's outer boundary.
# st_boundary() of the DR polygon includes both the ocean coast AND the Haiti
# land border — municipalities along that border would be wrongly classified
# as coastal without this correction.
haiti_gadm  <- gadm(country = "HTI", level = 0, path = "data/raw/satellite/")
haiti_sf    <- st_as_sf(haiti_gadm) %>% st_transform(32619)
haiti_buf   <- st_buffer(st_union(haiti_sf), dist = 2000)
dr_boundary <- st_union(dr_sf) %>% st_boundary()
coastline   <- st_difference(dr_boundary, haiti_buf)

coast_zone    <- st_buffer(coastline, dist = 500)
touches_coast <- lengths(st_intersects(dr_sf, coast_zone)) > 0

coastal_lookup <- dr_sf %>%
  st_drop_geometry() %>%
  mutate(
    muni_key     = clean_name(municipio),
    coastal_type = case_when(
      touches_coast ~ "coastal",
      TRUE          ~ "not_coastal"
    )
  ) %>%
  select(muni_key, coastal_type)

cat("Coastal type breakdown:\n")
print(table(coastal_lookup$coastal_type))
cat("\n")

# ── Step 3: Build monthly panel ───────────────────────────────
# Add muni_key and quarter to satellite data.
# Join quarterly ENCFT income by (muni_key, ANO, quarter).
# Result: each monthly satellite row carries the income from the
# ENCFT wave covering that month's quarter.

sat_monthly <- sat_coastal %>%
  mutate(
    muni_key = clean_name(municipio),
    quarter  = ceiling(month / 3)
  )

panel <- sat_monthly %>%
  left_join(
    encft_q_match %>% mutate(TRIMESTRE = as.integer(TRIMESTRE) %% 10),
    by = c("muni_key", "year" = "ANO", "quarter" = "TRIMESTRE")
  ) %>%
  left_join(sat_instr,    by = c("year", "month")) %>%
  left_join(coastal_lookup, by = "muni_key") %>%
  filter(!is.na(ID_MUNICIPIO))   # drop satellite rows with no ENCFT match

cat("Monthly panel before filtering:", nrow(panel), "rows\n")

# ── Step 4: Create analysis variables ─────────────────────────
panel <- panel %>%
  mutate(
    log_income    = log(ingreso_medio + 1),
    log_income_t1 = log(ingreso_T1 + 1),
    log_income_t2 = log(ingreso_T2 + 1),
    log_income_t3 = log(ingreso_T3 + 1),
    muni_fe          = factor(ID_MUNICIPIO),
    year_fe          = factor(year),
    month_fe         = factor(month),
    quarter_fe       = factor(quarter),
    # Year×month FE: absorbs every common national shock at monthly frequency
    # (seasonal Sargassum pattern × year, national economic events, etc.)
    year_month_fe    = factor(paste0(year, "_m", sprintf("%02d", month))),
    # Year×quarter FE: absorbs mechanical within-quarter income constancy;
    # this is the "deal with quarterly representativeness" control
    year_quarter_fe  = factor(paste0(year, "_q", quarter))
  ) %>%
  filter(!is.na(log_income), !is.na(afai_sargassum), !is.na(afai_ocean_coverage))

cat("Final monthly panel:", nrow(panel), "observations\n")
cat("Municipalities:     ", n_distinct(panel$ID_MUNICIPIO), "\n")
cat("Years:              ", paste(sort(unique(panel$year)), collapse = ", "), "\n")
cat("Months per muni-year (median):",
    median(panel %>% count(ID_MUNICIPIO, year) %>% pull(n)), "\n")
cat("Log income mean:    ", round(mean(panel$log_income, na.rm = TRUE), 3), "\n")
cat("\nMunicipalities by coastal type:\n")
print(table(panel$coastal_type, useNA = "ifany"))

# ── Step 5: Annual panel for robustness ───────────────────────
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

panel_annual <- encft_a_match %>%
  left_join(sat_annual,     by = c("muni_key", "ANO" = "year")) %>%
  left_join(instr_annual,   by = c("ANO" = "year")) %>%
  left_join(coastal_lookup, by = "muni_key") %>%
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
saveRDS(panel,        file.path(path_processed, "panel_monthly.rds"))
saveRDS(panel_annual, file.path(path_processed, "panel_analysis.rds"))
write_csv(panel,        file.path(path_processed, "panel_monthly.csv"))
write_csv(panel_annual, file.path(path_processed, "panel_analysis.csv"))
cat("\nSaved panel_monthly.rds (main) and panel_analysis.rds (annual robustness)\n")
