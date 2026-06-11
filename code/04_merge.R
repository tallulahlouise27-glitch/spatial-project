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
cat("Loading ENCFT household-month panel (main)...\n")
hogar_m <- readRDS(file.path(path_processed, "encft_hogar_monthly.rds"))

cat("Loading ENCFT quarterly municipality panel (robustness)...\n")
encft_q <- readRDS(file.path(path_processed, "encft_municipio_quarterly.rds"))

cat("Loading ENCFT annual panel (robustness)...\n")
encft_a <- readRDS(file.path(path_processed, "encft_municipio.rds"))

cat("Loading satellite coastal data (monthly)...\n")
sat_coastal <- readRDS(file.path(path_processed, "satellite_coastal.rds"))

cat("Loading satellite instrument data (monthly)...\n")
sat_instr <- readRDS(file.path(path_processed, "satellite_instrument.rds"))

cat("ENCFT household-month rows:", nrow(hogar_m), "\n")
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

hogar_m_match <- apply_crosswalk(hogar_m, year_col = "ANO")
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

# Compute minimum distance from each municipality boundary to the ocean coastline.
# Municipalities > 20km from the ocean cannot have legitimate ocean Sargassum values:
# their 20km satellite buffers only capture inland water bodies (Lago Enriquillo,
# Lago Azuei/Saumâtre) rather than ocean pixels. Their afai_sargassum values
# will be set to NA and excluded from all regressions.
ocean_dist_lookup <- dr_sf %>%
  mutate(
    muni_key     = clean_name(municipio),
    ocean_dist_m = as.numeric(st_distance(dr_sf, coastline))
  ) %>%
  st_drop_geometry() %>%
  select(muni_key, ocean_dist_m)

n_contaminated <- ocean_dist_lookup %>%
  filter(ocean_dist_m > 20000) %>%
  nrow()
cat("Municipalities > 20km from ocean (Sargassum values will be set NA):", n_contaminated, "\n")

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

# ── Step 3: Build household-level monthly panel ────────────────
# One row per household per month. Satellite treatment is merged in
# at the municipality-month level — all households in the same
# municipality-month receive the same Sargassum value.
# Standard errors must be clustered at municipality level because
# the treatment does not vary within municipality-months.

sat_monthly <- sat_coastal %>%
  mutate(
    muni_key = clean_name(municipio),
    quarter  = ceiling(month / 3)
  ) %>%
  select(muni_key, year, month, quarter,
         afai_sargassum, afai_coverage, afai_mean)

# ── IDW assignment for non-coastal municipalities ─────────────
# Each non-coastal municipality receives a weighted average of ALL
# coastal municipalities' Sargassum values, weighted by 1/distance^2.
# This avoids the hub-and-spoke problem of nearest-feature matching
# (where one coastal municipality dominates many inland ones) and gives
# each inland municipality a smooth, geographically sensible exposure value.

# Pepillo Salcedo is geographically coastal but its satellite data is
# contaminated (river plume, sensor clipping). Treat it as non-coastal
# so it receives IDW-weighted values from its neighbours instead.
contaminated_munis <- c("pepillo salcedo")

coastal_keys     <- coastal_lookup %>%
  filter(coastal_type == "coastal", !muni_key %in% contaminated_munis) %>%
  pull(muni_key)
non_coastal_keys <- coastal_lookup %>%
  filter(coastal_type == "not_coastal" | muni_key %in% contaminated_munis) %>%
  pull(muni_key)

coastal_geom <- dr_sf %>%
  mutate(muni_key = clean_name(municipio)) %>%
  filter(muni_key %in% coastal_keys) %>%
  st_transform(32619)

non_coastal_geom <- dr_sf %>%
  mutate(muni_key = clean_name(municipio)) %>%
  filter(muni_key %in% non_coastal_keys) %>%
  st_transform(32619)

# Distance matrix (metres): rows = non-coastal, cols = coastal
dist_matrix <- st_distance(
  st_centroid(non_coastal_geom),
  st_centroid(coastal_geom)
)
dist_mat <- matrix(as.numeric(dist_matrix),
                   nrow = nrow(non_coastal_geom),
                   ncol = nrow(coastal_geom))

# IDW weights: w = 1/d^2, normalised so each non-coastal row sums to 1
weight_mat <- 1 / dist_mat^2
weight_mat <- weight_mat / rowSums(weight_mat)

# Long-form weight table: one row per non-coastal × coastal pair
idw_weights <- expand_grid(
  nc_idx = seq_len(nrow(non_coastal_geom)),
  c_idx  = seq_len(nrow(coastal_geom))
) %>%
  mutate(
    muni_key       = non_coastal_geom$muni_key[nc_idx],
    source_coastal = coastal_geom$muni_key[c_idx],
    weight         = weight_mat[cbind(nc_idx, c_idx)]
  ) %>%
  select(muni_key, source_coastal, weight)

cat("IDW weights computed:", nrow(idw_weights), "non-coastal × coastal pairs\n")
cat("Non-coastal municipalities:", n_distinct(idw_weights$muni_key), "\n\n")

sat_coastal_vals <- sat_monthly %>% filter(muni_key %in% coastal_keys)

# For each non-coastal municipality × month, compute IDW-weighted Sargassum
sat_non_coastal_vals <- idw_weights %>%
  left_join(
    sat_coastal_vals %>% rename(source_coastal = muni_key),
    by = "source_coastal",
    relationship = "many-to-many"
  ) %>%
  group_by(muni_key, year, month, quarter) %>%
  summarise(
    afai_sargassum = weighted.mean(afai_sargassum, w = weight, na.rm = TRUE),
    afai_coverage  = weighted.mean(afai_coverage,  w = weight, na.rm = TRUE),
    afai_mean      = weighted.mean(afai_mean,      w = weight, na.rm = TRUE),
    .groups = "drop"
  )

sat_monthly <- bind_rows(sat_coastal_vals, sat_non_coastal_vals)
cat("Satellite rows after IDW (coastal + non-coastal):", nrow(sat_monthly), "\n\n")

# Build satellite + instrument lookup: one row per muni-month
sat_lookup <- sat_monthly %>%
  left_join(sat_instr,        by = c("year", "month")) %>%
  left_join(coastal_lookup,   by = "muni_key") %>%
  left_join(ocean_dist_lookup, by = "muni_key")

panel <- hogar_m_match %>%
  rename(month = MES, year = ANO) %>%
  mutate(TRIMESTRE = NULL) %>%
  inner_join(sat_lookup, by = c("muni_key", "year", "month"))

cat("Household-month panel before filtering:", nrow(panel), "rows\n")

# ── Step 3b: Quarterly robustness panel ───────────────────────
# Collapses Sargassum to quarterly averages and joins quarterly ENCFT.
# Use to verify that monthly results are not driven by within-quarter
# survey noise. Fixed effects in robustness: municipality + year×quarter.

sat_quarterly_raw <- sat_coastal %>%
  mutate(
    muni_key = clean_name(municipio),
    quarter  = ceiling(month / 3)
  ) %>%
  group_by(muni_key, year, quarter) %>%
  summarise(
    afai_sargassum = mean(afai_sargassum, na.rm = TRUE),
    afai_coverage  = mean(afai_coverage,  na.rm = TRUE),
    afai_mean      = mean(afai_mean,      na.rm = TRUE),
    .groups = "drop"
  )

# Apply same IDW weights to quarterly panel
sat_q_coastal_vals <- sat_quarterly_raw %>% filter(muni_key %in% coastal_keys)

sat_q_non_coastal_vals <- idw_weights %>%
  left_join(
    sat_q_coastal_vals %>% rename(source_coastal = muni_key),
    by = "source_coastal",
    relationship = "many-to-many"
  ) %>%
  group_by(muni_key, year, quarter) %>%
  summarise(
    afai_sargassum = weighted.mean(afai_sargassum, w = weight, na.rm = TRUE),
    afai_coverage  = weighted.mean(afai_coverage,  w = weight, na.rm = TRUE),
    afai_mean      = weighted.mean(afai_mean,      w = weight, na.rm = TRUE),
    .groups = "drop"
  )

sat_quarterly <- bind_rows(sat_q_coastal_vals, sat_q_non_coastal_vals)

instr_quarterly <- sat_instr %>%
  mutate(quarter = ceiling(month / 3)) %>%
  group_by(year, quarter) %>%
  summarise(
    afai_ocean_mean     = mean(afai_ocean_mean,     na.rm = TRUE),
    afai_ocean_coverage = mean(afai_ocean_coverage, na.rm = TRUE),
    .groups = "drop"
  )

panel_quarterly <- sat_quarterly %>%
  left_join(
    encft_q_match %>% mutate(TRIMESTRE = as.integer(TRIMESTRE) %% 10),
    by = c("muni_key", "year" = "ANO", "quarter" = "TRIMESTRE")
  ) %>%
  left_join(instr_quarterly,  by = c("year", "quarter")) %>%
  left_join(coastal_lookup,   by = "muni_key") %>%
  left_join(ocean_dist_lookup, by = "muni_key") %>%
  filter(!is.na(ID_MUNICIPIO))

cat("Quarterly robustness panel before filtering:", nrow(panel_quarterly), "rows\n")

# ── Step 4: Create analysis variables ─────────────────────────
add_common_vars <- function(df) {
  df %>%
    mutate(
      afai_sargassum  = if_else(ocean_dist_m > 20000, NA_real_, afai_sargassum),
      muni_fe         = factor(ID_MUNICIPIO),
      year_fe         = factor(as.integer(.data[["year"]])),
      quarter_fe      = factor(as.integer(.data[["quarter"]])),
      year_month_fe   = factor(paste0(as.integer(.data[["year"]]), "_m",
                                      sprintf("%02d", as.integer(.data[["month"]])))),
      year_quarter_fe = factor(paste0(as.integer(.data[["year"]]), "_q",
                                      as.integer(.data[["quarter"]])))
    )
}

# Household-level panel variables
panel <- panel %>%
  add_common_vars() %>%
  mutate(
    log_income    = log(ingreso_hogar + 1),
    log_income_pc = log(ingreso_pc    + 1)
  ) %>%
  filter(!is.na(log_income), !is.na(afai_sargassum), !is.na(afai_ocean_coverage))

# Municipality-level quarterly robustness panel variables
panel_quarterly <- panel_quarterly %>%
  mutate(month = NA_integer_) %>%
  add_common_vars() %>%
  mutate(
    log_income    = log(ingreso_medio + 1),
    log_income_t1 = log(ingreso_T1 + 1),
    log_income_t2 = log(ingreso_T2 + 1),
    log_income_t3 = log(ingreso_T3 + 1)
  ) %>%
  filter(!is.na(log_income), !is.na(afai_sargassum), !is.na(afai_ocean_coverage))

cat("Final household-month panel:", nrow(panel), "observations\n")
cat("Households:         ", n_distinct(panel$ID_HOGAR), "\n")
cat("Municipalities:     ", n_distinct(panel$ID_MUNICIPIO), "\n")
cat("Years:              ", paste(sort(unique(panel$year)), collapse = ", "), "\n")
cat("Log income mean:    ", round(mean(panel$log_income, na.rm = TRUE), 3), "\n")
cat("\nMunicipalities by coastal type:\n")
print(table(panel$coastal_type, useNA = "ifany"))

cat("\nQuarterly robustness panel:", nrow(panel_quarterly), "observations\n")
cat("Municipalities:     ", n_distinct(panel_quarterly$ID_MUNICIPIO), "\n")

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
saveRDS(panel,           file.path(path_processed, "panel_monthly.rds"))
saveRDS(panel_quarterly, file.path(path_processed, "panel_quarterly.rds"))
saveRDS(panel_annual,    file.path(path_processed, "panel_analysis.rds"))
write_csv(panel,           file.path(path_processed, "panel_monthly.csv"))
write_csv(panel_quarterly, file.path(path_processed, "panel_quarterly.csv"))
write_csv(panel_annual,    file.path(path_processed, "panel_analysis.csv"))
cat("\nSaved panel_monthly.rds (main), panel_quarterly.rds (robustness), panel_analysis.rds (annual robustness)\n")
