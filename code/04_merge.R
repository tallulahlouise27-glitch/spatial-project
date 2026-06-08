# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 4: Merge ENCFT survey data with satellite data
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)

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
    afai_mean_annual = mean(afai_mean,     na.rm = TRUE),
    afai_cov_annual  = mean(afai_coverage, na.rm = TRUE),
    # Peak season (May–Sep) when Sargassum is most intense in Caribbean
    afai_peak        = mean(afai_coverage[month %in% 5:9], na.rm = TRUE),
    .groups = "drop"
  )

instr_annual <- sat_instr %>%
  group_by(year) %>%
  summarise(
    afai_ocean_annual = mean(afai_ocean_mean,     na.rm = TRUE),
    afai_ocean_cov    = mean(afai_ocean_coverage, na.rm = TRUE),
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

# ── Step 3: Merge ─────────────────────────────────────────────
panel <- encft_match %>%
  left_join(sat_match, by = c("muni_key", "ANO" = "year")) %>%
  left_join(instr_annual, by = c("ANO" = "year"))

cat("Merge result:", nrow(panel), "rows\n")
cat("Municipalities matched with satellite data:",
    sum(!is.na(panel$chla_mean_annual)), "of", nrow(panel), "\n")

# Municipalities with NO satellite data in any year (truly unmatched)
truly_unmatched <- panel %>%
  group_by(DES_MUNICIPIO, DES_PROVINCIA) %>%
  summarise(any_sat = any(!is.na(chla_mean_annual)), .groups = "drop") %>%
  filter(!any_sat)

cat("\nMunicipalities with no satellite data at all:", nrow(truly_unmatched),
    "(likely landlocked — will be dropped from regression)\n")

# ── Step 4: Create analysis variables ─────────────────────────
panel <- panel %>%
  mutate(
    # Log income (add 1 to handle zeros)
    log_income    = log(ingreso_medio + 1),
    # Sargassum exposure: coverage fraction (0–1) needs no log transform
    # afai_cov_annual  = annual mean fraction of coastal pixels with Sargassum
    # afai_peak        = same, restricted to peak season May–Sep
    # afai_ocean_cov   = open-ocean coverage (instrument shift)
    # Year and municipality as factors for fixed effects
    year_fe  = factor(ANO),
    muni_fe  = factor(ID_MUNICIPIO)
  ) %>%
  # Keep only rows with all key variables present
  filter(!is.na(log_income), !is.na(afai_cov_annual), !is.na(afai_ocean_cov))

cat("\nFinal analysis panel:", nrow(panel), "observations\n")
cat("Municipalities:       ", n_distinct(panel$ID_MUNICIPIO), "\n")
cat("Years:                ", paste(sort(unique(panel$ANO)), collapse = ", "), "\n")
cat("Employment rate mean: ", round(mean(panel$tasa_empleo, na.rm=TRUE), 3), "\n")
cat("Log income mean:      ", round(mean(panel$log_income, na.rm=TRUE), 3), "\n")

# ── Step 5: Save ──────────────────────────────────────────────
saveRDS(panel, file.path(path_processed, "panel_analysis.rds"))
write_csv(panel, file.path(path_processed, "panel_analysis.csv"))
cat("\nSaved to data/processed/panel_analysis.rds\n")
