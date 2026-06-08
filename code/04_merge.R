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
    chla_mean_annual = mean(chla_coastal_mean, na.rm = TRUE),
    chla_max_annual  = max(chla_coastal_max,  na.rm = TRUE),
    # Peak season (May–Sep) when Sargassum is most intense in Caribbean
    chla_peak = mean(chla_coastal_mean[month %in% 5:9], na.rm = TRUE),
    .groups = "drop"
  )

instr_annual <- sat_instr %>%
  group_by(year) %>%
  summarise(
    chla_ocean_annual = mean(chla_ocean_mean, na.rm = TRUE),
    chla_ocean_peak   = mean(chla_ocean_mean[month %in% 5:9], na.rm = TRUE),
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

encft_match <- encft %>%
  mutate(
    ANO          = as.integer(ANO),
    muni_key     = clean_name(DES_MUNICIPIO),
    prov_key     = clean_name(DES_PROVINCIA)
  )

sat_match <- sat_annual %>%
  mutate(
    muni_key = clean_name(municipio),
    prov_key = clean_name(provincia)
  )

# ── Step 3: Merge ─────────────────────────────────────────────
panel <- encft_match %>%
  left_join(sat_match, by = c("muni_key", "prov_key", "ANO" = "year")) %>%
  left_join(instr_annual, by = c("ANO" = "year"))

cat("Merge result:", nrow(panel), "rows\n")
cat("Municipalities matched with satellite data:",
    sum(!is.na(panel$chla_mean_annual)), "of", nrow(panel), "\n")

# Flag unmatched municipalities for inspection
unmatched <- panel %>%
  filter(is.na(chla_mean_annual)) %>%
  distinct(DES_MUNICIPIO, DES_PROVINCIA) %>%
  arrange(DES_PROVINCIA, DES_MUNICIPIO)

if (nrow(unmatched) > 0) {
  cat("\nUnmatched municipalities (", nrow(unmatched), "):\n")
  print(unmatched)
}

# ── Step 4: Create analysis variables ─────────────────────────
panel <- panel %>%
  mutate(
    # Log income (add 1 to handle zeros)
    log_income    = log(ingreso_medio + 1),
    # Sargassum exposure variables (log + 1 to handle near-zero values)
    log_chla      = log(chla_mean_annual + 1),
    log_chla_peak = log(chla_peak + 1),
    log_chla_ocean= log(chla_ocean_annual + 1),
    # Year and municipality as factors for fixed effects
    year_fe  = factor(ANO),
    muni_fe  = factor(ID_MUNICIPIO)
  ) %>%
  # Keep only rows with all key variables present
  filter(!is.na(log_income), !is.na(log_chla), !is.na(log_chla_ocean))

cat("\nFinal analysis panel:", nrow(panel), "observations\n")
cat("Municipalities:       ", n_distinct(panel$ID_MUNICIPIO), "\n")
cat("Years:                ", paste(sort(unique(panel$ANO)), collapse = ", "), "\n")
cat("Employment rate mean: ", round(mean(panel$tasa_empleo, na.rm=TRUE), 3), "\n")
cat("Log income mean:      ", round(mean(panel$log_income, na.rm=TRUE), 3), "\n")

# ── Step 5: Save ──────────────────────────────────────────────
saveRDS(panel, file.path(path_processed, "panel_analysis.rds"))
write_csv(panel, file.path(path_processed, "panel_analysis.csv"))
cat("\nSaved to data/processed/panel_analysis.rds\n")
