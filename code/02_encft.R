# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 2: Load and clean ENCFT data (2016–2026)
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(readxl)

path_raw       <- "data/raw/"
path_processed <- "data/processed/"

# Variables to keep from the Miembros sheet
vars_keep <- c(
  "TRIMESTRE", "ANO", "MES",
  "ID_PROVINCIA", "DES_PROVINCIA", "ID_MUNICIPIO", "DES_MUNICIPIO", "ZONA",
  "ID_HOGAR", "MIEMBRO",
  "FACTOR_EXPANSION",
  "SEXO", "EDAD",
  "OCUPADO", "DESOCUPADO", "INACTIVO", "PEA",
  "INGRESO_ASALARIADO", "INGRESO_INDEPENDIENTES",
  "PESCA_NO_REMUN", "PESCA_NO_REMUN_MONTO",
  "RAMA_PRINCIPAL", "RAMA_PRINCIPAL_COD"
)

# Find the .xlsx file inside each year's folder
get_filepath <- function(year) {
  folder <- file.path(path_raw, paste0("Base ENCFT ", year))
  files  <- list.files(folder, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(files) == 0) stop("No xlsx found for year: ", year)
  files[1]
}

# Load one year's Miembros sheet, keep only needed columns
load_year <- function(year) {
  f <- get_filepath(year)
  message("Loading ", year, ": ", basename(f))
  df   <- read_excel(f, sheet = "Miembros")
  cols <- intersect(vars_keep, names(df))
  df   <- df[, cols]
  df$year_file <- year
  # Convert everything to character so all years stack cleanly
  # (column types vary across survey waves — we recode to numeric after stacking)
  df <- mutate(df, across(everything(), as.character))
  df
}

# Load all years and stack
years      <- 2016:2026
encft_raw  <- bind_rows(lapply(years, load_year))

message("\nRows loaded: ", nrow(encft_raw))

# ── Clean individual-level data ──────────────────────────────
encft <- encft_raw %>%
  mutate(
    ANO          = as.integer(ANO),
    TRIMESTRE    = as.integer(TRIMESTRE),
    ID_PROVINCIA = as.integer(ID_PROVINCIA),
    ID_MUNICIPIO = as.integer(ID_MUNICIPIO),
    ingreso_total = coalesce(as.numeric(INGRESO_ASALARIADO), 0) +
                    coalesce(as.numeric(INGRESO_INDEPENDIENTES), 0),
    es_pescador   = as.integer(!is.na(PESCA_NO_REMUN) & PESCA_NO_REMUN == 1)
  )

# ── Aggregate to household × quarter level ───────────────────
hogar <- encft %>%
  group_by(ID_HOGAR, ANO, TRIMESTRE,
           ID_PROVINCIA, DES_PROVINCIA,
           ID_MUNICIPIO, DES_MUNICIPIO,
           ZONA, FACTOR_EXPANSION) %>%
  summarise(
    n_miembros    = n(),
    ingreso_hogar = sum(ingreso_total, na.rm = TRUE),
    n_ocupados    = sum(as.integer(OCUPADO) == 1, na.rm = TRUE),
    tasa_empleo   = mean(as.integer(OCUPADO) == 1, na.rm = TRUE),
    hogar_pesca   = as.integer(any(es_pescador == 1, na.rm = TRUE)),
    .groups = "drop"
  )

# ── Aggregate to municipality × year level (for spatial merge) ──
municipio <- encft %>%
  group_by(ID_MUNICIPIO, DES_MUNICIPIO, ID_PROVINCIA, DES_PROVINCIA, ANO) %>%
  summarise(
    n_hogares          = n_distinct(ID_HOGAR),
    ingreso_medio      = weighted.mean(ingreso_total, w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    tasa_empleo        = weighted.mean(as.integer(OCUPADO) == 1,
                                       w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    share_pesca        = weighted.mean(es_pescador,
                                       w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    .groups = "drop"
  )

# ── Save ─────────────────────────────────────────────────────
saveRDS(hogar,     file.path(path_processed, "encft_hogar.rds"))
saveRDS(municipio, file.path(path_processed, "encft_municipio.rds"))
write_csv(municipio, file.path(path_processed, "encft_municipio.csv"))

# ── Summary ──────────────────────────────────────────────────
cat("\n============ ENCFT CLEAN SUMMARY ============\n")
cat("Household-quarter rows:  ", nrow(hogar), "\n")
cat("Municipality-year rows:  ", nrow(municipio), "\n")
cat("Years covered:           ", paste(sort(unique(municipio$ANO)), collapse = ", "), "\n")
cat("Provinces:               ", n_distinct(municipio$ID_PROVINCIA), "\n")
cat("Municipalities:          ", n_distinct(municipio$ID_MUNICIPIO), "\n")
cat("Avg income (municipality level, unweighted):\n")
print(summary(municipio$ingreso_medio))
cat("Avg employment rate:     ", round(mean(municipio$tasa_empleo, na.rm = TRUE), 3), "\n")
cat("Avg fishing share:       ", round(mean(municipio$share_pesca, na.rm = TRUE), 3), "\n")
cat("=============================================\n")
