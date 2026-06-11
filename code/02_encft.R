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
  "RAMA_PRINCIPAL", "RAMA_PRINCIPAL_COD",
  "GRUPO_SECTOR"
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
    MES          = as.integer(MES),
    ID_PROVINCIA = as.integer(ID_PROVINCIA),
    ID_MUNICIPIO = as.integer(ID_MUNICIPIO),
    ingreso_total = coalesce(as.numeric(INGRESO_ASALARIADO), 0) +
                    coalesce(as.numeric(INGRESO_INDEPENDIENTES), 0),
    # 311 = marine fishing, 312 = freshwater fishing (CIIU/ISIC codes)
    # PESCA_NO_REMUN only captures unpaid subsistence fishing (near-zero) — use occupation code instead
    es_pescador   = as.integer(as.numeric(RAMA_PRINCIPAL_COD) %in% c(311, 312)),
    # GRUPO_SECTOR is the ONE's official formal/informal classification.
    # "Informal" + "Servicio Doméstico" = informal; "Formal" = formal;
    # "Sin sector" = not employed (unemployed/inactive) — coded NA.
    es_informal   = case_when(
      GRUPO_SECTOR %in% c("Informal", "Servicio Doméstico") ~ 1L,
      GRUPO_SECTOR == "Formal"                                    ~ 0L,
      TRUE                                                        ~ NA_integer_
    )
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

# ── Assign national income tertiles by year ───────────────────
# Per-capita household income is the welfare measure for ranking.
# Tertile cutoffs are computed nationally within each year so that
# T1 = poorest third of DR households, T2 = middle, T3 = richest third.
# Zero-income households are included — they fall naturally into T1.

hogar <- hogar %>%
  mutate(ingreso_pc = ingreso_hogar / pmax(n_miembros, 1))

tertile_cuts <- hogar %>%
  group_by(ANO) %>%
  summarise(
    t1_cut = quantile(ingreso_pc, 1/3, na.rm = TRUE),
    t2_cut = quantile(ingreso_pc, 2/3, na.rm = TRUE),
    .groups = "drop"
  )

hogar <- hogar %>%
  left_join(tertile_cuts, by = "ANO") %>%
  mutate(
    tertile = case_when(
      ingreso_pc <= t1_cut ~ "T1",
      ingreso_pc <= t2_cut ~ "T2",
      TRUE                 ~ "T3"
    )
  )

cat("Tertile cutoffs by year (per-capita RD$/month):\n")
print(tertile_cuts)

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

# Merge tertile mean incomes into municipality panel
municipio_tertiles <- hogar %>%
  group_by(ID_MUNICIPIO, DES_MUNICIPIO, ID_PROVINCIA, DES_PROVINCIA, ANO, tertile) %>%
  summarise(
    ingreso_tertile = weighted.mean(ingreso_hogar, w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = tertile,
    values_from = ingreso_tertile,
    names_prefix = "ingreso_"
  )

municipio <- municipio %>%
  left_join(municipio_tertiles,
            by = c("ID_MUNICIPIO", "DES_MUNICIPIO", "ID_PROVINCIA", "DES_PROVINCIA", "ANO"))

# ── Aggregate to municipality × year × quarter ────────────────
# ENCFT is designed to be representative at the quarterly level.
# This aggregation preserves that structure and allows matching
# with quarterly satellite data (3-month average AFAI per quarter).

municipio_quarterly <- encft %>%
  group_by(ID_MUNICIPIO, DES_MUNICIPIO, ID_PROVINCIA, DES_PROVINCIA, ANO, TRIMESTRE) %>%
  summarise(
    n_hogares     = n_distinct(ID_HOGAR),
    ingreso_medio = weighted.mean(ingreso_total, w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    tasa_empleo   = weighted.mean(as.integer(OCUPADO) == 1,
                                  w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    share_pesca   = weighted.mean(es_pescador,
                                  w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    .groups = "drop"
  )

# Tertile income by quarter: use annual national tertile assignments
# (tertile membership is a property of each household's annual income
# relative to the national distribution, not quarter-specific)
municipio_tertiles_q <- hogar %>%
  group_by(ID_MUNICIPIO, DES_MUNICIPIO, ID_PROVINCIA, DES_PROVINCIA, ANO, TRIMESTRE, tertile) %>%
  summarise(
    ingreso_tertile = weighted.mean(ingreso_hogar, w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = tertile,
    values_from = ingreso_tertile,
    names_prefix = "ingreso_"
  )

municipio_quarterly <- municipio_quarterly %>%
  left_join(municipio_tertiles_q,
            by = c("ID_MUNICIPIO", "DES_MUNICIPIO", "ID_PROVINCIA", "DES_PROVINCIA",
                   "ANO", "TRIMESTRE"))

# ── Aggregate to municipality × year × month ──────────────────
# Uses MES to capture genuine within-quarter monthly income variation.
# Each household is surveyed once per quarter; MES records which month.
# Income reported reflects conditions at interview time, so monthly
# averages capture real within-quarter variation — not just repeated
# quarterly values. Monthly samples are smaller than quarterly, so
# estimates are noisier for small municipalities.

hogar_monthly <- encft %>%
  filter(!is.na(MES)) %>%
  group_by(ID_HOGAR, ANO, TRIMESTRE, MES,
           ID_PROVINCIA, DES_PROVINCIA,
           ID_MUNICIPIO, DES_MUNICIPIO,
           ZONA, FACTOR_EXPANSION) %>%
  summarise(
    n_miembros    = n(),
    ingreso_hogar = sum(ingreso_total, na.rm = TRUE),
    tasa_empleo   = mean(as.integer(OCUPADO) == 1, na.rm = TRUE),
    hogar_pesca   = as.integer(any(es_pescador == 1, na.rm = TRUE)),
    # Formal/informal sector classification at household level.
    # n_empleados: members with a known sector (formal OR informal).
    # share_informal: share of those who are informal/domestic workers.
    n_empleados   = sum(!is.na(es_informal)),
    n_informales  = sum(es_informal == 1, na.rm = TRUE),
    share_informal = if_else(
      n_empleados > 0,
      n_informales / n_empleados,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(ingreso_pc = ingreso_hogar / pmax(n_miembros, 1)) %>%
  left_join(tertile_cuts, by = "ANO") %>%
  mutate(
    tertile = case_when(
      ingreso_pc <= t1_cut ~ "T1",
      ingreso_pc <= t2_cut ~ "T2",
      TRUE                 ~ "T3"
    )
  )

municipio_monthly <- hogar_monthly %>%
  group_by(ID_MUNICIPIO, DES_MUNICIPIO, ID_PROVINCIA, DES_PROVINCIA, ANO, TRIMESTRE, MES) %>%
  summarise(
    n_hogares     = n_distinct(ID_HOGAR),
    ingreso_medio = weighted.mean(ingreso_hogar, w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    tasa_empleo   = weighted.mean(tasa_empleo,   w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    share_pesca   = weighted.mean(hogar_pesca,   w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    .groups = "drop"
  )

municipio_tertiles_m <- hogar_monthly %>%
  group_by(ID_MUNICIPIO, DES_MUNICIPIO, ID_PROVINCIA, DES_PROVINCIA,
           ANO, TRIMESTRE, MES, tertile) %>%
  summarise(
    ingreso_tertile = weighted.mean(ingreso_hogar, w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from   = tertile,
    values_from  = ingreso_tertile,
    names_prefix = "ingreso_"
  )

municipio_monthly <- municipio_monthly %>%
  left_join(municipio_tertiles_m,
            by = c("ID_MUNICIPIO", "DES_MUNICIPIO", "ID_PROVINCIA", "DES_PROVINCIA",
                   "ANO", "TRIMESTRE", "MES"))

cat("Monthly municipality rows:", nrow(municipio_monthly), "\n")
cat("Median households per municipality-month:",
    median(municipio_monthly$n_hogares, na.rm = TRUE), "\n")
cat("Min households per municipality-month:",
    min(municipio_monthly$n_hogares, na.rm = TRUE), "\n\n")

# ── Save ─────────────────────────────────────────────────────
saveRDS(hogar,               file.path(path_processed, "encft_hogar.rds"))
saveRDS(hogar_monthly,       file.path(path_processed, "encft_hogar_monthly.rds"))
saveRDS(municipio,           file.path(path_processed, "encft_municipio.rds"))
saveRDS(municipio_quarterly, file.path(path_processed, "encft_municipio_quarterly.rds"))
saveRDS(municipio_monthly,   file.path(path_processed, "encft_municipio_monthly.rds"))
write_csv(municipio,           file.path(path_processed, "encft_municipio.csv"))
write_csv(municipio_quarterly, file.path(path_processed, "encft_municipio_quarterly.csv"))
write_csv(municipio_monthly,   file.path(path_processed, "encft_municipio_monthly.csv"))

# ── Summary ──────────────────────────────────────────────────
cat("\n============ ENCFT CLEAN SUMMARY ============\n")
cat("Household-quarter rows:       ", nrow(hogar), "\n")
cat("Municipality-year rows:       ", nrow(municipio), "\n")
cat("Municipality-year-quarter rows:", nrow(municipio_quarterly), "\n")
cat("Years covered:           ", paste(sort(unique(municipio$ANO)), collapse = ", "), "\n")
cat("Provinces:               ", n_distinct(municipio$ID_PROVINCIA), "\n")
cat("Municipalities:          ", n_distinct(municipio$ID_MUNICIPIO), "\n")
cat("Avg income (municipality level, unweighted):\n")
print(summary(municipio$ingreso_medio))
cat("Avg employment rate:     ", round(mean(municipio$tasa_empleo, na.rm = TRUE), 3), "\n")
cat("Avg fishing share:       ", round(mean(municipio$share_pesca, na.rm = TRUE), 3), "\n")
cat("Avg income by tertile (municipality means):\n")
cat("  T1 (bottom third): RD$", round(mean(municipio$ingreso_T1, na.rm=TRUE)), "\n")
cat("  T2 (middle third): RD$", round(mean(municipio$ingreso_T2, na.rm=TRUE)), "\n")
cat("  T3 (top third):    RD$", round(mean(municipio$ingreso_T3, na.rm=TRUE)), "\n")
cat("=============================================\n")
