# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 9: Diagnostic figures for checking the analysis
#
# Produces figures/diagnostics/ with 7 checks:
#   1. Map: which municipalities fall in which income tertile?
#   2. Map: Sargassum values INCLUDING non-coastal municipalities
#   3. Time series: Sargassum by coastal vs non-coastal
#   4. Raw scatter: peak-season Sargassum vs peak-season income
#   5. Panel coverage: which municipality × year cells have data?
#   6. Income trends over time by tertile group
#   7. Coefficient comparison: original vs non-zero tertiles
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(sf)
library(geodata)
library(fixest)

path_processed <- "data/processed/"
path_raw_sat   <- "data/raw/satellite/"
path_diag      <- "figures/diagnostics/"
dir.create(path_diag, showWarnings = FALSE, recursive = TRUE)

clean_name <- function(x) {
  x <- tolower(x)
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  x <- gsub("[^a-z0-9 ]", "", x)
  x <- trimws(x)
  gsub("[ ]+", " ", x)
}

# ── Load data ─────────────────────────────────────────────────
cat("Loading data...\n")
panel <- readRDS(file.path(path_processed, "panel_monthly.rds"))
sat   <- readRDS(file.path(path_processed, "satellite_coastal.rds"))

dr_gadm <- gadm(country = "DOM", level = 2, path = path_raw_sat)
dr_sf   <- st_as_sf(dr_gadm) %>%
  select(municipio = NAME_2, provincia = NAME_1) %>%
  st_transform(4326) %>%
  mutate(muni_key = clean_name(municipio))

# Coastal classification (same as script 04)
dr_proj    <- st_transform(dr_sf, 32619)
haiti_gadm <- gadm(country = "HTI", level = 0, path = path_raw_sat)
haiti_sf   <- st_as_sf(haiti_gadm) %>% st_transform(32619)
haiti_buf  <- st_buffer(st_union(haiti_sf), dist = 2000)
coastline  <- st_difference(st_union(dr_proj) %>% st_boundary(), haiti_buf)
coast_zone <- st_buffer(coastline, dist = 500)
touches    <- lengths(st_intersects(dr_proj, coast_zone)) > 0

dr_sf <- dr_sf %>%
  mutate(coastal_type = if_else(touches, "Coastal", "Not coastal"))

cat("Data loaded.\n\n")

# ═══════════════════════════════════════════════════════════════
# FIGURE 1: Map — which income tertile does each municipality
#           typically fall in?
# ═══════════════════════════════════════════════════════════════
cat("Figure 1: Tertile assignment map...\n")

# panel$muni_key is already the GADM-crosswalked name (built in script 04).
# Join directly on muni_key — do NOT use clean_name(DES_MUNICIPIO),
# which uses raw ENCFT names that don't match GADM.
muni_tertile <- panel %>%
  group_by(muni_key) %>%
  summarise(
    mean_income = mean(ingreso_medio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    avg_tertile = case_when(
      ntile(mean_income, 3) == 1 ~ "Bottom third\n(typically poorest)",
      ntile(mean_income, 3) == 2 ~ "Middle third",
      ntile(mean_income, 3) == 3 ~ "Top third\n(typically richest)"
    )
  )

map1 <- dr_sf %>%
  left_join(muni_tertile, by = "muni_key") %>%
  mutate(avg_tertile = replace_na(avg_tertile, "No data (not in panel)"))

p1 <- ggplot(map1) +
  geom_sf(aes(fill = avg_tertile), colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      "Bottom third\n(typically poorest)" = "#d73027",
      "Middle third"                      = "#fee090",
      "Top third\n(typically richest)"    = "#1a9850",
      "No data (not in panel)"            = "grey85"
    ),
    name = NULL,
    na.value = "grey85"
  ) +
  labs(
    title    = "CHECK 1: Average income level by municipality",
    subtitle = "Does this match your expectations? Tourist zones and Santo Domingo should be green.\nRural and remote areas should be red."
  ) +
  theme_void() +
  theme(
    legend.position   = "bottom",
    legend.text       = element_text(size = 9),
    plot.title        = element_text(size = 12, face = "bold", margin = margin(b = 4)),
    plot.subtitle     = element_text(size = 9, colour = "grey30"),
    plot.margin       = margin(10, 10, 10, 10)
  )

ggsave(file.path(path_diag, "check1_income_tertile_map.png"),
       p1, width = 7, height = 5.5, dpi = 150)
cat("  Saved: check1_income_tertile_map.png\n")

# ═══════════════════════════════════════════════════════════════
# FIGURE 2: Map — which municipalities were excluded because their
#           Sargassum values come from inland lakes, not the ocean?
#           Shows ocean distance and exclusion status.
# ═══════════════════════════════════════════════════════════════
cat("Figure 2: Contamination exclusion map...\n")

# Compute ocean distance for all GADM municipalities
dr_proj2   <- st_transform(dr_sf, 32619)
coastline2 <- st_difference(st_union(dr_proj2) %>% st_boundary(), haiti_buf)

ocean_dist <- dr_sf %>%
  st_transform(32619) %>%
  mutate(ocean_dist_km = as.numeric(st_distance(
    st_transform(dr_sf, 32619), coastline2)) / 1000) %>%
  st_drop_geometry() %>%
  select(muni_key, ocean_dist_km)

panel_munis <- panel %>% distinct(muni_key)

map2_data <- dr_sf %>%
  left_join(ocean_dist, by = "muni_key") %>%
  mutate(
    in_panel   = muni_key %in% panel_munis$muni_key,
    excl_label = case_when(
      !in_panel & ocean_dist_km > 20 ~ "Excluded: >20km from ocean\n(inland lake contamination)",
      !in_panel                      ~ "Not in panel (no ENCFT data)",
      ocean_dist_km <= 0.5           ~ "In panel: coastal",
      TRUE                           ~ "In panel: non-coastal (within 20km)"
    )
  )

p2 <- ggplot(map2_data) +
  geom_sf(aes(fill = excl_label), colour = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = c(
      "In panel: coastal"                              = "#2166ac",
      "In panel: non-coastal (within 20km)"            = "#a6d4f9",
      "Excluded: >20km from ocean\n(inland lake contamination)" = "#d73027",
      "Not in panel (no ENCFT data)"                   = "grey85"
    ),
    name = NULL
  ) +
  labs(
    title    = "CHECK 2: Which municipalities were excluded for inland lake contamination?",
    subtitle = "Red = excluded from analysis because their 'Sargassum' values came from Lago Enriquillo\nor Lago Azuei/Saumâtre (inland saltwater lakes), not the ocean.\nBlue = coastal (in analysis). Light blue = non-coastal but within 20km of ocean (in analysis)."
  ) +
  theme_void() +
  theme(
    legend.position   = "bottom",
    legend.text       = element_text(size = 8),
    plot.title        = element_text(size = 11, face = "bold", margin = margin(b = 4)),
    plot.subtitle     = element_text(size = 8, colour = "grey30"),
    plot.margin       = margin(10, 10, 10, 10)
  )

ggsave(file.path(path_diag, "check2_contamination_exclusion_map.png"),
       p2, width = 8, height = 6, dpi = 150)
cat("  Saved: check2_contamination_exclusion_map.png\n")

# ═══════════════════════════════════════════════════════════════
# FIGURE 3: Time series — Sargassum by coastal vs non-coastal
#           Coastal should be higher. If not, the spatial
#           assignment may be wrong.
# ═══════════════════════════════════════════════════════════════
cat("Figure 3: Sargassum time series by coastal type...\n")

ts_coast <- panel %>%
  mutate(
    date         = as.Date(paste(year, sprintf("%02d", month), "01", sep = "-")),
    coastal_label = if_else(coastal_type == "coastal", "Coastal", "Not coastal")
  ) %>%
  group_by(date, coastal_label) %>%
  summarise(mean_sarg = mean(afai_sargassum, na.rm = TRUE), .groups = "drop")

p3 <- ggplot(ts_coast, aes(x = date, y = mean_sarg, colour = coastal_label)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(
    values = c("Coastal" = "#2166ac", "Not coastal" = "#d7301f"),
    name   = NULL
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  scale_y_continuous(labels = scales::label_scientific()) +
  labs(
    x        = NULL,
    y        = "Mean Sargassum intensity (excess AFAI)",
    title    = "CHECK 3: Sargassum intensity — coastal vs non-coastal over time",
    subtitle = "Coastal (blue) should generally be higher than not-coastal (red).\nIf the lines cross repeatedly or not-coastal is higher, the spatial assignment may be wrong."
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    legend.position = "bottom"
  )

ggsave(file.path(path_diag, "check3_sargassum_by_coast_type.png"),
       p3, width = 9, height = 4.5, dpi = 150)
cat("  Saved: check3_sargassum_by_coast_type.png\n")

# ═══════════════════════════════════════════════════════════════
# FIGURE 4: Raw scatter — peak-season Sargassum vs income
#           This shows the raw data before any fixed effects.
#           There should be a visible negative slope.
# ═══════════════════════════════════════════════════════════════
cat("Figure 4: Raw scatter of Sargassum vs income (peak season)...\n")

scatter_data <- panel %>%
  filter(month %in% 5:9) %>%
  group_by(ID_MUNICIPIO, DES_MUNICIPIO, coastal_type) %>%
  summarise(
    mean_sarg   = mean(afai_sargassum, na.rm = TRUE),
    mean_income = mean(log_income,     na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(coastal_label = if_else(coastal_type == "coastal", "Coastal", "Not coastal"))

p4 <- ggplot(scatter_data, aes(x = mean_sarg, y = mean_income, colour = coastal_label)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
  scale_colour_manual(
    values = c("Coastal" = "#2166ac", "Not coastal" = "#d7301f"),
    name   = NULL
  ) +
  scale_x_continuous(labels = scales::label_scientific()) +
  labs(
    x        = "Mean peak-season Sargassum intensity (excess AFAI)",
    y        = "Mean log household income",
    title    = "CHECK 4: Raw relationship between Sargassum and income (peak season, no fixed effects)",
    subtitle = "Each point is one municipality. Lines are simple linear fits.\nA downward slope means more Sargassum is associated with lower income — before any controls."
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    legend.position = "bottom"
  )

ggsave(file.path(path_diag, "check4_raw_scatter.png"),
       p4, width = 8, height = 5.5, dpi = 150)
cat("  Saved: check4_raw_scatter.png\n")

# ═══════════════════════════════════════════════════════════════
# FIGURE 5: Panel coverage heatmap
#           Which municipality × year cells have data?
#           Large gaps could indicate data problems.
# ═══════════════════════════════════════════════════════════════
cat("Figure 5: Panel coverage heatmap...\n")

coverage <- panel %>%
  group_by(DES_MUNICIPIO, year) %>%
  summarise(n_months = n(), .groups = "drop") %>%
  mutate(
    DES_MUNICIPIO = str_trunc(str_to_title(DES_MUNICIPIO), 22),
    coverage      = case_when(
      n_months == 12 ~ "Full (12 months)",
      n_months >= 9  ~ "Partial (9–11)",
      n_months >= 5  ~ "Partial (5–8)",
      TRUE           ~ "Sparse (<5)"
    )
  )

# Order municipalities by total coverage
muni_order <- coverage %>%
  group_by(DES_MUNICIPIO) %>%
  summarise(total = sum(n_months)) %>%
  arrange(desc(total)) %>%
  pull(DES_MUNICIPIO)

coverage <- coverage %>%
  mutate(DES_MUNICIPIO = factor(DES_MUNICIPIO, levels = rev(muni_order)))

p5 <- ggplot(coverage, aes(x = factor(year), y = DES_MUNICIPIO, fill = coverage)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c(
      "Full (12 months)" = "#1a9850",
      "Partial (9–11)"   = "#fee08b",
      "Partial (5–8)"    = "#f46d43",
      "Sparse (<5)"      = "#d73027"
    ),
    name = "Months of\ndata"
  ) +
  labs(
    x        = NULL,
    y        = NULL,
    title    = "CHECK 5: Data coverage — months of data per municipality per year",
    subtitle = "Green = full year of data. Red = sparse. Large red areas suggest data gaps worth investigating."
  ) +
  theme_minimal() +
  theme(
    axis.text.y   = element_text(size = 6),
    axis.text.x   = element_text(size = 9),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    legend.position = "right"
  )

ggsave(file.path(path_diag, "check5_panel_coverage.png"),
       p5, width = 9, height = 14, dpi = 150)
cat("  Saved: check5_panel_coverage.png\n")

# ═══════════════════════════════════════════════════════════════
# FIGURE 6: Income over time by tertile
#           Do the three groups move in a sensible way?
#           All should trend upward; T3 should be highest.
# ═══════════════════════════════════════════════════════════════
cat("Figure 6: Income trends by tertile over time...\n")

income_trends <- panel %>%
  group_by(year) %>%
  summarise(
    bottom = mean(ingreso_T1, na.rm = TRUE),
    middle = mean(ingreso_T2, na.rm = TRUE),
    top    = mean(ingreso_T3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(bottom, middle, top),
               names_to = "tertile", values_to = "mean_income") %>%
  mutate(tertile = factor(tertile,
                          levels = c("top", "middle", "bottom"),
                          labels = c("Top third", "Middle third", "Bottom third")))

p6 <- ggplot(income_trends, aes(x = year, y = mean_income, colour = tertile)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(
    values = c("Top third" = "#1a9850", "Middle third" = "#fe9929", "Bottom third" = "#d73027"),
    name   = NULL
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = 2017:2025) +
  labs(
    x        = NULL,
    y        = "Mean household income (RD$/month)",
    title    = "CHECK 6: Income over time by tertile group",
    subtitle = "All groups should trend upward. Top > Middle > Bottom at every year.\nA COVID dip around 2020 would be expected. Anything else may indicate a data problem."
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 12, face = "bold"),
    plot.subtitle   = element_text(size = 9, colour = "grey30"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(path_diag, "check6_income_trends_by_tertile.png"),
       p6, width = 9, height = 5, dpi = 150)
cat("  Saved: check6_income_trends_by_tertile.png\n")

# ═══════════════════════════════════════════════════════════════
# FIGURE 7: Coefficient comparison — original vs non-zero tertiles
#           Do the two sets of results move in the expected
#           direction? Non-zero T1 should be more significant.
# ═══════════════════════════════════════════════════════════════
cat("Figure 7: Coefficient comparison, original vs non-zero tertiles...\n")

# Build non-zero tertile variables for the panel
hogar_nz <- readRDS(file.path(path_processed, "encft_hogar.rds")) %>%
  filter(ingreso_pc > 0)

cuts_nz <- hogar_nz %>%
  group_by(ANO) %>%
  summarise(
    t1_cut_nz = quantile(ingreso_pc, 1/3, na.rm = TRUE),
    t2_cut_nz = quantile(ingreso_pc, 2/3, na.rm = TRUE),
    .groups = "drop"
  )

nz_q <- hogar_nz %>%
  left_join(cuts_nz, by = "ANO") %>%
  mutate(
    tertile_nz = case_when(
      ingreso_pc <= t1_cut_nz ~ "T1",
      ingreso_pc <= t2_cut_nz ~ "T2",
      TRUE                    ~ "T3"
    ),
    quarter = as.integer(TRIMESTRE) %% 10L
  ) %>%
  group_by(ID_MUNICIPIO, ANO, quarter, tertile_nz) %>%
  summarise(
    ingreso_nz = weighted.mean(ingreso_hogar,
                               w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = tertile_nz, values_from = ingreso_nz,
              names_glue = "ingreso_{tertile_nz}_nz")

panel_full <- panel %>%
  left_join(nz_q, by = c("ID_MUNICIPIO", "year" = "ANO", "quarter")) %>%
  mutate(
    log_income_t1_nz = log(ingreso_T1_nz + 1),
    log_income_t2_nz = log(ingreso_T2_nz + 1),
    log_income_t3_nz = log(ingreso_T3_nz + 1)
  )

panel_peak <- filter(panel_full, month %in% 5:9)

# Run models
models <- list(
  feols(log_income_t1    ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe),
  feols(log_income_t2    ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe),
  feols(log_income_t3    ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe),
  feols(log_income_t1_nz ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe),
  feols(log_income_t2_nz ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe),
  feols(log_income_t3_nz ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)
)

coef_df <- tibble(
  label   = c("Bottom (original)", "Middle (original)", "Top (original)",
               "Bottom (non-zero)", "Middle (non-zero)", "Top (non-zero)"),
  sample  = c(rep("Includes zero-income households", 3),
              rep("Excludes zero-income households", 3)),
  tertile = rep(c("Bottom", "Middle", "Top"), 2),
  est     = sapply(models, function(m) coef(m)["afai_sargassum"]),
  se      = sapply(models, function(m) se(m)["afai_sargassum"])
) %>%
  mutate(
    ci_lo  = est - 1.96 * se,
    ci_hi  = est + 1.96 * se,
    tertile = factor(tertile, levels = c("Bottom", "Middle", "Top")),
    sample  = factor(sample,
                     levels = c("Includes zero-income households",
                                "Excludes zero-income households"))
  )

p7 <- ggplot(coef_df, aes(x = est, y = tertile, colour = sample, shape = sample)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(
    values = c("Includes zero-income households" = "#d73027",
               "Excludes zero-income households" = "#1a9850"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("Includes zero-income households" = 16,
               "Excludes zero-income households" = 17),
    name = NULL
  ) +
  labs(
    x        = "Coefficient on Sargassum intensity (peak season)",
    y        = NULL,
    title    = "CHECK 7: Original vs non-zero tertile coefficients (peak season)",
    subtitle = "The bottom tertile coefficient should become more negative and more precise\n(tighter confidence interval) when zero-income households are excluded."
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 12, face = "bold"),
    plot.subtitle   = element_text(size = 9, colour = "grey30"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(path_diag, "check7_tertile_coefficient_comparison.png"),
       p7, width = 8, height = 4.5, dpi = 150)
cat("  Saved: check7_tertile_coefficient_comparison.png\n")

cat("\n====== DIAGNOSTICS COMPLETE ======\n")
cat("All 7 checks saved to figures/diagnostics/\n\n")
cat("What to look for in each figure:\n")
cat("  check1: Tourist zones & Santo Domingo should be green. Rural interior should be red.\n")
cat("  check2: Inland municipalities should have LOWER Sargassum than coastal ones.\n")
cat("  check3: Blue line (coastal) should generally be ABOVE red (non-coastal).\n")
cat("  check4: Downward slope = more Sargassum -> lower income. Should be negative.\n")
cat("  check5: Mostly green. Large red areas = data gaps worth investigating.\n")
cat("  check6: All lines upward over time. Top > Middle > Bottom at every year.\n")
cat("  check7: Green bottom point should be more negative/precise than red bottom point.\n")
