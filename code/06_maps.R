# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 6: Diagnostic and presentation maps
#
# Maps produced:
#   1. Coastal vs non-coastal municipality classification
#   2. Average Sargassum intensity per municipality (2017–2025)
#   3. Nearest-coastal crosswalk (which coastal neighbour each
#      inland municipality is assigned to)
#   4. Sargassum intensity by year (small multiples, spike vs quiet)
#   5. Peak vs off-peak season Sargassum intensity
#   6. Municipality inclusion in final panel vs excluded
#   7. Formal vs informal employment share by municipality
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(sf)
library(geodata)

.libPaths(c("/Users/t.pellissier-lloyd/R/library", .libPaths()))

path_processed <- "data/processed/"
path_figures   <- "figures/"
path_raw_sat   <- "data/raw/satellite/"
dir.create(path_figures, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────
cat("Loading data...\n")
sat       <- readRDS(file.path(path_processed, "satellite_coastal.rds"))
panel     <- readRDS(file.path(path_processed, "panel_monthly.rds"))
encft_m   <- readRDS(file.path(path_processed, "encft_hogar_monthly.rds"))

# ── Load DR municipality shapefile ────────────────────────────
dr_gadm <- gadm(country = "DOM", level = 2, path = path_raw_sat)
dr_sf   <- st_as_sf(dr_gadm) %>%
  select(municipio = NAME_2, provincia = NAME_1) %>%
  st_transform(4326)

clean_name <- function(x) {
  x %>% tolower() %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    gsub("[^a-z0-9 ]", "", .) %>%
    trimws() %>%
    gsub("\\s+", " ", .)
}

dr_sf <- dr_sf %>%
  mutate(muni_key = clean_name(municipio))

# ── Coastal classification ────────────────────────────────────
dr_sf_proj  <- st_transform(dr_sf, 32619)
haiti_gadm  <- gadm(country = "HTI", level = 0, path = path_raw_sat)
haiti_sf    <- st_as_sf(haiti_gadm) %>% st_transform(32619)
haiti_buf   <- st_buffer(st_union(haiti_sf), dist = 2000)
dr_boundary <- st_union(dr_sf_proj) %>% st_boundary()
coastline   <- st_difference(dr_boundary, haiti_buf)
coast_zone  <- st_buffer(coastline, dist = 500)
touches_coast <- lengths(st_intersects(dr_sf_proj, coast_zone)) > 0

coastal_lookup <- dr_sf %>%
  st_drop_geometry() %>%
  mutate(coastal_type = if_else(touches_coast, "Coastal", "Not coastal"))

dr_map <- dr_sf %>%
  left_join(coastal_lookup, by = "muni_key")

# ── IDW weights for map ───────────────────────────────────────
coastal_keys     <- coastal_lookup %>% filter(coastal_type == "Coastal")     %>% pull(muni_key)
non_coastal_keys <- coastal_lookup %>% filter(coastal_type == "Not coastal") %>% pull(muni_key)

coastal_geom <- dr_sf_proj %>%
  filter(muni_key %in% coastal_keys)

non_coastal_geom <- dr_sf_proj %>%
  filter(muni_key %in% non_coastal_keys)

dist_matrix <- st_distance(st_centroid(non_coastal_geom), st_centroid(coastal_geom))
dist_mat    <- matrix(as.numeric(dist_matrix),
                      nrow = nrow(non_coastal_geom), ncol = nrow(coastal_geom))
weight_mat  <- 1 / dist_mat^2
weight_mat  <- weight_mat / rowSums(weight_mat)

# For map: find the top-weighted coastal municipality per non-coastal one
top_weight <- apply(weight_mat, 1, which.max)
idw_dominant <- tibble(
  muni_key       = non_coastal_geom$muni_key,
  source_coastal = coastal_geom$muni_key[top_weight],
  top_weight     = apply(weight_mat, 1, max)
)

# ── Panel inclusion ───────────────────────────────────────────
in_panel <- panel %>%
  distinct(muni_key = clean_name(DES_MUNICIPIO)) %>%
  mutate(in_panel = TRUE)

# ── Average Sargassum per municipality ────────────────────────
sat_avg <- sat %>%
  mutate(muni_key = clean_name(municipio)) %>%
  group_by(muni_key) %>%
  summarise(
    afai_mean_overall = mean(afai_sargassum, na.rm = TRUE),
    afai_peak_mean    = mean(afai_sargassum[month %in% 5:9], na.rm = TRUE),
    afai_offpeak_mean = mean(afai_sargassum[!month %in% 5:9], na.rm = TRUE),
    .groups = "drop"
  )

# Annual Sargassum for small multiples
sat_annual <- sat %>%
  mutate(muni_key = clean_name(municipio)) %>%
  group_by(muni_key, year) %>%
  summarise(afai_annual = mean(afai_sargassum, na.rm = TRUE), .groups = "drop") %>%
  mutate(year_label = paste0(year,
    case_when(
      year %in% c(2018, 2021, 2022, 2023, 2025) ~ "\n(high)",
      year %in% c(2017, 2024)                    ~ "\n(low)",
      TRUE                                        ~ ""
    )
  ))

# Informal share per municipality
informal_share <- encft_m %>%
  filter(n_empleados > 0) %>%
  group_by(ID_MUNICIPIO, DES_MUNICIPIO) %>%
  summarise(
    mean_share_informal = weighted.mean(share_informal,
                                        w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(muni_key = clean_name(DES_MUNICIPIO))

# ── Join all layers to shapefile ──────────────────────────────
dr_map <- dr_map %>%
  left_join(sat_avg,       by = "muni_key") %>%
  left_join(in_panel,      by = "muni_key") %>%
  left_join(informal_share, by = "muni_key") %>%
  mutate(in_panel = replace_na(in_panel, FALSE))

# Assign nearest-coastal source for map lines
# ── Shared theme ──────────────────────────────────────────────
map_theme <- theme_void() +
  theme(
    plot.title    = element_text(size = 13, face = "bold", margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9,  colour = "grey40", margin = margin(b = 8)),
    legend.title  = element_text(size = 9),
    legend.text   = element_text(size = 8),
    plot.margin   = margin(10, 10, 10, 10)
  )

# ── Map 1: Coastal vs non-coastal classification ──────────────
cat("Producing Map 1: Coastal classification...\n")

p1 <- ggplot(dr_map) +
  geom_sf(aes(fill = coastal_type), colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c("Coastal" = "#1a6faf", "Not coastal" = "#d4a843"),
    name   = NULL
  ) +
  labs(
    title    = "Map 1: Coastal vs Non-Coastal Municipalities",
    subtitle = "Coastal = touches ocean; non-coastal = assigned nearest coastal neighbour's satellite values"
  ) +
  map_theme +
  theme(legend.position = "bottom")

ggsave(file.path(path_figures, "map1_coastal_classification.png"),
       p1, width = 8, height = 5, dpi = 150)
cat("  Saved map1_coastal_classification.png\n")

# ── Map 2: Average Sargassum intensity (coastal munis only) ──
cat("Producing Map 2: Average Sargassum intensity...\n")

dr_map_coast <- dr_map %>%
  mutate(afai_plot = if_else(coastal_type == "Coastal", afai_mean_overall, NA_real_))

p2 <- ggplot(dr_map_coast) +
  geom_sf(aes(fill = afai_plot), colour = "white", linewidth = 0.2) +
  scale_fill_gradient(
    low     = "#fff7bc",
    high    = "#8B0000",
    na.value = "grey85",
    name   = "Mean AFAI\n(2017–2025)",
    labels = scales::scientific
  ) +
  labs(
    title    = "Map 2: Mean Sargassum Intensity by Municipality (2017–2025)",
    subtitle = "Coastal municipalities only — direct satellite measurement. Grey = no satellite data (non-coastal)."
  ) +
  map_theme +
  theme(legend.position = "right")

ggsave(file.path(path_figures, "map2_sargassum_intensity.png"),
       p2, width = 8, height = 5, dpi = 150)
cat("  Saved map2_sargassum_intensity.png\n")

# ── Map 3: IDW colour-blended influence map ───────────────────
cat("Producing Map 3: IDW colour-blended influence...\n")

# Each coastal municipality gets a unique HCL colour, evenly spaced around
# the full hue wheel so they are as distinct from each other as possible.
# Each inland municipality gets a colour that is the IDW-weighted average
# of all coastal municipalities' RGB values — colour-mixing by distance.
# Two inland municipalities with similar colours draw from the same coastal sources.

# Order coastal municipalities by their angle around the island centre,
# so colours transition smoothly around the coastline.
n_coastal   <- nrow(coastal_geom)
island_centre <- st_centroid(st_union(coastal_geom))
coast_coords  <- st_coordinates(st_centroid(coastal_geom))
centre_coords <- st_coordinates(island_centre)
angles <- atan2(coast_coords[, 2] - centre_coords[2],
                coast_coords[, 1] - centre_coords[1])
order_idx  <- order(angles)

# Assign hues in geographic order — full rainbow around the coast
hues       <- seq(0, 360, length.out = n_coastal + 1)[seq_len(n_coastal)]
coast_cols <- rep(NA_character_, n_coastal)
coast_cols[order_idx] <- hcl(h = hues, c = 90, l = 55)

coast_rgb  <- col2rgb(coast_cols) / 255  # 3 × n_coastal

# IDW-weighted RGB blend for each non-coastal municipality
# weight_mat is n_non_coastal × n_coastal (rows sum to 1)
blended_rgb <- weight_mat %*% t(coast_rgb)  # n_non_coastal × 3

color_lookup <- bind_rows(
  tibble(muni_key  = coastal_geom$muni_key,
         fill_hex  = coast_cols),
  tibble(muni_key  = non_coastal_geom$muni_key,
         fill_hex  = rgb(blended_rgb[, 1], blended_rgb[, 2], blended_rgb[, 3]))
)

dr_map3 <- dr_sf %>%
  left_join(coastal_lookup, by = "muni_key") %>%
  left_join(color_lookup,   by = "muni_key") %>%
  mutate(fill_hex = replace_na(fill_hex, "#cccccc"))

p3 <- ggplot(dr_map3) +
  geom_sf(aes(fill = fill_hex), colour = "white", linewidth = 0.2) +
  geom_sf(data = filter(dr_map3, coastal_type == "Coastal"),
          fill = NA, colour = "black", linewidth = 0.7) +
  scale_fill_identity() +
  labs(
    title    = "Map 3: IDW Colour-Blended Coastal Influence",
    subtitle = paste0(
      "Each coastal municipality (bold border) has a unique colour.\n",
      "Inland municipalities show a blended colour weighted by distance to each coastal source — ",
      "similar colours mean similar satellite data."
    )
  ) +
  map_theme

ggsave(file.path(path_figures, "map3_idw_influence.png"),
       p3, width = 9, height = 5.5, dpi = 150)
cat("  Saved map3_idw_influence.png\n")

# ── Map 4: Panel inclusion ────────────────────────────────────
cat("Producing Map 4: Panel inclusion...\n")

p4 <- ggplot(dr_map) +
  geom_sf(aes(fill = in_panel), colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c("TRUE" = "#2ca25f", "FALSE" = "#d9d9d9"),
    labels = c("TRUE" = "In final panel", "FALSE" = "Excluded (no survey match)"),
    name   = NULL
  ) +
  labs(
    title    = "Map 4: Municipalities Included in Final Analysis Panel",
    subtitle = "Green = matched in both satellite and ENCFT survey data. Grey = no match found."
  ) +
  map_theme +
  theme(legend.position = "bottom")

ggsave(file.path(path_figures, "map4_panel_inclusion.png"),
       p4, width = 8, height = 5, dpi = 150)
cat("  Saved map4_panel_inclusion.png\n")

# ── Map 5: Peak vs off-peak Sargassum (side by side) ─────────
cat("Producing Map 5: Peak vs off-peak season Sargassum...\n")

dr_season <- dr_map %>%
  select(muni_key, geometry, coastal_type, afai_peak_mean, afai_offpeak_mean) %>%
  pivot_longer(
    cols      = c(afai_peak_mean, afai_offpeak_mean),
    names_to  = "season",
    values_to = "afai"
  ) %>%
  mutate(
    season = if_else(season == "afai_peak_mean",
                     "Peak season (May–Sep)", "Off-peak (Oct–Apr)"),
    afai   = if_else(coastal_type == "Coastal", afai, NA_real_)
  )

p5 <- ggplot(dr_season) +
  geom_sf(aes(fill = afai), colour = "white", linewidth = 0.15) +
  scale_fill_gradient(
    low = "#fff7bc", high = "#8B0000",
    na.value = "grey85",
    name   = "Mean AFAI",
    labels = scales::scientific
  ) +
  facet_wrap(~season) +
  labs(
    title    = "Map 5: Sargassum Intensity by Season",
    subtitle = "Peak season (May–Sep) should show higher values than off-peak — checks seasonal signal"
  ) +
  map_theme +
  theme(
    legend.position  = "right",
    strip.text       = element_text(size = 10, face = "bold"),
    panel.spacing    = unit(1, "lines")
  )

ggsave(file.path(path_figures, "map5_peak_offpeak.png"),
       p5, width = 12, height = 5, dpi = 150)
cat("  Saved map5_peak_offpeak.png\n")

# ── Map 6: Year-by-year small multiples ───────────────────────
cat("Producing Map 6: Annual Sargassum intensity (small multiples)...\n")

dr_annual <- dr_sf %>%
  left_join(sat_annual, by = "muni_key") %>%
  left_join(coastal_lookup, by = "muni_key") %>%
  filter(!is.na(year)) %>%
  mutate(afai_plot = if_else(coastal_type == "Coastal", afai_annual, NA_real_))

p6 <- ggplot(dr_annual) +
  geom_sf(aes(fill = afai_plot), colour = "white", linewidth = 0.1) +
  scale_fill_gradient(
    low = "#fff7bc", high = "#8B0000",
    na.value = "grey85",
    name   = "Mean AFAI",
    labels = scales::scientific
  ) +
  facet_wrap(~year_label, nrow = 3) +
  labs(
    title    = "Map 6: Annual Sargassum Intensity by Municipality",
    subtitle = "High = documented spike years (2018, 2021, 2022, 2023, 2025). Low = documented quiet years (2017, 2024)."
  ) +
  map_theme +
  theme(
    legend.position = "right",
    strip.text      = element_text(size = 7.5),
    panel.spacing   = unit(0.3, "lines")
  )

ggsave(file.path(path_figures, "map6_annual_sargassum.png"),
       p6, width = 12, height = 9, dpi = 150)
cat("  Saved map6_annual_sargassum.png\n")

# ── Map 7: Informal employment share ─────────────────────────
cat("Producing Map 7: Informal employment share...\n")

p7 <- ggplot(dr_map) +
  geom_sf(aes(fill = mean_share_informal), colour = "white", linewidth = 0.2) +
  scale_fill_gradient(
    low      = "#f7fbff",
    high     = "#08306b",
    na.value = "grey85",
    name     = "Share informal\n(0–1)",
    limits   = c(0, 1)
  ) +
  labs(
    title    = "Map 7: Informal Employment Share by Municipality",
    subtitle = "Average share of employed household members in informal or domestic service sector (GRUPO_SECTOR, ENCFT)"
  ) +
  map_theme +
  theme(legend.position = "right")

ggsave(file.path(path_figures, "map7_informal_share.png"),
       p7, width = 8, height = 5, dpi = 150)
cat("  Saved map7_informal_share.png\n")

# ── Copy all maps to Desktop ──────────────────────────────────
cat("\nCopying maps to Desktop...\n")
map_files <- list.files(path_figures, pattern = "^map[0-9].*\\.png$", full.names = TRUE)
for (f in map_files) {
  dest <- file.path("~/Desktop", basename(f))
  file.copy(f, dest, overwrite = TRUE)
  cat("  Copied:", basename(f), "\n")
}

cat("\n====== ALL MAPS DONE ======\n")
cat("Maps saved to figures/ and Desktop\n")
