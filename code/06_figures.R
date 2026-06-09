# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 6: Descriptive figures
#
# Produces three figures:
#   1. map_coastal_type.png      — coastal vs not-coastal classification
#   2. map_sargassum_intensity.png — mean Sargassum exposure by municipality
#   3. ts_sargassum.png          — monthly Sargassum intensity 2017–2025
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(sf)
library(geodata)

path_processed <- "data/processed/"
path_raw_sat   <- "data/raw/satellite/"
path_figures   <- "figures/"
dir.create(path_figures, showWarnings = FALSE)

clean_name <- function(x) {
  x <- tolower(x)
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  x <- gsub("[^a-z0-9 ]", "", x)
  x <- trimws(x)
  gsub("[ ]+", " ", x)
}

# ── Load data ─────────────────────────────────────────────────
cat("Loading satellite coastal data...\n")
sat <- readRDS(file.path(path_processed, "satellite_coastal.rds"))

cat("Loading GADM boundaries...\n")
dr_gadm <- gadm(country = "DOM", level = 2, path = path_raw_sat)
dr_sf   <- st_as_sf(dr_gadm) %>%
  select(municipio = NAME_2, provincia = NAME_1) %>%
  st_transform(4326)

# ── Coastal classification (matches script 04 logic exactly) ──
dr_sf_proj <- st_transform(dr_sf, 32619)
coastline  <- st_union(dr_sf_proj) %>% st_boundary()
coast_zone <- st_buffer(coastline, dist = 500)
touches    <- lengths(st_intersects(dr_sf_proj, coast_zone)) > 0

dr_sf <- dr_sf %>%
  mutate(
    muni_key     = clean_name(municipio),
    coastal_type = if_else(touches, "Coastal", "Not coastal")
  )

# ── Municipality-level mean Sargassum ─────────────────────────
muni_mean <- sat %>%
  mutate(muni_key = clean_name(municipio)) %>%
  group_by(muni_key) %>%
  summarise(mean_sargassum = mean(afai_sargassum, na.rm = TRUE), .groups = "drop")

dr_sf <- dr_sf %>%
  left_join(muni_mean, by = "muni_key") %>%
  mutate(mean_sargassum = if_else(coastal_type == "Coastal", mean_sargassum, NA_real_))

cat("Municipalities with Sargassum data:", sum(!is.na(dr_sf$mean_sargassum)), "\n\n")

# ── Figure 1: Coastal classification map ─────────────────────
cat("Building Figure 1: coastal classification map...\n")

p1 <- ggplot(dr_sf) +
  geom_sf(aes(fill = coastal_type), colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c("Coastal" = "#2166ac", "Not coastal" = "#d1e5f0"),
    name   = NULL
  ) +
  labs(
    title    = "Municipal coastal classification",
    subtitle = "Coastal = direct shoreline contact; Not coastal = all others"
  ) +
  theme_void() +
  theme(
    legend.position      = "bottom",
    legend.text          = element_text(size = 10),
    plot.title           = element_text(size = 13, face = "bold", margin = margin(b = 4)),
    plot.subtitle        = element_text(size = 9, colour = "grey40"),
    plot.margin          = margin(10, 10, 10, 10)
  )

ggsave(file.path(path_figures, "map_coastal_type.png"),
       p1, width = 7, height = 5, dpi = 150)
cat("  Saved: figures/map_coastal_type.png\n")

# ── Figure 2: Mean Sargassum intensity map ────────────────────
cat("Building Figure 2: Sargassum intensity map...\n")

p2 <- ggplot(dr_sf) +
  geom_sf(aes(fill = mean_sargassum), colour = "white", linewidth = 0.2) +
  scale_fill_distiller(
    palette  = "YlOrRd",
    direction = 1,
    na.value = "grey92",
    name     = "Mean excess\nAFAI",
    labels   = scales::label_scientific()
  ) +
  labs(
    title    = "Mean Sargassum intensity by municipality, 2017–2025",
    subtitle = "Average monthly excess AFAI within 20 km of coastline. Grey = no coastal pixels."
  ) +
  theme_void() +
  theme(
    legend.position      = "right",
    legend.title         = element_text(size = 9),
    legend.text          = element_text(size = 8),
    plot.title           = element_text(size = 13, face = "bold", margin = margin(b = 4)),
    plot.subtitle        = element_text(size = 9, colour = "grey40"),
    plot.margin          = margin(10, 10, 10, 10)
  )

ggsave(file.path(path_figures, "map_sargassum_intensity.png"),
       p2, width = 7, height = 5, dpi = 150)
cat("  Saved: figures/map_sargassum_intensity.png\n")

# ── Figure 3: Monthly Sargassum time series 2017–2025 ─────────
cat("Building Figure 3: time series...\n")

ts_data <- sat %>%
  filter(!is.na(afai_sargassum)) %>%
  group_by(year, month) %>%
  summarise(mean_sargassum = mean(afai_sargassum, na.rm = TRUE), .groups = "drop") %>%
  mutate(date = as.Date(paste(year, sprintf("%02d", month), "01", sep = "-")))

# Shaded bands for documented high-Sargassum years
spike_years <- c(2018, 2021, 2022, 2023, 2025)
shading <- tibble(
  xmin = as.Date(paste0(spike_years, "-01-01")),
  xmax = as.Date(paste0(spike_years, "-12-31")),
  ymin = -Inf,
  ymax =  Inf
)

p3 <- ggplot(ts_data, aes(x = date, y = mean_sargassum)) +
  geom_rect(
    data        = shading,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill        = "#fee08b",
    alpha       = 0.45
  ) +
  geom_line(colour = "#2166ac", linewidth = 0.7) +
  geom_point(size = 1.4, colour = "#2166ac") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand      = expansion(add = 15)
  ) +
  scale_y_continuous(labels = scales::label_scientific()) +
  labs(
    x        = NULL,
    y        = "Mean Sargassum intensity (excess AFAI)",
    title    = "Monthly Sargassum intensity offshore Dominican Republic, 2017–2025",
    subtitle = "Mean across all coastal municipalities. Yellow shading = documented high-Sargassum years (2018, 2021–2023, 2025)."
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 9)
  )

ggsave(file.path(path_figures, "ts_sargassum.png"),
       p3, width = 9, height = 4, dpi = 150)
cat("  Saved: figures/ts_sargassum.png\n")

cat("\n====== FIGURES COMPLETE ======\n")
cat("Output in figures/:\n")
cat("  map_coastal_type.png\n")
cat("  map_sargassum_intensity.png\n")
cat("  ts_sargassum.png\n")
