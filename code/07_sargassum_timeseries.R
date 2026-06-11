# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 7: Average Sargassum exposure time series (2017–2025)
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)

path_processed <- "data/processed/"
path_figures   <- "figures/"
dir.create(path_figures, showWarnings = FALSE)

sat <- readRDS(file.path(path_processed, "satellite_coastal.rds"))

# Monthly average across all coastal municipalities
monthly_avg <- sat %>%
  group_by(year, month) %>%
  summarise(
    afai_sargassum = mean(afai_sargassum, na.rm = TRUE),
    afai_coverage  = mean(afai_coverage,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(date = as.Date(sprintf("%d-%02d-01", year, month)))

# 3-month rolling average for smoothed trend line
monthly_avg <- monthly_avg %>%
  arrange(date) %>%
  mutate(rolling3 = stats::filter(afai_sargassum, rep(1/3, 3), sides = 2) %>% as.numeric())

# Shaded peak season bands (May–Sep each year)
peak_bands <- monthly_avg %>%
  filter(month >= 5, month <= 9) %>%
  group_by(year) %>%
  summarise(xmin = min(date), xmax = max(date), .groups = "drop")

p <- ggplot(monthly_avg, aes(x = date)) +
  # Peak season shading
  geom_rect(
    data = peak_bands,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "#f4c542", alpha = 0.15
  ) +
  # Monthly values
  geom_line(aes(y = afai_sargassum), colour = "#8B9BAB", linewidth = 0.5) +
  geom_point(aes(y = afai_sargassum), colour = "#8B9BAB", size = 1.2) +
  # 3-month rolling average
  geom_line(aes(y = rolling3), colour = "#1B6CA8", linewidth = 1.1, na.rm = TRUE) +
  # Year boundaries
  geom_vline(
    xintercept = as.Date(paste0(2018:2025, "-01-01")),
    colour = "grey70", linewidth = 0.3, linetype = "dashed"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = 0.01)
  ) +
  scale_y_continuous(
    labels = scales::scientific,
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  labs(
    title    = "Average Sargassum Exposure — Dominican Republic (2017–2025)",
    subtitle = "Mean AFAI across all coastal municipalities · Blue line = 3-month rolling average · Yellow bands = peak season (May–Sep)",
    x        = NULL,
    y        = "Mean AFAI (Alternative Floating Algae Index)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x      = element_text(angle = 0, hjust = 0.5)
  )

out_file <- file.path(path_figures, "sargassum_timeseries.png")
ggsave(out_file, p, width = 12, height = 5, dpi = 150, bg = "white")
cat("Saved:", out_file, "\n")

desktop_dir <- "~/Desktop/sargassum_maps"
dir.create(desktop_dir, showWarnings = FALSE)
file.copy(out_file, file.path(desktop_dir, "sargassum_timeseries.png"), overwrite = TRUE)
cat("Copied to Desktop.\n")
