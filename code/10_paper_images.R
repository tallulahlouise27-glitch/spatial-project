# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 10: Export publication-quality images
#
# Saves everything to figures/paper_images/:
#   Figures (PNG, 300 dpi):
#     fig1_sargassum_map.png         — mean Sargassum by municipality
#     fig2_main_results.png          — main coefficient plot
#     fig3_tertile_results.png       — tertile coefficient plot
#     fig4_coastal_results.png       — coastal subgroup plot
#     fig5_income_map.png            — income level by municipality
#     fig6_income_trends.png         — income over time by tertile
#   Tables (PNG via flextable):
#     table1_main.png
#     table2_tertile.png
#     table3_coastal.png
#     table4_nz_tertile.png
#     tableA1_fe_progression.png
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(sf)
library(geodata)
library(fixest)
library(flextable)
library(magick)

path_processed <- "data/processed/"
path_raw_sat   <- "data/raw/satellite/"
out_dir        <- "figures/paper_images/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

clean_name <- function(x) {
  x <- tolower(x)
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  x <- gsub("[^a-z0-9 ]", "", x)
  trimws(gsub("[ ]+", " ", x))
}

# ── Theme for all figures ─────────────────────────────────────
theme_paper <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9, colour = "grey40", margin = margin(b = 8)),
      plot.caption     = element_text(size = 8, colour = "grey50", hjust = 0),
      legend.position  = "bottom",
      legend.text      = element_text(size = 9),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(12, 12, 12, 12)
    )
}

save_fig <- function(p, name, w = 7, h = 5) {
  path <- file.path(out_dir, name)
  ggsave(path, p, width = w, height = h, dpi = 300, bg = "white")
  cat("  Saved:", name, "\n")
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

dr_proj    <- st_transform(dr_sf, 32619)
haiti_gadm <- gadm(country = "HTI", level = 0, path = path_raw_sat)
haiti_sf   <- st_as_sf(haiti_gadm) %>% st_transform(32619)
haiti_buf  <- st_buffer(st_union(haiti_sf), dist = 2000)
coastline  <- st_difference(st_union(dr_proj) %>% st_boundary(), haiti_buf)
coast_zone <- st_buffer(coastline, dist = 500)
touches    <- lengths(st_intersects(dr_proj, coast_zone)) > 0
dr_sf      <- dr_sf %>% mutate(coastal = if_else(touches, "Coastal", "Not coastal"))

# ── Re-run models ─────────────────────────────────────────────
cat("Running models...\n")
panel_peak <- filter(panel, month %in% 5:9)

m_full   <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster = ~muni_fe)
m_peak   <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)
m_t1     <- feols(log_income_t1 ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster = ~muni_fe)
m_t2     <- feols(log_income_t2 ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster = ~muni_fe)
m_t3     <- feols(log_income_t3 ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster = ~muni_fe)
m_t1p    <- feols(log_income_t1 ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)
m_t2p    <- feols(log_income_t2 ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)
m_t3p    <- feols(log_income_t3 ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)
m_coast  <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, filter(panel_peak, coastal_type == "coastal"),     cluster = ~muni_fe)
m_inland <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, filter(panel_peak, coastal_type == "not_coastal"), cluster = ~muni_fe)

# FE progression
m_none  <- feols(log_income ~ afai_sargassum, panel_peak, cluster = ~muni_fe)
m_muni  <- feols(log_income ~ afai_sargassum | muni_fe, panel_peak, cluster = ~muni_fe)
m_time  <- feols(log_income ~ afai_sargassum | year_month_fe, panel_peak, cluster = ~muni_fe)
m_both  <- m_peak

extract_coef <- function(m, label, season) {
  b <- coef(m)["afai_sargassum"]
  s <- se(m)["afai_sargassum"]
  p <- pvalue(m)["afai_sargassum"]
  tibble(label = label, season = season, est = b, se = s, pval = p,
         ci_lo = b - 1.96*s, ci_hi = b + 1.96*s,
         sig = cut(p, c(-Inf,0.001,0.01,0.05,0.1,Inf),
                   labels = c("p < 0.001","p < 0.01","p < 0.05","p < 0.10","n.s.")))
}

# ════════════════════════════════════════════════════════════════
# FIGURE 1: Sargassum intensity map
# ════════════════════════════════════════════════════════════════
cat("\nFigure 1: Sargassum map...\n")

muni_sarg <- sat %>%
  mutate(muni_key = clean_name(municipio)) %>%
  group_by(muni_key) %>%
  summarise(
    mean_sarg      = mean(afai_sargassum, na.rm = TRUE),
    peak_sarg      = mean(afai_sargassum[month %in% 5:9], na.rm = TRUE),
    .groups = "drop"
  )

map_sarg <- dr_sf %>%
  left_join(muni_sarg, by = "muni_key") %>%
  filter(!is.na(mean_sarg), muni_key %in% panel$muni_key)

p1 <- ggplot() +
  geom_sf(data = dr_sf, fill = "grey92", colour = "white", linewidth = 0.15) +
  geom_sf(data = map_sarg, aes(fill = peak_sarg), colour = "white", linewidth = 0.15) +
  scale_fill_distiller(
    palette   = "YlOrRd",
    direction = 1,
    name      = "Mean excess\nAFAI (May–Sep)",
    labels    = scales::label_scientific(),
    na.value  = "grey92"
  ) +
  labs(
    title    = "Figure 1: Peak-season Sargassum intensity by municipality",
    subtitle = "Mean excess AFAI (May–September, 2017–2025). Grey = municipalities not in analysis sample.",
    caption  = "Source: NOAA AOML AFAI satellite data. Sample: 86 municipalities."
  ) +
  theme_void() +
  theme(
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 11, margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 8.5, colour = "grey40"),
    plot.caption     = element_text(size = 8, colour = "grey50", hjust = 0),
    plot.margin      = margin(10, 10, 10, 10)
  )

save_fig(p1, "fig1_sargassum_map.png", w = 7.5, h = 5.5)

# ════════════════════════════════════════════════════════════════
# FIGURE 2: Main results coefficient plot
# ════════════════════════════════════════════════════════════════
cat("Figure 2: Main results plot...\n")

coef2 <- bind_rows(
  extract_coef(m_full, "Full year",    "Full year"),
  extract_coef(m_peak, "Peak season\n(May–Sep)", "Peak season")
) %>%
  mutate(label = factor(label, levels = c("Peak season\n(May–Sep)", "Full year")))

p2 <- ggplot(coef2, aes(x = est, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.12, linewidth = 0.7) +
  geom_point(aes(colour = sig), size = 4) +
  scale_colour_manual(
    values = c("p < 0.001"="#2166ac","p < 0.01"="#4393c3","p < 0.05"="#74add1",
               "p < 0.10"="#fee090","n.s."="#d9d9d9"),
    name = "Significance", drop = FALSE
  ) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x       = "Coefficient on Sargassum intensity",
    y       = NULL,
    title   = "Figure 2: Effect of Sargassum on household income",
    subtitle = "OLS with municipality and year×month fixed effects. Error bars = 95% CI.\nDependent variable: log household income. N = 86 municipalities.",
    caption = "Standard errors clustered by municipality."
  ) +
  theme_paper()

save_fig(p2, "fig2_main_results.png", w = 7, h = 4)

# ════════════════════════════════════════════════════════════════
# FIGURE 3: Tertile coefficient plot
# ════════════════════════════════════════════════════════════════
cat("Figure 3: Tertile results plot...\n")

coef3 <- bind_rows(
  extract_coef(m_t1,  "Bottom",  "Full year"),
  extract_coef(m_t2,  "Middle",  "Full year"),
  extract_coef(m_t3,  "Top",     "Full year"),
  extract_coef(m_t1p, "Bottom",  "Peak season (May–Sep)"),
  extract_coef(m_t2p, "Middle",  "Peak season (May–Sep)"),
  extract_coef(m_t3p, "Top",     "Peak season (May–Sep)")
) %>%
  mutate(
    label  = factor(label,  levels = c("Bottom","Middle","Top")),
    season = factor(season, levels = c("Full year","Peak season (May–Sep)"))
  )

p3 <- ggplot(coef3, aes(x = est, y = label, colour = sig, shape = season)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, linewidth = 0.7,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.55)) +
  scale_colour_manual(
    values = c("p < 0.001"="#2166ac","p < 0.01"="#4393c3","p < 0.05"="#74add1",
               "p < 0.10"="#fee090","n.s."="#d9d9d9"),
    name = "Significance", drop = FALSE
  ) +
  scale_shape_manual(values = c("Full year" = 16, "Peak season (May–Sep)" = 17),
                     name = NULL) +
  labs(
    x       = "Coefficient on Sargassum intensity",
    y       = "Income tertile",
    title   = "Figure 3: Effect of Sargassum by income tertile",
    subtitle = "Separate regressions for bottom, middle, and top income thirds.\nCircles = full year. Triangles = peak season (May–Sep).",
    caption = "OLS with municipality and year×month fixed effects. SEs clustered by municipality."
  ) +
  theme_paper()

save_fig(p3, "fig3_tertile_results.png", w = 7.5, h = 4.5)

# ════════════════════════════════════════════════════════════════
# FIGURE 4: Coastal vs non-coastal coefficient plot
# ════════════════════════════════════════════════════════════════
cat("Figure 4: Coastal subgroup plot...\n")

coef4 <- bind_rows(
  extract_coef(m_coast,  "Coastal",     "Peak season"),
  extract_coef(m_inland, "Not coastal", "Peak season")
) %>%
  mutate(label = factor(label, levels = c("Not coastal","Coastal")))

p4 <- ggplot(coef4, aes(x = est, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.12, linewidth = 0.7) +
  geom_point(aes(colour = sig), size = 4) +
  scale_colour_manual(
    values = c("p < 0.001"="#2166ac","p < 0.01"="#4393c3","p < 0.05"="#74add1",
               "p < 0.10"="#fee090","n.s."="#d9d9d9"),
    name = "Significance", drop = FALSE
  ) +
  labs(
    x       = "Coefficient on Sargassum intensity (peak season)",
    y       = NULL,
    title   = "Figure 4: Peak-season effect by coastal proximity",
    subtitle = "Coastal = municipalities whose boundary touches the ocean shoreline (45 of 86).\nNot coastal = all others, within 20km of coast (41 of 86).",
    caption = "OLS with municipality and year×month fixed effects. SEs clustered by municipality."
  ) +
  theme_paper()

save_fig(p4, "fig4_coastal_results.png", w = 7, h = 4)

# ════════════════════════════════════════════════════════════════
# FIGURE 5: Income map
# ════════════════════════════════════════════════════════════════
cat("Figure 5: Income map...\n")

muni_inc <- panel %>%
  group_by(muni_key) %>%
  summarise(mean_income = mean(ingreso_medio, na.rm = TRUE), .groups = "drop") %>%
  mutate(tertile = case_when(
    ntile(mean_income, 3) == 1 ~ "Bottom third",
    ntile(mean_income, 3) == 2 ~ "Middle third",
    ntile(mean_income, 3) == 3 ~ "Top third"
  ))

map_inc <- dr_sf %>%
  left_join(muni_inc, by = "muni_key") %>%
  mutate(tertile = replace_na(tertile, "Not in sample"))

p5 <- ggplot(map_inc) +
  geom_sf(aes(fill = tertile), colour = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = c("Bottom third" = "#d73027", "Middle third" = "#fee090",
               "Top third" = "#1a9850", "Not in sample" = "grey88"),
    name = NULL
  ) +
  labs(
    title    = "Figure 5: Average household income level by municipality",
    subtitle = "Municipalities ranked by mean household income and divided into thirds.\nGrey = not in analysis sample.",
    caption  = "Source: ENCFT household survey, 2017–2025."
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 9),
    plot.title      = element_text(face = "bold", size = 11, margin = margin(b = 4)),
    plot.subtitle   = element_text(size = 8.5, colour = "grey40"),
    plot.caption    = element_text(size = 8, colour = "grey50", hjust = 0),
    plot.margin     = margin(10, 10, 10, 10)
  )

save_fig(p5, "fig5_income_map.png", w = 7.5, h = 5.5)

# ════════════════════════════════════════════════════════════════
# FIGURE 6: Income trends over time
# ════════════════════════════════════════════════════════════════
cat("Figure 6: Income trends...\n")

income_trends <- panel %>%
  group_by(year) %>%
  summarise(
    Bottom = mean(ingreso_T1, na.rm = TRUE),
    Middle = mean(ingreso_T2, na.rm = TRUE),
    Top    = mean(ingreso_T3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-year, names_to = "tertile", values_to = "income") %>%
  mutate(tertile = factor(tertile, levels = c("Top","Middle","Bottom")))

p6 <- ggplot(income_trends, aes(x = year, y = income, colour = tertile)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_colour_manual(
    values = c("Top" = "#1a9850", "Middle" = "#fe9929", "Bottom" = "#d73027"),
    name   = "Income group"
  ) +
  scale_y_continuous(labels = scales::comma_format(prefix = "RD$")) +
  scale_x_continuous(breaks = 2017:2025) +
  labs(
    x       = NULL,
    y       = "Mean household income (RD$/month)",
    title   = "Figure 6: Household income over time by income tertile",
    subtitle = "Annual means for bottom, middle, and top thirds of the national income distribution.",
    caption = "Source: ENCFT household survey. COVID-19 dip visible in 2020."
  ) +
  theme_paper()

save_fig(p6, "fig6_income_trends.png", w = 8.5, h = 5)

# ════════════════════════════════════════════════════════════════
# TABLES — use flextable to render as PNG images
# ════════════════════════════════════════════════════════════════
cat("\nBuilding table images...\n")

ft_theme <- function(ft) {
  ft %>%
    theme_vanilla() %>%
    fontsize(size = 9, part = "all") %>%
    fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>%
    align(align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "all") %>%
    set_table_properties(width = 1, layout = "autofit") %>%
    padding(padding = 4, part = "all") %>%
    border_outer(part = "all", border = fp_border_default(width = 1.5)) %>%
    border_inner_h(part = "header", border = fp_border_default(width = 1)) %>%
    hline_bottom(part = "header", border = fp_border_default(width = 1.5))
}

fmt_coef <- function(b, s, p) {
  stars <- case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*",
                     p < 0.1  ~ "†",   TRUE       ~ "")
  sprintf("%.2f%s\n(%.2f)", b, stars, s)
}

make_reg_ft <- function(models, col_labels, row_label = "Sargassum intensity",
                        note = "") {
  rows <- map(models, function(m) {
    b <- coef(m)["afai_sargassum"]
    s <- se(m)["afai_sargassum"]
    p <- pvalue(m)["afai_sargassum"]
    stars <- case_when(p<0.001~"***",p<0.01~"**",p<0.05~"*",p<0.1~"†",TRUE~"")
    list(
      coef  = sprintf("%.2f%s", b, stars),
      se    = sprintf("(%.2f)", s),
      r2    = sprintf("%.3f", r2(m, type="r2")),
      obs   = formatC(nobs(m), format="d", big.mark=",")
    )
  })

  df <- data.frame(
    ` ` = c(row_label, "", "R²", "Observations"),
    check.names = FALSE
  )
  for (i in seq_along(models)) {
    df[[col_labels[i]]] <- c(rows[[i]]$coef, rows[[i]]$se,
                              rows[[i]]$r2,  rows[[i]]$obs)
  }

  ft <- flextable(df) %>%
    ft_theme() %>%
    italic(i = 2, part = "body") %>%
    color(i = 2, color = "grey40", part = "body") %>%
    hline(i = 2, border = fp_border_default(width = 0.5, color = "grey80")) %>%
    hline(i = 3, border = fp_border_default(width = 1, color = "black"))

  if (nchar(note) > 0) {
    ft <- add_footer_lines(ft, note) %>%
      fontsize(size = 8, part = "footer") %>%
      color(color = "grey40", part = "footer") %>%
      italic(part = "footer")
  }
  ft
}

save_ft <- function(ft, name, w = 6, h = NULL) {
  path <- file.path(out_dir, name)
  save_as_image(ft, path = path, zoom = 3, expand = 10)
  cat("  Saved:", name, "\n")
}

note_main <- "Note: *** p<0.001  ** p<0.01  * p<0.05  † p<0.10. Standard errors (in parentheses) clustered by municipality.\nAll models include municipality and year×month fixed effects."

# Table 1: Main results
ft1 <- make_reg_ft(
  list(m_full, m_peak),
  c("Full Year", "Peak Season\n(May–Sep)"),
  note = note_main
)
ft1 <- set_caption(ft1, "Table 1: Effect of Sargassum on household income — main results")
save_ft(ft1, "table1_main.png")

# Table 2: Tertile results
ft2 <- make_reg_ft(
  list(m_t1, m_t2, m_t3, m_t1p, m_t2p, m_t3p),
  c("Bottom","Middle","Top","Bottom (Peak)","Middle (Peak)","Top (Peak)"),
  note = paste(note_main, "\nIncome tertiles defined nationally within each year.")
)
ft2 <- set_caption(ft2, "Table 2: Effect of Sargassum by income tertile")
save_ft(ft2, "table2_tertile.png", w = 8)

# Table 3: Coastal subgroup (peak season)
ft3 <- make_reg_ft(
  list(m_coast, m_inland),
  c("Coastal", "Not Coastal"),
  note = paste(note_main, "\nPeak season (May–Sep) only. Coastal = boundary touches ocean shoreline.")
)
ft3 <- set_caption(ft3, "Table 3: Peak-season effect by coastal proximity")
save_ft(ft3, "table3_coastal.png")

# Table 4: FE progression
ft4 <- make_reg_ft(
  list(m_none, m_muni, m_time, m_both),
  c("No FEs", "Municipality\nFE only", "Year×Month\nFE only", "Both FEs"),
  note = paste(note_main, "\nPeak season (May–Sep) only.")
)
ft4 <- set_caption(ft4, "Table A1: Fixed effects progression — peak season")
save_ft(ft4, "tableA1_fe_progression.png")

cat("\n====== ALL IMAGES SAVED ======\n")
cat("Location: figures/paper_images/\n")
cat("\nFigures:\n")
figs <- list.files(out_dir, pattern="^fig.*\\.png")
for(f in figs) cat(" ", f, "\n")
cat("\nTables:\n")
tbls <- list.files(out_dir, pattern="^table.*\\.png")
for(f in tbls) cat(" ", f, "\n")
