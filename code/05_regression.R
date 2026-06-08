# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 5: Panel IV regression
#
# Model:
#   log(income)_it = β × afai_coverage_it + α_i + γ_t + ε_it
#
#   Where:
#     i = municipality, t = year
#     α_i = municipality fixed effect
#     γ_t = year fixed effect
#     afai_coverage = fraction of coastal pixels with positive AFAI
#                     (direct Sargassum surface coverage measure)
#
# Instrument: Bartik shift-share
#   z_bartik = open_ocean_afai_cov_t × baseline_afai_cov_i
#   "shift"  = open-ocean Sargassum abundance (exogenous Atlantic supply)
#   "share"  = municipality's pre-determined coastal exposure
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(fixest)
library(stargazer)

path_processed <- "data/processed/"
path_results   <- "results/"
path_figures   <- "figures/"
dir.create(path_results, showWarnings = FALSE)
dir.create(path_figures, showWarnings = FALSE)

# ── Load analysis panel ───────────────────────────────────────
panel <- readRDS(file.path(path_processed, "panel_analysis.rds"))
cat("Panel loaded:", nrow(panel), "observations\n\n")

# ── Construct Bartik shift-share instrument ───────────────────
# open_ocean_afai_cov varies only by year → absorbed by year FE alone.
# Interacting with municipality baseline creates municipality×year variation.
# Exclusion: Atlantic Sargassum supply affects DR income only through
# local coastal inundation, scaled by each municipality's exposure.

panel <- panel %>%
  group_by(ID_MUNICIPIO) %>%
  mutate(
    baseline_afai = mean(afai_cov_annual[ANO == min(ANO)], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    z_bartik = afai_ocean_cov * baseline_afai
  )

cat("Instrument (z_bartik) summary:\n")
print(summary(panel$z_bartik))
cat("\n")

# ── Model 1: OLS ──────────────────────────────────────────────
ols_income <- feols(
  log_income ~ afai_cov_annual | muni_fe + year_fe,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 2: IV — main causal estimate ────────────────────────
iv_income <- feols(
  log_income ~ 1 | muni_fe + year_fe | afai_cov_annual ~ z_bartik,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 3: Effect on employment rate ────────────────────────
iv_employ <- feols(
  tasa_empleo ~ 1 | muni_fe + year_fe | afai_cov_annual ~ z_bartik,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 4: Peak season exposure (May–Sep) ───────────────────
iv_peak <- feols(
  log_income ~ 1 | muni_fe + year_fe | afai_peak ~ z_bartik,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 5: Fishing municipalities only ──────────────────────
fishing_muni <- panel %>%
  filter(share_pesca > median(share_pesca, na.rm = TRUE))

iv_fishing <- feols(
  log_income ~ 1 | muni_fe + year_fe | afai_cov_annual ~ z_bartik,
  data    = fishing_muni,
  cluster = ~muni_fe
)

# ── Print results ─────────────────────────────────────────────
cat("====== REGRESSION RESULTS ======\n\n")

cat("--- First stage (instrument strength) ---\n")
first_stage <- feols(
  afai_cov_annual ~ z_bartik | muni_fe + year_fe,
  data    = panel,
  cluster = ~muni_fe
)
print(summary(first_stage))

cat("\n--- Main results ---\n")
etable(
  ols_income, iv_income, iv_employ, iv_peak, iv_fishing,
  headers = c("OLS Income", "IV Income", "IV Employ", "IV Peak", "IV Fishing"),
  se.below = TRUE
)

# ── Instrument validity checks ────────────────────────────────
cat("\n====== INSTRUMENT CHECKS ======\n")
fs <- fitstat(iv_income, "ivf")
cat("First-stage F-statistic:", round(fs$ivf$stat, 2), "\n")
cat("Rule of thumb: F > 10 = strong instrument\n\n")

# ── Export tables ─────────────────────────────────────────────
sink(file.path(path_results, "regression_main.txt"))
etable(
  ols_income, iv_income, iv_employ,
  headers  = c("OLS Income", "IV Income", "IV Employment"),
  se.below = TRUE,
  title    = "Effect of Sargassum Exposure on Household Welfare, DR 2017-2025"
)
sink()

sink(file.path(path_results, "regression_main_latex.tex"))
etable(
  ols_income, iv_income, iv_employ,
  headers  = c("OLS Income", "IV Income", "IV Employment"),
  se.below = TRUE,
  tex      = TRUE,
  title    = "Effect of Sargassum Exposure on Household Welfare, DR 2017--2025"
)
sink()

# ── Coefficient plot ──────────────────────────────────────────
coef_data <- data.frame(
  model    = c("OLS", "IV (Main)", "IV (Peak\nSeason)", "IV (Fishing\nMunicipalities)"),
  estimate = c(coef(ols_income)["afai_cov_annual"],
               coef(iv_income)["fit_afai_cov_annual"],
               coef(iv_peak)["fit_afai_peak"],
               coef(iv_fishing)["fit_afai_cov_annual"]),
  se       = c(se(ols_income)["afai_cov_annual"],
               se(iv_income)["fit_afai_cov_annual"],
               se(iv_peak)["fit_afai_peak"],
               se(iv_fishing)["fit_afai_cov_annual"])
) %>%
  mutate(
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se,
    model = factor(model, levels = rev(model))
  )

p <- ggplot(coef_data, aes(x = estimate, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2) +
  geom_point(size = 3) +
  labs(
    x = "Coefficient on Sargassum coastal coverage (fraction of pixels)",
    y = NULL,
    title = "Effect of Sargassum Exposure on Log Household Income",
    subtitle = "Municipality + year fixed effects. 95% confidence intervals."
  ) +
  theme_minimal()

ggsave(file.path(path_figures, "coef_plot.png"), p, width = 7, height = 4, dpi = 150)
cat("Coefficient plot saved to figures/coef_plot.png\n")

cat("\n====== DONE ======\n")
cat("Results saved to results/regression_main.txt\n")
cat("LaTeX table saved to results/regression_main_latex.tex\n")
