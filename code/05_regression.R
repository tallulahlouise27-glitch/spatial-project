# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 5: Panel IV regression
#
# Model:
#   log(income)_it = β × log(coastal_chla)_it + α_i + γ_t + ε_it
#
#   Where:
#     i = municipality, t = year
#     α_i = municipality fixed effect (absorbs time-invariant local factors)
#     γ_t = year fixed effect (absorbs nationwide shocks)
#     Instrument for log(coastal_chla): log(open_ocean_chla) — offshore
#       Sargassum abundance, exogenous to DR local economic conditions
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(fixest)
library(stargazer)

path_processed <- "data/processed/"
path_results   <- "results/"
path_figures   <- "figures/"

# ── Load analysis panel ───────────────────────────────────────
panel <- readRDS(file.path(path_processed, "panel_analysis.rds"))
cat("Panel loaded:", nrow(panel), "observations\n\n")

# ── Model 1: OLS with municipality + year fixed effects ───────
# Baseline: does coastal chlorophyll correlate with lower income?
ols_income <- feols(
  log_income ~ log_chla | muni_fe + year_fe,
  data    = panel,
  cluster = ~muni_fe   # cluster standard errors at municipality level
)

# ── Model 2: IV regression (2SLS) ─────────────────────────────
# Instrument: open-ocean chlorophyll (log_chla_ocean)
# This is the main causal estimate
iv_income <- feols(
  log_income ~ 1 | muni_fe + year_fe | log_chla ~ log_chla_ocean,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 3: Effect on employment rate ────────────────────────
iv_employ <- feols(
  tasa_empleo ~ 1 | muni_fe + year_fe | log_chla ~ log_chla_ocean,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 4: Peak season exposure (May–Sep) ───────────────────
iv_peak <- feols(
  log_income ~ 1 | muni_fe + year_fe | log_chla_peak ~ log_chla_ocean,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 5: Fishing municipalities only ──────────────────────
fishing_muni <- panel %>%
  filter(share_pesca > median(share_pesca, na.rm = TRUE))

iv_fishing <- feols(
  log_income ~ 1 | muni_fe + year_fe | log_chla ~ log_chla_ocean,
  data    = fishing_muni,
  cluster = ~muni_fe
)

# ── Print results ─────────────────────────────────────────────
cat("====== REGRESSION RESULTS ======\n\n")

cat("--- First stage (instrument strength) ---\n")
first_stage <- feols(
  log_chla ~ log_chla_ocean | muni_fe + year_fe,
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

# First-stage F-statistic (should be > 10 for strong instrument)
fs <- fitstat(iv_income, "ivf")
cat("First-stage F-statistic:", round(fs$ivf$stat, 2), "\n")
cat("Rule of thumb: F > 10 = strong instrument\n\n")

# ── Export tables ─────────────────────────────────────────────
# Main results table
sink(file.path(path_results, "regression_main.txt"))
etable(
  ols_income, iv_income, iv_employ,
  headers  = c("OLS Income", "IV Income", "IV Employment"),
  se.below = TRUE,
  title    = "Effect of Sargassum Exposure on Household Welfare, DR 2016-2025"
)
sink()

# LaTeX table for dissertation
sink(file.path(path_results, "regression_main_latex.tex"))
etable(
  ols_income, iv_income, iv_employ,
  headers  = c("OLS Income", "IV Income", "IV Employment"),
  se.below = TRUE,
  tex      = TRUE,
  title    = "Effect of Sargassum Exposure on Household Welfare, DR 2016--2025"
)
sink()

# ── Coefficient plot ──────────────────────────────────────────
coef_data <- data.frame(
  model    = c("OLS", "IV (Main)", "IV (Peak\nSeason)", "IV (Fishing\nMunicipalities)"),
  estimate = c(coef(ols_income)["log_chla"],
               coef(iv_income)["fit_log_chla"],
               coef(iv_peak)["fit_log_chla_peak"],
               coef(iv_fishing)["fit_log_chla"]),
  se       = c(se(ols_income)["log_chla"],
               se(iv_income)["fit_log_chla"],
               se(iv_peak)["fit_log_chla_peak"],
               se(iv_fishing)["fit_log_chla"])
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
    x = "Coefficient on log(coastal chlorophyll)",
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
