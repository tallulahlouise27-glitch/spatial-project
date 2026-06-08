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

# ── Models 6–8: By income tertile ────────────────────────────
# Tests whether Sargassum hits poorer households harder.
# Expected: negative coefficient largest (most negative) for T1,
# shrinking toward zero for T3.

iv_t1 <- feols(
  log_income_t1 ~ 1 | muni_fe + year_fe | afai_cov_annual ~ z_bartik,
  data    = panel,
  cluster = ~muni_fe
)

iv_t2 <- feols(
  log_income_t2 ~ 1 | muni_fe + year_fe | afai_cov_annual ~ z_bartik,
  data    = panel,
  cluster = ~muni_fe
)

iv_t3 <- feols(
  log_income_t3 ~ 1 | muni_fe + year_fe | afai_cov_annual ~ z_bartik,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 9: Stacked IV with tertile interactions ─────────────
# Single regression that formally tests whether the Sargassum effect
# differs across income tertiles.
#
# Method: reshape panel to long format (3 rows per municipality-year,
# one per tertile), then instrument AFAI and its interactions with
# T2/T3 dummies using the Bartik instrument and its interactions.
#
# Coefficients:
#   fit_afai_cov_annual = effect on T1 (base group, poorest third)
#   fit_afai_x_t2       = how much SMALLER the effect is for T2 vs T1
#   fit_afai_x_t3       = how much SMALLER the effect is for T3 vs T1
#
# If Sargassum hits the poor hardest, we expect:
#   fit_afai_cov_annual < 0  (negative effect on T1)
#   fit_afai_x_t2 > 0        (T2 less affected than T1)
#   fit_afai_x_t3 > 0        (T3 even less affected)

panel_long <- panel %>%
  select(ID_MUNICIPIO, ANO, muni_fe, year_fe,
         afai_cov_annual, z_bartik,
         log_income_t1, log_income_t2, log_income_t3) %>%
  pivot_longer(
    cols      = c(log_income_t1, log_income_t2, log_income_t3),
    names_to  = "tertile",
    values_to = "log_income_tertile"
  ) %>%
  mutate(
    tertile = factor(tertile,
                     levels = c("log_income_t1", "log_income_t2", "log_income_t3"),
                     labels = c("T1", "T2", "T3")),
    is_t2 = as.integer(tertile == "T2"),
    is_t3 = as.integer(tertile == "T3"),
    # Interactions: AFAI × tertile dummy (T1 is base group)
    afai_x_t2     = afai_cov_annual * is_t2,
    afai_x_t3     = afai_cov_annual * is_t3,
    # Instruments: Bartik × tertile dummy
    z_bartik_x_t2 = z_bartik * is_t2,
    z_bartik_x_t3 = z_bartik * is_t3,
    # Municipality × tertile and year × tertile fixed effects
    muni_tertile_fe = interaction(muni_fe, tertile),
    year_tertile_fe = interaction(year_fe, tertile)
  ) %>%
  filter(!is.na(log_income_tertile))

# Three endogenous variables, three instruments (exactly identified)
iv_stacked <- feols(
  log_income_tertile ~ 1 |
    muni_tertile_fe + year_tertile_fe |
    afai_cov_annual + afai_x_t2 + afai_x_t3 ~
    z_bartik        + z_bartik_x_t2 + z_bartik_x_t3,
  data    = panel_long,
  cluster = ~ID_MUNICIPIO   # cluster at municipality, not municipality×tertile
)

# ── Models 10–12: By coastal proximity ────────────────────────
# Tests whether the Sargassum effect is stronger for municipalities
# with direct coastline vs near-coastal vs inland.

# Check sample sizes before running — small groups produce unreliable results
MIN_MUNICIPALITIES <- 10
coastal_counts <- panel %>%
  group_by(coastal_type) %>%
  summarise(n_munis = n_distinct(ID_MUNICIPIO), n_obs = n(), .groups = "drop")

cat("\nSample size by coastal type:\n")
print(coastal_counts)

small_groups <- coastal_counts %>% filter(n_munis < MIN_MUNICIPALITIES)
if (nrow(small_groups) > 0) {
  cat("\nWARNING: the following groups have fewer than", MIN_MUNICIPALITIES,
      "municipalities — IV estimates will be unreliable:\n")
  print(small_groups$coastal_type)
}

safe_iv <- function(type, data) {
  sub <- filter(data, coastal_type == type)
  n   <- n_distinct(sub$ID_MUNICIPIO)
  if (n < MIN_MUNICIPALITIES) {
    cat("SKIPPING", type, "(only", n, "municipalities)\n")
    return(NULL)
  }
  feols(log_income ~ 1 | muni_fe + year_fe | afai_cov_annual ~ z_bartik,
        data = sub, cluster = ~muni_fe)
}

iv_coastal      <- safe_iv("coastal",      panel)
iv_near_coastal <- safe_iv("near_coastal", panel)
iv_inland       <- safe_iv("inland",       panel)

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

cat("\n--- Heterogeneity by income tertile (separate regressions) ---\n")
etable(
  iv_t1, iv_t2, iv_t3,
  headers  = c("T1 (Bottom)", "T2 (Middle)", "T3 (Top)"),
  se.below = TRUE
)

cat("\n--- Stacked regression with tertile interactions ---\n")
cat("Base group = T1. Interaction coefficients show difference from T1.\n")
print(summary(iv_stacked))
cat("\nWald test: are T2 and T3 effects significantly different from T1?\n")
tryCatch(
  print(wald(iv_stacked, "afai_x_t2|afai_x_t3")),
  error = function(e) cat("Wald test not available:", conditionMessage(e), "\n")
)

cat("\n--- Heterogeneity by coastal proximity ---\n")
het_models   <- Filter(Negate(is.null), list(iv_coastal, iv_near_coastal, iv_inland))
het_headers  <- c("Coastal", "Near-Coastal", "Inland")[!sapply(list(iv_coastal, iv_near_coastal, iv_inland), is.null)]
if (length(het_models) > 0) {
  do.call(etable, c(het_models, list(headers = het_headers, se.below = TRUE)))
} else {
  cat("No coastal-type subgroups had sufficient municipalities to estimate.\n")
}

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

sink(file.path(path_results, "regression_tertile_latex.tex"))
etable(
  iv_t1, iv_t2, iv_t3,
  headers  = c("T1 (Bottom Third)", "T2 (Middle Third)", "T3 (Top Third)"),
  se.below = TRUE,
  tex      = TRUE,
  title    = "Effect of Sargassum Exposure by Income Tertile, DR 2017--2025"
)
sink()

sink(file.path(path_results, "regression_stacked_latex.tex"))
etable(
  iv_stacked,
  se.below = TRUE,
  tex      = TRUE,
  title    = "Stacked IV: Differential Sargassum Effects by Income Tertile, DR 2017--2025"
)
sink()

if (length(het_models) > 0) {
  sink(file.path(path_results, "regression_heterogeneity_latex.tex"))
  do.call(etable, c(het_models, list(
    headers = het_headers,
    se.below = TRUE,
    tex   = TRUE,
    title = "Heterogeneity by Coastal Proximity, DR 2017--2025"
  )))
  sink()
}

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

# ── Tertile coefficient plot ──────────────────────────────────
# Shows the three income-tertile IV estimates side by side.
# If Sargassum hits the poor hardest, the T1 bar should be the
# most negative and T3 closest to zero.
tertile_coef <- data.frame(
  group    = c("T1 (Bottom third)", "T2 (Middle third)", "T3 (Top third)"),
  estimate = c(coef(iv_t1)["fit_afai_cov_annual"],
               coef(iv_t2)["fit_afai_cov_annual"],
               coef(iv_t3)["fit_afai_cov_annual"]),
  se       = c(se(iv_t1)["fit_afai_cov_annual"],
               se(iv_t2)["fit_afai_cov_annual"],
               se(iv_t3)["fit_afai_cov_annual"])
) %>%
  mutate(
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se,
    group = factor(group, levels = rev(group))
  )

p_tertile <- ggplot(tertile_coef, aes(x = estimate, y = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2) +
  geom_point(size = 3) +
  labs(
    x = "Coefficient on Sargassum coastal coverage",
    y = NULL,
    title = "Sargassum Effect by Income Tertile",
    subtitle = "IV estimates. Municipality × tertile and year × tertile fixed effects. 95% CIs."
  ) +
  theme_minimal()

ggsave(file.path(path_figures, "coef_plot_tertile.png"),
       p_tertile, width = 7, height = 3, dpi = 150)
cat("Tertile coefficient plot saved to figures/coef_plot_tertile.png\n")

cat("\n====== DONE ======\n")
cat("Results saved to results/regression_main.txt\n")
cat("LaTeX table saved to results/regression_main_latex.tex\n")
