# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 5: Panel IV regression
#
# Model:
#   log(income)_it = β × afai_sargassum_it + α_i + γ_t + ε_it
#
#   Where:
#     i = municipality, t = year
#     α_i = municipality fixed effect
#     γ_t = year fixed effect
#     afai_sargassum = mean(pmax(AFAI − 0.001, 0)) per municipality-year
#                     captures both spatial coverage AND mat density
#
# Instrument: Bartik shift-share
#   z_bartik = afai_ocean_cov_peak_t × baseline_sargassum_i
#   "shift"  = peak-season (May–Sep) open-ocean Sargassum coverage
#              (exogenous Atlantic supply shock, varies year to year)
#   "share"  = municipality's full-sample mean Sargassum exposure
#              (pre-determined geographic exposure, varies cross-sectionally)
#
# Why peak season? The annual average dilutes the Sargassum signal with
# winter months that have near-zero coverage. May–Sep captures the period
# when the Atlantic Sargassum belt is active and heading toward the Caribbean.
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
# Monthly panel: municipality × year × month.
# Income from quarterly ENCFT matched to the correct quarter for each month.
panel <- readRDS(file.path(path_processed, "panel_monthly.rds"))

# ── Construct Bartik shift-share instrument ───────────────────
# Monthly Bartik: z_bartik = afai_ocean_coverage_mt × baseline_sargassum_i
#
# afai_ocean_coverage varies by month × year (Atlantic Sargassum signal).
# baseline_sargassum is each municipality's full-sample mean exposure.
# Their product creates municipality × month × year variation.
# year×month FEs absorb the common ocean signal; identification comes from
# the differential response of high-exposure vs low-exposure municipalities.
#
# Exclusion restriction: Atlantic open-ocean Sargassum abundance affects
# DR household income only through local coastal inundation.

panel <- panel %>%
  group_by(ID_MUNICIPIO) %>%
  mutate(
    baseline_sargassum = mean(afai_sargassum, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    z_bartik       = afai_ocean_coverage * baseline_sargassum,
    is_not_coastal = as.integer(coastal_type == "not_coastal"),
    # Interaction terms for coastal heterogeneity test
    afai_x_not_coastal    = afai_sargassum * is_not_coastal,
    z_bartik_x_not_coastal = z_bartik      * is_not_coastal
  )

cat("Baseline sargassum by municipality (full-sample mean):\n")
print(summary(panel$baseline_sargassum))
cat("\nInstrument (z_bartik) summary:\n")
print(summary(panel$z_bartik))
cat("\n")

# ── Model 1: OLS ──────────────────────────────────────────────
ols_income <- feols(
  log_income ~ afai_sargassum | muni_fe + year_month_fe,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 2: IV — main causal estimate ────────────────────────
iv_income <- feols(
  log_income ~ 1 | muni_fe + year_month_fe | afai_sargassum ~ z_bartik,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 3: Effect on employment rate ────────────────────────
iv_employ <- feols(
  tasa_empleo ~ 1 | muni_fe + year_month_fe | afai_sargassum ~ z_bartik,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 4: Peak season (Q2 and Q3 only — Apr–Sep) ──────────
# Restricts to the two quarters when Sargassum is most intense.
iv_peak <- feols(
  log_income ~ 1 | muni_fe + year_month_fe | afai_sargassum ~ z_bartik,
  data    = filter(panel, month %in% 5:9),
  cluster = ~muni_fe
)


# ── Models 6–8: By income tertile ────────────────────────────
# Tests whether Sargassum hits poorer households harder.
# Expected: negative coefficient largest (most negative) for T1,
# shrinking toward zero for T3.

iv_t1 <- feols(
  log_income_t1 ~ 1 | muni_fe + year_month_fe | afai_sargassum ~ z_bartik,
  data    = panel,
  cluster = ~muni_fe
)

iv_t2 <- feols(
  log_income_t2 ~ 1 | muni_fe + year_month_fe | afai_sargassum ~ z_bartik,
  data    = panel,
  cluster = ~muni_fe
)

iv_t3 <- feols(
  log_income_t3 ~ 1 | muni_fe + year_month_fe | afai_sargassum ~ z_bartik,
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
#   fit_afai_sargassum = effect on T1 (base group, poorest third)
#   fit_afai_x_t2       = how much SMALLER the effect is for T2 vs T1
#   fit_afai_x_t3       = how much SMALLER the effect is for T3 vs T1
#
# If Sargassum hits the poor hardest, we expect:
#   fit_afai_sargassum < 0  (negative effect on T1)
#   fit_afai_x_t2 > 0        (T2 less affected than T1)
#   fit_afai_x_t3 > 0        (T3 even less affected)

panel_long <- panel %>%
  select(ID_MUNICIPIO, year, month, muni_fe, year_month_fe,
         afai_sargassum, z_bartik,
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
    afai_x_t2     = afai_sargassum * is_t2,
    afai_x_t3     = afai_sargassum * is_t3,
    z_bartik_x_t2 = z_bartik * is_t2,
    z_bartik_x_t3 = z_bartik * is_t3,
    muni_tertile_fe         = interaction(muni_fe, tertile),
    year_month_tertile_fe = interaction(year_month_fe, tertile)
  ) %>%
  filter(!is.na(log_income_tertile))

# Three endogenous variables, three instruments (exactly identified)
iv_stacked <- feols(
  log_income_tertile ~ 1 |
    muni_tertile_fe + year_month_tertile_fe |
    afai_sargassum + afai_x_t2 + afai_x_t3 ~
    z_bartik        + z_bartik_x_t2 + z_bartik_x_t3,
  data    = panel_long,
  cluster = ~ID_MUNICIPIO   # cluster at municipality, not municipality×tertile
)

# ── Models 10–11: By coastal proximity ───────────────────────
# Two categories: coastal (direct shoreline) vs not_coastal (all others).
# Tests whether the Sargassum effect is concentrated in coastal municipalities.

MIN_MUNICIPALITIES <- 10
coastal_counts <- panel %>%
  group_by(coastal_type) %>%
  summarise(n_munis = n_distinct(ID_MUNICIPIO), n_obs = n(), .groups = "drop")

cat("\nSample size by coastal type:\n")
print(coastal_counts)

safe_iv <- function(type, data) {
  sub <- filter(data, coastal_type == type)
  n   <- n_distinct(sub$ID_MUNICIPIO)
  if (n < MIN_MUNICIPALITIES) {
    cat("SKIPPING", type, "(only", n, "municipalities)\n")
    return(NULL)
  }
  feols(log_income ~ 1 | muni_fe + year_month_fe | afai_sargassum ~ z_bartik,
        data = sub, cluster = ~muni_fe)
}

iv_coastal     <- safe_iv("coastal",     panel)
iv_not_coastal <- safe_iv("not_coastal", panel)

# Formal interaction test: does the effect differ between coastal and not-coastal?
# Base group = coastal. Coefficient on afai_x_not_coastal = difference in effect.
# Two endogenous variables, two instruments (exactly identified).
iv_coastal_interact <- feols(
  log_income ~ 1 | muni_fe + year_month_fe |
    afai_sargassum + afai_x_not_coastal ~
    z_bartik       + z_bartik_x_not_coastal,
  data    = panel,
  cluster = ~muni_fe
)

# ── Print results ─────────────────────────────────────────────
cat("====== REGRESSION RESULTS ======\n\n")

cat("--- First stage (instrument strength) ---\n")
first_stage <- feols(
  afai_sargassum ~ z_bartik | muni_fe + year_month_fe,
  data    = panel,
  cluster = ~muni_fe
)
print(summary(first_stage))

cat("\n--- Main results ---\n")
etable(
  ols_income, iv_income, iv_employ, iv_peak,
  headers = c("OLS Income", "IV Income", "IV Employ", "IV Peak Season"),
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
het_models  <- Filter(Negate(is.null), list(iv_coastal, iv_not_coastal))
het_headers <- c("Coastal", "Not Coastal")[
  !c(is.null(iv_coastal), is.null(iv_not_coastal))
]
if (length(het_models) == 0) {
  cat("No coastal-type subgroups had sufficient municipalities to estimate.\n")
} else if (length(het_models) == 1) {
  etable(het_models[[1]], headers = het_headers, se.below = TRUE)
} else if (length(het_models) == 2) {
  etable(het_models[[1]], het_models[[2]], headers = het_headers, se.below = TRUE)
}

cat("\n--- Formal coastal interaction test ---\n")
cat("Base group = Coastal. Coefficient on afai_x_not_coastal = difference in effect.\n")
print(summary(iv_coastal_interact))
cat("\nWald test: is the not-coastal effect significantly different from coastal?\n")
tryCatch(
  print(wald(iv_coastal_interact, "afai_x_not_coastal")),
  error = function(e) cat("Wald test not available:", conditionMessage(e), "\n")
)
} else {
  etable(het_models[[1]], het_models[[2]], het_models[[3]], headers = het_headers, se.below = TRUE)
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

if (length(het_models) >= 1) {
  sink(file.path(path_results, "regression_heterogeneity_latex.tex"))
  do.call(etable, c(het_models,
                    list(headers = het_headers, se.below = TRUE, tex = TRUE,
                         title = "Heterogeneity by Coastal Proximity, DR 2017--2025")))
  sink()
}

sink(file.path(path_results, "regression_coastal_interact_latex.tex"))
etable(iv_coastal_interact,
       se.below = TRUE, tex = TRUE,
       title = "Formal Test: Differential Sargassum Effect by Coastal Proximity, DR 2017--2025")
sink()

# ── Coefficient plot ──────────────────────────────────────────
coef_data <- data.frame(
  model    = c("OLS", "IV (Main)", "IV (Peak\nSeason)"),
  estimate = c(coef(ols_income)["afai_sargassum"],
               coef(iv_income)["fit_afai_sargassum"],
               coef(iv_peak)["fit_afai_sargassum"]),
  se       = c(se(ols_income)["afai_sargassum"],
               se(iv_income)["fit_afai_sargassum"],
               se(iv_peak)["fit_afai_sargassum"])
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
    x = "Coefficient on Sargassum intensity (mean excess AFAI above 0.001)",
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
  estimate = c(coef(iv_t1)["fit_afai_sargassum"],
               coef(iv_t2)["fit_afai_sargassum"],
               coef(iv_t3)["fit_afai_sargassum"]),
  se       = c(se(iv_t1)["fit_afai_sargassum"],
               se(iv_t2)["fit_afai_sargassum"],
               se(iv_t3)["fit_afai_sargassum"])
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
    x = "Coefficient on Sargassum intensity (mean excess AFAI above 0.001)",
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
