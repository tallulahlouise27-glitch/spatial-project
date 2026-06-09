# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 5: Panel regressions with two-way fixed effects
#
# Model:
#   log(income)_it = β × afai_sargassum_it + α_i + γ_t + ε_it
#
#   Where:
#     i = municipality, t = year × month
#     α_i = municipality fixed effect (absorbs time-invariant geography,
#           baseline wealth, distance from coast, institutions)
#     γ_t = year × month fixed effect (absorbs all common national shocks
#           at monthly frequency — seasonality, economic cycles, etc.)
#
# Identification: within-municipality, within-time-period variation in
# Sargassum intensity. OLS with two-way FEs — strongly suggestive of
# causality but not fully causal due to remaining time-varying confounders.
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(fixest)

path_processed <- "data/processed/"
path_results   <- "results/"
path_figures   <- "figures/"
dir.create(path_results, showWarnings = FALSE)
dir.create(path_figures, showWarnings = FALSE)

# ── Load analysis panel ───────────────────────────────────────
panel <- readRDS(file.path(path_processed, "panel_monthly.rds"))

# ── Create interaction terms for heterogeneity tests ─────────
panel <- panel %>%
  mutate(
    is_not_coastal     = as.integer(coastal_type == "not_coastal"),
    afai_x_not_coastal = afai_sargassum * is_not_coastal
  )

cat("Panel summary:\n")
cat("Observations:", nrow(panel), "\n")
cat("Municipalities:", n_distinct(panel$ID_MUNICIPIO), "\n")
cat("Years:", paste(sort(unique(panel$year)), collapse = ", "), "\n\n")

# ── Model 1: Main — effect on log income ─────────────────────
ols_income <- feols(
  log_income ~ afai_sargassum | muni_fe + year_month_fe,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 2: Effect on employment rate ───────────────────────
ols_employ <- feols(
  tasa_empleo ~ afai_sargassum | muni_fe + year_month_fe,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 3: Peak season only (May–Sep) ──────────────────────
# Restricts to months when Sargassum is most intense.
ols_peak <- feols(
  log_income ~ afai_sargassum | muni_fe + year_month_fe,
  data    = filter(panel, month %in% 5:9),
  cluster = ~muni_fe
)

# ── Models 4–6: By income tertile ────────────────────────────
# Tests whether Sargassum hits poorer households harder.
# Expected: most negative for T1, smaller for T3.

ols_t1 <- feols(
  log_income_t1 ~ afai_sargassum | muni_fe + year_month_fe,
  data    = panel,
  cluster = ~muni_fe
)

ols_t2 <- feols(
  log_income_t2 ~ afai_sargassum | muni_fe + year_month_fe,
  data    = panel,
  cluster = ~muni_fe
)

ols_t3 <- feols(
  log_income_t3 ~ afai_sargassum | muni_fe + year_month_fe,
  data    = panel,
  cluster = ~muni_fe
)

# ── Model 7: Stacked regression with tertile interactions ─────
# Single regression that formally tests whether the Sargassum effect
# differs across income tertiles.
#
# Coefficients:
#   afai_sargassum = effect on T1 (base group, poorest third)
#   afai_x_t2      = difference in effect for T2 vs T1
#   afai_x_t3      = difference in effect for T3 vs T1

panel_long <- panel %>%
  select(ID_MUNICIPIO, year, month, muni_fe, year_month_fe,
         afai_sargassum,
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
    is_t2          = as.integer(tertile == "T2"),
    is_t3          = as.integer(tertile == "T3"),
    afai_x_t2      = afai_sargassum * is_t2,
    afai_x_t3      = afai_sargassum * is_t3,
    muni_tertile_fe       = interaction(muni_fe, tertile),
    year_month_tertile_fe = interaction(year_month_fe, tertile)
  ) %>%
  filter(!is.na(log_income_tertile))

ols_stacked <- feols(
  log_income_tertile ~ afai_sargassum + afai_x_t2 + afai_x_t3 |
    muni_tertile_fe + year_month_tertile_fe,
  data    = panel_long,
  cluster = ~ID_MUNICIPIO
)

# ── Models 8–9: By coastal proximity ─────────────────────────
MIN_MUNICIPALITIES <- 10

coastal_counts <- panel %>%
  group_by(coastal_type) %>%
  summarise(n_munis = n_distinct(ID_MUNICIPIO), n_obs = n(), .groups = "drop")

cat("Sample size by coastal type:\n")
print(coastal_counts)
cat("\n")

safe_ols <- function(type, data) {
  sub <- filter(data, coastal_type == type)
  n   <- n_distinct(sub$ID_MUNICIPIO)
  if (n < MIN_MUNICIPALITIES) {
    cat("SKIPPING", type, "(only", n, "municipalities)\n")
    return(NULL)
  }
  feols(log_income ~ afai_sargassum | muni_fe + year_month_fe,
        data = sub, cluster = ~muni_fe)
}

ols_coastal     <- safe_ols("coastal",     panel)
ols_not_coastal <- safe_ols("not_coastal", panel)

# ── Model 10: Formal coastal interaction test ─────────────────
# Base group = coastal. Coefficient on afai_x_not_coastal is the
# difference in the Sargassum effect for not-coastal municipalities.
ols_coastal_interact <- feols(
  log_income ~ afai_sargassum + afai_x_not_coastal | muni_fe + year_month_fe,
  data    = panel,
  cluster = ~muni_fe
)

# ── Models 11–12: Peak season by coastal proximity ───────────
panel_peak <- filter(panel, month %in% 5:9)

safe_ols_peak <- function(type, data) {
  sub <- filter(data, coastal_type == type)
  n   <- n_distinct(sub$ID_MUNICIPIO)
  if (n < MIN_MUNICIPALITIES) {
    cat("SKIPPING", type, "(only", n, "municipalities)\n")
    return(NULL)
  }
  feols(log_income ~ afai_sargassum | muni_fe + year_month_fe,
        data = sub, cluster = ~muni_fe)
}

ols_peak_coastal     <- safe_ols_peak("coastal",     panel_peak)
ols_peak_not_coastal <- safe_ols_peak("not_coastal", panel_peak)

# ── Models 13–15: Peak season by income tertile ───────────────
ols_peak_t1 <- feols(
  log_income_t1 ~ afai_sargassum | muni_fe + year_month_fe,
  data    = panel_peak,
  cluster = ~muni_fe
)

ols_peak_t2 <- feols(
  log_income_t2 ~ afai_sargassum | muni_fe + year_month_fe,
  data    = panel_peak,
  cluster = ~muni_fe
)

ols_peak_t3 <- feols(
  log_income_t3 ~ afai_sargassum | muni_fe + year_month_fe,
  data    = panel_peak,
  cluster = ~muni_fe
)

# ── Print results ─────────────────────────────────────────────
cat("====== REGRESSION RESULTS ======\n\n")

cat("--- Main results ---\n")
etable(
  ols_income, ols_employ, ols_peak,
  headers  = c("Income", "Employment", "Peak Season\n(May-Sep)"),
  se.below = TRUE
)

cat("\n--- Heterogeneity by income tertile (separate regressions) ---\n")
etable(
  ols_t1, ols_t2, ols_t3,
  headers  = c("T1 (Bottom)", "T2 (Middle)", "T3 (Top)"),
  se.below = TRUE
)

cat("\n--- Stacked regression with tertile interactions ---\n")
cat("Base group = T1. Interaction coefficients show difference from T1.\n")
print(summary(ols_stacked))
cat("\nWald test: are T2 and T3 effects significantly different from T1?\n")
tryCatch(
  print(wald(ols_stacked, "afai_x_t2|afai_x_t3")),
  error = function(e) cat("Wald test not available:", conditionMessage(e), "\n")
)

cat("\n--- Heterogeneity by coastal proximity ---\n")
het_models  <- Filter(Negate(is.null), list(ols_coastal, ols_not_coastal))
het_headers <- c("Coastal", "Not Coastal")[
  !c(is.null(ols_coastal), is.null(ols_not_coastal))
]
if (length(het_models) == 0) {
  cat("No coastal-type subgroups had sufficient municipalities to estimate.\n")
} else if (length(het_models) == 1) {
  etable(het_models[[1]], headers = het_headers, se.below = TRUE)
} else {
  etable(het_models[[1]], het_models[[2]], headers = het_headers, se.below = TRUE)
}

cat("\n--- Peak season by coastal proximity ---\n")
peak_coast_models  <- Filter(Negate(is.null), list(ols_peak_coastal, ols_peak_not_coastal))
peak_coast_headers <- c("Coastal", "Not Coastal")[
  !c(is.null(ols_peak_coastal), is.null(ols_peak_not_coastal))
]
if (length(peak_coast_models) == 2) {
  etable(peak_coast_models[[1]], peak_coast_models[[2]],
         headers = peak_coast_headers, se.below = TRUE)
}

cat("\n--- Peak season by income tertile ---\n")
etable(
  ols_peak_t1, ols_peak_t2, ols_peak_t3,
  headers  = c("T1 (Bottom)", "T2 (Middle)", "T3 (Top)"),
  se.below = TRUE
)

cat("\n--- Formal coastal interaction test ---\n")
cat("Base group = Coastal. Coefficient on afai_x_not_coastal = difference in effect.\n")
print(summary(ols_coastal_interact))
cat("\nWald test: is the not-coastal effect significantly different from coastal?\n")
tryCatch(
  print(wald(ols_coastal_interact, "afai_x_not_coastal")),
  error = function(e) cat("Wald test not available:", conditionMessage(e), "\n")
)

# ── Export tables ─────────────────────────────────────────────
sink(file.path(path_results, "regression_main.txt"))
etable(
  ols_income, ols_employ, ols_peak,
  headers  = c("Income", "Employment", "Peak Season (May-Sep)"),
  se.below = TRUE,
  title    = "Effect of Sargassum Exposure on Household Welfare, DR 2017-2025"
)
sink()

sink(file.path(path_results, "regression_main_latex.tex"))
etable(
  ols_income, ols_employ, ols_peak,
  headers  = c("Income", "Employment", "Peak Season (May--Sep)"),
  se.below = TRUE,
  tex      = TRUE,
  title    = "Effect of Sargassum Exposure on Household Welfare, DR 2017--2025"
)
sink()

sink(file.path(path_results, "regression_tertile_latex.tex"))
etable(
  ols_t1, ols_t2, ols_t3,
  headers  = c("T1 (Bottom Third)", "T2 (Middle Third)", "T3 (Top Third)"),
  se.below = TRUE,
  tex      = TRUE,
  title    = "Effect of Sargassum Exposure by Income Tertile, DR 2017--2025"
)
sink()

sink(file.path(path_results, "regression_stacked_latex.tex"))
etable(
  ols_stacked,
  se.below = TRUE,
  tex      = TRUE,
  title    = "Stacked Regression: Differential Sargassum Effects by Income Tertile, DR 2017--2025"
)
sink()

if (length(het_models) >= 1) {
  sink(file.path(path_results, "regression_heterogeneity_latex.tex"))
  do.call(etable, c(het_models,
                    list(headers = het_headers, se.below = TRUE, tex = TRUE,
                         title = "Heterogeneity by Coastal Proximity, DR 2017--2025")))
  sink()
}

if (length(peak_coast_models) == 2) {
  sink(file.path(path_results, "regression_peak_coastal_latex.tex"))
  etable(peak_coast_models[[1]], peak_coast_models[[2]],
         headers = peak_coast_headers, se.below = TRUE, tex = TRUE,
         title = "Peak Season Sargassum Effect by Coastal Proximity, DR 2017--2025")
  sink()
}

sink(file.path(path_results, "regression_peak_tertile_latex.tex"))
etable(
  ols_peak_t1, ols_peak_t2, ols_peak_t3,
  headers  = c("T1 (Bottom Third)", "T2 (Middle Third)", "T3 (Top Third)"),
  se.below = TRUE, tex = TRUE,
  title    = "Peak Season Sargassum Effect by Income Tertile, DR 2017--2025"
)
sink()

sink(file.path(path_results, "regression_coastal_interact_latex.tex"))
etable(ols_coastal_interact,
       se.below = TRUE, tex = TRUE,
       title = "Formal Test: Differential Sargassum Effect by Coastal Proximity, DR 2017--2025")
sink()

# ── Coefficient plot: main models ────────────────────────────
coef_data <- data.frame(
  model    = c("Income\n(full year)", "Employment\n(full year)", "Income\n(peak season)"),
  estimate = c(coef(ols_income)["afai_sargassum"],
               coef(ols_employ)["afai_sargassum"],
               coef(ols_peak)["afai_sargassum"]),
  se       = c(se(ols_income)["afai_sargassum"],
               se(ols_employ)["afai_sargassum"],
               se(ols_peak)["afai_sargassum"])
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
    x        = "Coefficient on Sargassum intensity (mean excess AFAI above 0.001)",
    y        = NULL,
    title    = "Effect of Sargassum Exposure on Household Welfare",
    subtitle = "Municipality + year×month fixed effects. Clustered SEs. 95% confidence intervals."
  ) +
  theme_minimal()

ggsave(file.path(path_figures, "coef_plot.png"), p, width = 7, height = 4, dpi = 150)
cat("Coefficient plot saved to figures/coef_plot.png\n")

# ── Coefficient plot: by income tertile ──────────────────────
tertile_coef <- data.frame(
  group    = c("T1 (Bottom third)", "T2 (Middle third)", "T3 (Top third)"),
  estimate = c(coef(ols_t1)["afai_sargassum"],
               coef(ols_t2)["afai_sargassum"],
               coef(ols_t3)["afai_sargassum"]),
  se       = c(se(ols_t1)["afai_sargassum"],
               se(ols_t2)["afai_sargassum"],
               se(ols_t3)["afai_sargassum"])
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
    x        = "Coefficient on Sargassum intensity (mean excess AFAI above 0.001)",
    y        = NULL,
    title    = "Sargassum Effect by Income Tertile",
    subtitle = "Municipality + year×month fixed effects. Clustered SEs. 95% CIs."
  ) +
  theme_minimal()

ggsave(file.path(path_figures, "coef_plot_tertile.png"),
       p_tertile, width = 7, height = 3, dpi = 150)
cat("Tertile coefficient plot saved to figures/coef_plot_tertile.png\n")

cat("\n====== DONE ======\n")
cat("Results saved to results/\n")
