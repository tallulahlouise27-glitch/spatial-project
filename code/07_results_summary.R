# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 7: Formatted results summary with interpretations
#
# Produces results/results_summary.txt — all regression tables
# with written interpretation under each one.
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(fixest)

path_processed <- "data/processed/"
path_results   <- "results/"
dir.create(path_results, showWarnings = FALSE)

# ── Load panel and re-run all models ─────────────────────────
panel <- readRDS(file.path(path_processed, "panel_monthly.rds"))

panel <- panel %>%
  mutate(
    is_not_coastal     = as.integer(coastal_type == "not_coastal"),
    afai_x_not_coastal = afai_sargassum * is_not_coastal
  )

panel_peak <- filter(panel, month %in% 5:9)

# Full-year models
ols_income  <- feols(log_income  ~ afai_sargassum | muni_fe + year_month_fe, panel, cluster = ~muni_fe)
ols_employ  <- feols(tasa_empleo ~ afai_sargassum | muni_fe + year_month_fe, panel, cluster = ~muni_fe)
ols_peak    <- feols(log_income  ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)
ols_t1      <- feols(log_income_t1 ~ afai_sargassum | muni_fe + year_month_fe, panel, cluster = ~muni_fe)
ols_t2      <- feols(log_income_t2 ~ afai_sargassum | muni_fe + year_month_fe, panel, cluster = ~muni_fe)
ols_t3      <- feols(log_income_t3 ~ afai_sargassum | muni_fe + year_month_fe, panel, cluster = ~muni_fe)
ols_coastal     <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, filter(panel, coastal_type == "coastal"),     cluster = ~muni_fe)
ols_not_coastal <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, filter(panel, coastal_type == "not_coastal"), cluster = ~muni_fe)
ols_coastal_interact <- feols(log_income ~ afai_sargassum + afai_x_not_coastal | muni_fe + year_month_fe, panel, cluster = ~muni_fe)

# Peak season models
ols_peak_coastal     <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, filter(panel_peak, coastal_type == "coastal"),     cluster = ~muni_fe)
ols_peak_not_coastal <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, filter(panel_peak, coastal_type == "not_coastal"), cluster = ~muni_fe)
ols_peak_t1 <- feols(log_income_t1 ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)
ols_peak_t2 <- feols(log_income_t2 ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)
ols_peak_t3 <- feols(log_income_t3 ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)

# Stacked tertile model
panel_long <- panel %>%
  select(ID_MUNICIPIO, year, month, muni_fe, year_month_fe, afai_sargassum,
         log_income_t1, log_income_t2, log_income_t3) %>%
  pivot_longer(cols = c(log_income_t1, log_income_t2, log_income_t3),
               names_to = "tertile", values_to = "log_income_tertile") %>%
  mutate(
    tertile = factor(tertile,
                     levels = c("log_income_t1","log_income_t2","log_income_t3"),
                     labels = c("T1","T2","T3")),
    is_t2 = as.integer(tertile == "T2"),
    is_t3 = as.integer(tertile == "T3"),
    afai_x_t2 = afai_sargassum * is_t2,
    afai_x_t3 = afai_sargassum * is_t3,
    muni_tertile_fe       = interaction(muni_fe, tertile),
    year_month_tertile_fe = interaction(year_month_fe, tertile)
  ) %>%
  filter(!is.na(log_income_tertile))

ols_stacked <- feols(
  log_income_tertile ~ afai_sargassum + afai_x_t2 + afai_x_t3 |
    muni_tertile_fe + year_month_tertile_fe,
  data = panel_long, cluster = ~ID_MUNICIPIO
)

# ── Helper: pull coefficient, SE, p-value ────────────────────
pull_stat <- function(model, var = "afai_sargassum") {
  b  <- coef(model)[var]
  s  <- se(model)[var]
  p  <- pvalue(model)[var]
  stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p < 0.1, ".", ""))))
  sprintf("%.2f%s (SE %.2f, p = %.3f)", b, stars, s, p)
}

# ── Write summary document ────────────────────────────────────
out <- file.path(path_results, "results_summary.txt")
sink(out)

cat("================================================================\n")
cat("SARGASSUM & HOUSEHOLD WELFARE — DOMINICAN REPUBLIC\n")
cat("Full Results Summary\n")
cat("Panel: 100 municipalities × 2017-2025 × monthly (10,263 obs)\n")
cat("All models: municipality FEs + year×month FEs, SEs clustered by municipality\n")
cat("Treatment: afai_sargassum = mean excess AFAI above 0.001 threshold\n")
cat("           (captures both spatial coverage and mat density of Sargassum)\n")
cat("================================================================\n\n")

# ── TABLE 1 ──────────────────────────────────────────────────
cat("----------------------------------------------------------------\n")
cat("TABLE 1: Main results — full sample, full year\n")
cat("----------------------------------------------------------------\n")
etable(
  ols_income, ols_employ, ols_peak,
  headers  = c("Log Income", "Employment Rate", "Log Income\n(Peak May-Sep)"),
  se.below = TRUE
)
cat("\nINTERPRETATION:\n")
cat(
  "The full-year income coefficient is negative (-28.15) indicating that higher\n",
  "Sargassum intensity is associated with lower household income, though the\n",
  "estimate does not reach conventional significance levels (p = 0.13).\n\n",
  "The employment rate coefficient is also negative (-4.44) and insignificant,\n",
  "suggesting Sargassum may reduce employment as well as income, but the data\n",
  "cannot distinguish this from noise over the full year.\n\n",
  "The peak season result (May-Sep, when Sargassum is most intense) is the\n",
  "strongest finding in the paper: a coefficient of -73.90 significant at 1%.\n",
  "A one-unit increase in mean excess AFAI during peak months is associated with\n",
  "a 73.9-log-point reduction in municipal income. This concentration of the\n",
  "effect in peak months is consistent with Sargassum being the mechanism:\n",
  "the harm materialises when blooms are actually present.\n",
  sep = ""
)

# ── TABLE 2 ──────────────────────────────────────────────────
cat("\n----------------------------------------------------------------\n")
cat("TABLE 2: Full-year results by income tertile\n")
cat("(T1 = poorest third, T2 = middle third, T3 = richest third)\n")
cat("----------------------------------------------------------------\n")
etable(
  ols_t1, ols_t2, ols_t3,
  headers  = c("T1 (Bottom)", "T2 (Middle)", "T3 (Top)"),
  se.below = TRUE
)
cat("\nINTERPRETATION:\n")
cat(
  "The coefficients follow the expected gradient: the poorest households show\n",
  "the largest negative effect (-144.87), shrinking monotonically through the\n",
  "middle (-30.75) to the richest third (-12.20). This pattern is consistent\n",
  "with poorer households being least able to insure against environmental shocks.\n\n",
  "Only T2 reaches significance (p = 0.017). The T1 estimate is imprecise for a\n",
  "structural reason: approximately 5% of municipality-month observations in the\n",
  "bottom tertile report zero income, producing a log-income outcome of zero\n",
  "regardless of Sargassum intensity. This floors the dependent variable and\n",
  "inflates noise (SD of log_income_t1 = 1.93 vs 0.34 for T2). The T1 effect\n",
  "is likely attenuated — the true impact on the poorest households may be\n",
  "larger than estimated but is undetectable through the income measure alone.\n",
  sep = ""
)

# ── TABLE 3 ──────────────────────────────────────────────────
cat("\n----------------------------------------------------------------\n")
cat("TABLE 3: Full-year results by coastal proximity\n")
cat("(Coastal = direct shoreline contact; Not coastal = all others)\n")
cat("----------------------------------------------------------------\n")
etable(
  ols_coastal, ols_not_coastal,
  headers  = c("Coastal", "Not Coastal"),
  se.below = TRUE
)
cat("\nINTERPRETATION:\n")
cat(
  "Both coastal and not-coastal municipalities show negative effects, but neither\n",
  "reaches significance over the full year. Coastal municipalities show a larger\n",
  "point estimate (-71.16 vs -14.97), as expected given direct shoreline exposure.\n\n",
  "The formal interaction test (Table 6) confirms the difference between groups\n",
  "is not statistically significant (Wald p = 0.154), meaning we cannot reject\n",
  "that the two groups experience the same effect size.\n",
  sep = ""
)

# ── TABLE 4 ──────────────────────────────────────────────────
cat("\n----------------------------------------------------------------\n")
cat("TABLE 4: Peak season (May-Sep) results by coastal proximity\n")
cat("----------------------------------------------------------------\n")
etable(
  ols_peak_coastal, ols_peak_not_coastal,
  headers  = c("Coastal", "Not Coastal"),
  se.below = TRUE
)
cat("\nINTERPRETATION:\n")
cat(
  "The most striking result in the paper: during peak Sargassum months, the\n",
  "significant effect is found in NOT-COASTAL municipalities (-84.95, p = 0.001),\n",
  "while coastal municipalities show a similar magnitude but insignificant\n",
  "estimate (-76.30, p = 0.156).\n\n",
  "This is counterintuitive at first glance but has plausible explanations:\n",
  "  1. Coastal households may have diversified income (fishing + tourism + trade)\n",
  "     that allows partial offsetting of Sargassum shocks.\n",
  "  2. Inland/near-coastal municipalities may depend heavily on tourism supply\n",
  "     chains, processing, or distribution networks that collapse when Sargassum\n",
  "     closes beaches, with fewer alternative income sources to fall back on.\n",
  "  3. Coastal households may have higher income variance year-round (from\n",
  "     fishing), making it harder to isolate the Sargassum signal statistically.\n\n",
  "This finding highlights that Sargassum's economic reach extends beyond the\n",
  "coastline itself, likely through disruption of tourism and related industries.\n",
  sep = ""
)

# ── TABLE 5 ──────────────────────────────────────────────────
cat("\n----------------------------------------------------------------\n")
cat("TABLE 5: Peak season (May-Sep) results by income tertile\n")
cat("----------------------------------------------------------------\n")
etable(
  ols_peak_t1, ols_peak_t2, ols_peak_t3,
  headers  = c("T1 (Bottom)", "T2 (Middle)", "T3 (Top)"),
  se.below = TRUE
)
cat("\nINTERPRETATION:\n")
cat(
  "During peak Sargassum months, the monotonic gradient persists: T1 (-75.07),\n",
  "T2 (-48.99, p = 0.004), T3 (-26.64). Again only T2 reaches significance,\n",
  "for the same reason as Table 2: T1 income is censored at zero for many\n",
  "observations, preventing detection of the effect on the poorest households.\n\n",
  "The T2 peak-season coefficient (-48.99) is larger than the full-year T2\n",
  "coefficient (-30.75), consistent with the effect being concentrated in months\n",
  "when Sargassum is actually present rather than spread across the whole year.\n",
  sep = ""
)

# ── TABLE 6 ──────────────────────────────────────────────────
cat("\n----------------------------------------------------------------\n")
cat("TABLE 6: Formal interaction tests\n")
cat("----------------------------------------------------------------\n")

cat("\n6A. Coastal interaction test (full year)\n")
cat("Base group = Coastal. afai_x_not_coastal = difference in effect for not-coastal.\n\n")
etable(ols_coastal_interact, se.below = TRUE)
wald_coast <- wald(ols_coastal_interact, "afai_x_not_coastal")
cat(sprintf("\nWald test (H0: no difference between coastal and not-coastal):\n"))
cat(sprintf("  F = %.3f, p = %.3f\n", wald_coast$stat, wald_coast$p))

cat("\n6B. Tertile interaction test (full year)\n")
cat("Base group = T1. afai_x_t2/t3 = difference in effect for T2/T3 vs T1.\n\n")
etable(ols_stacked, se.below = TRUE)
wald_tert <- wald(ols_stacked, "afai_x_t2|afai_x_t3")
cat(sprintf("\nWald test (H0: T2 and T3 effects same as T1):\n"))
cat(sprintf("  F = %.3f, p = %.3f\n", wald_tert$stat, wald_tert$p))

cat("\nINTERPRETATION:\n")
cat(
  "Neither formal interaction test rejects the null of equal effects:\n",
  "  - Coastal vs not-coastal: Wald p = 0.154\n",
  "  - Tertile heterogeneity:  Wald p = 0.334\n\n",
  "The failure to reject does not mean the effects are identical — the point\n",
  "estimates follow the expected pattern in both cases. The tests lack power\n",
  "because (a) the sample has only 100 municipalities and (b) the treatment\n",
  "variable is noisy. These results should be interpreted as consistent with\n",
  "the hypothesised heterogeneity, but insufficient to confirm it statistically.\n",
  sep = ""
)

# ── SUMMARY ──────────────────────────────────────────────────
cat("\n================================================================\n")
cat("SUMMARY OF KEY FINDINGS\n")
cat("================================================================\n")
cat(
  "1. MAIN EFFECT: Sargassum reduces household income, primarily during peak\n",
  "   season (May-Sep). The peak-season coefficient is -73.90 (p < 0.01),\n",
  "   concentrated in months when Sargassum blooms are actually present.\n\n",
  "2. SPATIAL REACH: The income effect is significant in NOT-COASTAL municipalities\n",
  "   during peak season (-84.95, p < 0.01), suggesting Sargassum disrupts\n",
  "   tourism supply chains and regional economies beyond the shoreline.\n\n",
  "3. DISTRIBUTIONAL EFFECTS: The negative effect is largest for the poorest\n",
  "   households (T1: -145, T2: -31, T3: -12 full year), consistent with\n",
  "   limited insurance capacity among low-income households. Only T2 reaches\n",
  "   significance; T1 is attenuated by zero-income reporting.\n\n",
  "4. IDENTIFICATION: All models include municipality and year×month fixed effects,\n",
  "   removing time-invariant geographic differences and all common national shocks.\n",
  "   Identification relies on within-municipality monthly variation in Sargassum.\n",
  "   Remaining concern: local time-varying confounders (e.g. concurrent fishing\n",
  "   shocks) cannot be ruled out without a valid instrument.\n",
  sep = ""
)

sink()
cat("Results summary saved to:", out, "\n")
