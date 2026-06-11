# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 8: Export all results tables to a single HTML document
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(fixest)
library(knitr)

path_processed <- "data/processed/"
path_results   <- "results/"
dir.create(path_results, showWarnings = FALSE)

# ── Load panels ───────────────────────────────────────────────
panel           <- readRDS(file.path(path_processed, "panel_monthly.rds"))
panel_quarterly <- readRDS(file.path(path_processed, "panel_quarterly.rds"))

# ── Derived variables ─────────────────────────────────────────
panel <- panel %>%
  mutate(
    is_not_coastal       = as.integer(coastal_type == "not_coastal"),
    afai_x_not_coastal   = afai_sargassum * is_not_coastal,
    tertile_factor       = factor(tertile, levels = c("T1","T2","T3")),
    is_t2                = as.integer(tertile == "T2"),
    is_t3                = as.integer(tertile == "T3"),
    afai_x_t2            = afai_sargassum * is_t2,
    afai_x_t3            = afai_sargassum * is_t3,
    muni_tertile_fe      = interaction(muni_fe, tertile_factor),
    year_month_tertile_fe= interaction(year_month_fe, tertile_factor),
    is_informal_hh       = as.integer(n_empleados > 0 & share_informal > 0.5),
    afai_x_informal      = afai_sargassum * is_informal_hh
  )

cuts_nz <- panel %>%
  filter(ingreso_pc > 0) %>%
  group_by(year) %>%
  summarise(t1_cut_nz = quantile(ingreso_pc, 1/3, na.rm=TRUE),
            t2_cut_nz = quantile(ingreso_pc, 2/3, na.rm=TRUE), .groups="drop")

panel <- panel %>%
  left_join(cuts_nz, by="year") %>%
  mutate(
    tertile_nz    = case_when(
      ingreso_pc == 0         ~ NA_character_,
      ingreso_pc <= t1_cut_nz ~ "T1",
      ingreso_pc <= t2_cut_nz ~ "T2",
      TRUE                    ~ "T3"
    ),
    log_income_nz = if_else(ingreso_pc > 0, log_income, NA_real_)
  )

panel_quarterly <- panel_quarterly %>%
  mutate(is_not_coastal = as.integer(coastal_type == "not_coastal"),
         afai_x_nc      = afai_sargassum * is_not_coastal)

panel_peak <- filter(panel, month %in% 5:9)
pq_peak    <- filter(panel_quarterly, quarter %in% 2:3)

# Positive-income-only quarterly tertile means (re-aggregated from household level)
tertile_nz_q <- panel %>%
  filter(ingreso_pc > 0) %>%
  group_by(muni_fe, year, quarter, tertile) %>%
  summarise(log_income_nz = mean(log_income, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = tertile, values_from = log_income_nz,
              names_prefix = "log_income_nz_")

panel_quarterly <- panel_quarterly %>%
  left_join(tertile_nz_q, by = c("muni_fe", "year", "quarter"))

pq_peak <- filter(panel_quarterly, quarter %in% 2:3)

# ── Run all models ────────────────────────────────────────────
cat("Running regressions...\n")

# Quarterly — main
ols_q_main    <- feols(log_income  ~ afai_sargassum | muni_fe + year_quarter_fe, panel_quarterly, cluster=~muni_fe)
ols_q_peak    <- feols(log_income  ~ afai_sargassum | muni_fe + year_quarter_fe, pq_peak,         cluster=~muni_fe)
ols_q_employ  <- feols(tasa_empleo ~ afai_sargassum | muni_fe + year_quarter_fe, panel_quarterly, cluster=~muni_fe)

# Quarterly — income tertiles
ols_q_t1      <- feols(log_income_t1 ~ afai_sargassum | muni_fe + year_quarter_fe, panel_quarterly, cluster=~muni_fe)
ols_q_t2      <- feols(log_income_t2 ~ afai_sargassum | muni_fe + year_quarter_fe, panel_quarterly, cluster=~muni_fe)
ols_q_t3      <- feols(log_income_t3 ~ afai_sargassum | muni_fe + year_quarter_fe, panel_quarterly, cluster=~muni_fe)
ols_q_peak_t1 <- feols(log_income_t1 ~ afai_sargassum | muni_fe + year_quarter_fe, pq_peak,         cluster=~muni_fe)
ols_q_peak_t2 <- feols(log_income_t2 ~ afai_sargassum | muni_fe + year_quarter_fe, pq_peak,         cluster=~muni_fe)
ols_q_peak_t3 <- feols(log_income_t3 ~ afai_sargassum | muni_fe + year_quarter_fe, pq_peak,         cluster=~muni_fe)

# Quarterly — positive-income-only tertiles, peak season
ols_q_pk_t1_nz <- feols(log_income_nz_T1 ~ afai_sargassum | muni_fe + year_quarter_fe, pq_peak, cluster=~muni_fe)
ols_q_pk_t2_nz <- feols(log_income_nz_T2 ~ afai_sargassum | muni_fe + year_quarter_fe, pq_peak, cluster=~muni_fe)
ols_q_pk_t3_nz <- feols(log_income_nz_T3 ~ afai_sargassum | muni_fe + year_quarter_fe, pq_peak, cluster=~muni_fe)

# Quarterly — coastal proximity
ols_q_coastal     <- feols(log_income ~ afai_sargassum | muni_fe + year_quarter_fe, filter(panel_quarterly, coastal_type=="coastal"),     cluster=~muni_fe)
ols_q_not_coastal <- feols(log_income ~ afai_sargassum | muni_fe + year_quarter_fe, filter(panel_quarterly, coastal_type=="not_coastal"), cluster=~muni_fe)
ols_q_pk_coastal  <- feols(log_income ~ afai_sargassum | muni_fe + year_quarter_fe, filter(pq_peak, coastal_type=="coastal"),             cluster=~muni_fe)
ols_q_pk_nc       <- feols(log_income ~ afai_sargassum | muni_fe + year_quarter_fe, filter(pq_peak, coastal_type=="not_coastal"),         cluster=~muni_fe)
ols_q_coastal_int <- feols(log_income ~ afai_sargassum + afai_x_nc | muni_fe + year_quarter_fe, panel_quarterly, cluster=~muni_fe)
ols_q_pk_coast_int<- feols(log_income ~ afai_sargassum + afai_x_nc | muni_fe + year_quarter_fe, pq_peak,         cluster=~muni_fe)

# Monthly — robustness
ols_income_pc  <- feols(log_income_pc ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster=~muni_fe)
ols_employ     <- feols(tasa_empleo   ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster=~muni_fe)
ols_peak       <- feols(log_income_pc ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster=~muni_fe)

# FE progression (for Table 0)
ols_income     <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, panel, cluster=~muni_fe)

cat("All models estimated.\n")

# ── Helpers to build clean HTML tables from model objects ──────

stars <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
}

fmt_coef <- function(model, var) {
  b <- coef(model)[var]
  p <- pvalue(model)[var]
  if (is.na(b)) return("")
  paste0(formatC(b, format="f", digits=4), stars(p))
}

fmt_se <- function(model, var) {
  s <- se(model)[var]
  if (is.na(s)) return("")
  paste0("(", formatC(s, format="f", digits=4), ")")
}

# Build a standard 2-panel table (coefs + fit stats)
# vars: named list, e.g. list("Sargassum intensity" = "afai_sargassum")
reg_table_html <- function(models, col_names, vars,
                           fe_line = "Municipality + Year × Month",
                           title = NULL, note = NULL) {

  # --- coefficient rows ---
  coef_rows <- lapply(names(vars), function(label) {
    var <- vars[[label]]
    coef_vals <- sapply(models, fmt_coef, var=var)
    se_vals   <- sapply(models, fmt_se,   var=var)
    list(label=label, coef=coef_vals, se=se_vals)
  })

  n_obs <- sapply(models, function(m) formatC(nobs(m), format="d", big.mark=","))
  wr2   <- sapply(models, function(m) {
    v <- r2(m)["wr2"]
    if (is.na(v)) "—" else formatC(v, format="f", digits=4)
  })

  dep_vars <- sapply(models, function(m) {
    dv <- deparse(m$fml[[2]])
    switch(dv,
           "log_income"    = "Log income",
           "log_income_pc" = "Log income p.c.",
           "log_income_t1"    = "Log income (T1)",
           "log_income_t2"    = "Log income (T2)",
           "log_income_t3"    = "Log income (T3)",
           "log_income_nz_T1" = "Log income >0 (T1)",
           "log_income_nz_T2" = "Log income >0 (T2)",
           "log_income_nz_T3" = "Log income >0 (T3)",
           "log_income_nz" = "Log income (excl. zero)",
           "tasa_empleo"   = "Employment rate",
           dv)
  })

  n_cols <- length(models)
  th_inner <- paste0("<th>", col_names, "</th>", collapse="")
  thead <- paste0('<thead><tr><th></th>', th_inner, '</tr></thead>')

  tbody_rows <- character(0)

  # Dep var row — only shown when columns have different dependent variables
  if (length(unique(dep_vars)) > 1) {
    tbody_rows <- c(tbody_rows,
      paste0('<tr class="dep-row"><td>Dependent variable</td>',
             paste0('<td>', dep_vars, '</td>', collapse=""), '</tr>'))
  }

  # Coefficient rows
  for (cr in coef_rows) {
    tbody_rows <- c(tbody_rows,
      paste0('<tr class="coef-row"><td>', cr$label, '</td>',
             paste0('<td>', cr$coef, '</td>', collapse=""), '</tr>'),
      paste0('<tr class="se-row"><td></td>',
             paste0('<td>', cr$se, '</td>', collapse=""), '</tr>'))
  }

  # Fit stats
  tbody_rows <- c(tbody_rows,
    '<tr class="stat-divider"><td colspan="100%"></td></tr>',
    paste0('<tr class="stat-row"><td>Observations</td>',
           paste0('<td>', n_obs, '</td>', collapse=""), '</tr>'),
    paste0('<tr class="stat-row"><td>Within R²</td>',
           paste0('<td>', wr2, '</td>', collapse=""), '</tr>'))

  tbody <- paste0('<tbody>', paste(tbody_rows, collapse="\n"), '</tbody>')

  tbl <- paste0('<table class="results-table">', thead, tbody, '</table>')

  out <- character(0)
  if (!is.null(title)) out <- c(out, paste0('<h3 class="table-title">', title, '</h3>'))
  out <- c(out, tbl)
  if (!is.null(note))  out <- c(out, paste0('<p class="table-note">', note, '</p>'))
  paste(out, collapse="\n")
}

note_fe   <- "<em>Notes:</em> Municipality and year &times; month fixed effects. Standard errors clustered at municipality level in parentheses. *** p&lt;0.001, ** p&lt;0.01, * p&lt;0.05."
note_q    <- "<em>Notes:</em> Quarterly panel. Municipality and year &times; quarter fixed effects. Standard errors clustered at municipality level. *** p&lt;0.001, ** p&lt;0.01, * p&lt;0.05."
note_int  <- "<em>Notes:</em> Municipality and year &times; month fixed effects throughout. Interaction coefficients show the difference in the Sargassum effect relative to the base group. SEs clustered at municipality level. *** p&lt;0.001, ** p&lt;0.01, * p&lt;0.05."

sat_var <- list("Sargassum intensity (AFAI)" = "afai_sargassum")

# ── Fixed-effects progression table ──────────────────────────
fe_prog_html <- function(models, col_names, var='afai_sargassum', title=NULL, note=NULL) {
  fe_labels <- c("Municipality FE", "Year × Month FE")

  has_muni <- sapply(models, function(m) 'muni_fe'       %in% m$fixef_vars)
  has_time <- sapply(models, function(m) 'year_month_fe' %in% m$fixef_vars)

  n_obs <- sapply(models, function(m) formatC(nobs(m), format='d', big.mark=','))
  r2_val <- sapply(models, function(m) {
    wr <- r2(m)['wr2']
    if (is.na(wr)) formatC(r2(m)['r2'], format='f', digits=4)
    else           formatC(wr,           format='f', digits=4)
  })
  r2_label <- sapply(models, function(m) if (is.na(r2(m)['wr2'])) 'R²' else 'Within R²')
  # use a single label for the row; flag if mixed
  r2_row_label <- if (length(unique(r2_label)) == 1) r2_label[1] else 'R² / Within R²'

  th_inner <- paste0('<th>', col_names, '</th>', collapse='')
  thead    <- paste0('<thead><tr><th></th>', th_inner, '</tr></thead>')

  rows <- character(0)

  # Coefficient + SE
  rows <- c(rows,
    paste0('<tr class="coef-row"><td>Sargassum intensity (AFAI)</td>',
           paste0('<td>', sapply(models, fmt_coef, var=var), '</td>', collapse=''), '</tr>'),
    paste0('<tr class="se-row"><td></td>',
           paste0('<td>', sapply(models, fmt_se, var=var), '</td>', collapse=''), '</tr>'))

  # FE indicators
  rows <- c(rows, '<tr class="stat-divider"><td colspan="100%"></td></tr>')
  rows <- c(rows,
    paste0('<tr class="fe-row"><td>Municipality FE</td>',
           paste0('<td>', ifelse(has_muni, 'Yes', 'No'), '</td>', collapse=''), '</tr>'),
    paste0('<tr class="fe-row"><td>Year &times; Month FE</td>',
           paste0('<td>', ifelse(has_time, 'Yes', 'No'), '</td>', collapse=''), '</tr>'))

  # Fit stats
  rows <- c(rows, '<tr class="stat-divider"><td colspan="100%"></td></tr>',
    paste0('<tr class="stat-row"><td>Observations</td>',
           paste0('<td>', n_obs, '</td>', collapse=''), '</tr>'),
    paste0('<tr class="stat-row"><td>', r2_row_label, '</td>',
           paste0('<td>', r2_val, '</td>', collapse=''), '</tr>'))

  tbody <- paste0('<tbody>', paste(rows, collapse='\n'), '</tbody>')
  tbl   <- paste0('<table class="results-table">', thead, tbody, '</table>')

  out <- character(0)
  if (!is.null(title)) out <- c(out, paste0('<h3 class="table-title">', title, '</h3>'))
  out <- c(out, tbl)
  if (!is.null(note))  out <- c(out, paste0('<p class="table-note">', note, '</p>'))
  paste(out, collapse='\n')
}

ols_no_fe    <- feols(log_income ~ afai_sargassum,              panel, cluster=~muni_fe)
ols_muni_fe  <- feols(log_income ~ afai_sargassum | muni_fe,   panel, cluster=~muni_fe)
ols_time_fe  <- feols(log_income ~ afai_sargassum | year_month_fe, panel, cluster=~muni_fe)

html_fe_prog <- fe_prog_html(
  models    = list(ols_no_fe, ols_muni_fe, ols_time_fe, ols_income),
  col_names = c("No FE", "Municipality FE", "Year × Month FE", "Both FE"),
  title     = "Table 0. Effect of Sargassum on Log Household Income: Fixed Effects Progression",
  note      = "Standard errors clustered at municipality level in parentheses. R² reported for specification (1); within R² for specifications (2)–(4). *** p&lt;0.001, ** p&lt;0.01, * p&lt;0.05."
)

# ── Build all table HTML blocks ───────────────────────────────

# Section 1: Main quarterly results
html_t1 <- reg_table_html(
  models    = list(ols_q_main, ols_q_peak, ols_q_employ),
  col_names = c("Per-capita Income<br>Full Year", "Per-capita Income<br>Peak Season (Q2–Q3)", "Employment Rate<br>Full Year"),
  vars      = sat_var,
  fe_line   = "Municipality + Year × Quarter",
  title     = "Table 1. Main Results: Effect of Sargassum Exposure on Household Welfare",
  note      = note_q
)

# Section 2: Quarterly tertile heterogeneity
html_t2 <- reg_table_html(
  models    = list(ols_q_t1, ols_q_t2, ols_q_t3),
  col_names = c("Bottom Third", "Middle Third", "Top Third"),
  vars      = sat_var,
  fe_line   = "Municipality + Year × Quarter",
  title     = "Table 2. Heterogeneity by Income Tertile — Full Year",
  note      = paste(note_q, "Dependent variable is mean log per-capita income for households in each tertile within the municipality-quarter. Tertiles assigned annually from the national income distribution.")
)

html_t3 <- reg_table_html(
  models    = list(ols_q_peak_t1, ols_q_peak_t2, ols_q_peak_t3),
  col_names = c("Bottom Third", "Middle Third", "Top Third"),
  vars      = sat_var,
  fe_line   = "Municipality + Year × Quarter",
  title     = "Table 3. Heterogeneity by Income Tertile — Peak Season (Q2–Q3, April–September)",
  note      = note_q
)

html_t3b <- reg_table_html(
  models    = list(ols_q_pk_t1_nz, ols_q_pk_t2_nz, ols_q_pk_t3_nz),
  col_names = c("Bottom Third", "Middle Third", "Top Third"),
  vars      = sat_var,
  fe_line   = "Municipality + Year × Quarter",
  title     = "Table 3b. Heterogeneity by Income Tertile — Peak Season, Positive-Income Households Only",
  note      = paste(note_q, "Dependent variable is mean log per-capita income averaged over households with strictly positive reported income (ingreso_pc &gt; 0) within each municipality-quarter-tertile cell. Tertile assignments use the full distribution (including zero-income households).")
)

# Section 3: Quarterly coastal proximity
html_t4 <- reg_table_html(
  models    = list(ols_q_coastal, ols_q_not_coastal, ols_q_pk_coastal, ols_q_pk_nc),
  col_names = c("Coastal<br>Full Year", "Non-Coastal<br>Full Year",
                "Coastal<br>Peak Season", "Non-Coastal<br>Peak Season"),
  vars      = sat_var,
  fe_line   = "Municipality + Year × Quarter",
  title     = "Table 4. Heterogeneity by Coastal Proximity",
  note      = paste(note_q, "Coastal: municipalities with a shoreline. Non-coastal: inland municipalities assigned Sargassum exposure via inverse-distance weighting from coastal neighbours.")
)

html_t5 <- reg_table_html(
  models    = list(ols_q_coastal_int, ols_q_pk_coast_int),
  col_names = c("Full Year", "Peak Season (Q2–Q3)"),
  vars      = list("Sargassum intensity (AFAI)"                    = "afai_sargassum",
                   "&nbsp;&nbsp;&nbsp;× Non-coastal (interaction)" = "afai_x_nc"),
  fe_line   = "Municipality + Year × Quarter",
  title     = "Table 5. Interaction Test: Differential Sargassum Effect by Coastal Proximity",
  note      = paste(note_q, "Base group = coastal municipalities.")
)

# Section 4: Monthly household-level robustness
html_t6 <- reg_table_html(
  models    = list(ols_income_pc, ols_employ, ols_peak),
  col_names = c("Per-capita Income<br>Full Year", "Employment Rate<br>Full Year", "Per-capita Income<br>Peak Season (May–Sep)"),
  vars      = sat_var,
  title     = "Table 6. Robustness: Household-Level Monthly Panel",
  note      = paste(note_fe, "Unit of observation is the household-month. Income averaged within municipality-month cells; treatment (AFAI) varies at municipality-month level. Larger standard errors reflect idiosyncratic household income noise.")
)

# ── CSS ───────────────────────────────────────────────────────
css <- '
<style>
  body {
    font-family: "Times New Roman", Times, serif;
    font-size: 16px;
    line-height: 1.55;
    color: #111;
    max-width: 1000px;
    margin: 40px auto;
    padding: 0 32px 60px;
  }
  h1  { font-size: 24px; margin-bottom: 4px; }
  .doc-subtitle { font-size: 15px; color: #444; margin-top: 0; margin-bottom: 28px; }
  h2  { font-size: 18px; margin-top: 48px; margin-bottom: 8px;
        border-bottom: 1.5px solid #555; padding-bottom: 3px; }
  h3.table-title {
    font-size: 16px; font-weight: bold;
    margin-top: 36px; margin-bottom: 4px;
  }
  table.results-table {
    border-collapse: collapse;
    width: 100%;
    font-size: 15px;
    margin-top: 4px;
  }
  table.results-table thead tr th {
    border-top: 2px solid #111;
    border-bottom: 1px solid #111;
    padding: 6px 12px;
    text-align: center;
    font-weight: bold;
    background: #fff;
  }
  table.results-table thead tr th:first-child { text-align: left; }
  table.results-table td { padding: 4px 12px; text-align: center; }
  table.results-table td:first-child { text-align: left; }
  tr.dep-row td  { border-bottom: 1px solid #bbb; font-style: italic; padding-bottom: 5px; }
  tr.coef-row td { padding-top: 5px; font-weight: bold; }
  tr.se-row td   { padding-bottom: 3px; color: #333; font-size: 13px; }
  tr.fe-header td {
    font-variant: small-caps;
    letter-spacing: 0.04em;
    color: #555;
    font-size: 13px;
    padding-top: 8px;
    border-top: 1px solid #bbb;
  }
  tr.fe-row td   { color: #333; }
  tr.stat-divider td { border-top: 1px solid #bbb; padding: 0; height: 4px; }
  tr.stat-row:last-child td { border-bottom: 2px solid #111; padding-bottom: 6px; }
  p.table-note {
    font-size: 13px; color: #444; margin-top: 5px; margin-bottom: 32px;
  }
  hr.section-break {
    border: none; border-top: 1px dashed #ccc; margin: 48px 0 0;
  }
</style>
'

# ── Assemble full HTML document ───────────────────────────────
body_html <- paste0('
<h1>Sargassum and Household Welfare &mdash; Dominican Republic</h1>
<p class="doc-subtitle">
  Regression Results | Two-Way Fixed Effects (Municipality + Year &times; Quarter)<br>
  Survey: ENCFT 2017&ndash;2025 &nbsp;|&nbsp; Satellite: NOAA AOML AFAI 7-day composite<br>
  Treatment variable: mean AFAI for the municipality&ndash;quarter (afai_sargassum).<br>
  Pepillo Salcedo excluded from coastal data (Río Masacre river plume contamination);<br>
  assigned inverse-distance-weighted values from coastal neighbours in all regressions.
</p>

<h2>Section 1 &mdash; Main Results (Quarterly Panel)</h2>
', html_fe_prog, html_t1, '

<hr class="section-break">
<h2>Section 2 &mdash; Heterogeneity by Income Tertile</h2>
', html_t2, html_t3, html_t3b, '

<hr class="section-break">
<h2>Section 3 &mdash; Heterogeneity by Coastal Proximity</h2>
', html_t4, html_t5, '

<hr class="section-break">
<h2>Section 4 &mdash; Robustness: Household-Level Monthly Panel</h2>
', html_t6
)

full_html <- paste0(
  '<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">',
  '<title>Sargassum DR &mdash; Results Tables</title>',
  css,
  '</head><body>',
  body_html,
  '</body></html>'
)

out_html <- file.path(path_results, "results_tables.html")
writeLines(full_html, out_html)
cat("Saved:", out_html, "\n")

desktop_copy <- path.expand("~/Desktop/sargassum_results_tables.html")
file.copy(out_html, desktop_copy, overwrite=TRUE)
cat("Copied to Desktop.\n")
system(paste("open -a Safari", shQuote(desktop_copy)))
cat("Opened in Safari.\n")
