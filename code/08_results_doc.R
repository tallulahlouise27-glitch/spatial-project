# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 8: Export all results tables to a single HTML document
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(fixest)

path_processed <- "data/processed/"
path_results   <- "results/"
dir.create(path_results, showWarnings = FALSE)

# ── Load quarterly panel ──────────────────────────────────────
panel_quarterly <- readRDS(file.path(path_processed, "panel_quarterly.rds"))

panel_quarterly <- panel_quarterly %>%
  mutate(
    is_not_coastal    = as.integer(coastal_type == "not_coastal"),
    afai_sargassum_1k = afai_sargassum * 1000,   # rescale to milli-AFAI units
    afai_x_nc_1k      = afai_sargassum_1k * is_not_coastal
  )

pq_peak <- filter(panel_quarterly, quarter %in% 2:3)

# ── Run all models ────────────────────────────────────────────
cat("Running regressions...\n")

# FE progression
ols_q_no_fe   <- feols(log_income ~ afai_sargassum,                        panel_quarterly, cluster=~muni_fe)
ols_q_muni_fe <- feols(log_income ~ afai_sargassum_1k | muni_fe,              panel_quarterly, cluster=~muni_fe)
ols_q_time_fe <- feols(log_income ~ afai_sargassum_1k | year_quarter_fe,      panel_quarterly, cluster=~muni_fe)
ols_q_main    <- feols(log_income ~ afai_sargassum_1k | muni_fe + year_quarter_fe, panel_quarterly, cluster=~muni_fe)

# Main
ols_q_peak    <- feols(log_income  ~ afai_sargassum_1k | muni_fe + year_quarter_fe, pq_peak,         cluster=~muni_fe)
ols_q_employ  <- feols(tasa_empleo ~ afai_sargassum_1k | muni_fe + year_quarter_fe, panel_quarterly, cluster=~muni_fe)

# Income tertiles
ols_q_t1      <- feols(log_income_t1 ~ afai_sargassum_1k | muni_fe + year_quarter_fe, panel_quarterly, cluster=~muni_fe)
ols_q_t2      <- feols(log_income_t2 ~ afai_sargassum_1k | muni_fe + year_quarter_fe, panel_quarterly, cluster=~muni_fe)
ols_q_t3      <- feols(log_income_t3 ~ afai_sargassum_1k | muni_fe + year_quarter_fe, panel_quarterly, cluster=~muni_fe)
ols_q_peak_t1 <- feols(log_income_t1 ~ afai_sargassum_1k | muni_fe + year_quarter_fe, pq_peak,         cluster=~muni_fe)
ols_q_peak_t2 <- feols(log_income_t2 ~ afai_sargassum_1k | muni_fe + year_quarter_fe, pq_peak,         cluster=~muni_fe)
ols_q_peak_t3 <- feols(log_income_t3 ~ afai_sargassum_1k | muni_fe + year_quarter_fe, pq_peak,         cluster=~muni_fe)

# Coastal proximity
ols_q_coastal     <- feols(log_income ~ afai_sargassum_1k | muni_fe + year_quarter_fe, filter(panel_quarterly, coastal_type=="coastal"),     cluster=~muni_fe)
ols_q_not_coastal <- feols(log_income ~ afai_sargassum_1k | muni_fe + year_quarter_fe, filter(panel_quarterly, coastal_type=="not_coastal"), cluster=~muni_fe)
ols_q_pk_coastal  <- feols(log_income ~ afai_sargassum_1k | muni_fe + year_quarter_fe, filter(pq_peak, coastal_type=="coastal"),             cluster=~muni_fe)
ols_q_pk_nc       <- feols(log_income ~ afai_sargassum_1k | muni_fe + year_quarter_fe, filter(pq_peak, coastal_type=="not_coastal"),         cluster=~muni_fe)
ols_q_coastal_int <- feols(log_income ~ afai_sargassum_1k + afai_x_nc_1k | muni_fe + year_quarter_fe, panel_quarterly, cluster=~muni_fe)
ols_q_pk_coast_int<- feols(log_income ~ afai_sargassum_1k + afai_x_nc_1k | muni_fe + year_quarter_fe, pq_peak,         cluster=~muni_fe)

cat("All models estimated.\n")

# ── Helpers ───────────────────────────────────────────────────

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

reg_table_html <- function(models, col_names, vars,
                           title = NULL, note = NULL) {
  coef_rows <- lapply(names(vars), function(label) {
    var <- vars[[label]]
    list(label  = label,
         coef   = sapply(models, fmt_coef, var=var),
         se     = sapply(models, fmt_se,   var=var))
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
      "log_income_t1" = "Log income (T1)",
      "log_income_t2" = "Log income (T2)",
      "log_income_t3" = "Log income (T3)",
      "tasa_empleo"   = "Employment rate",
      dv)
  })

  th_inner <- paste0("<th>", col_names, "</th>", collapse="")
  thead    <- paste0('<thead><tr><th></th>', th_inner, '</tr></thead>')

  rows <- character(0)

  if (length(unique(dep_vars)) > 1)
    rows <- c(rows, paste0('<tr class="dep-row"><td>Dependent variable</td>',
                           paste0('<td>', dep_vars, '</td>', collapse=""), '</tr>'))

  for (cr in coef_rows) {
    rows <- c(rows,
      paste0('<tr class="coef-row"><td>', cr$label, '</td>',
             paste0('<td>', cr$coef, '</td>', collapse=""), '</tr>'),
      paste0('<tr class="se-row"><td></td>',
             paste0('<td>', cr$se, '</td>', collapse=""), '</tr>'))
  }

  rows <- c(rows,
    '<tr class="stat-divider"><td colspan="100%"></td></tr>',
    paste0('<tr class="stat-row"><td>Observations</td>',
           paste0('<td>', n_obs, '</td>', collapse=""), '</tr>'),
    paste0('<tr class="stat-row"><td>Within R²</td>',
           paste0('<td>', wr2, '</td>', collapse=""), '</tr>'))

  tbody <- paste0('<tbody>', paste(rows, collapse="\n"), '</tbody>')
  tbl   <- paste0('<table class="results-table">', thead, tbody, '</table>')

  out <- character(0)
  if (!is.null(title)) out <- c(out, paste0('<h3 class="table-title">', title, '</h3>'))
  out <- c(out, tbl)
  if (!is.null(note))  out <- c(out, paste0('<p class="table-note">', note, '</p>'))
  paste(out, collapse="\n")
}

note_q   <- "<em>Notes:</em> Municipality and year &times; quarter fixed effects. Standard errors clustered at municipality level in parentheses. Treatment variable is AFAI &times; 10&sup3; (milli-AFAI); a one-unit change corresponds to a 0.001 increase in raw AFAI. *** p&lt;0.001, ** p&lt;0.01, * p&lt;0.05."
note_int <- "<em>Notes:</em> Municipality and year &times; quarter fixed effects throughout. Treatment variable is AFAI &times; 10&sup3;. Interaction coefficient shows the difference in the Sargassum effect for non-coastal relative to coastal municipalities. SEs clustered at municipality level. *** p&lt;0.001, ** p&lt;0.01, * p&lt;0.05."
sat_var  <- list("Sargassum intensity (AFAI × 10³)" = "afai_sargassum_1k")

# ── Fixed-effects progression table ──────────────────────────
fe_prog_html <- function(models, col_names, var='afai_sargassum_1k', title=NULL, note=NULL) {
  has_muni <- sapply(models, function(m) 'muni_fe'        %in% m$fixef_vars)
  has_time <- sapply(models, function(m) 'year_quarter_fe' %in% m$fixef_vars)

  n_obs <- sapply(models, function(m) formatC(nobs(m), format='d', big.mark=','))
  r2_val <- sapply(models, function(m) {
    wr <- r2(m)['wr2']
    if (is.na(wr)) formatC(r2(m)['r2'], format='f', digits=4)
    else           formatC(wr,           format='f', digits=4)
  })
  r2_row_label <- if (all(!is.na(sapply(models, function(m) r2(m)['wr2'])))) 'Within R²' else 'R² / Within R²'

  th_inner <- paste0('<th>', col_names, '</th>', collapse='')
  thead    <- paste0('<thead><tr><th></th>', th_inner, '</tr></thead>')

  rows <- c(
    paste0('<tr class="coef-row"><td>Sargassum intensity (AFAI &times; 10&sup3;)</td>',
           paste0('<td>', sapply(models, fmt_coef, var=var), '</td>', collapse=''), '</tr>'),
    paste0('<tr class="se-row"><td></td>',
           paste0('<td>', sapply(models, fmt_se, var=var), '</td>', collapse=''), '</tr>'),
    '<tr class="stat-divider"><td colspan="100%"></td></tr>',
    paste0('<tr class="fe-row"><td>Municipality FE</td>',
           paste0('<td>', ifelse(has_muni, 'Yes', 'No'), '</td>', collapse=''), '</tr>'),
    paste0('<tr class="fe-row"><td>Year &times; Quarter FE</td>',
           paste0('<td>', ifelse(has_time, 'Yes', 'No'), '</td>', collapse=''), '</tr>'),
    '<tr class="stat-divider"><td colspan="100%"></td></tr>',
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

html_fe_prog <- fe_prog_html(
  models    = list(ols_q_no_fe, ols_q_muni_fe, ols_q_time_fe, ols_q_main),
  col_names = c("No FE", "Municipality FE", "Year &times; Quarter FE", "Both FE"),
  title     = "Table 0. Fixed Effects Progression: Sargassum Effect on Log Household Income",
  note      = "Standard errors clustered at municipality level in parentheses. R&sup2; reported for specification (1); within R&sup2; for (2)&ndash;(4). *** p&lt;0.001, ** p&lt;0.01, * p&lt;0.05."
)

# ── Build tables ──────────────────────────────────────────────

html_t1 <- reg_table_html(
  models    = list(ols_q_main, ols_q_peak, ols_q_employ),
  col_names = c("Income<br>Full Year", "Income<br>Peak Season (Q2&ndash;Q3)", "Employment Rate<br>Full Year"),
  vars      = sat_var,
  title     = "Table 1. Main Results: Effect of Sargassum Exposure on Household Welfare",
  note      = note_q
)

html_t2 <- reg_table_html(
  models    = list(ols_q_t1, ols_q_t2, ols_q_t3),
  col_names = c("Bottom Third", "Middle Third", "Top Third"),
  vars      = sat_var,
  title     = "Table 2. Heterogeneity by Income Tertile &mdash; Full Year",
  note      = paste(note_q, "Dependent variable is mean log household income for households in each tertile within the municipality-quarter cell. Tertiles assigned annually from the national per-capita income distribution.")
)

html_t3 <- reg_table_html(
  models    = list(ols_q_peak_t1, ols_q_peak_t2, ols_q_peak_t3),
  col_names = c("Bottom Third", "Middle Third", "Top Third"),
  vars      = sat_var,
  title     = "Table 3. Heterogeneity by Income Tertile &mdash; Peak Season (Q2&ndash;Q3)",
  note      = note_q
)

html_t4 <- reg_table_html(
  models    = list(ols_q_coastal, ols_q_not_coastal, ols_q_pk_coastal, ols_q_pk_nc),
  col_names = c("Coastal<br>Full Year", "Non-Coastal<br>Full Year",
                "Coastal<br>Peak Season", "Non-Coastal<br>Peak Season"),
  vars      = sat_var,
  title     = "Table 4. Heterogeneity by Coastal Proximity",
  note      = paste(note_q, "Coastal: municipalities with a shoreline. Non-coastal: inland municipalities assigned Sargassum exposure via inverse-distance weighting from coastal neighbours.")
)

html_t5 <- reg_table_html(
  models    = list(ols_q_coastal_int, ols_q_pk_coast_int),
  col_names = c("Full Year", "Peak Season (Q2&ndash;Q3)"),
  vars      = list("Sargassum intensity (AFAI)"                    = "afai_sargassum",
                   "&nbsp;&nbsp;&nbsp;&times; Non-coastal"         = "afai_x_nc_1k"),
  title     = "Table 5. Interaction Test: Coastal vs Non-Coastal Sargassum Effect",
  note      = paste(note_int, "Base group = coastal municipalities.")
)

# ── CSS ───────────────────────────────────────────────────────
css <- '
<style>
  body {
    font-family: "Times New Roman", Times, serif;
    font-size: 19px;
    line-height: 1.55;
    color: #111;
    max-width: 1000px;
    margin: 40px auto;
    padding: 0 32px 60px;
  }
  h1  { font-size: 28px; margin-bottom: 4px; }
  .doc-subtitle { font-size: 18px; color: #444; margin-top: 0; margin-bottom: 28px; }
  h2  { font-size: 22px; margin-top: 48px; margin-bottom: 8px;
        border-bottom: 1.5px solid #555; padding-bottom: 3px; }
  h3.table-title {
    font-size: 19px; font-weight: bold;
    margin-top: 36px; margin-bottom: 4px;
  }
  table.results-table {
    border-collapse: collapse;
    width: 100%;
    font-size: 18px;
    margin-top: 4px;
  }
  table.results-table thead tr th {
    border-top: 2px solid #111;
    border-bottom: 1px solid #111;
    padding: 7px 14px;
    text-align: center;
    font-weight: bold;
    background: #fff;
  }
  table.results-table thead tr th:first-child { text-align: left; }
  table.results-table td { padding: 5px 14px; text-align: center; }
  table.results-table td:first-child { text-align: left; }
  tr.dep-row td  { border-bottom: 1px solid #bbb; font-style: italic; padding-bottom: 5px; }
  tr.coef-row td { padding-top: 5px; font-weight: bold; }
  tr.se-row td   { padding-bottom: 3px; color: #333; font-size: 16px; }
  tr.fe-row td   { color: #333; }
  tr.stat-divider td { border-top: 1px solid #bbb; padding: 0; height: 4px; }
  tr.stat-row:last-child td { border-bottom: 2px solid #111; padding-bottom: 6px; }
  p.table-note {
    font-size: 16px; color: #444; margin-top: 5px; margin-bottom: 32px;
  }
  hr.section-break {
    border: none; border-top: 1px dashed #ccc; margin: 48px 0 0;
  }
</style>
'

# ── Assemble document ─────────────────────────────────────────
body_html <- paste0('
<h1>Sargassum and Household Welfare &mdash; Dominican Republic</h1>
<p class="doc-subtitle">
  Regression Results | Municipality + Year &times; Quarter Fixed Effects<br>
  Survey: ENCFT 2017&ndash;2025 &nbsp;|&nbsp; Satellite: NOAA AOML AFAI 7-day composite<br>
  Unit of observation: municipality &times; quarter (N = 3,095).<br>
  Pepillo Salcedo excluded from coastal satellite data (Río Masacre river plume);<br>
  assigned IDW-weighted values from coastal neighbours in all regressions.
</p>

<h2>Section 1 &mdash; Main Results</h2>
', html_fe_prog, html_t1, '

<hr class="section-break">
<h2>Section 2 &mdash; Heterogeneity by Income Tertile</h2>
', html_t2, html_t3, '

<hr class="section-break">
<h2>Section 3 &mdash; Heterogeneity by Coastal Proximity</h2>
', html_t4, html_t5
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
