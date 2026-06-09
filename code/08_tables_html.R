# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 8: Publication-quality HTML tables for Word / Google Docs
#
# Output: results/tables_for_paper.html
# Usage:  open the HTML file in a browser, select a table,
#         copy and paste directly into Word or Google Docs.
# ============================================================

options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(tidyverse)
library(fixest)

path_processed <- "data/processed/"
path_results   <- "results/"

# ── Load panel and rebuild all models ────────────────────────
panel <- readRDS(file.path(path_processed, "panel_monthly.rds"))

panel <- panel %>%
  mutate(
    is_not_coastal     = as.integer(coastal_type == "not_coastal"),
    afai_x_not_coastal = afai_sargassum * is_not_coastal
  )

# Non-zero tertile income variables
hogar_raw <- readRDS(file.path(path_processed, "encft_hogar.rds"))
hogar_nz  <- hogar_raw %>% filter(ingreso_pc > 0)

cuts_nz <- hogar_nz %>%
  group_by(ANO) %>%
  summarise(
    t1_cut_nz = quantile(ingreso_pc, 1/3, na.rm = TRUE),
    t2_cut_nz = quantile(ingreso_pc, 2/3, na.rm = TRUE),
    .groups = "drop"
  )

nz_q <- hogar_nz %>%
  left_join(cuts_nz, by = "ANO") %>%
  mutate(
    tertile_nz = case_when(
      ingreso_pc <= t1_cut_nz ~ "T1",
      ingreso_pc <= t2_cut_nz ~ "T2",
      TRUE                    ~ "T3"
    ),
    quarter = as.integer(TRIMESTRE) %% 10L
  ) %>%
  group_by(ID_MUNICIPIO, ANO, quarter, tertile_nz) %>%
  summarise(
    ingreso_nz = weighted.mean(ingreso_hogar,
                               w = as.numeric(FACTOR_EXPANSION), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = tertile_nz,
    values_from = ingreso_nz,
    names_glue  = "ingreso_{tertile_nz}_nz"
  )

panel <- panel %>%
  left_join(nz_q, by = c("ID_MUNICIPIO", "year" = "ANO", "quarter")) %>%
  mutate(
    log_income_t1_nz = log(ingreso_T1_nz + 1),
    log_income_t2_nz = log(ingreso_T2_nz + 1),
    log_income_t3_nz = log(ingreso_T3_nz + 1)
  )

panel_peak <- filter(panel, month %in% 5:9)

cat("Running models...\n")

ols_income  <- feols(log_income    ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster = ~muni_fe)
ols_employ  <- feols(tasa_empleo   ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster = ~muni_fe)
ols_peak    <- feols(log_income    ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)

ols_t1 <- feols(log_income_t1 ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster = ~muni_fe)
ols_t2 <- feols(log_income_t2 ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster = ~muni_fe)
ols_t3 <- feols(log_income_t3 ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster = ~muni_fe)

ols_t1_nz      <- feols(log_income_t1_nz ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster = ~muni_fe)
ols_t2_nz      <- feols(log_income_t2_nz ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster = ~muni_fe)
ols_t3_nz      <- feols(log_income_t3_nz ~ afai_sargassum | muni_fe + year_month_fe, panel,      cluster = ~muni_fe)

ols_peak_t1 <- feols(log_income_t1 ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)
ols_peak_t2 <- feols(log_income_t2 ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)
ols_peak_t3 <- feols(log_income_t3 ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)

ols_peak_t1_nz <- feols(log_income_t1_nz ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)
ols_peak_t2_nz <- feols(log_income_t2_nz ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)
ols_peak_t3_nz <- feols(log_income_t3_nz ~ afai_sargassum | muni_fe + year_month_fe, panel_peak, cluster = ~muni_fe)

ols_coastal          <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, filter(panel,      coastal_type == "coastal"),     cluster = ~muni_fe)
ols_not_coastal      <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, filter(panel,      coastal_type == "not_coastal"), cluster = ~muni_fe)
ols_peak_coastal     <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, filter(panel_peak, coastal_type == "coastal"),     cluster = ~muni_fe)
ols_peak_not_coastal <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, filter(panel_peak, coastal_type == "not_coastal"), cluster = ~muni_fe)

# Fixed effects progression models (for appendix)
fe_none     <- feols(log_income ~ afai_sargassum,                          panel, cluster = ~muni_fe)
fe_muni     <- feols(log_income ~ afai_sargassum | muni_fe,                panel, cluster = ~muni_fe)
fe_time     <- feols(log_income ~ afai_sargassum | year_month_fe,          panel, cluster = ~muni_fe)
fe_both     <- feols(log_income ~ afai_sargassum | muni_fe + year_month_fe, panel, cluster = ~muni_fe)

cat("Models done. Building HTML tables...\n")

# ── Helper functions ──────────────────────────────────────────

stars <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}

# Extract one cell: coefficient with stars, SE below in parentheses
cell <- function(model, var = "afai_sargassum") {
  b <- coef(model)[var]
  s <- se(model)[var]
  p <- pvalue(model)[var]
  if (is.na(b)) return(list(coef = "&mdash;", se = ""))
  st <- stars(p)
  list(
    coef = sprintf("%.3f%s", b, st),
    se   = sprintf("(%.3f)", s)
  )
}

fmt_n   <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_r2  <- function(x) sprintf("%.3f", x)

# Build one full HTML regression table
# models: named list of feols objects
# col_labels: character vector, one per model
# depvar_labels: character vector (dependent variable names), one per model
# title: table title string
# notes: footer note string
make_reg_table <- function(models, col_labels, depvar_labels, title, notes,
                           var_label = "Sargassum intensity",
                           panel_headers = NULL) {
  ncols <- length(models)
  col_nos <- paste0("(", seq_len(ncols), ")")

  # Build rows
  coef_row <- sapply(models, cell)
  coefs    <- sapply(coef_row["coef", ], identity)
  ses      <- sapply(coef_row["se",   ], identity)

  obs_row  <- sapply(models, function(m) fmt_n(nobs(m)))
  r2_row   <- sapply(models, function(m) fmt_r2(r2(m, type = "r2")))

  td <- function(x, cls = "num") sprintf('<td class="%s">%s</td>', cls, x)
  th <- function(x, cls = "hdr") sprintf('<th class="%s">%s</th>', cls, x)

  header_labs  <- paste(sapply(col_labels, function(x) th(x, "hdr")), collapse = "\n          ")

  # Build combined header block: panel spanning row (if any) + column labels
  if (!is.null(panel_headers)) {
    spans <- paste(sapply(panel_headers, function(h) {
      sprintf('<th colspan="%d" class="panel-hdr">%s</th>', h$span, h$label)
    }), collapse = "\n          ")
    header_block <- sprintf(
      '<tr class="topline panel-hdr-row">
      <th class="rowlabel"></th>
      %s
    </tr>
    <tr>
      <th class="rowlabel"></th>
      %s
    </tr>', spans, header_labs)
  } else {
    header_block <- sprintf(
      '<tr class="topline">
      <th class="rowlabel"></th>
      %s
    </tr>', header_labs)
  }

  header_depvars <- paste(sapply(depvar_labels, function(x) td(x, "depvar")), collapse = "\n          ")

  coef_cells  <- paste(sapply(coefs, td), collapse = "\n          ")
  se_cells    <- paste(sapply(ses,   function(x) td(x, "se")), collapse = "\n          ")
  obs_cells   <- paste(sapply(obs_row, td), collapse = "\n          ")
  r2_cells    <- paste(sapply(r2_row,  td), collapse = "\n          ")

  sprintf('
<div class="tbl-wrap">
<p class="tbl-title">%s</p>
<table>
  <colgroup>
    <col class="rowlabel">%s
  </colgroup>
  <thead>
    %s
    <tr class="midline">
      <td colspan="%d"></td>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td class="rowlabel">%s</td>
      %s
    </tr>
    <tr>
      <td class="rowlabel"></td>
      %s
    </tr>
    <tr class="stat-sep">
      <td class="rowlabel">Observations</td>
      %s
    </tr>
    <tr>
      <td class="rowlabel">R&sup2;</td>
      %s
    </tr>
  </tbody>
  <tfoot>
    <tr class="bottomline"><td colspan="%d"></td></tr>
    <tr>
      <td colspan="%d" class="notes">%s</td>
    </tr>
  </tfoot>
</table>
</div>
',
    title,
    paste(rep('<col class="datacol">', ncols), collapse = ""),
    header_block,
    ncols + 1,
    var_label,
    coef_cells,
    se_cells,
    obs_cells,
    r2_cells,
    ncols + 1,
    ncols + 1,
    notes
  )
}

# ── CSS and page shell ────────────────────────────────────────
css <- '
<style>
  body {
    font-family: "Times New Roman", Times, serif;
    font-size: 12pt;
    color: #000;
    max-width: 800px;
    margin: 40px auto;
    padding: 0 20px;
    line-height: 1.4;
  }
  h1 { font-size: 14pt; font-weight: bold; margin-top: 0; }
  h2 { font-size: 12pt; font-weight: normal; font-style: italic; margin-bottom: 4px; }
  p.intro { font-size: 11pt; margin-bottom: 24px; color: #333; }

  .tbl-wrap { margin: 32px 0 48px 0; page-break-inside: avoid; }
  .tbl-title {
    font-size: 12pt;
    font-weight: bold;
    margin-bottom: 6px;
  }

  table {
    border-collapse: collapse;
    width: 100%;
    font-size: 11pt;
  }
  col.rowlabel { width: 200px; }
  col.datacol  { }

  th, td { padding: 3px 10px; text-align: center; border: none; }
  th.rowlabel, td.rowlabel { text-align: left; padding-left: 0; }
  th.hdr { font-weight: normal; font-style: italic; }

  td.depvar { font-style: italic; font-size: 10pt; color: #333; }
  td.se { font-size: 10.5pt; color: #222; }
  td.num { }
  td.notes {
    text-align: left;
    font-size: 10pt;
    padding-top: 6px;
    padding-left: 0;
    color: #333;
  }

  tr.topline td, tr.topline th { border-top: 2px solid #000; padding-top: 6px; }
  tr.midline td { border-top: 1px solid #000; padding: 2px 0; height: 4px; }
  tr.bottomline td { border-top: 1.5px solid #000; padding: 0; height: 2px; }
  tr.fe-sep td { border-top: 1px solid #ccc; padding-top: 5px; }
  tr.stat-sep td { border-top: 1px solid #ccc; padding-top: 5px; }
  th.panel-hdr {
    font-weight: normal;
    font-style: italic;
    border-bottom: 1px solid #000;
    padding-bottom: 2px;
  }
  tr.panel-hdr-row th.rowlabel { border-bottom: none; }
  tr.depvar-row td { font-size: 10pt; padding-bottom: 2px; }

  hr.page-break { border: none; border-top: 1px dashed #ccc; margin: 60px 0; }
  .instructions {
    background: #f5f5f5;
    border-left: 3px solid #aaa;
    padding: 10px 14px;
    font-size: 10.5pt;
    margin-bottom: 32px;
    font-family: Arial, sans-serif;
  }
</style>
'

notes_main <- paste0(
  "<em>Notes:</em> All models estimated by OLS with municipality and year&times;month fixed effects. ",
  "Standard errors clustered by municipality in parentheses. ",
  "Treatment variable is mean excess AFAI above the 0.001 threshold (Sargassum intensity). ",
  "Peak season = May&ndash;September. ",
  "Within R&sup2; measures variation explained after removing fixed effects. ",
  "* p&lt;0.05, ** p&lt;0.01, *** p&lt;0.001."
)

notes_tertile <- paste0(
  "<em>Notes:</em> All models estimated by OLS with municipality and year&times;month fixed effects. ",
  "Standard errors clustered by municipality in parentheses. ",
  "Tertiles defined nationally by year using per-capita household income. ",
  "T1 = bottom third, T2 = middle third, T3 = top third. ",
  "Dependent variable is log mean income of each tertile group within the municipality&ndash;month cell. ",
  "* p&lt;0.05, ** p&lt;0.01, *** p&lt;0.001."
)

notes_tertile_nz <- paste0(
  "<em>Notes:</em> All models estimated by OLS with municipality and year&times;month fixed effects. ",
  "Standard errors clustered by municipality in parentheses. ",
  "Tertiles defined nationally by year after <strong>excluding zero-income households</strong> (~17.6% of the sample). ",
  "T1 = bottom third, T2 = middle third, T3 = top third of positive-income households. ",
  "* p&lt;0.05, ** p&lt;0.01, *** p&lt;0.001."
)

notes_coastal <- paste0(
  "<em>Notes:</em> All models estimated by OLS with municipality and year&times;month fixed effects. ",
  "Standard errors clustered by municipality in parentheses. ",
  "Coastal municipalities have direct ocean shoreline contact (48 municipalities); ",
  "not-coastal municipalities are all others (52 municipalities). ",
  "The Haiti land border is excluded from the coastline classification. ",
  "* p&lt;0.05, ** p&lt;0.01, *** p&lt;0.001."
)

# ── Build each table ──────────────────────────────────────────

tbl1 <- make_reg_table(
  models        = list(ols_income, ols_employ, ols_peak),
  col_labels    = c("Log Income", "Employment Rate", "Log Income (Peak)"),
  depvar_labels = c("log(income)", "employment rate", "log(income) — peak"),
  title         = "Table 1. Effect of Sargassum Exposure on Household Welfare",
  notes         = notes_main
)

tbl2 <- make_reg_table(
  models      = list(ols_t1, ols_t2, ols_t3,
                     ols_peak_t1, ols_peak_t2, ols_peak_t3),
  col_labels  = c("Bottom", "Middle", "Top",
                  "Bottom", "Middle", "Top"),
  depvar_labels = c("log income", "log income", "log income",
                    "log income", "log income", "log income"),
  title         = "Table 2. Sargassum Effect by Income Tertile (Full Year and Peak Season)",
  notes         = notes_tertile,
  panel_headers = list(
    list(label = "Full Year",             span = 3),
    list(label = "Peak Season (May–Sep)", span = 3)
  )
)

tbl3 <- make_reg_table(
  models        = list(ols_coastal, ols_not_coastal,
                       ols_peak_coastal, ols_peak_not_coastal),
  col_labels    = c("Coastal", "Not Coastal", "Coastal (Peak)", "Not Coastal (Peak)"),
  depvar_labels = c("log income", "log income", "log income", "log income"),
  title         = "Table 3. Sargassum Effect by Coastal Proximity (Full Year and Peak Season)",
  notes         = notes_coastal
)

tbl4 <- make_reg_table(
  models      = list(ols_t1_nz, ols_t2_nz, ols_t3_nz,
                     ols_peak_t1_nz, ols_peak_t2_nz, ols_peak_t3_nz),
  col_labels  = c("Bottom", "Middle", "Top",
                  "Bottom", "Middle", "Top"),
  depvar_labels = c("log income", "log income", "log income",
                    "log income", "log income", "log income"),
  title         = "Table 4. Sargassum Effect by Income Tertile — Positive-Income Households Only",
  notes         = notes_tertile_nz,
  panel_headers = list(
    list(label = "Full Year",             span = 3),
    list(label = "Peak Season (May–Sep)", span = 3)
  )
)

# ── Tertile cutoffs table (non-regression) ────────────────────
hogar    <- readRDS(file.path(path_processed, "encft_hogar.rds"))
cuts_all <- hogar %>%
  group_by(ANO) %>%
  summarise(t1 = quantile(ingreso_pc, 1/3, na.rm = TRUE),
            t2 = quantile(ingreso_pc, 2/3, na.rm = TRUE), .groups = "drop")

cuts_nz2 <- hogar %>%
  filter(ingreso_pc > 0) %>%
  group_by(ANO) %>%
  summarise(t1_nz = quantile(ingreso_pc, 1/3, na.rm = TRUE),
            t2_nz = quantile(ingreso_pc, 2/3, na.rm = TRUE), .groups = "drop")

tbl_cuts <- cuts_all %>%
  left_join(cuts_nz2, by = "ANO") %>%
  filter(ANO >= 2017, ANO <= 2025)

cuts_rows <- paste(apply(tbl_cuts, 1, function(r) {
  sprintf(
    '<tr><td class="rowlabel">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td></tr>',
    r["ANO"],
    formatC(as.numeric(r["t1"]),    format="f", digits=0, big.mark=","),
    formatC(as.numeric(r["t2"]),    format="f", digits=0, big.mark=","),
    formatC(as.numeric(r["t1_nz"]), format="f", digits=0, big.mark=","),
    formatC(as.numeric(r["t2_nz"]), format="f", digits=0, big.mark=",")
  )
}), collapse = "\n")

tbl_appendix <- sprintf('
<div class="tbl-wrap">
<p class="tbl-title">Table A1. National Income Tertile Cutoffs by Year</p>
<table>
  <colgroup>
    <col class="rowlabel">
    <col class="datacol"><col class="datacol">
    <col class="datacol"><col class="datacol">
  </colgroup>
  <thead>
    <tr class="topline">
      <th class="rowlabel"></th>
      <th colspan="2" class="hdr">All households</th>
      <th colspan="2" class="hdr">Positive-income households only</th>
    </tr>
    <tr>
      <th class="rowlabel">Year</th>
      <th class="hdr">Bottom/Middle boundary</th>
      <th class="hdr">Middle/Top boundary</th>
      <th class="hdr">Bottom/Middle boundary</th>
      <th class="hdr">Middle/Top boundary</th>
    </tr>
    <tr class="midline"><td colspan="5"></td></tr>
  </thead>
  <tbody>
    %s
  </tbody>
  <tfoot>
    <tr class="bottomline"><td colspan="5"></td></tr>
    <tr>
      <td colspan="5" class="notes">
        <em>Notes:</em> Units are RD$ per capita per month (household total income divided by number of members).
        Cutoffs are national percentiles of the per-capita income distribution within each year.
        Positive-income sample excludes the approximately 17.6%% of surveyed households reporting zero income.
      </td>
    </tr>
  </tfoot>
</table>
</div>
', cuts_rows)

# ── FE progression table (Table A2) ──────────────────────────
fe_models <- list(fe_none, fe_muni, fe_time, fe_both)
fe_coefs  <- sapply(fe_models, function(m) {
  b <- coef(m)["afai_sargassum"]
  s <- se(m)["afai_sargassum"]
  p <- pvalue(m)["afai_sargassum"]
  sprintf("%.3f%s", b, stars(p))
})
fe_ses <- sapply(fe_models, function(m) {
  sprintf("(%.3f)", se(m)["afai_sargassum"])
})
fe_obs <- sapply(fe_models, function(m) fmt_n(nobs(m)))
fe_r2  <- sapply(fe_models, function(m) fmt_r2(r2(m, type = "r2")))

fe_yn <- function(vals) paste(sapply(vals, function(v)
  sprintf('<td class="num">%s</td>', v)), collapse = "\n          ")

tbl_fe_prog <- sprintf('
<div class="tbl-wrap">
<p class="tbl-title">Table A2. Adding Fixed Effects Progressively</p>
<table>
  <colgroup>
    <col class="rowlabel">
    <col class="datacol"><col class="datacol">
    <col class="datacol"><col class="datacol">
  </colgroup>
  <thead>
    <tr class="topline">
      <th class="rowlabel"></th>
      <th class="hdr">No fixed effects</th>
      <th class="hdr">Municipality FE</th>
      <th class="hdr">Year&times;Month FE</th>
      <th class="hdr">Both FEs</th>
    </tr>
    <tr class="midline"><td colspan="5"></td></tr>
  </thead>
  <tbody>
    <tr>
      <td class="rowlabel">Sargassum intensity</td>
      %s
    </tr>
    <tr>
      <td class="rowlabel"></td>
      %s
    </tr>
    <tr class="fe-sep">
      <td class="rowlabel">Municipality FE</td>
      <td class="num">No</td><td class="num">Yes</td>
      <td class="num">No</td><td class="num">Yes</td>
    </tr>
    <tr>
      <td class="rowlabel">Year&times;Month FE</td>
      <td class="num">No</td><td class="num">No</td>
      <td class="num">Yes</td><td class="num">Yes</td>
    </tr>
    <tr class="stat-sep">
      <td class="rowlabel">Observations</td>
      %s
    </tr>
    <tr>
      <td class="rowlabel">R&sup2;</td>
      %s
    </tr>
  </tbody>
  <tfoot>
    <tr class="bottomline"><td colspan="5"></td></tr>
    <tr>
      <td colspan="5" class="notes">
        <em>Notes:</em> Dependent variable is log household income. Each column adds fixed effects
        progressively to the baseline OLS regression. Standard errors clustered by municipality
        in parentheses. The two-way fixed effects specification in column (4) is the main
        specification used throughout the paper.
        * p&lt;0.05, ** p&lt;0.01, *** p&lt;0.001.
      </td>
    </tr>
  </tfoot>
</table>
</div>
',
  fe_yn(fe_coefs),
  fe_yn(fe_ses),
  fe_yn(fe_obs),
  fe_yn(fe_r2)
)

# ── Assemble HTML document ────────────────────────────────────
html <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Regression Tables — Sargassum DR</title>
%s
</head>
<body>
<h1>Sargassum &amp; Household Welfare in the Dominican Republic</h1>
<h2>Regression Tables for Paper</h2>

<div class="instructions">
  <strong>How to use:</strong> Open this file in any web browser (Safari, Chrome, etc.).
  Click just before a table title, drag to select the whole table, then copy (&amp;#8984;C on Mac)
  and paste directly into Word or Google Docs. The table formatting will be preserved.
  To select all tables at once: Edit &rarr; Select All, then paste.
</div>

<p class="intro">
  Panel: 86 municipalities &times; 2017&ndash;2025 &times; monthly (9,089 observations).<br>
  All models include municipality fixed effects and year&times;month fixed effects,
  with standard errors clustered by municipality.
</p>

%s
<hr class="page-break">
<p style="font-style:italic; font-size:11pt; margin-bottom:8px;">
  <strong>Column headers for Tables 2 and 4:</strong>
  Columns (1)&ndash;(3) report full-year results; columns (4)&ndash;(6) report peak-season (May&ndash;Sep) results.
</p>
<p style="font-style:italic; font-size:11pt; margin-bottom:8px;">
  <strong>Column headers for Table 3:</strong>
  Columns (1)&ndash;(2) report full-year results; columns (3)&ndash;(4) report peak-season results.
</p>
<hr class="page-break">

%s

<hr class="page-break">
<p style="font-size:11pt; font-weight:bold;">Appendix</p>
%s

</body>
</html>',
  css,
  paste(tbl1, tbl2, tbl3, sep = "\n"),
  tbl4,
  paste(tbl_fe_prog, tbl_appendix, sep = "\n")
)

out <- file.path(path_results, "tables_for_paper.html")
writeLines(html, out)
cat("Saved:", out, "\n")
cat("Open this file in a browser, then copy-paste each table into Word or Google Docs.\n")
