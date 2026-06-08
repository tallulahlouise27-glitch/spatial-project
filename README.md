# Sargassum & Household Welfare — Dominican Republic

**Research question:** What is the impact of Sargassum algal blooms on household economic conditions in the Dominican Republic?

## Data sources
- **ENCFT** (Encuesta Nacional de Condiciones de Vida y Trabajo): Dominican Republic household survey, 2016–2026. Stored locally in `data/raw/` and backed up on OneDrive.
- **NASA Earthdata**: Satellite data on Sargassum coastal exposure. To be downloaded and added to `data/raw/satellite/`.

## Method
Panel regression with household and year fixed effects, using an instrumental variable (IV) for Sargassum coastal exposure.

## Folder structure
```
spatial-project/
├── code/          R scripts (numbered in order)
├── data/
│   ├── raw/       Original data — not pushed to GitHub
│   └── processed/ Cleaned data — not pushed to GitHub
├── figures/       Output plots
└── results/       Regression output tables
```

## Scripts
| Script | Purpose |
|--------|---------|
| `code/01_setup.R` | Load packages and set file paths |
| `code/02_encft.R` | Load and clean ENCFT survey data |
| `code/03_satellite.R` | Download and process NASA Sargassum data |
| `code/04_merge.R` | Merge survey and satellite data |
| `code/05_regression.R` | Panel IV regression and results |
