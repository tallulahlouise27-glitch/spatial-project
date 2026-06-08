# ============================================================
# Sargassum & Household Welfare — Dominican Republic
# Script 1: Setup and package loading
# ============================================================

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Install packages if not already installed
packages <- c(
  "tidyverse",   # data cleaning and manipulation
  "sf",          # spatial data (shapefiles etc.)
  "terra",       # raster/satellite data
  "rerddap",     # download satellite data from ERDDAP servers
  "geodata",     # download GADM administrative boundaries
  "plm",         # panel data regression
  "fixest",      # fast fixed effects regression
  "AER",         # instrumental variable regression (ivreg)
  "stargazer"    # export regression tables
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================
# File paths — update these if your folder moves
# ============================================================

# Path to your raw ENCFT data
path_encft <- "data/raw/"

# Path to AFAI satellite data (NOAA AOML / USF Optical Oceanography Lab)
path_satellite <- "data/raw/satellite/"

# Where to save cleaned data
path_processed <- "data/processed/"

# Where to save results and figures
path_results <- "results/"
path_figures <- "figures/"
