# =============================================================================
# dependencies.R — Install all packages required by the MorPhiC Bulk RNA-seq app
#
# Usage:  Rscript dependencies.R
# =============================================================================

install.packages(c(
  "shiny",
  "bslib",
  "bsicons",
  "DT",
  "plotly",
  "duckdb",
  "DBI"
))
