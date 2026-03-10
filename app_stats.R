# =============================================================================
# app_stats.R
# Reads IMPC viability/phenotype data, OMIM phenotypes, and Orphanet disorder
# data from EGEx-db.duckdb, appends HGNC gene symbols, and writes cleaned
# tables into morphic_bulkRNA3.duckdb for use by the Shiny app.
#
# Tables written to morphic_bulkRNA3 (main schema):
#   impc_viability   — one row per gene:  symbol | GeneID | impc_viability
#   impc_phenotypes  — long format:       symbol | GeneID | impc_phenotypes | impc_zygosity
#   omim_phenotypes  — long format:       symbol | GeneID | omim_phenotype  | <other cols>
#   orphanet_data    — long format:       symbol | GeneID | disorder_name   | <other cols>
#
# Run from: /Users/gabrielm/Desktop/MorPhiC/MorPhiC-Analysis/
# =============================================================================

library(DBI)
library(duckdb)
library(dplyr)

# ── Paths ─────────────────────────────────────────────────────────────────────
EGEX_PATH    <- "./EGEx-db/EGEx-db.duckdb"
MORPHIC_PATH <- "./bulkRNA-seq/bulkRNA-seq-analysis/DuckDB/morphic_bulkRNA3.duckdb"

# ── Connect to EGEx-db (read-only) ────────────────────────────────────────────
message("Connecting to EGEx-db …")
egex <- dbConnect(duckdb(), EGEX_PATH, read_only = TRUE)

cat("Tables available in EGEx-db:\n")
print(dbListTables(egex))

# Gene ID → HGNC symbol lookup (shared by all tables)
genes_tbl <- dbReadTable(egex, "genes") |>
  select(GeneID, symbol) |>
  filter(!is.na(symbol), nzchar(trimws(symbol)))
message(sprintf("  Gene lookup loaded: %d entries", nrow(genes_tbl)))

# =============================================================================
# IMPC Viability
#   Source cols: GeneID + impc_viability
#   After clean: symbol | GeneID | impc_viability
# =============================================================================
message("\nProcessing IMPC_viability …")

impc_viability_raw <- dbReadTable(egex, "IMPC_viability")
cat(sprintf("  Raw rows: %d  |  cols: %s\n",
            nrow(impc_viability_raw), paste(names(impc_viability_raw), collapse = ", ")))

impc_viability <- impc_viability_raw |>
  left_join(genes_tbl, by = "GeneID") |>
  filter(!is.na(symbol), symbol != "WT", !is.na(impc_viability)) |>
  relocate(symbol, .before = everything()) |>
  distinct()

message(sprintf("  Cleaned rows: %d", nrow(impc_viability)))

# =============================================================================
# IMPC Phenotypes (long format — one row per gene × phenotype × zygosity)
#   Source cols: GeneID | impc_phenotypes | impc_zygosity
#   After clean: symbol | GeneID | impc_phenotypes | impc_zygosity
# =============================================================================
message("\nProcessing IMPC_phenotypes …")

impc_phenotypes_raw <- dbReadTable(egex, "IMPC_phenotypes")
cat(sprintf("  Raw rows: %d  |  cols: %s\n",
            nrow(impc_phenotypes_raw), paste(names(impc_phenotypes_raw), collapse = ", ")))

impc_phenotypes <- impc_phenotypes_raw |>
  left_join(genes_tbl, by = "GeneID") |>
  filter(!is.na(symbol), symbol != "WT", !is.na(impc_phenotypes)) |>
  relocate(symbol, .before = everything()) |>
  distinct()

message(sprintf("  Cleaned rows: %d  (%d unique genes)",
                nrow(impc_phenotypes), n_distinct(impc_phenotypes$symbol)))

# =============================================================================
# OMIM Phenotypes (long format — one row per gene × phenotype)
#   Source cols: GeneID | phenotypes | <any other cols>
#   After clean: symbol | GeneID | omim_phenotype | <other cols>
#   Note: "phenotypes" renamed to "omim_phenotype" to avoid ambiguity
# =============================================================================
message("\nProcessing OMIM_phenotypes …")

omim_raw <- dbReadTable(egex, "OMIM_phenotypes")
cat(sprintf("  Raw rows: %d  |  cols: %s\n",
            nrow(omim_raw), paste(names(omim_raw), collapse = ", ")))

omim_phenotypes <- omim_raw |>
  left_join(genes_tbl, by = "GeneID") |>
  filter(!is.na(symbol), !is.na(phenotypes)) |>
  rename(omim_phenotype = phenotypes) |>
  relocate(symbol, .before = everything()) |>
  distinct()

message(sprintf("  Cleaned rows: %d  (%d unique genes)",
                nrow(omim_phenotypes), n_distinct(omim_phenotypes$symbol)))

# =============================================================================
# Orphanet Data (long format — one row per gene × disorder)
#   Source cols: GeneID | disorder_name | <any other cols>
#   After clean: symbol | GeneID | disorder_name | <other cols>
# =============================================================================
message("\nProcessing OrphaData …")

orpha_raw <- dbReadTable(egex, "OrphaData")
cat(sprintf("  Raw rows: %d  |  cols: %s\n",
            nrow(orpha_raw), paste(names(orpha_raw), collapse = ", ")))

orphanet_data <- orpha_raw |>
  left_join(genes_tbl, by = "GeneID") |>
  filter(!is.na(symbol), !is.na(disorder_name)) |>
  relocate(symbol, .before = everything()) |>
  distinct()

message(sprintf("  Cleaned rows: %d  (%d unique genes)",
                nrow(orphanet_data), n_distinct(orphanet_data$symbol)))

# ── Close EGEx-db ─────────────────────────────────────────────────────────────
dbDisconnect(egex, shutdown = TRUE)
message("\nEGEx-db closed.")

# =============================================================================
# Write to morphic_bulkRNA3.duckdb
# =============================================================================
message("\nWriting to morphic_bulkRNA3 …")
morphic <- dbConnect(duckdb(), MORPHIC_PATH)

dbWriteTable(morphic, "impc_viability",  impc_viability,  overwrite = TRUE)
message("  ✓ impc_viability")

dbWriteTable(morphic, "impc_phenotypes", impc_phenotypes, overwrite = TRUE)
message("  ✓ impc_phenotypes")

dbWriteTable(morphic, "omim_phenotypes", omim_phenotypes, overwrite = TRUE)
message("  ✓ omim_phenotypes")

dbWriteTable(morphic, "orphanet_data",   orphanet_data,   overwrite = TRUE)
message("  ✓ orphanet_data")

dbDisconnect(morphic, shutdown = TRUE)
message("morphic_bulkRNA3 closed.")

# =============================================================================
# Verification — reconnect read-only and print row counts + column names
# =============================================================================
message("\n── Verification ──────────────────────────────────────────────────────")
morphic <- dbConnect(duckdb(), MORPHIC_PATH, read_only = TRUE)

for (tbl in c("impc_viability", "impc_phenotypes", "omim_phenotypes", "orphanet_data")) {
  n    <- dbGetQuery(morphic, sprintf('SELECT COUNT(*) AS n FROM main."%s"', tbl))$n
  cols <- paste(dbListFields(morphic, tbl), collapse = ", ")
  message(sprintf("  %-22s  %6d rows   [%s]", tbl, n, cols))
}

# Quick sanity: show head of each table
message("\nSample rows per table:")
for (tbl in c("impc_viability", "impc_phenotypes", "omim_phenotypes", "orphanet_data")) {
  message(sprintf("\n  -- %s --", tbl))
  print(dbGetQuery(morphic, sprintf('SELECT * FROM main."%s" LIMIT 3', tbl)))
}

dbDisconnect(morphic, shutdown = TRUE)
message("\nDone.")
