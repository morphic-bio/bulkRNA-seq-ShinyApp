# =============================================================================
# make-db.R — Build the app DB (morphic_bulkRNA3.duckdb) from scratch
#
# Run from repo root:  Rscript RNAseq_app/make-db.R
#
# Pipeline:
#   1. Split multi-KO-strategy tables into per-strategy tables
#   2. Add DEG column, keep only: gene_ID, hgnc_symbol, log2FoldChange, DEG
#   3. Add study_info from bulk-study-info2.xlsx
#   4. Add reference gene-set tables from EGEx-db (Reactome, GO, PANTHER)
#   5. Add phenotype / disease tables from EGEx-db (IMPC, HPO, OMIM, Orphanet)
#   6. Copy final DB to RNAseq_app/ for the Shiny app
#
# Sources:
#   - Raw DESeq2 results:  bulkRNA-seq/bulkRNA-seq-analysis/DuckDB/morphic_bulkRNA.duckdb
#   - Study metadata:      bulkRNA-seq/bulk-study-info2.xlsx
#   - Gene annotations:    EGEx-db/EGEx-db.duckdb
# =============================================================================

library(DBI); library(duckdb)

RAW_DB   <- "./bulkRNA-seq/bulkRNA-seq-analysis/DuckDB/morphic_bulkRNA.duckdb"
TMP_DB   <- "./bulkRNA-seq/bulkRNA-seq-analysis/DuckDB/morphic_bulkRNA2.duckdb"
BUILD_DB <- "./bulkRNA-seq/bulkRNA-seq-analysis/DuckDB/morphic_bulkRNA3.duckdb"
EGEX_DB  <- "./EGEx-db/EGEx-db.duckdb"
APP_DB   <- "./RNAseq_app/morphic_bulkRNA3.duckdb"

cat("========================================\n")
cat("make-db.R — Building app DB\n")
cat("========================================\n\n")

# =============================================================================
# 1. Split multi-KO-strategy tables
# =============================================================================
cat("── Step 1: Split multi-strategy tables ──\n")

base   <- c("gene_ID", "hgnc_symbol")
strats <- c("CE", "KO", "PTC")

unlink(TMP_DB)
con <- dbConnect(duckdb(), TMP_DB)
dbExecute(con, sprintf("ATTACH '%s' AS src", RAW_DB))

m <- dbGetQuery(con, "
  SELECT table_name, COUNT(*) n
  FROM information_schema.columns
  WHERE table_catalog='src' AND table_schema='main'
  GROUP BY table_name
")

# Single-strategy tables (≤9 cols)
for (t in m$table_name[m$n <= 9]) {
  src_cols <- dbGetQuery(con, sprintf(
    "SELECT column_name FROM information_schema.columns
     WHERE table_catalog='src' AND table_schema='main' AND table_name='%s'",
    gsub("'", "''", t)))$column_name
  lfc <- src_cols[grepl("log2FoldChange$", src_cols)][1]
  pv  <- src_cols[grepl("pvalue$", src_cols)][1]
  keep <- c(sprintf('"%s"', base), sprintf('"%s"', lfc), sprintf('"%s"', pv))
  dbExecute(con, sprintf('CREATE TABLE "%s" AS SELECT %s FROM src.main."%s"',
                         t, paste(keep, collapse = ", "), t))
}

# Multi-strategy tables (>9 cols) — split per strategy
for (t in m$table_name[m$n > 9]) {
  gene_suffix <- sub(".*_", "", t)
  for (s in strats) {
    lfc_col <- sprintf("%s_%s_log2FoldChange", gene_suffix, s)
    pv_col  <- sprintf("%s_%s_pvalue", gene_suffix, s)
    sel <- paste(c(sprintf('"%s"', base),
                   sprintf('"%s"', lfc_col),
                   sprintf('"%s"', pv_col)), collapse = ", ")
    dbExecute(con, sprintf('CREATE TABLE "%s_%s" AS SELECT %s FROM src.main."%s"',
                           t, s, sel, t))
  }
}

dbExecute(con, "DETACH src")
n_tmp <- length(dbListTables(con))
dbDisconnect(con, shutdown = TRUE)
cat("  Created", n_tmp, "tables in temp DB\n\n")

# =============================================================================
# 2. Add DEG column, keep only gene_ID + hgnc_symbol + log2FC + DEG
# =============================================================================
cat("── Step 2: Add DEG column, drop pvalue ──\n")

unlink(BUILD_DB)
con <- dbConnect(duckdb(), BUILD_DB)
dbExecute(con, sprintf("ATTACH '%s' AS src", TMP_DB))

tabs <- dbGetQuery(con, "
  SELECT table_name
  FROM information_schema.tables
  WHERE table_catalog='src' AND table_schema='main' AND table_type='BASE TABLE'
  ORDER BY table_name
")$table_name

n_ok <- 0L
for (t in tabs) {
  cols <- dbGetQuery(con, sprintf("
    SELECT column_name
    FROM information_schema.columns
    WHERE table_catalog='src' AND table_schema='main' AND table_name='%s'
  ", gsub("'", "''", t)))$column_name
  lfc <- cols[grepl("log2FoldChange$", cols)][1]
  pv  <- cols[grepl("pvalue$", cols)][1]
  if (is.na(lfc) || is.na(pv)) next
  dbExecute(con, sprintf(
    'CREATE TABLE main."%s" AS
     SELECT "gene_ID", "hgnc_symbol", "%s",
       CASE
         WHEN "%s" < 0.05 AND "%s" >= 1  THEN \'up\'
         WHEN "%s" < 0.05 AND "%s" <= -1 THEN \'down\'
         ELSE \'no\'
       END AS DEG
     FROM src.main."%s"',
    t, lfc, pv, lfc, pv, lfc, t
  ))
  n_ok <- n_ok + 1L
}

dbExecute(con, "DETACH src")
dbDisconnect(con, shutdown = TRUE)
cat("  Created", n_ok, "assay tables (4 columns each)\n\n")

# Clean up temp DB
unlink(TMP_DB)

# =============================================================================
# 3. Add study_info
# =============================================================================
cat("── Step 3: Add study_info ──\n")

study_info <- openxlsx::read.xlsx("./bulkRNA-seq/bulk-study-info3.xlsx")
# names(study_info)[names(study_info) == "ASCL2"] <- "Gene"

con <- dbConnect(duckdb(), BUILD_DB)
dbWriteTable(con, "study_info", study_info, overwrite = TRUE)
dbDisconnect(con, shutdown = TRUE)
cat("  study_info:", nrow(study_info), "rows,", ncol(study_info), "columns\n")
cat("  Unique KO genes:", length(unique(na.omit(study_info$Gene))), "\n\n")

# =============================================================================
# 4. Add reference gene-set tables from EGEx-db
# =============================================================================
cat("── Step 4: Add reference gene-set tables ──\n")

con <- dbConnect(duckdb(), BUILD_DB)
dbExecute(con, sprintf("ATTACH '%s' AS egex (READ_ONLY)", EGEX_DB))

ref_specs <- list(
  ref_reactome = "
    SELECT DISTINCT r.path_name AS category, g.symbol
    FROM egex.Reactome r
    JOIN egex.genes g ON r.GeneID = g.GeneID
    WHERE r.path_name IS NOT NULL AND g.symbol IS NOT NULL",
  ref_go_bp = "
    SELECT DISTINCT bp.go_term_BP AS category, g.symbol
    FROM egex.GO_BP bp
    JOIN egex.genes g ON bp.GeneID = g.GeneID
    WHERE bp.go_term_BP IS NOT NULL AND g.symbol IS NOT NULL",
  ref_go_mf = "
    SELECT DISTINCT mf.go_term_MF AS category, g.symbol
    FROM egex.GO_MF mf
    JOIN egex.genes g ON mf.GeneID = g.GeneID
    WHERE mf.go_term_MF IS NOT NULL AND g.symbol IS NOT NULL",
  ref_panther_class = "
    SELECT DISTINCT p.CLASS_TERM AS category, g.symbol
    FROM egex.PANTHERdb p
    JOIN egex.genes g ON p.GeneID = g.GeneID
    WHERE p.CLASS_TERM IS NOT NULL AND g.symbol IS NOT NULL",
  ref_panther_family = "
    SELECT DISTINCT p.FAMILY_TERM AS category, g.symbol
    FROM egex.PANTHERdb p
    JOIN egex.genes g ON p.GeneID = g.GeneID
    WHERE p.FAMILY_TERM IS NOT NULL AND g.symbol IS NOT NULL"
)

for (tbl in names(ref_specs)) {
  dbExecute(con, sprintf("DROP TABLE IF EXISTS %s", tbl))
  dbExecute(con, sprintf("CREATE TABLE %s AS %s", tbl, ref_specs[[tbl]]))
  n <- dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM %s", tbl))$n
  nc <- dbGetQuery(con, sprintf("SELECT COUNT(DISTINCT category) AS n FROM %s", tbl))$n
  cat(sprintf("  %s: %s rows, %s categories\n", tbl, format(n, big.mark = ","), format(nc, big.mark = ",")))
}

dbExecute(con, "DETACH egex")
dbDisconnect(con, shutdown = TRUE)
cat("\n")

# =============================================================================
# 5. Add phenotype / disease tables from EGEx-db
# =============================================================================
cat("── Step 5: Add phenotype / disease tables ──\n")

con <- dbConnect(duckdb(), BUILD_DB)
dbExecute(con, sprintf("ATTACH '%s' AS egex (READ_ONLY)", EGEX_DB))

pheno_specs <- list(
  impc_phenotypes = "
    SELECT DISTINCT g.symbol, p.GeneID, p.impc_phenotypes, p.impc_zygosity
    FROM egex.IMPC_phenotypes p
    JOIN egex.genes g ON p.GeneID = g.GeneID
    WHERE g.symbol IS NOT NULL AND p.impc_phenotypes IS NOT NULL",
  hpo_phenotypes = "
    SELECT DISTINCT g.symbol, p.GeneID, p.hpo_id, p.hpo_name
    FROM egex.HPO_phenotypes p
    JOIN egex.genes g ON p.GeneID = g.GeneID
    WHERE g.symbol IS NOT NULL AND p.hpo_name IS NOT NULL",
  omim_phenotypes = "
    SELECT DISTINCT g.symbol, p.GeneID, p.phenotypes AS omim_phenotype
    FROM egex.OMIM_phenotypes p
    JOIN egex.genes g ON p.GeneID = g.GeneID
    WHERE g.symbol IS NOT NULL AND p.phenotypes IS NOT NULL",
  orphanet_data = "
    SELECT DISTINCT g.symbol, p.GeneID, p.disorder_name, p.orpha_code, p.association_type
    FROM egex.OrphaData p
    JOIN egex.genes g ON p.GeneID = g.GeneID
    WHERE g.symbol IS NOT NULL AND p.disorder_name IS NOT NULL"
)

for (tbl in names(pheno_specs)) {
  dbExecute(con, sprintf('DROP TABLE IF EXISTS "%s"', tbl))
  dbExecute(con, sprintf('CREATE TABLE "%s" AS %s', tbl, pheno_specs[[tbl]]))
  n <- dbGetQuery(con, sprintf('SELECT COUNT(*) AS n FROM "%s"', tbl))$n
  cat(sprintf("  %s: %s rows\n", tbl, format(n, big.mark = ",")))
}

dbExecute(con, "DETACH egex")
dbDisconnect(con, shutdown = TRUE)
cat("\n")

# =============================================================================
# 6. Copy to app directory
# =============================================================================
cat("── Step 6: Copy to app directory ──\n")

file.copy(BUILD_DB, APP_DB, overwrite = TRUE)
cat("  Copied to:", APP_DB, "\n\n")

# =============================================================================
# Final summary
# =============================================================================
cat("========================================\n")
cat("=== Final DB Summary ===\n")
cat("========================================\n")

con <- dbConnect(duckdb(), APP_DB, read_only = TRUE)
all_tabs <- dbListTables(con)
meta_tabs <- c("study_info", "ref_reactome", "ref_go_bp", "ref_go_mf",
               "ref_panther_class", "ref_panther_family",
               "impc_phenotypes", "hpo_phenotypes", "omim_phenotypes", "orphanet_data")
assay_tabs <- setdiff(all_tabs, meta_tabs)

cat("Assay tables:", length(assay_tabs), "\n")
cat("Metadata tables:", length(intersect(all_tabs, meta_tabs)), "\n")

# Show sample assay columns
t1 <- assay_tabs[1]
cat("\nSample assay table:", t1, "\n")
cols <- dbGetQuery(con, sprintf("PRAGMA table_info('%s')", gsub("'", "''", t1)))
cat("Columns:", paste(cols$name, collapse = ", "), "\n")

dbDisconnect(con, shutdown = TRUE)

size_mb <- round(file.size(APP_DB) / 1024 / 1024, 1)
cat("\nApp DB size:", size_mb, "MB\n")
cat("Location:", normalizePath(APP_DB), "\n")
cat("\nDone!\n")
