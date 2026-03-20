# =============================================================================
# helpers.R — shared constants and functions used across modules
# Source this before any module that needs .HP_*, load_deg_genes, etc.
# =============================================================================

library(DT)

# ── Annotation source maps (used by studies_overview + deg_annotations) ───────

.HP_SOURCES <- c(
  "Reactome Pathways"       = "reactome",
  "GO Biological Process"   = "go_bp",
  "GO Molecular Function"   = "go_mf",
  "PANTHER Protein Class"   = "panther_class",
  "PANTHER Protein Family"  = "panther_family",
  "PanelApp Panels"         = "panelapp"
)
.HP_GRP <- c(
  reactome       = "path_name",
  go_bp          = "go_term_BP",
  go_mf          = "go_term_MF",
  panther_class  = "CLASS_TERM",
  panther_family = "FAMILY_TERM",
  panelapp       = "panel_name"
)
.HP_N <- c(
  reactome       = "pathway_genes_n",
  go_bp          = "biological_process_genes_n",
  go_mf          = "molecular_function_genes_n",
  panther_class  = "class_genes_n",
  panther_family = "family_genes_n",
  panelapp       = "panel_genes_n"
)
.HP_DEG_SRC <- c(
  reactome       = "Reactome",
  go_bp          = "GO_BP",
  go_mf          = "GO_MF",
  panther_class  = "PANTHER",
  panther_family = "PANTHER",
  panelapp       = "PanelApp"
)
.HP_DEG_CATTYPE <- c(
  reactome       = NA_character_,
  go_bp          = NA_character_,
  go_mf          = NA_character_,
  panther_class  = "protein_class",
  panther_family = "protein_family",
  panelapp       = NA_character_
)
.HP_REF_TBL <- c(
  reactome       = "ref_reactome",
  go_bp          = "ref_go_bp",
  go_mf          = "ref_go_mf",
  panther_class  = "ref_panther_class",
  panther_family = "ref_panther_family"
)
.HP_CLRS <- list(up = "#E63946", down = "#457B9D", all = "#6A4C93")

# ── Metadata badge colour helpers ─────────────────────────────────────────────

.BADGE_PALETTE <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#86BCB6", "#8CD17D", "#D4A6C8", "#D7B5A6", "#A0CBE8"
)

#' Render a coloured rounded-pill badge for a metadata value.
#' @param val   Single value to display
#' @param pool  All values in that category (used to assign a stable colour)
.meta_badge <- function(val, pool) {
  if (is.null(val) || is.na(val) || !nzchar(as.character(val))) return(NULL)
  vals <- sort(unique(na.omit(as.character(pool))))
  idx  <- match(as.character(val), vals)
  if (is.na(idx)) idx <- 1L
  bg   <- .BADGE_PALETTE[((idx - 1L) %% length(.BADGE_PALETTE)) + 1L]
  rgb_v <- col2rgb(bg)
  lum   <- (0.299 * rgb_v[1] + 0.587 * rgb_v[2] + 0.114 * rgb_v[3]) / 255
  tc    <- if (lum > 0.5) "#212529" else "#ffffff"
  tags$span(
    class = "badge rounded-pill",
    style = sprintf("background-color:%s; color:%s; font-size:0.68rem; font-weight:500;", bg, tc),
    as.character(val)
  )
}

# ── Gene chip colour helpers ──────────────────────────────────────────────────

.hp_chip_colour <- function(is_up, lfc_abs) {
  lvl <- if (!is.na(lfc_abs) && lfc_abs >= 2) 3L else
         if (!is.na(lfc_abs) && lfc_abs >= 1) 2L else 1L
  if (is_up) c("#F1948A", "#E74C3C", "#C0392B")[lvl]
  else       c("#85C1E9", "#5DADE2", "#2C7BB6")[lvl]
}
.hp_chip_style <- function(bg) {
  paste0("background-color:", bg, "; color:white; font-size:0.78rem; font-weight:600;",
         " padding:0.28em 0.6em; border-radius:4px; display:inline-block; margin:2px 3px 2px 0;")
}

# ── DEG symbol SQL fallback ───────────────────────────────────────────────────

.display_sym_sql <-
  "CASE WHEN hgnc_symbol IS NULL OR hgnc_symbol = ''
        THEN gene_ID ELSE hgnc_symbol END"

# ── logFC cell formatter (R side) ─────────────────────────────────────────────

.fmt_dir_cell <- function(is_up, is_dn, log2fc = NA_real_) {
  if (!is_up && !is_dn) return("")
  if (is_up  && is_dn)  return("\u2191\u2193")
  dir <- if (is_up) "\u2191" else "\u2193"
  if (is.na(log2fc)) return(dir)
  paste0(dir, " ", ifelse(log2fc >= 0, "+", ""), round(log2fc, 2))
}

# ── logFC cell formatter (JS DT rowCallback) ──────────────────────────────────

.LFC_ROW_CB <- DT::JS(
  "function(row, data, index) {",
  "  var UP='\u2191', DN='\u2193', UPDN='\u2191\u2193';",
  "  $(row).find('td').each(function() {",
  "    var v=$(this).text().trim(); if(!v) return;",
  "    var both=v.substring(0,2)===UPDN, isUp=!both&&v.charAt(0)===UP, isDn=!both&&v.charAt(0)===DN;",
  "    if(!isUp&&!isDn&&!both) return;",
  "    var bg,col;",
  "    if(both){bg='#f3e8fd';col='#7B2D8B';}",
  "    else{",
  "      var lfc=v.indexOf(' ')>-1?Math.abs(parseFloat(v.split(' ')[1])):0;",
  "      var l=lfc>=2?2:(lfc>=1?1:0);",
  "      if(isUp){bg=['#fde8ea','#f8c1c4','#f4989b'][l];col=['#E74C3C','#C0392B','#922B21'][l];}",
  "      else    {bg=['#ddeef7','#b8d8ef','#92c4e8'][l];col=['#457B9D','#2C6E8A','#1A5276'][l];}",
  "    }",
  "    $(this).css({'background-color':bg,'color':col,'font-weight':'bold','text-align':'center'});",
  "  });",
  "}"
)

# ── Download dropdown nav_item helper ─────────────────────────────────────────

.dl_dropdown_nav <- function(ns, pfx) {
  nav_item(
    dropdownButton(
      downloadButton(ns(paste0(pfx, "_dl_csv")),
                     tagList(bsicons::bs_icon("table", size = "0.8rem"),
                             " Download Table (CSV)"),
                     class = "btn btn-outline-secondary btn-sm w-100 mb-2"),
      downloadButton(ns(paste0(pfx, "_dl_png")),
                     tagList(bsicons::bs_icon("card-list", size = "0.8rem"),
                             " Download Plot (PNG)"),
                     class = "btn btn-outline-secondary btn-sm w-100"),
      circle   = FALSE,
      status   = "outline-secondary",
      size     = "sm",
      icon     = bsicons::bs_icon("download", size = "0.8rem"),
      label    = "Download",
      width    = "220px",
      inputId  = ns(paste0(pfx, "_dl_dd"))
    )
  )
}

# ── Info tooltip helper ─────────────────────────────────────────────────────

.info_tip <- function(text) {
  tooltip(
    bsicons::bs_icon("info-circle-fill", size = "0.85rem", class = "text-muted ms-1"),
    text
  )
}

# ── Shared DT header styling ───────────────────────────────────────────────

.dt_header_js <- DT::JS("
  function(settings, json) {
    $(this.api().table().header()).css({
      'background-color': '#f1f3f5',
      'color': '#495057',
      'font-size': '0.82rem',
      'font-weight': '600',
      'letter-spacing': '0.01em',
      'border-bottom': '2px solid #dee2e6'
    });
    $(this.api().table().node()).css('font-size', '0.85rem');
  }
")

# ── Empty state placeholder ────────────────────────────────────────────────

.empty_state <- function(icon, title, subtitle = NULL) {
  div(class = "text-center py-5", style = "color: #6c757d;",
      bsicons::bs_icon(icon, size = "2.5rem"),
      tags$h6(class = "mt-3 fw-semibold", title),
      if (!is.null(subtitle)) tags$p(class = "text-muted small", subtitle))
}

# ── Data loader ───────────────────────────────────────────────────────────────

load_deg_genes <- function(con, assay_names) {
  rows <- lapply(assay_names, function(tbl) {
    col_names  <- tryCatch(dbListFields(con, tbl), error = function(e) character(0))
    lfc_col    <- col_names[grepl("log2FoldChange$", col_names)]
    lfc_clause <- if (length(lfc_col) > 0)
      sprintf(', "%s" AS log2fc', lfc_col[1])
    else
      ", CAST(NULL AS DOUBLE) AS log2fc"
    res <- tryCatch(
      dbGetQuery(con, paste0(
        'SELECT ', .display_sym_sql, ' AS symbol, DEG', lfc_clause,
        ' FROM "', tbl, '"',
        " WHERE DEG IN ('up','down')"
      )),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res$assay <- tbl
    res
  })
  do.call(rbind, rows)
}

# ── On-the-fly overlap computation ──────────────────────────────────────────

#' Compute DEG overlap with a reference gene-set table.
#' @param con     DBI connection to app DuckDB
#' @param assay   Assay table name (character)
#' @param ref_tbl Reference table name (e.g. "ref_reactome")
#' @param dir     Direction: "all", "up", or "down"
#' @param topn    Max rows to return
#' @param maxsize Maximum category size (filter out broad terms)
#' @param min_lfc Minimum absolute log2FC threshold (NULL = no filter)
#' @param max_lfc Maximum absolute log2FC threshold (NULL = no filter)
#' @return data.frame with: category, n_pathway, n_overlap, pct, genes_up, genes_down
compute_deg_overlap <- function(con, assay, ref_tbl, dir, topn, maxsize,
                                order_col = "pct",
                                min_lfc = NULL, max_lfc = NULL) {
  # Build optional log2FC filter clause
  lfc_clause <- ""
  if (!is.null(min_lfc) || !is.null(max_lfc)) {
    lfc_col <- tryCatch({
      cols <- dbListFields(con, assay)
      cols[grepl("log2FoldChange$", cols)][1]
    }, error = function(e) NA_character_)
    if (!is.na(lfc_col)) {
      parts <- character(0)
      if (!is.null(min_lfc) && is.finite(min_lfc))
        parts <- c(parts, sprintf('"%s" >= %s', lfc_col, min_lfc))
      if (!is.null(max_lfc) && is.finite(max_lfc))
        parts <- c(parts, sprintf('"%s" <= %s', lfc_col, max_lfc))
      if (length(parts) > 0)
        lfc_clause <- paste0(" AND ", paste(parts, collapse = " AND "))
    }
  }
  n_col <- paste0("n_", dir)
  order_expr <- if (identical(order_col, "n_overlap"))
    sprintf("o.%s DESC", n_col) else "pct DESC"
  sql <- sprintf('
    WITH degs AS (
      SELECT COALESCE(NULLIF(TRIM(hgnc_symbol), \'\'), gene_ID) AS symbol, DEG
      FROM main."%s"
      WHERE DEG IN (\'up\', \'down\')%s
    ),
    cat_sizes AS (
      SELECT category, COUNT(DISTINCT symbol) AS n_pathway
      FROM %s
      GROUP BY category
    ),
    ovlp AS (
      SELECT r.category,
             STRING_AGG(DISTINCT CASE WHEN d.DEG = \'up\'   THEN d.symbol END, \'; \') AS genes_up,
             STRING_AGG(DISTINCT CASE WHEN d.DEG = \'down\' THEN d.symbol END, \'; \') AS genes_down,
             COUNT(DISTINCT CASE WHEN d.DEG = \'up\'   THEN d.symbol END) AS n_up,
             COUNT(DISTINCT CASE WHEN d.DEG = \'down\' THEN d.symbol END) AS n_down,
             COUNT(DISTINCT d.symbol) AS n_all
      FROM %s r
      INNER JOIN degs d ON r.symbol = d.symbol
      GROUP BY r.category
    )
    SELECT o.category,
           cs.n_pathway,
           o.%s AS n_overlap,
           ROUND(100.0 * o.%s / cs.n_pathway, 2) AS pct,
           o.genes_up,
           o.genes_down
    FROM ovlp o
    JOIN cat_sizes cs ON o.category = cs.category
    WHERE cs.n_pathway <= %d
    ORDER BY %s
    LIMIT %d',
    assay, lfc_clause, ref_tbl, ref_tbl, n_col, n_col, maxsize, order_expr, topn)
  tryCatch(dbGetQuery(con, sql), error = function(e) NULL)
}

#' Compute KO gene overlap with a reference gene-set table.
#' @param con     DBI connection to app DuckDB
#' @param ref_tbl Reference table name (e.g. "ref_reactome")
#' @return data.frame with: category, n_pathway, overlap_n, pct, genes
compute_ko_overlap <- function(con, ref_tbl) {
  sql <- sprintf('
    WITH ko_genes AS (
      SELECT DISTINCT Gene AS symbol
      FROM study_info
      WHERE Gene IS NOT NULL
    ),
    cat_sizes AS (
      SELECT category, COUNT(DISTINCT symbol) AS n_pathway
      FROM %s
      GROUP BY category
    ),
    ovlp AS (
      SELECT r.category,
             COUNT(DISTINCT k.symbol) AS overlap_n,
             STRING_AGG(DISTINCT k.symbol, \'; \' ORDER BY k.symbol) AS genes
      FROM %s r
      INNER JOIN ko_genes k ON r.symbol = k.symbol
      GROUP BY r.category
    )
    SELECT o.category,
           cs.n_pathway,
           o.overlap_n,
           ROUND(100.0 * o.overlap_n / (SELECT COUNT(*) FROM ko_genes), 2) AS pct,
           o.genes
    FROM ovlp o
    JOIN cat_sizes cs ON o.category = cs.category
    WHERE cs.n_pathway >= 1
    ORDER BY pct DESC',
    ref_tbl, ref_tbl)
  tryCatch(dbGetQuery(con, sql), error = function(e) NULL)
}
