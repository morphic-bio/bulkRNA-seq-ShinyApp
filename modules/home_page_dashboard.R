# =============================================================================
# Module: home_page_dashboard
# Description: Landing overview page.
#
# Layout (top → bottom):
#   1. Study Breakdown          — grouped bar (assays & KO genes by group)
#   2. N Assays per KO Gene     — horizontal bar, separate card
#   3. Assay Overview Table     — flat DT, all columns, no row interaction
#   4. KO Gene Pathway Overlap  — morphic_overlap bar + table (KO mode only)
#   ── DEG Analysis divider ──
#   5. Assay Comparison         — filter+picker, metadata header, tabs:
#                                  Gene Matrix | UpSet (up/down)
#   6. DEG Pathway Overlap      — filter+picker, bar|table tabs → clicked row
#                                  opens coloured chip panel (up/down separated)
#
# Data contract
#   con_r        — reactive: open read-only DuckDB connection
#   study_info_r — reactive: study_info data.frame
#
# Helpers re-used from global scope (sourced via deg_summary_dashboard.R):
#   .display_sym_sql, .LFC_ROW_CB, .fmt_dir_cell, load_deg_genes
# =============================================================================

library(plotly)
library(DT)
library(UpSetR)

# ── Module-local constants ─────────────────────────────────────────────────────

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
.HP_CLRS <- list(up = "#E63946", down = "#457B9D", all = "#6A4C93")

# Chip colour: 3-level scale per direction (no "both" — one assay = one direction)
.hp_chip_colour <- function(is_up, lfc_abs) {
  lvl <- if (!is.na(lfc_abs) && lfc_abs >= 2) 3L else if (!is.na(lfc_abs) && lfc_abs >= 1) 2L else 1L
  if (is_up) c("#F1948A", "#E74C3C", "#C0392B")[lvl]
  else       c("#85C1E9", "#5DADE2", "#2C7BB6")[lvl]
}
.hp_chip_style <- function(bg) {
  paste0("background-color:", bg, "; color:white; font-size:0.78rem; font-weight:600;",
         " padding:0.28em 0.6em; border-radius:4px; display:inline-block; margin:2px 3px 2px 0;")
}

# ── UI ────────────────────────────────────────────────────────────────────────

home_page_dashboardUI <- function(id) {
  ns <- NS(id)

  tagList(

    # ── 1. Study Breakdown ────────────────────────────────────────────────────
    card(
      card_header(
        class = "d-flex justify-content-between align-items-center",
        span(tagList(bsicons::bs_icon("bar-chart-fill"), " Study Breakdown")),
        div(class = "d-flex align-items-center gap-2",
            tags$label("Group by:", class = "mb-0 small fw-semibold"),
            selectInput(ns("group_by"), NULL,
                        choices  = c("Model System" = "Model_system",
                                     "KO Strategy"  = "KO_strat",
                                     "Cell Line"    = "Cell_Line",
                                     "DPC"          = "DPC"),
                        selected = "Model_system", width = "150px"))
      ),
      plotlyOutput(ns("breakdown_bar"), height = "360px"),
      full_screen = TRUE
    ),

    # ── 2. N Assays per KO Gene ───────────────────────────────────────────────
    card(
      card_header(tagList(bsicons::bs_icon("bar-chart-fill"), " N Assays per KO Gene")),
      plotlyOutput(ns("assays_per_gene"), height = "420px"),
      full_screen = TRUE
    ),

    # ── 3. Assay Overview Table ───────────────────────────────────────────────
    card(
      card_header(tagList(bsicons::bs_icon("table"), " Assay Overview")),
      DTOutput(ns("assay_tbl")),
      full_screen = TRUE
    ),

    # ── 4. KO Gene Pathway Overlap ────────────────────────────────────────────
    card(
      card_header(tagList(bsicons::bs_icon("intersect"), " KO Gene Pathway Overlap")),
      layout_columns(
        col_widths = c(3, 9),
        card(card_body(
          selectInput(ns("ko_source"), "Annotation Source",
                      choices = .HP_SOURCES, selected = "reactome"),
          sliderInput(ns("ko_topn"),    "Top N categories",  5, 50, 20, step = 5),
          sliderInput(ns("ko_minsize"), "Min category size", 1, 500, 10)
        )),
        navset_tab(
          nav_panel("Bar Chart", plotlyOutput(ns("ko_bar"), height = "450px")),
          nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl"))))
        )
      ),
      full_screen = TRUE
    ),

    # ── DEG Analysis divider ──────────────────────────────────────────────────
    div(
      class = "my-4",
      tags$hr(style = "border-top: 2px solid #adb5bd;"),
      div(class = "text-center", style = "margin-top: -1.1em;",
          tags$span(
            class = "px-3 fw-semibold",
            style = paste0("background:#fff; color:#6c757d; font-size:0.75rem;",
                           " text-transform:uppercase; letter-spacing:0.12em;"),
            bsicons::bs_icon("scissors"), " DEG Analysis"
          ))
    ),

    # ── 5. Assay Comparison ───────────────────────────────────────────────────
    card(
      card_header(tagList(bsicons::bs_icon("grid-3x3-gap-fill"), " Assay Comparison")),
      layout_columns(
        col_widths = c(3, 9),

        # Left: filter + picker
        card(card_body(
          tags$p(class = "fw-semibold small mb-1", "Filter assay list:"),
          selectInput(ns("cmp_filter_gene"),  "KO Gene",
                      choices = NULL, multiple = TRUE, width = "100%"),
          selectInput(ns("cmp_filter_model"), "Model System",
                      choices = NULL, multiple = TRUE, width = "100%"),
          selectInput(ns("cmp_filter_kos"),   "KO Strategy",
                      choices = NULL, multiple = TRUE, width = "100%"),
          tags$hr(class = "my-2"),
          selectizeInput(ns("cmp_assays"), "Select assays to compare:",
                         choices = NULL, selected = NULL, multiple = TRUE,
                         options = list(placeholder = "Choose assays\u2026",
                                        plugins     = list("remove_button"),
                                        maxItems    = NULL),
                         width = "100%"),
          selectInput(ns("cmp_dir"), "Direction",
                      choices  = c("All" = "all", "Up only" = "up",
                                   "Down only" = "down"),
                      selected = "all"),
          tags$p(class = "text-muted small mt-1", "Select 2 or more assays.")
        )),

        # Right: metadata + tabbed results
        card(card_body(
          uiOutput(ns("cmp_meta_above")),
          navset_tab(
            nav_panel("Gene Matrix",
                      div(class = "pt-2", DTOutput(ns("cmp_matrix_tbl")))),
            nav_panel("UpSet Plot",
                      div(class = "pt-2",
                          div(class = "d-flex align-items-center gap-3 mb-2",
                              tags$label("Direction:", class = "mb-0 small fw-semibold"),
                              radioButtons(ns("cmp_upset_dir"), NULL,
                                           choices  = c("\u2191 Up" = "up",
                                                        "\u2193 Down" = "down",
                                                        "\u2191\u2193 All DEGs" = "all"),
                                           selected = "up", inline = TRUE)),
                          plotOutput(ns("cmp_upset"), height = "460px"))),
            nav_panel("Common DEGs",
                      div(class = "pt-2", uiOutput(ns("cmp_consensus"))))
          )
        ))
      ),
      full_screen = TRUE
    ),

    # ── 6. DEG Pathway Overlap ────────────────────────────────────────────────
    card(
      card_header(tagList(bsicons::bs_icon("diagram-3"), " DEG Pathway Overlap")),
      layout_columns(
        col_widths = c(3, 9),

        # Left: filter + picker + controls
        card(card_body(
          tags$p(class = "fw-semibold small mb-1", "Filter assay list:"),
          selectInput(ns("deg_filter_gene"),  "KO Gene",
                      choices = NULL, multiple = TRUE, width = "100%"),
          selectInput(ns("deg_filter_model"), "Model System",
                      choices = NULL, multiple = TRUE, width = "100%"),
          selectInput(ns("deg_filter_kos"),   "KO Strategy",
                      choices = NULL, multiple = TRUE, width = "100%"),
          tags$hr(class = "my-2"),
          selectizeInput(ns("deg_assay"), "Select assay:",
                         choices = NULL, selected = NULL, multiple = FALSE,
                         options = list(placeholder = "Choose an assay\u2026"),
                         width = "100%"),
          tags$hr(class = "my-2"),
          selectInput(ns("deg_source"), "Annotation Source",
                      choices = .HP_SOURCES, selected = "reactome"),
          selectInput(ns("deg_dir"), "Direction",
                      choices  = c("Up & Down" = "updn",
                                   "Up only"    = "up",
                                   "Down only"  = "down"),
                      selected = "updn"),
          sliderInput(ns("deg_topn"),    "Top N pathways",    5, 50, 20, step = 5),
          sliderInput(ns("deg_minsize"), "Min category size", 1, 500, 10)
        )),

        # Right: bar | table tabs, then gene chip panel below
        card(
          navset_tab(
            nav_panel("Bar Chart",
                      plotlyOutput(ns("deg_bar"), height = "420px")),
            nav_panel("Pathway Table",
                      div(class = "p-2",
                          tags$p(class = "text-muted small mb-1",
                                 bsicons::bs_icon("table"),
                                 " Click a row to see genes for that pathway."),
                          DTOutput(ns("deg_pathway_tbl"))))
          ),
          uiOutput(ns("deg_gene_panel"))
        )
      ),
      full_screen = TRUE
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

home_page_dashboardServer <- function(id, con_r, study_info_r) {
  moduleServer(id, function(input, output, session) {

    # =========================================================================
    # 1. Study Breakdown
    # =========================================================================

    output$breakdown_bar <- renderPlotly({
      si      <- study_info_r()
      grp_col <- input$group_by
      req(grp_col %in% colnames(si))
      agg_a <- tapply(si$Assay, si[[grp_col]], function(x) length(unique(x)))
      agg_g <- tapply(si$Gene,  si[[grp_col]], function(x) length(unique(na.omit(x))))
      df <- data.frame(group    = names(agg_a),
                       n_assays = as.integer(agg_a),
                       n_genes  = as.integer(agg_g[names(agg_a)]),
                       stringsAsFactors = FALSE)
      df <- df[order(-df$n_assays), ]
      df$group <- factor(df$group, levels = rev(df$group))
      grp_lbl <- switch(grp_col, Model_system="Model System",
                        KO_strat="KO Strategy", Cell_Line="Cell Line", DPC="DPC")
      plot_ly(df) |>
        add_bars(x = ~n_assays, y = ~group, name = "Assays", orientation = "h",
                 marker = list(color = "#4E79A7"),
                 hovertemplate = "<b>%{y}</b><br>Assays: %{x}<extra></extra>") |>
        add_bars(x = ~n_genes, y = ~group, name = "KO Genes", orientation = "h",
                 marker = list(color = "#F28E2B"),
                 hovertemplate = "<b>%{y}</b><br>KO Genes: %{x}<extra></extra>") |>
        layout(barmode = "group",
               xaxis = list(title = "Count"),
               yaxis = list(title = grp_lbl, automargin = TRUE),
               legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.18),
               margin = list(l = 5, r = 5, t = 35, b = 5),
               title  = list(text = paste("Assays & KO Genes by", grp_lbl),
                             font = list(size = 13))) |>
        config(displayModeBar = FALSE)
    })

    # =========================================================================
    # 2. N Assays per KO Gene
    # =========================================================================

    output$assays_per_gene <- renderPlotly({
      si  <- study_info_r()
      agg <- sort(tapply(si$Assay, si$Gene, function(x) length(unique(x))),
                  decreasing = FALSE)
      df  <- data.frame(gene     = names(agg), n_assays = as.integer(agg),
                        stringsAsFactors = FALSE)
      df$gene <- factor(df$gene, levels = df$gene)
      plot_ly(df, x = ~n_assays, y = ~gene, type = "bar", orientation = "h",
              marker = list(color = "#59A14F"),
              hovertemplate = "<b>%{y}</b><br>Assays: %{x}<extra></extra>") |>
        layout(xaxis  = list(title = "N Assays", dtick = 1),
               yaxis  = list(title = "", automargin = TRUE, tickfont = list(size = 10)),
               margin = list(l = 5, r = 5, t = 10, b = 5)) |>
        config(displayModeBar = FALSE)
    })

    # =========================================================================
    # 3. Assay Overview Table — flat, no row interaction
    # =========================================================================

    output$assay_tbl <- DT::renderDT({
      si   <- study_info_r()
      want <- c("Gene", "Model_system", "KO_strat", "Cell_Line", "DPC",
                "Differentation_time_point", "Condition_levels", "Assay")
      cols <- want[want %in% colnames(si)]
      tbl  <- si[!duplicated(si$Assay), cols, drop = FALSE]
      tbl  <- tbl[order(tbl$Gene, tbl$Model_system, tbl$Assay), ]
      pretty <- c(Gene = "KO Gene", Model_system = "Model System",
                  KO_strat = "KO Strategy", Cell_Line = "Cell Line", DPC = "DPC",
                  Differentation_time_point = "Diff. Time Point",
                  Condition_levels = "Condition Levels", Assay = "Assay")
      colnames(tbl) <- pretty[colnames(tbl)]
      DT::datatable(tbl, rownames = FALSE, selection = "none", filter = "top",
                    extensions = "Buttons",
                    options = list(pageLength = 25, scrollX = TRUE,
                                   scrollY = "500px", scrollCollapse = TRUE,
                                   dom = "Bfrtip", buttons = list("csv", "excel"),
                                   columnDefs = list(list(className = "dt-left",
                                                          targets = "_all"))),
                    class = "compact row-border hover")
    })

    # =========================================================================
    # 4. KO Gene Pathway Overlap
    # =========================================================================

    ko_data <- reactive({
      src     <- input$ko_source; top_n <- input$ko_topn; minsize <- input$ko_minsize
      req(src, top_n, minsize)
      cat_col <- .HP_GRP[[src]]; n_col <- .HP_N[[src]]; con <- con_r()
      tryCatch(
        dbGetQuery(con, sprintf(
          'SELECT "%s" AS category, "%s" AS n_pathway,
                  overlap_n, pct_morphic_overlap AS pct, overlap_genes AS genes
           FROM morphic_overlap."%s" WHERE "%s" >= %d ORDER BY pct DESC LIMIT %d',
          cat_col, n_col, src, n_col, minsize, top_n)),
        error = function(e) NULL)
    })

    output$ko_bar <- renderPlotly({
      df <- ko_data(); req(df, nrow(df) > 0)
      df$pct <- round(df$pct, 1)
      df$category <- ifelse(nchar(df$category) > 55,
                            paste0(substr(df$category, 1, 52), "\u2026"), df$category)
      df$category <- factor(df$category, levels = rev(df$category))
      src_label <- names(.HP_SOURCES)[.HP_SOURCES == input$ko_source]
      plot_ly(df, x = ~pct, y = ~category, type = "bar", orientation = "h",
              marker = list(color = "#B07AA1"), customdata = ~genes,
              hovertemplate = paste0("<b>%{y}</b><br>Overlap: %{x}%<br>",
                                    "Genes: %{customdata}<extra></extra>")) |>
        layout(xaxis = list(title = "% MorPhiC KO Genes in Category", range = c(0, 105)),
               yaxis = list(title = "", automargin = TRUE),
               margin = list(l = 5, r = 5, t = 35, b = 5),
               title = list(text = paste0("Top ", nrow(df), " \u2014 ", src_label),
                            font = list(size = 13))) |>
        config(displayModeBar = FALSE)
    })

    output$ko_tbl <- DT::renderDT({
      df <- ko_data(); req(df, nrow(df) > 0)
      df$pct <- round(df$pct, 1); df$n_pathway <- as.integer(df$n_pathway)
      df$overlap_n <- as.integer(df$overlap_n)
      # Build gene → model-systems lookup from study_info
      si          <- study_info_r()
      gene_models <- tapply(si$Model_system, si$Gene,
                            function(x) paste(sort(unique(na.omit(x))), collapse = ", "))
      df$model_systems <- vapply(df$genes, function(gs_str) {
        if (is.null(gs_str) || is.na(gs_str) || !nzchar(trimws(gs_str))) return("")
        gs <- trimws(strsplit(gs_str, ";")[[1]]); gs <- gs[nzchar(gs)]
        ms <- sort(unique(na.omit(unlist(lapply(gs, function(g) {
          m <- gene_models[g]
          if (length(m) == 0 || is.na(m)) return(NA_character_)
          strsplit(as.character(m), ", ")[[1]]
        })))))
        paste(ms, collapse = ", ")
      }, character(1))
      tbl <- df[, c("category", "n_pathway", "overlap_n", "pct", "model_systems", "genes")]
      names(tbl) <- c("Category", "Category Size", "KO Overlap N", "% Overlap",
                      "Model Systems", "Overlapping Genes")
      DT::datatable(tbl, rownames = FALSE, filter = "top", extensions = "Buttons",
                    options = list(pageLength = 25, scrollX = TRUE, scrollY = "440px",
                                   scrollCollapse = TRUE, dom = "Bfrtip",
                                   buttons = list("csv", "excel")),
                    class = "compact row-border hover")
    })

    # =========================================================================
    # 5. Assay Comparison
    # =========================================================================

    # ── Populate filter choices ───────────────────────────────────────────────
    observe({
      si <- study_info_r()
      updateSelectInput(session, "cmp_filter_gene",
                        choices  = sort(unique(na.omit(si$Gene))),
                        selected = character(0))
      updateSelectInput(session, "cmp_filter_model",
                        choices  = sort(unique(na.omit(si$Model_system))),
                        selected = character(0))
      updateSelectInput(session, "cmp_filter_kos",
                        choices  = sort(unique(na.omit(si$KO_strat))),
                        selected = character(0))
    })

    # ── Filtered assay list ───────────────────────────────────────────────────
    filtered_cmp_assays <- reactive({
      si      <- study_info_r()
      si_uniq <- si[!duplicated(si$Assay), ]
      if (length(input$cmp_filter_gene)  > 0)
        si_uniq <- si_uniq[si_uniq$Gene         %in% input$cmp_filter_gene,  ]
      if (length(input$cmp_filter_model) > 0)
        si_uniq <- si_uniq[si_uniq$Model_system %in% input$cmp_filter_model, ]
      if (length(input$cmp_filter_kos)   > 0)
        si_uniq <- si_uniq[si_uniq$KO_strat     %in% input$cmp_filter_kos,   ]
      sort(si_uniq$Assay)
    })

    observeEvent(filtered_cmp_assays(), {
      ch   <- filtered_cmp_assays()
      kept <- intersect(isolate(input$cmp_assays), ch)
      updateSelectizeInput(session, "cmp_assays",
                           choices = ch, selected = kept, server = TRUE)
    })

    # ── DEG gene data (lazy, selected assays only) ────────────────────────────
    cmp_genes_data <- reactive({
      assays <- input$cmp_assays
      req(assays, length(assays) >= 2)
      con <- con_r()
      withProgress(message = "Loading DEG data\u2026", value = 0.5, {
        rows <- load_deg_genes(con, assays)
      })
      if (is.null(rows)) return(NULL)
      meta_cols <- intersect(c("Assay", "Gene", "Model_system", "KO_strat", "DPC",
                               "Cell_Line", "Differentation_time_point",
                               "Condition_levels"),
                             colnames(study_info_r()))
      meta <- study_info_r()[!duplicated(study_info_r()$Assay), meta_cols]
      merge(rows, meta, by.x = "assay", by.y = "Assay", all.x = TRUE)
    })

    # ── Per-assay gene sets ───────────────────────────────────────────────────
    cmp_sets <- reactive({
      genes  <- cmp_genes_data(); req(genes)
      assays <- input$cmp_assays; req(assays, length(assays) >= 2)
      dir    <- input$cmp_dir
      sub    <- genes[genes$assay %in% assays, ]
      if (dir == "up")   sub <- sub[sub$DEG == "up",   ]
      if (dir == "down") sub <- sub[sub$DEG == "down", ]
      lapply(setNames(assays, assays),
             function(a) unique(na.omit(sub$symbol[sub$assay == a])))
    })

    # ── Metadata header (all selected assays) ─────────────────────────────────
    output$cmp_meta_above <- renderUI({
      assays <- input$cmp_assays; req(assays, length(assays) >= 2)
      si <- study_info_r()

      bdg <- function(val, bg_col) {
        v <- as.character(val)
        if (is.null(val) || is.na(val) || !nzchar(v)) return(NULL)
        tags$span(class = "badge me-1", style = paste0("background:", bg_col,
                  "; color:white; font-size:0.72rem;"), v)
      }

      cards <- lapply(assays, function(a) {
        r <- si[si$Assay == a, ][1, ]
        extra <- tagList(
          if ("Differentation_time_point" %in% names(r) &&
              !is.na(r$Differentation_time_point) &&
              nzchar(as.character(r$Differentation_time_point)))
            bdg(r$Differentation_time_point, "#6f42c1"),
          if ("Condition_levels" %in% names(r) &&
              !is.na(r$Condition_levels) && nzchar(as.character(r$Condition_levels)))
            bdg(r$Condition_levels, "#d63384"),
          if ("Replicate" %in% names(r) && !is.na(r$Replicate))
            bdg(r$Replicate, "#6c757d")
        )
        div(class = "border rounded p-2 mb-1", style = "background:#f8f9fa; flex:1; min-width:180px;",
            tags$div(class = "small fw-semibold text-break mb-1",
                     style = "font-size:0.68rem; word-break:break-all; color:#495057;", a),
            tags$div(
              bdg(r$Gene,         "#198754"),
              bdg(r$Model_system, "#0d6efd"),
              bdg(r$KO_strat,     "#0dcaf0"),
              bdg(r$DPC,          "#e67e22"),
              bdg(r$Cell_Line,    "#6c757d"),
              extra
            ))
      })
      tagList(div(class = "d-flex flex-wrap gap-2 mb-2", cards), tags$hr(class = "my-2"))
    })

    # ── Full gene matrix ──────────────────────────────────────────────────────
    output$cmp_matrix_tbl <- DT::renderDT({
      sets     <- cmp_sets(); req(sets, length(sets) >= 2)
      genes_df <- cmp_genes_data(); assays <- names(sets); dir <- input$cmp_dir
      sub <- genes_df[genes_df$assay %in% assays, ]
      if (dir == "up")   sub <- sub[sub$DEG == "up",   ]
      if (dir == "down") sub <- sub[sub$DEG == "down", ]
      all_syms <- sort(unique(na.omit(sub$symbol)))
      if (length(all_syms) == 0)
        return(DT::datatable(data.frame(Message = "No DEGs found for selected assays."),
                             rownames = FALSE, options = list(dom = "t")))
      mat <- data.frame(Gene = all_syms, stringsAsFactors = FALSE)
      for (a in assays) {
        a_sub <- sub[sub$assay == a, ]
        is_up <- all_syms %in% a_sub$symbol[a_sub$DEG == "up"]
        is_dn <- all_syms %in% a_sub$symbol[a_sub$DEG == "down"]
        lbl   <- ifelse(nchar(a) > 40, paste0("\u2026", substr(a, nchar(a)-39, nchar(a))), a)
        lfc_v <- if ("log2fc" %in% names(a_sub))
          setNames(a_sub$log2fc, a_sub$symbol) else setNames(numeric(0), character(0))
        mat[[lbl]] <- mapply(function(g, up, dn) {
          lfc <- lfc_v[g]; .fmt_dir_cell(up, dn, if (length(lfc) == 0) NA_real_ else lfc)
        }, all_syms, is_up, is_dn, SIMPLIFY = TRUE, USE.NAMES = FALSE)
      }
      ko_cols <- names(mat)[-1]
      mat$N_Assays <- rowSums(mat[, ko_cols, drop = FALSE] != "")
      mat <- mat[, c("Gene", "N_Assays", ko_cols)]
      mat <- mat[order(-mat$N_Assays), ]
      DT::datatable(mat, rownames = FALSE, filter = "top", extensions = "Buttons",
                    options = list(pageLength = 25, scrollX = TRUE,
                                   scrollY = "460px", scrollCollapse = TRUE,
                                   dom = "Bfrtip", buttons = list("csv", "excel"),
                                   fixedColumns = list(leftColumns = 2),
                                   columnDefs = list(
                                     list(className = "dt-center",
                                          targets = seq_along(ko_cols) + 1),
                                     list(width = "90px",
                                          targets = seq_along(ko_cols) + 1)),
                                   rowCallback = .LFC_ROW_CB),
                    class = "compact cell-border hover") |>
        DT::formatStyle("N_Assays",
                        background = DT::styleColorBar(c(0, length(assays)), "#d4e6f1"),
                        backgroundSize = "100% 80%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "center")
    })

    # ── UpSet plot ────────────────────────────────────────────────────────────
    output$cmp_upset <- renderPlot({
      genes_df <- cmp_genes_data(); req(genes_df)
      assays   <- input$cmp_assays; req(assays, length(assays) >= 2)
      dir      <- input$cmp_upset_dir

      sub <- if (dir == "all") {
        genes_df[genes_df$assay %in% assays & genes_df$DEG %in% c("up", "down"), ]
      } else {
        genes_df[genes_df$assay %in% assays & genes_df$DEG == dir, ]
      }
      if (nrow(sub) == 0 || length(unique(na.omit(sub$symbol))) == 0) {
        dir_lbl <- switch(dir, up = "up-regulated", down = "down-regulated", all = "")
        plot.new()
        text(0.5, 0.5, paste0("No ", if (nzchar(dir_lbl)) paste0(dir_lbl, " ") else "",
                               "DEGs found for selected assays."),
             cex = 1.2, col = "#6c757d")
        return(invisible())
      }

      all_syms <- sort(unique(na.omit(sub$symbol)))
      # Build binary membership matrix (assay × gene)
      mat <- as.data.frame(lapply(setNames(assays, assays), function(a) {
        as.integer(all_syms %in% sub$symbol[sub$assay == a])
      }))
      rownames(mat) <- all_syms
      # Sanitise column names for UpSetR (no special chars or leading digits)
      safe <- gsub("[^A-Za-z0-9_]", "_", names(mat))
      safe <- sub("^([0-9])", "X\\1", safe)
      safe <- make.unique(safe, sep = "_")
      names(mat) <- safe

      clr <- switch(dir, up = .HP_CLRS$up, down = .HP_CLRS$down, all = .HP_CLRS$all)
      UpSetR::upset(
        mat,
        sets           = rev(names(mat)),
        order.by       = "freq",
        main.bar.color = clr,
        sets.bar.color = clr,
        text.scale     = c(1.2, 1, 1, 1, 1, 0.9),
        keep.order     = TRUE,
        mb.ratio       = c(0.60, 0.40)
      )
    })

    # ── Consensus DEGs: intersection across all selected assays ──────────────
    output$cmp_consensus <- renderUI({
      genes_df <- cmp_genes_data(); req(genes_df)
      assays   <- input$cmp_assays; req(assays, length(assays) >= 2)
      n        <- length(assays)

      up_per <- lapply(assays, function(a)
        unique(na.omit(genes_df$symbol[genes_df$assay == a & genes_df$DEG == "up"])))
      dn_per <- lapply(assays, function(a)
        unique(na.omit(genes_df$symbol[genes_df$assay == a & genes_df$DEG == "down"])))

      up_all <- sort(Reduce(intersect, up_per))
      dn_all <- sort(Reduce(intersect, dn_per))

      chip <- function(g, bg) {
        tags$span(
          style = paste0("background-color:", bg, "; color:white; font-size:0.78rem;",
                         " font-weight:600; padding:0.28em 0.6em; border-radius:4px;",
                         " display:inline-block; margin:2px 3px 2px 0;"), g)
      }

      make_section <- function(genes, label, bg, border_col) {
        n_g <- length(genes)
        hdr <- tags$div(class = "fw-bold mb-1", style = paste0("color:", border_col, ";"),
                        label,
                        if (n_g > 0) tags$span(class = "fw-normal text-muted ms-1",
                                               paste0("(", n_g, " gene",
                                                      if (n_g == 1) "" else "s", ")")))
        body <- if (n_g == 0)
          tags$p(class = "text-muted small mb-0", "No genes shared across all assays.")
        else
          div(class = "d-flex flex-wrap", lapply(genes, chip, bg = bg))
        div(class = "mb-3 p-2 rounded",
            style = paste0("border-left: 4px solid ", border_col, "; background:#fafafa;"),
            hdr, body)
      }

      tagList(
        tags$p(class = "text-muted small mb-3",
               paste0("Genes regulated in the same direction across all ",
                      n, " selected assays.")),
        make_section(up_all, "\u2191 Up in all assays", "#E74C3C", "#C0392B"),
        make_section(dn_all, "\u2193 Down in all assays", "#5DADE2", "#2C7BB6")
      )
    })

    # =========================================================================
    # 6. DEG Pathway Overlap
    # =========================================================================

    # ── Populate DEG filter choices ───────────────────────────────────────────
    observe({
      si <- study_info_r()
      updateSelectInput(session, "deg_filter_gene",
                        choices  = sort(unique(na.omit(si$Gene))),
                        selected = character(0))
      updateSelectInput(session, "deg_filter_model",
                        choices  = sort(unique(na.omit(si$Model_system))),
                        selected = character(0))
      updateSelectInput(session, "deg_filter_kos",
                        choices  = sort(unique(na.omit(si$KO_strat))),
                        selected = character(0))
    })

    # ── Filtered assay list ───────────────────────────────────────────────────
    filtered_deg_assays <- reactive({
      si      <- study_info_r()
      si_uniq <- si[!duplicated(si$Assay), ]
      if (length(input$deg_filter_gene)  > 0)
        si_uniq <- si_uniq[si_uniq$Gene         %in% input$deg_filter_gene,  ]
      if (length(input$deg_filter_model) > 0)
        si_uniq <- si_uniq[si_uniq$Model_system %in% input$deg_filter_model, ]
      if (length(input$deg_filter_kos)   > 0)
        si_uniq <- si_uniq[si_uniq$KO_strat     %in% input$deg_filter_kos,   ]
      sort(si_uniq$Assay)
    })

    observeEvent(filtered_deg_assays(), {
      ch   <- filtered_deg_assays()
      kept <- if (isolate(input$deg_assay) %in% ch) isolate(input$deg_assay) else NULL
      updateSelectizeInput(session, "deg_assay",
                           choices = ch, selected = kept, server = TRUE)
    })

    # ── Load DEG overlap for chosen assay ─────────────────────────────────────
    deg_overlap_data <- reactive({
      assay   <- input$deg_assay; src <- input$deg_source; dir <- input$deg_dir
      top_n   <- input$deg_topn;  minsize <- input$deg_minsize
      req(assay, src, dir, top_n, minsize)
      db_src   <- .HP_DEG_SRC[[src]]; cat_type <- .HP_DEG_CATTYPE[[src]]
      n_col    <- paste0("n_", dir); pct_col <- paste0("pct_", dir)
      where_type <- if (!is.na(cat_type))
        sprintf(" AND category_type = '%s'", cat_type) else ""
      con <- con_r()
      tryCatch(
        dbGetQuery(con, sprintf(
          'SELECT category_name AS category, category_genes_n AS n_pathway,
                  "%s" AS n_overlap, "%s" AS pct,
                  overlap_up AS genes_up, overlap_down AS genes_down
           FROM deg_overlap."%s"
           WHERE source = \'%s\'%s AND category_genes_n >= %d
           ORDER BY "%s" DESC LIMIT %d',
          n_col, pct_col, assay, db_src, where_type, minsize, pct_col, top_n)),
        error = function(e) NULL)
    })

    output$deg_bar <- renderPlotly({
      df <- deg_overlap_data(); req(df, nrow(df) > 0)
      df$pct <- round(df$pct, 1)
      df$category <- ifelse(nchar(df$category) > 55,
                            paste0(substr(df$category, 1, 52), "\u2026"), df$category)
      df$category <- factor(df$category, levels = rev(df$category))
      dir_clr   <- switch(input$deg_dir, updn=.HP_CLRS$all, up=.HP_CLRS$up, down=.HP_CLRS$down)
      src_label <- names(.HP_SOURCES)[.HP_SOURCES == input$deg_source]
      plot_ly(df, x = ~pct, y = ~category, type = "bar", orientation = "h",
              marker = list(color = dir_clr),
              hovertemplate = "<b>%{y}</b><br>Overlap: %{x}%<extra></extra>") |>
        layout(xaxis = list(title = "% DEGs in Category"),
               yaxis = list(title = "", automargin = TRUE),
               margin = list(l = 5, r = 5, t = 35, b = 5),
               title = list(text = paste0("DEG Overlap \u2014 ", src_label),
                            font = list(size = 13))) |>
        config(displayModeBar = FALSE)
    })

    output$deg_pathway_tbl <- DT::renderDT({
      df <- deg_overlap_data(); req(df, nrow(df) > 0)
      df$pct       <- round(df$pct, 1)
      df$n_overlap <- as.integer(df$n_overlap)
      df$n_pathway <- as.integer(df$n_pathway)
      count_genes  <- function(x) {
        if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(0L)
        length(trimws(strsplit(x, ";")[[1]]))
      }
      df$n_up   <- vapply(df$genes_up,   count_genes, integer(1))
      df$n_down <- vapply(df$genes_down, count_genes, integer(1))
      tbl       <- df[, c("category", "n_pathway", "n_overlap", "pct", "n_up", "n_down")]
      names(tbl) <- c("Category", "Cat. Size", "N Overlap", "% Overlap", "N Up", "N Down")
      DT::datatable(tbl, rownames = FALSE, filter = "top", selection = "single",
                    options = list(pageLength = 15, scrollX = TRUE,
                                   scrollY = "300px", scrollCollapse = TRUE,
                                   dom = "frtip"),
                    class = "compact row-border hover") |>
        DT::formatStyle("% Overlap",
                        background = DT::styleColorBar(c(0, max(df$pct, na.rm = TRUE)), "#b3cde3"),
                        backgroundSize = "100% 60%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "center")
    })

    # ── Track selected pathway row ────────────────────────────────────────────
    selected_pathway_row <- reactiveVal(NULL)

    observeEvent(input$deg_pathway_tbl_rows_selected, {
      idx <- input$deg_pathway_tbl_rows_selected
      df  <- deg_overlap_data()
      if (length(idx) == 1 && !is.null(df) && idx <= nrow(df))
        selected_pathway_row(df[idx, ])
      else
        selected_pathway_row(NULL)
    })
    observeEvent(list(input$deg_assay, input$deg_source, input$deg_dir,
                      input$deg_topn, input$deg_minsize), {
      selected_pathway_row(NULL)
    }, ignoreInit = TRUE)

    # ── Fetch logFC for genes in selected pathway + grey gene list ────────────
    deg_pathway_logfc <- reactive({
      row   <- selected_pathway_row(); assay <- input$deg_assay
      src   <- input$deg_source
      req(row, assay, src)

      parse_genes <- function(x) {
        if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(character(0))
        v <- trimws(strsplit(x, ";")[[1]]); v[nzchar(v)]
      }
      up_genes  <- parse_genes(row$genes_up[1])
      dn_genes  <- parse_genes(row$genes_down[1])
      all_genes <- unique(c(up_genes, dn_genes))
      if (length(all_genes) == 0) return(NULL)

      con       <- con_r()
      col_names <- tryCatch(dbListFields(con, assay), error = function(e) character(0))
      lfc_col   <- col_names[grepl("log2FoldChange$", col_names)]

      lfc_map <- NULL
      if (length(lfc_col) > 0) {
        gene_sql <- paste0("'", paste(gsub("'", "''", all_genes), collapse="','"), "'")
        lfc_df   <- tryCatch(
          dbGetQuery(con, sprintf(
            'SELECT hgnc_symbol, "%s" AS log2fc FROM main."%s"
             WHERE hgnc_symbol IN (%s)',
            lfc_col[1], assay, gene_sql)),
          error = function(e) NULL)
        if (!is.null(lfc_df) && nrow(lfc_df) > 0)
          lfc_map <- setNames(lfc_df$log2fc, lfc_df$hgnc_symbol)
      }

      # Try to get all pathway genes from morphic_overlap for grey chip display
      grey_genes <- character(0)
      tryCatch({
        cat_col  <- .HP_GRP[[src]]
        cat_name <- row$category[1]
        mo_cols  <- dbGetQuery(con, sprintf(
          "SELECT column_name FROM information_schema.columns
           WHERE table_schema = 'morphic_overlap' AND table_name = '%s'",
          src))$column_name
        # Find columns that look like gene-name lists (not counts, not the KO overlap list)
        gene_list_cols <- mo_cols[
          grepl("gene", mo_cols, ignore.case = TRUE) &
          !grepl("_n$|^overlap_genes$|^overlap_n$", mo_cols) &
          mo_cols != .HP_N[[src]]
        ]
        if (length(gene_list_cols) > 0) {
          safe_cat <- gsub("'", "''", cat_name)
          pg_df <- dbGetQuery(con, sprintf(
            'SELECT "%s" AS gene_list FROM morphic_overlap."%s"
             WHERE "%s" = \'%s\' LIMIT 1',
            gene_list_cols[1], src, cat_col, safe_cat))
          if (nrow(pg_df) > 0 && !is.na(pg_df$gene_list[1]) &&
              nzchar(trimws(as.character(pg_df$gene_list[1]))))
            grey_genes <- setdiff(parse_genes(as.character(pg_df$gene_list[1])), all_genes)
        }
      }, error = function(e) NULL)

      list(up = up_genes, down = dn_genes, lfc_map = lfc_map, grey = grey_genes)
    })

    # ── Gene chip panel: Up | Down sections, sorted by |logFC| ───────────────
    output$deg_gene_panel <- renderUI({
      row <- selected_pathway_row()
      if (is.null(row)) return(NULL)

      info <- deg_pathway_logfc()
      req(info)

      up_genes   <- info$up; dn_genes <- info$down; lfc_map <- info$lfc_map
      grey_genes <- if (!is.null(info$grey)) info$grey else character(0)

      make_chips <- function(genes, is_up) {
        if (length(genes) == 0) return(NULL)
        # Attach logFC and sort by magnitude desc
        lfc_vals <- if (!is.null(lfc_map)) lfc_map[genes] else rep(NA_real_, length(genes))
        lfc_vals[is.na(lfc_vals)] <- 0
        ord   <- if (is_up) order(-lfc_vals) else order(lfc_vals)   # up: high→low; dn: low→high
        genes <- genes[ord]; lfc_vals <- lfc_vals[ord]

        lapply(seq_along(genes), function(i) {
          g       <- genes[i]
          lfc     <- lfc_vals[i]
          lfc_abs <- abs(lfc)
          bg  <- .hp_chip_colour(is_up, lfc_abs)
          arrow   <- if (is_up) "\u2191" else "\u2193"
          lfc_txt <- if (!is.na(lfc) && lfc != 0)
            paste0(" (", ifelse(lfc >= 0, "+", ""), round(lfc, 2), ")") else ""
          tags$span(style = .hp_chip_style(bg), paste0(g, " ", arrow, lfc_txt))
        })
      }

      up_chips <- make_chips(up_genes, TRUE)
      dn_chips <- make_chips(dn_genes, FALSE)
      n_up     <- length(up_genes); n_dn <- length(dn_genes)
      n_total  <- as.integer(row$n_pathway[1])

      # ── Grey chips: remaining pathway genes (not DEGs) ────────────────────
      grey_section <- if (length(grey_genes) > 0) {
        grey_chips <- lapply(sort(grey_genes), function(g) {
          tags$span(
            style = paste0("background-color:#adb5bd; color:white; font-size:0.75rem;",
                           " padding:0.2em 0.55em; border-radius:4px;",
                           " display:inline-block; margin:2px 3px 2px 0;"),
            g)
        })
        div(
          class = "mb-3 p-2 rounded",
          style = "border-left: 4px solid #adb5bd; background:#fafafa;",
          tags$div(class = "fw-bold mb-1", style = "color:#6c757d;",
                   "Other pathway genes",
                   tags$span(class = "fw-normal text-muted ms-1",
                             paste0("(", length(grey_genes), " gene",
                                    if (length(grey_genes) == 1L) "" else "s", ")"))),
          div(class = "d-flex flex-wrap", grey_chips)
        )
      } else NULL

      section <- function(label, n, chips, border_col) {
        if (n == 0) return(NULL)
        div(
          class = "mb-3 p-2 rounded",
          style = paste0("border-left: 4px solid ", border_col, "; background:#fafafa;"),
          tags$div(class = "fw-bold mb-1", style = paste0("color:", border_col, ";"),
                   label, tags$span(class = "fw-normal text-muted ms-1",
                                    paste0("(", n, " gene", if (n == 1) "" else "s",
                                           " — sorted by |logFC|)"))),
          div(class = "d-flex flex-wrap", chips)
        )
      }

      legend <- div(
        class = "text-muted small mb-2",
        "Colour intensity = |logFC|:",
        tags$span(style = .hp_chip_style("#F1948A"), "light < 1"),
        tags$span(style = .hp_chip_style("#E74C3C"), "medium 1\u20132"),
        tags$span(style = .hp_chip_style("#C0392B"), "dark \u22652"),
        " | ",
        tags$span(style = .hp_chip_style("#85C1E9"), "light < 1"),
        tags$span(style = .hp_chip_style("#5DADE2"), "medium 1\u20132"),
        tags$span(style = .hp_chip_style("#2C7BB6"), "dark \u22652")
      )

      tagList(
        tags$hr(),
        div(class = "px-2 pb-2",
            tags$h6(class = "fw-bold mb-1",
                    paste0(row$category[1], " \u2014 ",
                           n_up + n_dn, " DEGs / ", n_total, " total in pathway")),
            legend,
            section("\u2191 Up-regulated",   n_up, up_chips, "#C0392B"),
            section("\u2193 Down-regulated", n_dn, dn_chips, "#2C7BB6"),
            grey_section
        )
      )
    })

    # =========================================================================
    # Suspend all outputs until Home tab is visited
    # =========================================================================
    invisible(lapply(
      c("breakdown_bar", "assays_per_gene",
        "assay_tbl",
        "ko_bar", "ko_tbl",
        "cmp_meta_above", "cmp_matrix_tbl", "cmp_upset", "cmp_consensus",
        "deg_bar", "deg_pathway_tbl", "deg_gene_panel"),
      function(oid) outputOptions(output, oid, suspendWhenHidden = TRUE)
    ))

  })
}
