# =============================================================================
# Module: studies_overview
# Page 1 — Perturbed Genes Overview
#
# Sections:
#   1. Top card with tabs: Studies Overview | Study Breakdown | N Assays/Gene
#   2. Perturbation Coverage & Annotation — one card, tabs:
#        Reactome | GO BP | GO MF | PANTHER Class
#        Phenotype Associations (IMPC, HPO)
#        Disease Associations   (OMIM, Orphanet)
#
# Requires helpers.R to be sourced first (.HP_SOURCES, .HP_GRP, .HP_N).
# Data contract:
#   con_r        — reactive: open read-only DuckDB connection
#   study_info_r — reactive: study_info data.frame
# =============================================================================

library(plotly)
library(DT)
library(dplyr)

# ── UI ────────────────────────────────────────────────────────────────────────

studies_overviewUI <- function(id) {
  ns <- NS(id)

  tagList(

    # ── 1. Studies Overview / Study Breakdown / N Assays per Gene (tabs) ───
    card(
      card_header(
        tagList(bsicons::bs_icon("table"), " Studies Overview")
      ),
      navset_tab(

        # ── Tab: Studies Overview table ──────────────────────────────────
        nav_panel(
          tagList("Studies Overview",
                  .info_tip("Browse all MorPhiC assays grouped by KO gene. Use the gene filter to narrow results.")),
          div(
            class = "px-3 pt-3 pb-2",

            # ── Filter bar ──────────────────────────────────────────────
            div(
              class = "d-flex align-items-end gap-3 mb-3 pb-3",
              style = "border-bottom: 1px solid #e9ecef;",

              # Gene filter
              div(
                style = "flex: 0 0 380px;",
                tags$label(
                  class = "form-label small fw-semibold text-secondary mb-1",
                  style = "letter-spacing: 0.02em;",
                  bsicons::bs_icon("funnel"), " Filter by KO Gene"
                ),
                selectizeInput(ns("tbl_gene_filter"), NULL,
                               choices = NULL, multiple = TRUE, width = "100%",
                               options = list(
                                 placeholder = "Type to search genes\u2026",
                                 plugins     = list("remove_button")
                               ))
              ),

              # Summary counts
              uiOutput(ns("tbl_summary"), inline = TRUE)
            ),

            # ── Table ───────────────────────────────────────────────────
            DTOutput(ns("assay_tbl"))
          )
        ),

        # ── Tab: Study Breakdown ─────────────────────────────────────────
        nav_panel(
          tagList("Study Breakdown",
                  .info_tip("Grouped bar chart showing N assays and N unique KO genes by selected attribute.")),
          div(
            class = "p-3",
            div(
              class = "d-flex align-items-center gap-2 pb-2",
              tags$label("Group by:", class = "mb-0 small fw-semibold"),
              selectInput(ns("group_by"), NULL,
                          choices  = c("Model System" = "Model_system",
                                       "KO Strategy"  = "KO_strat",
                                       "Cell Line"    = "Cell_Line",
                                       "DPC"          = "DPC"),
                          selected = "Model_system", width = "160px")
            ),
            plotlyOutput(ns("breakdown_bar"), height = "420px")
          )
        ),

        # ── Tab: N Assays per Gene ───────────────────────────────────────
        nav_panel(
          tagList("N Assays per Gene",
                  .info_tip("Number of assays available per KO gene across all studies.")),
          div(class = "p-3",
              plotlyOutput(ns("assays_per_gene"), height = "480px")))
      ),
      full_screen = TRUE
    ),

    # ── 2. Perturbation Coverage & Annotation ─────────────────────────────
    card(
      card_header(tagList(bsicons::bs_icon("intersect"),
                          " Perturbation Coverage & Annotation",
                          .info_tip("Percentage of MorPhiC KO genes found in pathway, phenotype, and disease annotation databases."))),
      navset_tab(

        # ── Reactome ─────────────────────────────────────────────────────────
        nav_panel(
          tagList("Reactome",
                  .info_tip("Overlap between MorPhiC KO genes and Reactome pathway gene sets.")),
          navset_tab(
            nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_reactome"),      height = "400px")),
            nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_reactome"))))
          )
        ),

        # ── GO Biological Process ─────────────────────────────────────────────
        nav_panel(
          tagList("GO Biological Process",
                  .info_tip("Overlap between MorPhiC KO genes and GO Biological Process terms.")),
          navset_tab(
            nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_go_bp"),         height = "400px")),
            nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_go_bp"))))
          )
        ),

        # ── GO Molecular Function ─────────────────────────────────────────────
        nav_panel(
          tagList("GO Molecular Function",
                  .info_tip("Overlap between MorPhiC KO genes and GO Molecular Function terms.")),
          navset_tab(
            nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_go_mf"),         height = "400px")),
            nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_go_mf"))))
          )
        ),

        # ── PANTHER Protein Class ─────────────────────────────────────────────
        nav_panel(
          tagList("PANTHER Protein Class",
                  .info_tip("Overlap between MorPhiC KO genes and PANTHER protein class annotations.")),
          navset_tab(
            nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_panther_class"), height = "400px")),
            nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_panther_class"))))
          )
        ),

        # ── Phenotype Associations ────────────────────────────────────────────
        nav_panel(
          tagList("Phenotype Associations",
                  .info_tip("Phenotype annotations for MorPhiC KO genes from IMPC mouse and HPO human databases.")),
          navset_tab(
            nav_panel(
              "IMPC Mouse Phenotypes",
              navset_tab(
                nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_impc"), height = "480px")),
                nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_impc"))))
              )
            ),
            nav_panel(
              "HPO Human Phenotypes",
              navset_tab(
                nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_hpo"),  height = "480px")),
                nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_hpo"))))
              )
            )
          )
        ),

        # ── Disease Associations ──────────────────────────────────────────────
        nav_panel(
          tagList("Disease Associations",
                  .info_tip("Disease annotations for MorPhiC KO genes from OMIM and Orphanet databases.")),
          navset_tab(
            nav_panel(
              "OMIM",
              navset_tab(
                nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_omim"),     height = "480px")),
                nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_omim"))))
              )
            ),
            nav_panel(
              "Orphanet",
              navset_tab(
                nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_orphanet"), height = "480px")),
                nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_orphanet"))))
              )
            )
          )
        )

      ), # end outer navset_tab
      full_screen = TRUE
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

studies_overviewServer <- function(id, con_r, study_info_r) {
  moduleServer(id, function(input, output, session) {

    # =========================================================================
    # 1. Studies Overview Table (RowGroup by Gene)
    # =========================================================================

    # Populate gene filter choices
    observe({
      si <- study_info_r()
      genes <- sort(unique(na.omit(si$Gene)))
      updateSelectizeInput(session, "tbl_gene_filter",
                           choices = genes, selected = character(0))
    })

    # Reactive: filtered table data (shared by table + summary)
    tbl_data_r <- reactive({
      si   <- study_info_r()
      want <- c("Gene", "Study_title", "Comparison", "Model_system",
                "Cell_Line", "Differentation_time_point", "Condition_details")
      cols <- want[want %in% colnames(si)]
      tbl  <- si[!duplicated(si$Assay), cols, drop = FALSE]

      sel_genes <- input$tbl_gene_filter
      if (length(sel_genes) > 0)
        tbl <- tbl[tbl$Gene %in% sel_genes, ]

      tbl[order(tbl$Gene, tbl$Model_system), ]
    })

    # Summary badges
    output$tbl_summary <- renderUI({
      tbl <- tbl_data_r()
      n_genes  <- length(unique(na.omit(tbl$Gene)))
      n_assays <- nrow(tbl)
      div(
        class = "d-flex gap-2 align-items-center mb-2",
        style = "padding-bottom: 2px;",
        tags$span(
          class = "badge rounded-pill",
          style = "background: #4E79A7; font-size: 0.78rem; font-weight: 500; padding: 6px 12px;",
          paste(n_genes, "genes")
        ),
        tags$span(
          class = "badge rounded-pill",
          style = "background: #6c757d; font-size: 0.78rem; font-weight: 500; padding: 6px 12px;",
          paste(n_assays, "assays")
        )
      )
    })

    output$assay_tbl <- DT::renderDT({
      tbl <- tbl_data_r()

      gene_idx <- 0L

      pretty <- c(Gene = "KO Gene", Study_title = "Study Title",
                  Comparison = "Comparison", Model_system = "Model System",
                  Cell_Line = "Cell Line",
                  Differentation_time_point = "Diff. Time Point",
                  Condition_details = "Condition Details")
      colnames(tbl) <- pretty[colnames(tbl)]

      # Truncate long study titles for display
      tbl[["Study Title"]] <- ifelse(
        nchar(tbl[["Study Title"]]) > 80,
        paste0(substr(tbl[["Study Title"]], 1, 77), "\u2026"),
        tbl[["Study Title"]]
      )

      DT::datatable(
        tbl,
        rownames   = FALSE,
        selection  = "none",
        extensions = "RowGroup",
        options    = list(
          pageLength     = 167,
          scrollX        = TRUE,
          scrollY        = "520px",
          scrollCollapse = TRUE,
          paging         = FALSE,
          dom            = "rt",
          rowGroup       = list(dataSrc = gene_idx),
          columnDefs     = list(
            list(visible = FALSE, targets = gene_idx),
            list(className = "dt-left", targets = "_all"),
            list(width = "220px", targets = 1),
            list(width = "130px", targets = 2)
          ),
          initComplete = DT::JS("
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
        ),
        class = "row-border hover nowrap"
      ) |>
        DT::formatStyle(
          "Study Title",
          fontSize = "0.82em", color = "#495057"
        ) |>
        DT::formatStyle(
          "Condition Details",
          fontSize = "0.8em", color = "#868e96", fontStyle = "italic"
        ) |>
        DT::formatStyle(
          "Comparison",
          fontWeight = "500", color = "#343a40"
        ) |>
        DT::formatStyle(
          "Model System",
          color = "#0d6efd", fontWeight = "500", fontSize = "0.85em"
        )
    })

    # =========================================================================
    # 2. Study Breakdown
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
      grp_lbl <- switch(grp_col, Model_system = "Model System",
                        KO_strat = "KO Strategy", Cell_Line = "Cell Line", DPC = "DPC")
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
    # 3. N Assays per KO Gene
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
    # 4. Perturbation Coverage & Annotation
    # =========================================================================

    # ── 4a. Pathway / GO / PANTHER  (on-the-fly via ref tables) ──────────────

    # Compute KO gene overlap on-the-fly via reference tables
    .ko_load <- function(src) {
      ref_tbl <- .HP_REF_TBL[[src]]
      reactive({
        con <- con_r()
        compute_ko_overlap(con, ref_tbl)
      })
    }

    ko_data_reactome      <- .ko_load("reactome")
    ko_data_go_bp         <- .ko_load("go_bp")
    ko_data_go_mf         <- .ko_load("go_mf")
    ko_data_panther_class <- .ko_load("panther_class")

    # Bar chart: top 10 rows, % KO genes in category
    .ko_render_bar <- function(df_r, src_label, bar_colour = "#B07AA1") {
      renderPlotly({
        df <- df_r(); req(df, nrow(df) > 0)
        top <- head(df, 10)
        top$pct <- round(top$pct, 1)
        top$category <- ifelse(nchar(top$category) > 55,
                               paste0(substr(top$category, 1, 52), "\u2026"),
                               top$category)
        top$category <- factor(top$category, levels = rev(top$category))
        plot_ly(top,
                x             = ~pct,
                y             = ~category,
                type          = "bar",
                orientation   = "h",
                marker        = list(color = bar_colour),
                customdata    = ~n_pathway,
                hovertemplate = paste0(
                  "<b>%{y}</b><br>",
                  "Overlap: %{x}%<br>",
                  "Category size: %{customdata} genes<extra></extra>")) |>
          layout(
            xaxis  = list(title = "% MorPhiC KO Genes in Category", range = c(0, 105)),
            yaxis  = list(title = "", automargin = TRUE),
            margin = list(l = 5, r = 5, t = 35, b = 5),
            title  = list(text = paste0("Top 10 \u2014 ", src_label),
                          font = list(size = 13))) |>
          config(displayModeBar = FALSE)
      })
    }

    # Table: all rows
    .ko_render_tbl <- function(df_r) {
      DT::renderDT({
        df <- df_r(); req(df, nrow(df) > 0)
        df$pct        <- round(df$pct, 1)
        df$n_pathway  <- as.integer(df$n_pathway)
        df$overlap_n  <- as.integer(df$overlap_n)
        tbl <- df[, c("category", "n_pathway", "overlap_n", "pct", "genes")]
        names(tbl) <- c("Category", "Category Size", "KO Overlap N", "% Overlap",
                        "Overlapping Genes")
        DT::datatable(tbl, rownames = FALSE, filter = "top",
                      extensions = "Buttons",
                      options = list(pageLength = 25, scrollX = TRUE,
                                     scrollY = "420px", scrollCollapse = TRUE,
                                     dom = "Bfrtip", buttons = list("csv", "excel")),
                      class = "compact row-border hover")
      })
    }

    # Wire up pathway sources
    output$ko_bar_reactome     <- .ko_render_bar(ko_data_reactome,     "Reactome Pathways",       "#4E79A7")
    output$ko_tbl_reactome     <- .ko_render_tbl(ko_data_reactome)
    output$ko_bar_go_bp        <- .ko_render_bar(ko_data_go_bp,        "GO Biological Process",   "#59A14F")
    output$ko_tbl_go_bp        <- .ko_render_tbl(ko_data_go_bp)
    output$ko_bar_go_mf        <- .ko_render_bar(ko_data_go_mf,        "GO Molecular Function",   "#F28E2B")
    output$ko_tbl_go_mf        <- .ko_render_tbl(ko_data_go_mf)
    output$ko_bar_panther_class <- .ko_render_bar(ko_data_panther_class, "PANTHER Protein Class",  "#B07AA1")
    output$ko_tbl_panther_class <- .ko_render_tbl(ko_data_panther_class)

    # ── 4b. Phenotype / Disease annotation tables (term-centric) ─────────────

    # Query distinct gene-annotation pairs from a flat annotation table.
    # Returns data.frame with columns: gene, annotation
    .pheno_load <- function(tbl_name, ph_col) {
      reactive({
        con <- con_r()
        tryCatch(
          dbGetQuery(con, sprintf(
            'SELECT DISTINCT s.Gene AS gene, p."%s" AS annotation
             FROM "%s" p
             INNER JOIN (SELECT DISTINCT Gene FROM study_info WHERE Gene IS NOT NULL) s
               ON p.symbol = s.Gene
             WHERE p."%s" IS NOT NULL',
            ph_col, tbl_name, ph_col)),
          error = function(e) NULL)
      })
    }

    ko_data_impc     <- .pheno_load("impc_phenotypes", "impc_phenotypes")
    ko_data_hpo      <- .pheno_load("hpo_phenotypes",  "hpo_name")
    ko_data_omim     <- .pheno_load("omim_phenotypes", "omim_phenotype")
    ko_data_orphanet <- .pheno_load("orphanet_data",   "disorder_name")

    # Bar chart: top 10 annotation terms by N KO genes (term-centric)
    .ko_render_pheno_bar <- function(df_r, src_label, bar_colour) {
      renderPlotly({
        df <- df_r(); req(df, nrow(df) > 0)
        agg <- df |>
          dplyr::group_by(annotation) |>
          dplyr::summarise(n_genes = dplyr::n_distinct(gene), .groups = "drop") |>
          dplyr::arrange(dplyr::desc(n_genes)) |>
          head(10)
        # Truncate long term names for display
        agg$annotation <- ifelse(nchar(agg$annotation) > 55,
                                 paste0(substr(agg$annotation, 1, 52), "\u2026"),
                                 agg$annotation)
        agg <- agg[order(agg$n_genes), ]
        agg$annotation <- factor(agg$annotation, levels = agg$annotation)
        plot_ly(agg, x = ~n_genes, y = ~annotation, type = "bar", orientation = "h",
                marker = list(color = bar_colour),
                hovertemplate = paste0(
                  "<b>%{y}</b><br>",
                  "KO Genes: %{x}<extra></extra>")) |>
          layout(
            xaxis  = list(title = "N KO Genes", dtick = 1),
            yaxis  = list(title = "", automargin = TRUE),
            margin = list(l = 5, r = 5, t = 35, b = 5),
            title  = list(text = paste0("Top 10 \u2014 ", src_label),
                          font = list(size = 13))) |>
          config(displayModeBar = FALSE)
      })
    }

    # Table: all annotation terms with N KO genes and gene list (term-centric)
    .ko_render_pheno_tbl <- function(df_r, ann_col_label) {
      DT::renderDT({
        df <- df_r(); req(df, nrow(df) > 0)
        tbl <- df |>
          dplyr::group_by(annotation) |>
          dplyr::summarise(
            n_genes  = dplyr::n_distinct(gene),
            gene_list = paste(sort(unique(gene)), collapse = "; "),
            .groups = "drop") |>
          dplyr::arrange(dplyr::desc(n_genes)) |>
          as.data.frame()
        names(tbl) <- c(ann_col_label, "N KO Genes", "KO Genes")
        DT::datatable(tbl, rownames = FALSE, filter = "top",
                      extensions = "Buttons",
                      options = list(pageLength = 25, scrollX = TRUE,
                                     scrollY = "420px", scrollCollapse = TRUE,
                                     dom = "Bfrtip", buttons = list("csv", "excel")),
                      class = "compact row-border hover")
      })
    }

    # Wire up phenotype sources
    output$ko_bar_impc     <- .ko_render_pheno_bar(ko_data_impc,     "IMPC Phenotype",   "#E63946")
    output$ko_tbl_impc     <- .ko_render_pheno_tbl(ko_data_impc,     "IMPC Phenotype")
    output$ko_bar_hpo      <- .ko_render_pheno_bar(ko_data_hpo,      "HPO Term",         "#D62728")
    output$ko_tbl_hpo      <- .ko_render_pheno_tbl(ko_data_hpo,      "HPO Term")

    # Wire up disease sources
    output$ko_bar_omim     <- .ko_render_pheno_bar(ko_data_omim,     "OMIM Disease",     "#457B9D")
    output$ko_tbl_omim     <- .ko_render_pheno_tbl(ko_data_omim,     "OMIM Phenotype")
    output$ko_bar_orphanet <- .ko_render_pheno_bar(ko_data_orphanet, "Orphanet Disorder","#1D6FA4")
    output$ko_tbl_orphanet <- .ko_render_pheno_tbl(ko_data_orphanet, "Disorder")

    # ── Suspend all outputs until visible ─────────────────────────────────────
    invisible(lapply(
      c("assay_tbl", "tbl_summary", "breakdown_bar", "assays_per_gene",
        "ko_bar_reactome",     "ko_tbl_reactome",
        "ko_bar_go_bp",        "ko_tbl_go_bp",
        "ko_bar_go_mf",        "ko_tbl_go_mf",
        "ko_bar_panther_class","ko_tbl_panther_class",
        "ko_bar_impc",         "ko_tbl_impc",
        "ko_bar_hpo",          "ko_tbl_hpo",
        "ko_bar_omim",         "ko_tbl_omim",
        "ko_bar_orphanet",     "ko_tbl_orphanet"),
      function(oid) outputOptions(output, oid, suspendWhenHidden = TRUE)
    ))

  })
}
