# =============================================================================
# Module: studies_overview
# Page 1 — Perturbed Genes Overview
#
# Sections:
#   1. Top card: Studies Overview table + Charts modal (4 pie charts)
#   2. Functional Annotations — Reactome | GO BP | GO MF | PANTHER
#   3. Disease & Phenotype Associations — IMPC | HPO | OMIM | Orphanet
#
# Requires helpers.R to be sourced first (.HP_REF_TBL, .dt_header_js, etc.).
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

    # ── 1. Studies Overview ──────────────────────────────────────────────
    card(
      card_header(
        class = "d-flex align-items-center gap-2 flex-wrap",
        tags$span(class = "fw-semibold",
                  tagList(bsicons::bs_icon("table"), " Studies Overview")),

        # Filters dropdown
        shinyWidgets::dropdownButton(
          inputId = ns("tbl_filters"),
          label   = "Filters",
          icon    = tagList(bsicons::bs_icon("funnel", size = "0.75rem"),
                            uiOutput(ns("tbl_filter_badge"), inline = TRUE)),
          circle  = FALSE, status = "default", size = "sm",
          width   = "520px", inline = TRUE,

          # Clear button
          actionButton(ns("tbl_clear_filters"), "Clear all filters",
                       class = "btn btn-outline-secondary btn-sm w-100 mb-2"),

          # Two-column grid
          div(
            style = "display: grid; grid-template-columns: 1fr 1fr; gap: 0 12px;",
            selectizeInput(ns("tbl_filter_gene"), "KO Gene",
                           choices = NULL, multiple = TRUE, width = "100%",
                           options = list(placeholder = "All", plugins = list("remove_button"))),
            selectizeInput(ns("tbl_filter_study"), "Study Title",
                           choices = NULL, multiple = TRUE, width = "100%",
                           options = list(placeholder = "All", plugins = list("remove_button"))),
            selectizeInput(ns("tbl_filter_comparison"), "Comparison",
                           choices = NULL, multiple = TRUE, width = "100%",
                           options = list(placeholder = "All", plugins = list("remove_button"))),
            selectizeInput(ns("tbl_filter_model"), "Model System",
                           choices = NULL, multiple = TRUE, width = "100%",
                           options = list(placeholder = "All", plugins = list("remove_button"))),
            selectizeInput(ns("tbl_filter_cell_line"), "Cell Line",
                           choices = NULL, multiple = TRUE, width = "100%",
                           options = list(placeholder = "All", plugins = list("remove_button"))),
            selectizeInput(ns("tbl_filter_dpc"), "DPC",
                           choices = NULL, multiple = TRUE, width = "100%",
                           options = list(placeholder = "All", plugins = list("remove_button"))),
            selectizeInput(ns("tbl_filter_diff"), "Diff. Time Point",
                           choices = NULL, multiple = TRUE, width = "100%",
                           options = list(placeholder = "All", plugins = list("remove_button"))),
            selectizeInput(ns("tbl_filter_condition"), "Condition Details",
                           choices = NULL, multiple = TRUE, width = "100%",
                           options = list(placeholder = "All", plugins = list("remove_button")))
          )
        ),

        # Summary counts
        uiOutput(ns("tbl_summary"), inline = TRUE),

        # Spacer
        div(style = "flex: 1;"),

        # Charts button (opens modal)
        actionButton(ns("show_charts"),
                     tagList(bsicons::bs_icon("bar-chart-fill", size = "0.75rem"), " Charts"),
                     class = "btn btn-outline-secondary btn-sm")
      ),
      card_body(
        class = "px-3 pt-2 pb-2",
        style = "min-height: 1200px;",

        # ── Table ───────────────────────────────────────────────────
        DTOutput(ns("assay_tbl"))
      ),
      full_screen = TRUE
    ),

    # ── 2. Gene Annotations ─────────────────────────────────────────────────
    navset_card_tab(
      title = tagList(bsicons::bs_icon("intersect"),
                      " Gene Annotations",
                      .info_tip("Pathway, protein, phenotype and disease annotations for MorPhiC KO genes.")),
      full_screen = TRUE,

      # ── Functional Annotations ─────────────────────────────────────────
      nav_menu(
        title = "Functional Annotations",
        icon  = bsicons::bs_icon("intersect"),

        nav_panel(
          tagList("Pathways",
                  .info_tip("Overlap between MorPhiC KO genes and Reactome pathway gene sets.")),
          div(class = "py-3 px-1",
              card(navset_tab(
                nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_reactome"), height = "400px")),
                nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_reactome"))))
              )))
        ),

        nav_panel(
          tagList("Biological Processes",
                  .info_tip("Overlap between MorPhiC KO genes and GO Biological Process terms.")),
          div(class = "py-3 px-1",
              card(navset_tab(
                nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_go_bp"), height = "400px")),
                nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_go_bp"))))
              )))
        ),

        nav_panel(
          tagList("Molecular Functions",
                  .info_tip("Overlap between MorPhiC KO genes and GO Molecular Function terms.")),
          div(class = "py-3 px-1",
              card(navset_tab(
                nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_go_mf"), height = "400px")),
                nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_go_mf"))))
              )))
        ),

        nav_panel(
          tagList("Protein Class",
                  .info_tip("Overlap between MorPhiC KO genes and PANTHER protein class annotations.")),
          div(class = "py-3 px-1",
              card(navset_tab(
                nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_panther_class"), height = "400px")),
                nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_panther_class"))))
              )))
        )

      ), # end nav_menu "Functional Annotations"

      # ── Disease & Phenotype Associations ───────────────────────────────
      nav_menu(
        title = "Disease & Phenotypes",
        icon  = bsicons::bs_icon("card-list"),

        nav_panel(
          tagList("Mouse Phenotypes",
                  .info_tip("Mouse phenotype annotations from the International Mouse Phenotyping Consortium.")),
          div(class = "py-3 px-1",
              card(navset_tab(
                nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_impc"), height = "480px")),
                nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_impc"))))
              )))
        ),

        nav_panel(
          tagList("Human Phenotypes",
                  .info_tip("Human phenotype annotations from the Human Phenotype Ontology.")),
          div(class = "py-3 px-1",
              card(navset_tab(
                nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_hpo"), height = "480px")),
                nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_hpo"))))
              )))
        ),

        nav_panel(
          tagList("Disease Associations",
                  .info_tip("Disease annotations from OMIM and Orphanet databases.")),
          div(class = "py-3 px-1",
            card(navset_tab(
              nav_panel("Bar Chart", plotlyOutput(ns("ko_bar_disease"), height = "480px")),
              nav_panel("Table",     div(class = "p-2", DTOutput(ns("ko_tbl_disease")))),
              nav_spacer(),
              nav_item(
                shinyWidgets::radioGroupButtons(
                  inputId  = ns("disease_source"),
                  label    = NULL,
                  choices  = c("OMIM" = "omim", "Orphanet" = "orphanet"),
                  selected = "omim",
                  size     = "sm",
                  status   = "outline-primary"
                )
              )
            ))
          )
        )


      ) # end nav_menu "Disease & Phenotype"
    ) # end navset_card_tab
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

studies_overviewServer <- function(id, con_r, study_info_r) {
  moduleServer(id, function(input, output, session) {

    # =========================================================================
    # 1. Studies Overview Table (RowGroup by Gene)
    # =========================================================================

    # Populate filter choices
    observe({
      si <- study_info_r()
      .upd <- function(id, col) {
        ch <- sort(unique(na.omit(si[[col]])))
        updateSelectizeInput(session, id, choices = ch, selected = character(0))
      }
      .upd("tbl_filter_gene",       "Gene")
      .upd("tbl_filter_study",      "Study_title")
      .upd("tbl_filter_comparison", "Comparison")
      .upd("tbl_filter_model",      "Model_system")
      .upd("tbl_filter_cell_line",  "Cell_Line")
      .upd("tbl_filter_dpc",        "DPC")
      .upd("tbl_filter_diff",       "Differentation_time_point")
      .upd("tbl_filter_condition",  "Condition_details")
    })

    # ── Active filter count badge ─────────────────────────────────────────────
    output$tbl_filter_badge <- renderUI({
      n <- sum(
        length(input$tbl_filter_gene)       > 0,
        length(input$tbl_filter_study)      > 0,
        length(input$tbl_filter_comparison) > 0,
        length(input$tbl_filter_model)      > 0,
        length(input$tbl_filter_cell_line)  > 0,
        length(input$tbl_filter_dpc)        > 0,
        length(input$tbl_filter_diff)       > 0,
        length(input$tbl_filter_condition)  > 0
      )
      if (n > 0) {
        tags$span(class = "badge rounded-pill bg-primary",
                  style = "font-size:0.65rem;", n)
      }
    })

    # ── Clear all filters ────────────────────────────────────────────────────
    observeEvent(input$tbl_clear_filters, {
      updateSelectizeInput(session, "tbl_filter_gene",       selected = character(0))
      updateSelectizeInput(session, "tbl_filter_study",      selected = character(0))
      updateSelectizeInput(session, "tbl_filter_comparison", selected = character(0))
      updateSelectizeInput(session, "tbl_filter_model",      selected = character(0))
      updateSelectizeInput(session, "tbl_filter_cell_line",  selected = character(0))
      updateSelectizeInput(session, "tbl_filter_dpc",        selected = character(0))
      updateSelectizeInput(session, "tbl_filter_diff",       selected = character(0))
      updateSelectizeInput(session, "tbl_filter_condition",  selected = character(0))
    })

    # Reactive: filtered table data (shared by table + summary)
    tbl_data_r <- reactive({
      si   <- study_info_r()
      want <- c("Gene", "Study_title", "Comparison", "Model_system",
                "Cell_Line", "DPC", "Differentation_time_point", "Condition_details")
      cols <- want[want %in% colnames(si)]
      tbl  <- si[!duplicated(si$Assay), cols, drop = FALSE]

      # Apply filters
      .filt <- function(tbl, col, vals) {
        if (length(vals) > 0 && col %in% colnames(tbl))
          tbl[tbl[[col]] %in% vals, , drop = FALSE]
        else tbl
      }
      tbl <- .filt(tbl, "Gene",                      input$tbl_filter_gene)
      tbl <- .filt(tbl, "Study_title",               input$tbl_filter_study)
      tbl <- .filt(tbl, "Comparison",                input$tbl_filter_comparison)
      tbl <- .filt(tbl, "Model_system",              input$tbl_filter_model)
      tbl <- .filt(tbl, "Cell_Line",                 input$tbl_filter_cell_line)
      tbl <- .filt(tbl, "DPC",                       input$tbl_filter_dpc)
      tbl <- .filt(tbl, "Differentation_time_point", input$tbl_filter_diff)
      tbl <- .filt(tbl, "Condition_details",         input$tbl_filter_condition)

      tbl[order(tbl$Gene, tbl$Model_system), ]
    })

    # Summary badges
    output$tbl_summary <- renderUI({
      tbl <- tbl_data_r()
      n_genes  <- length(unique(na.omit(tbl$Gene)))
      n_models <- length(unique(na.omit(tbl$Model_system)))
      n_assays <- nrow(tbl)
      div(
        class = "d-flex gap-2 align-items-center",
        tags$span(
          class = "badge rounded-pill",
          style = "background: #4E79A7; font-size: 0.78rem; font-weight: 500; padding: 6px 12px;",
          paste(n_genes, "genes")
        ),
        tags$span(
          class = "badge rounded-pill",
          style = "background: #59A14F; font-size: 0.78rem; font-weight: 500; padding: 6px 12px;",
          paste(n_models, "model systems")
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
                  Cell_Line = "Cell Line", DPC = "DPC",
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
          pageLength     = 750,
          scrollX        = TRUE,
          paging         = FALSE,
          dom            = "rt",
          rowGroup       = list(dataSrc = gene_idx),
          columnDefs     = list(
            list(visible = FALSE, targets = gene_idx),
            list(className = "dt-left", targets = "_all"),
            list(width = "220px", targets = 1),
            list(width = "130px", targets = 2)
          ),
          initComplete = .dt_header_js
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
          fontWeight = "500", fontSize = "0.85em"
        )
    })

    # =========================================================================
    # 2. Charts Modal — 4 pie charts with assay/gene toggle
    # =========================================================================

    # Helper: render a pie chart for a given grouping column
    .render_pie <- function(col_name, title_label) {
      renderPlotly({
        si <- study_info_r()
        use_assays <- if (!is.null(input$chart_metric)) input$chart_metric == "assays" else TRUE
        if (use_assays) {
          agg <- tapply(si$Assay, si[[col_name]], function(x) length(unique(x)))
        } else {
          agg <- tapply(si$Gene, si[[col_name]], function(x) length(unique(na.omit(x))))
        }
        df <- data.frame(label = names(agg), value = as.integer(agg), stringsAsFactors = FALSE)
        df <- df[order(-df$value), ]
        clrs <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#76B7B2",
                  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")
        plot_ly(df, labels = ~label, values = ~value, type = "pie",
                textinfo = "label+value",
                hovertemplate = "<b>%{label}</b><br>Count: %{value}<br>%{percent}<extra></extra>",
                marker = list(colors = clrs[seq_len(nrow(df))])) |>
          layout(title = list(text = title_label, font = list(size = 14)),
                 showlegend = FALSE,
                 margin = list(l = 10, r = 10, t = 40, b = 10)) |>
          config(displayModeBar = FALSE)
      })
    }

    output$pie_model_system <- .render_pie("Model_system", "Model System")
    output$pie_ko_strat     <- .render_pie("KO_strat",     "KO Strategy")
    output$pie_cell_line    <- .render_pie("Cell_Line",     "Cell Line")
    output$pie_dpc          <- .render_pie("DPC",           "DPC")

    # ── N Assays per Gene bar chart ───────────────────────────────────────
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

    # ── Charts Modal ──────────────────────────────────────────────────────
    observeEvent(input$show_charts, {
      ns <- session$ns
      showModal(modalDialog(
        title     = NULL,
        size      = "xl",
        easyClose = TRUE,
        div(
          class = "d-flex justify-content-center pb-3",
          shinyWidgets::radioGroupButtons(
            inputId  = ns("chart_metric"),
            label    = NULL,
            choices  = c("N Assays" = "assays", "N Genes" = "genes"),
            selected = if (!is.null(isolate(input$chart_metric))) isolate(input$chart_metric) else "assays",
            size     = "sm",
            status   = "outline-primary"
          )
        ),
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          card(card_body(plotlyOutput(ns("pie_model_system"), height = "300px")), full_screen = TRUE),
          card(card_body(plotlyOutput(ns("pie_ko_strat"),     height = "300px")), full_screen = TRUE),
          card(card_body(plotlyOutput(ns("pie_cell_line"),     height = "300px")), full_screen = TRUE),
          card(card_body(plotlyOutput(ns("pie_dpc"),           height = "300px")), full_screen = TRUE)
        ),
        tags$hr(),
        tags$h6(class = "fw-semibold mt-2", "N Assays per Gene"),
        plotlyOutput(ns("assays_per_gene"), height = "480px"),
        footer = modalButton("Close")
      ))
    })

    # =========================================================================
    # 4. Functional Annotations (Reactome, GO, PANTHER)
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

    # Bar chart: top 10 rows by N KO genes annotated
    .ko_render_bar <- function(df_r, src_label, bar_colour = "#B07AA1") {
      renderPlotly({
        df <- df_r(); req(df, nrow(df) > 0)
        df <- df[order(-df$overlap_n), ]
        top <- head(df, 10)
        top$category <- ifelse(nchar(top$category) > 55,
                               paste0(substr(top$category, 1, 52), "\u2026"),
                               top$category)
        top$category <- factor(top$category, levels = rev(top$category))
        plot_ly(top,
                x             = ~overlap_n,
                y             = ~category,
                type          = "bar",
                orientation   = "h",
                marker        = list(color = bar_colour),
                customdata    = ~n_pathway,
                hovertemplate = paste0(
                  "<b>%{y}</b><br>",
                  "KO Genes: %{x}<br>",
                  "Category size: %{customdata} genes<extra></extra>")) |>
          layout(
            xaxis  = list(title = "N KO Genes"),
            yaxis  = list(title = "", automargin = TRUE),
            margin = list(l = 5, r = 5, t = 35, b = 5),
            title  = list(text = paste0("Top 10 \u2014 ", src_label),
                          font = list(size = 13))) |>
          config(displayModeBar = FALSE)
      })
    }

    # Table: all rows — Category, Category Size, N Genes, Genes
    .ko_render_tbl <- function(df_r) {
      DT::renderDT({
        df <- df_r(); req(df, nrow(df) > 0)
        df <- df[order(-df$overlap_n), ]
        df$n_pathway  <- as.integer(df$n_pathway)
        df$overlap_n  <- as.integer(df$overlap_n)
        tbl <- df[, c("category", "n_pathway", "overlap_n", "genes")]
        names(tbl) <- c("Category", "Category Size", "N Genes", "Genes")
        DT::datatable(tbl, rownames = FALSE,
                      options = list(pageLength = 25, scrollX = TRUE,
                                     scrollY = "420px", scrollCollapse = TRUE,
                                     dom = "tip",
                                     initComplete = .dt_header_js),
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

    # =========================================================================
    # 5. Disease & Phenotype Associations (IMPC, HPO, OMIM, Orphanet)
    # =========================================================================

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

    # Table: all annotation terms — Category, N Genes, Genes
    .ko_render_pheno_tbl <- function(df_r, ann_col_label) {
      DT::renderDT({
        df <- df_r(); req(df, nrow(df) > 0)
        tbl <- df |>
          dplyr::group_by(annotation) |>
          dplyr::summarise(
            n_genes   = dplyr::n_distinct(gene),
            gene_list = paste(sort(unique(gene)), collapse = "; "),
            .groups = "drop") |>
          dplyr::arrange(dplyr::desc(n_genes)) |>
          as.data.frame()
        names(tbl) <- c(ann_col_label, "N Genes", "Genes")
        DT::datatable(tbl, rownames = FALSE,
                      options = list(pageLength = 25, scrollX = TRUE,
                                     scrollY = "420px", scrollCollapse = TRUE,
                                     dom = "tip",
                                     initComplete = .dt_header_js),
                      class = "compact row-border hover")
      })
    }

    # Wire up phenotype sources
    output$ko_bar_impc     <- .ko_render_pheno_bar(ko_data_impc,     "IMPC Phenotype",   "#E63946")
    output$ko_tbl_impc     <- .ko_render_pheno_tbl(ko_data_impc,     "IMPC Phenotype")
    output$ko_bar_hpo      <- .ko_render_pheno_bar(ko_data_hpo,      "HPO Term",         "#D62728")
    output$ko_tbl_hpo      <- .ko_render_pheno_tbl(ko_data_hpo,      "HPO Term")

    # Wire up disease sources (reactive based on selector)
    ko_data_disease <- reactive({
      src <- input$disease_source
      if (is.null(src) || src == "omim") ko_data_omim() else ko_data_orphanet()
    })

    output$ko_bar_disease <- renderPlotly({
      df <- ko_data_disease(); req(df, nrow(df) > 0)
      src <- if (!is.null(input$disease_source)) input$disease_source else "omim"
      src_label  <- if (src == "omim") "OMIM Disease" else "Orphanet Disorder"
      bar_colour <- if (src == "omim") "#457B9D" else "#1D6FA4"
      agg <- df |>
        dplyr::group_by(annotation) |>
        dplyr::summarise(n_genes = dplyr::n_distinct(gene), .groups = "drop") |>
        dplyr::arrange(dplyr::desc(n_genes)) |>
        head(10)
      agg$annotation <- ifelse(nchar(agg$annotation) > 55,
                               paste0(substr(agg$annotation, 1, 52), "\u2026"),
                               agg$annotation)
      agg <- agg[order(agg$n_genes), ]
      agg$annotation <- factor(agg$annotation, levels = agg$annotation)
      plot_ly(agg, x = ~n_genes, y = ~annotation, type = "bar", orientation = "h",
              marker = list(color = bar_colour),
              hovertemplate = paste0("<b>%{y}</b><br>", "KO Genes: %{x}<extra></extra>")) |>
        layout(
          xaxis  = list(title = "N KO Genes", dtick = 1),
          yaxis  = list(title = "", automargin = TRUE),
          margin = list(l = 5, r = 5, t = 35, b = 5),
          title  = list(text = paste0("Top 10 \u2014 ", src_label),
                        font = list(size = 13))) |>
        config(displayModeBar = FALSE)
    })

    output$ko_tbl_disease <- DT::renderDT({
      df <- ko_data_disease(); req(df, nrow(df) > 0)
      src <- if (!is.null(input$disease_source)) input$disease_source else "omim"
      ann_col_label <- if (src == "omim") "OMIM Phenotype" else "Disorder"
      tbl <- df |>
        dplyr::group_by(annotation) |>
        dplyr::summarise(
          n_genes   = dplyr::n_distinct(gene),
          gene_list = paste(sort(unique(gene)), collapse = "; "),
          .groups = "drop") |>
        dplyr::arrange(dplyr::desc(n_genes)) |>
        as.data.frame()
      names(tbl) <- c(ann_col_label, "N Genes", "Genes")
      DT::datatable(tbl, rownames = FALSE,
                    options = list(pageLength = 25, scrollX = TRUE,
                                   scrollY = "420px", scrollCollapse = TRUE,
                                   dom = "tip",
                                   initComplete = .dt_header_js),
                    class = "compact row-border hover")
    })

    # ── Suspend all outputs until visible ─────────────────────────────────────
    invisible(lapply(
      c("assay_tbl", "tbl_summary", "tbl_filter_badge",
        "pie_model_system", "pie_ko_strat", "pie_cell_line", "pie_dpc", "assays_per_gene",
        "ko_bar_reactome",     "ko_tbl_reactome",
        "ko_bar_go_bp",        "ko_tbl_go_bp",
        "ko_bar_go_mf",        "ko_tbl_go_mf",
        "ko_bar_panther_class","ko_tbl_panther_class",
        "ko_bar_impc",         "ko_tbl_impc",
        "ko_bar_hpo",          "ko_tbl_hpo",
        "ko_bar_disease",      "ko_tbl_disease"),
      function(oid) outputOptions(output, oid, suspendWhenHidden = TRUE)
    ))

  })
}
