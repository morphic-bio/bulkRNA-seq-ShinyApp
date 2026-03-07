# =============================================================================
# Module: deg_annotations
# Page 3 — Assay Analysis
#
# One card, shared assay picker, tabs grouped by:
#   Functional Annotations:
#     1. Pathway Overlap       — Reactome (bar + table + gene chip)
#     2. GO Overlap            — GO BP + GO MF (bar + table + gene chip each)
#     3. Protein Class         — PANTHER Class (bar + table + gene chip)
#   Disease & Phenotype:
#     4. Phenotype Annotations — IMPC + HPO (sub-tabs: Top / Coverage / Table)
#     5. Disease Annotations   — OMIM + Orphanet (sub-tabs: Top / Coverage / Table)
#
# Requires helpers.R sourced first.
# =============================================================================

library(plotly)
library(DT)
library(dplyr)

# ── Phenotype annotation constants ────────────────────────────────────────────

.PB_TABLE <- c(
  impc     = "impc_phenotypes",
  hpo      = "hpo_phenotypes",
  omim     = "omim_phenotypes",
  orphanet = "orphanet_data"
)
.PB_COL <- c(
  impc     = "impc_phenotypes",
  hpo      = "hpo_name",
  omim     = "omim_phenotype",
  orphanet = "disorder_name"
)
.PB_LABEL <- c(
  impc     = "Phenotype",
  hpo      = "HPO Term",
  omim     = "OMIM Phenotype",
  orphanet = "Disorder"
)
.PB_CLRS     <- c(up = "#E63946", down = "#457B9D")
.PB_ANN_CLRS <- c("Annotated" = "#6AB187", "No annotation" = "#CBD5D9")

# ── UI helper: left-column control panel for DEG overlap tabs ─────────────────

.deg_ctrl <- function(ns, pfx, topn_default = 20L, min_default = 10L) {
  div(
    class = "d-flex align-items-end flex-wrap gap-3 mb-2",
    div(style = "min-width:140px;",
        selectInput(ns(paste0(pfx, "_dir")), "Direction",
                    choices  = c("Up & Down" = "updn", "Up only" = "up", "Down only" = "down"),
                    selected = "updn", width = "100%")),
    div(style = "min-width:180px; flex:1; max-width:280px;",
        sliderInput(ns(paste0(pfx, "_topn")), "Top N pathways",
                    5, 50, topn_default, step = 5, width = "100%")),
    div(style = "min-width:180px; flex:1; max-width:280px;",
        sliderInput(ns(paste0(pfx, "_minsize")), "Min category size",
                    1, 500, min_default, width = "100%"))
  )
}

# ── UI ────────────────────────────────────────────────────────────────────────

deg_annotationsUI <- function(id) {
  ns <- NS(id)

  card(
    card_header(tagList(
      bsicons::bs_icon("diagram-3"), " Assay Analysis",
      .info_tip("Explore DEG annotations for a single assay: pathway overlaps, GO terms, protein classes, phenotype & disease associations.")
    )),
    full_screen = TRUE,

    # ── Assay picker (left) + filters (right) ────────────────────────────
    div(
      class = "border-bottom pb-2 mb-0 px-2 pt-2",
      layout_columns(
        col_widths = c(4, 8),

        # Left: assay picker
        div(
          div(class = "d-flex align-items-center justify-content-between mb-1",
              div(class = "d-flex align-items-center gap-2",
                  tags$span(class = "fw-semibold small", "Select assay"),
                  actionButton(ns("deg_clear_assay"),
                               tagList(bsicons::bs_icon("x-lg", size = "0.6rem"), " Clear selection"),
                               class = "btn btn-link btn-sm p-0",
                               style = "font-size: 0.7rem; text-decoration: none; color: #212529;")),
              uiOutput(ns("deg_match_count"), inline = TRUE)),
          selectizeInput(ns("deg_assay"), NULL,
                         choices = NULL, selected = NULL, multiple = FALSE,
                         options = list(placeholder = "Choose an assay\u2026"),
                         width = "100%"),
          uiOutput(ns("deg_no_assays"))
        ),

        # Right: collapsible filters
        div(
          div(class = "d-flex align-items-center gap-2 mb-1",
              tags$button(
                type = "button",
                class = "btn btn-outline-secondary btn-sm d-flex align-items-center gap-1",
                `data-bs-toggle` = "collapse",
                `data-bs-target` = paste0("#", ns("deg_filter_panel")),
                `aria-expanded` = "false",
                `aria-controls` = ns("deg_filter_panel"),
                bsicons::bs_icon("funnel", size = "0.75rem"),
                "Filters",
                uiOutput(ns("deg_filter_badge"), inline = TRUE)
              ),
              actionButton(ns("deg_clear_filters"),
                           tagList(bsicons::bs_icon("x-lg", size = "0.6rem"), " Clear"),
                           class = "btn btn-link btn-sm p-0",
                           style = "font-size: 0.7rem; text-decoration: none; color: #212529;")
          ),
          div(id = ns("deg_filter_panel"), class = "collapse",
              div(class = "pt-2",
                  layout_columns(
                    col_widths = c(4, 4, 4),
                    selectInput(ns("deg_filter_gene"),  "KO Gene",
                                choices = NULL, multiple = TRUE, width = "100%"),
                    selectInput(ns("deg_filter_model"), "Model System",
                                choices = NULL, multiple = TRUE, width = "100%"),
                    selectInput(ns("deg_filter_comp"),  "Comparison",
                                choices = NULL, multiple = TRUE, width = "100%")
                  )
              )
          )
        )
      )
    ),

    # ── Empty state when no assay selected ────────────────────────────────
    uiOutput(ns("deg_no_selection")),

    # ── Assay metadata card (shown when assay is selected) ─────────────────
    uiOutput(ns("deg_meta_above")),

    # ── Tabs grouped into Functional / Disease & Phenotype ──────────────────
    navset_tab(

      # ── Functional Annotations ─────────────────────────────────────────────
      nav_menu(
        title = "Functional Annotations",
        icon  = bsicons::bs_icon("intersect"),

      # ── Tab 1: Pathway Overlap (Reactome) ──────────────────────────────────
      nav_panel(
        tagList("Pathway Overlap",
                .info_tip("Overlap between assay DEGs and Reactome pathway gene sets.")),
        div(
          class = "py-3 px-1",
          uiOutput(ns("gs_ro")),
          .deg_ctrl(ns, "ro"),
          card(
            navset_tab(
              nav_panel("Bar Chart",
                        plotlyOutput(ns("ro_bar"), height = "420px")),
              nav_panel("Table",
                        div(class = "p-2",
                            tags$p(class = "text-muted small mb-1",
                                   bsicons::bs_icon("table"),
                                   " Click a row to see genes for that pathway."),
                            DTOutput(ns("ro_tbl"))))
            ),
            uiOutput(ns("ro_gene_panel"))
          )
        )
      ),

      # ── Tab 2: Biological Processes (GO BP) ────────────────────────────────
      nav_panel(
        tagList("Biological Processes",
                .info_tip("Overlap between assay DEGs and Gene Ontology Biological Process terms.")),
        div(
          class = "py-3 px-1",
          uiOutput(ns("gs_gobp")),
          .deg_ctrl(ns, "gobp"),
          card(
            navset_tab(
              nav_panel("Bar Chart",
                        plotlyOutput(ns("gobp_bar"), height = "420px")),
              nav_panel("Table",
                        div(class = "p-2",
                            tags$p(class = "text-muted small mb-1",
                                   bsicons::bs_icon("table"),
                                   " Click a row to see genes for that GO term."),
                            DTOutput(ns("gobp_tbl"))))
            ),
            uiOutput(ns("gobp_gene_panel"))
          )
        )
      ),

      # ── Tab 3: Molecular Functions (GO MF) ─────────────────────────────────
      nav_panel(
        tagList("Molecular Functions",
                .info_tip("Overlap between assay DEGs and Gene Ontology Molecular Function terms.")),
        div(
          class = "py-3 px-1",
          uiOutput(ns("gs_gomf")),
          .deg_ctrl(ns, "gomf"),
          card(
            navset_tab(
              nav_panel("Bar Chart",
                        plotlyOutput(ns("gomf_bar"), height = "420px")),
              nav_panel("Table",
                        div(class = "p-2",
                            tags$p(class = "text-muted small mb-1",
                                   bsicons::bs_icon("table"),
                                   " Click a row to see genes for that GO term."),
                            DTOutput(ns("gomf_tbl"))))
            ),
            uiOutput(ns("gomf_gene_panel"))
          )
        )
      ),

      # ── Tab 4: Protein Class Overlap (PANTHER Class) ───────────────────────
      nav_panel(
        tagList("Protein Class",
                .info_tip("Overlap between assay DEGs and PANTHER protein class annotations.")),
        div(
          class = "py-3 px-1",
          uiOutput(ns("gs_pc")),
          .deg_ctrl(ns, "pc"),
          card(
            navset_tab(
              nav_panel("Bar Chart",
                        plotlyOutput(ns("pc_bar"), height = "420px")),
              nav_panel("Table",
                        div(class = "p-2",
                            tags$p(class = "text-muted small mb-1",
                                   bsicons::bs_icon("table"),
                                   " Click a row to see genes for that class."),
                            DTOutput(ns("pc_tbl"))))
            ),
            uiOutput(ns("pc_gene_panel"))
          )
        )
      )

      ), # end nav_menu "Functional Annotations"

      # ── Disease & Phenotype Associations ───────────────────────────────────
      nav_menu(
        title = "Disease & Phenotype",
        icon  = bsicons::bs_icon("card-list"),

      # ── Tab 5: Mouse Phenotypes (IMPC) ─────────────────────────────────────
      nav_panel(
        tagList("Mouse Phenotypes",
                .info_tip("IMPC mouse phenotype annotations for the assay's DEGs.")),
        div(
          class = "py-3 px-1",
          uiOutput(ns("gs_ph")),
          div(
            class = "d-flex align-items-end flex-wrap gap-3 mb-2",
            div(style = "min-width:180px;",
                checkboxGroupInput(ns("ph_zyg"), "IMPC Zygosity",
                                   choices  = c("Homozygote" = "homozygote",
                                                "Heterozygote" = "heterozygote"),
                                   selected = c("homozygote", "heterozygote"),
                                   inline = TRUE)),
            div(style = "min-width:180px; flex:1; max-width:280px;",
                sliderInput(ns("ph_topn"), "Top N phenotypes", 5, 50, 20,
                            step = 5, width = "100%")),
            div(style = "min-width:200px;",
                radioButtons(ns("ph_dir"), "Direction",
                             choices  = c("Both" = "both",
                                          "\u2191 Up only" = "up",
                                          "\u2193 Down only" = "down"),
                             selected = "both", inline = TRUE))
          ),
          navset_tab(
            nav_panel(
              tagList("Top Phenotypes",
                      .info_tip("Most frequently annotated phenotypes among the assay's DEGs.")),
              card(full_screen = TRUE, class = "mt-2",
                   card_body(plotlyOutput(ns("ph_top"), height = "480px")))
            ),
            nav_panel(
              tagList("Annotation Coverage",
                      .info_tip("Proportion of up/down-regulated DEGs that have annotations in IMPC.")),
              card(class = "mt-2",
                   card_body(plotlyOutput(ns("ph_coverage"), height = "280px")))
            ),
            nav_panel(
              "Table",
              card(class = "mt-2",
                   card_body(
                     tags$p(class = "text-muted small mb-1",
                            bsicons::bs_icon("table"),
                            " Click a row to see up/down-regulated genes for that phenotype."),
                     DTOutput(ns("ph_tbl")))),
              uiOutput(ns("ph_gene_panel"))
            )
          )
        )
      ),

      # ── Tab 6: Human Phenotypes (HPO) ──────────────────────────────────────
      nav_panel(
        tagList("Human Phenotypes",
                .info_tip("HPO human phenotype annotations for the assay's DEGs.")),
        div(
          class = "py-3 px-1",
          uiOutput(ns("gs_hpo")),
          div(
            class = "d-flex align-items-end flex-wrap gap-3 mb-2",
            div(style = "min-width:180px; flex:1; max-width:280px;",
                sliderInput(ns("hpo_topn"), "Top N phenotypes", 5, 50, 20,
                            step = 5, width = "100%")),
            div(style = "min-width:200px;",
                radioButtons(ns("hpo_dir"), "Direction",
                             choices  = c("Both" = "both",
                                          "\u2191 Up only" = "up",
                                          "\u2193 Down only" = "down"),
                             selected = "both", inline = TRUE))
          ),
          navset_tab(
            nav_panel(
              tagList("Top Phenotypes",
                      .info_tip("Most frequently annotated HPO terms among the assay's DEGs.")),
              card(full_screen = TRUE, class = "mt-2",
                   card_body(plotlyOutput(ns("hpo_top"), height = "480px")))
            ),
            nav_panel(
              tagList("Annotation Coverage",
                      .info_tip("Proportion of up/down-regulated DEGs that have annotations in HPO.")),
              card(class = "mt-2",
                   card_body(plotlyOutput(ns("hpo_coverage"), height = "280px")))
            ),
            nav_panel(
              "Table",
              card(class = "mt-2",
                   card_body(
                     tags$p(class = "text-muted small mb-1",
                            bsicons::bs_icon("table"),
                            " Click a row to see up/down-regulated genes for that phenotype."),
                     DTOutput(ns("hpo_tbl")))),
              uiOutput(ns("hpo_gene_panel"))
            )
          )
        )
      ),

      # ── Tab 7: Disease Associations (OMIM / Orphanet) ──────────────────────
      nav_panel(
        tagList("Disease Associations",
                .info_tip("Disease annotations for DEGs from OMIM and Orphanet databases.")),
        div(
          class = "py-3 px-1",
          uiOutput(ns("gs_di")),
          div(
            class = "d-flex align-items-end flex-wrap gap-3 mb-2",
            div(style = "min-width:180px;",
                radioButtons(ns("di_ds"), "Dataset",
                             choices  = c("OMIM" = "omim", "Orphanet" = "orphanet"),
                             selected = "omim", inline = TRUE)),
            div(style = "min-width:180px; flex:1; max-width:280px;",
                sliderInput(ns("di_topn"), "Top N shown", 5, 50, 20,
                            step = 5, width = "100%")),
            div(style = "min-width:200px;",
                radioButtons(ns("di_dir"), "Direction",
                             choices  = c("Both" = "both",
                                          "\u2191 Up only" = "up",
                                          "\u2193 Down only" = "down"),
                             selected = "both", inline = TRUE))
          ),
          navset_tab(
            nav_panel(
              tagList("Top Annotations",
                      .info_tip("Most frequently annotated disease terms among the assay's DEGs.")),
              card(full_screen = TRUE, class = "mt-2",
                   card_body(plotlyOutput(ns("di_top"), height = "480px")))
            ),
            nav_panel(
              tagList("Annotation Coverage",
                      .info_tip("Proportion of up/down-regulated DEGs that have annotations in the selected database.")),
              card(class = "mt-2",
                   card_body(plotlyOutput(ns("di_coverage"), height = "280px")))
            ),
            nav_panel(
              "Table",
              card(class = "mt-2",
                   card_body(
                     tags$p(class = "text-muted small mb-1",
                            bsicons::bs_icon("table"),
                            " Click a row to see up/down-regulated genes for that disease."),
                     DTOutput(ns("di_tbl")))),
              uiOutput(ns("di_gene_panel"))
            )
          )
        )
      )
      ) # end nav_menu "Disease & Phenotype"
    ) # end navset_tab
  )   # end card
}

# ── Server ────────────────────────────────────────────────────────────────────

deg_annotationsServer <- function(id, con_r, study_info_r) {
  moduleServer(id, function(input, output, session) {

    # =========================================================================
    # Shared: assay filter + picker
    # =========================================================================

    observe({
      si <- study_info_r()
      updateSelectInput(session, "deg_filter_gene",
                        choices  = sort(unique(na.omit(si$Gene))),
                        selected = character(0))
      updateSelectInput(session, "deg_filter_model",
                        choices  = sort(unique(na.omit(si$Model_system))),
                        selected = character(0))
      updateSelectInput(session, "deg_filter_comp",
                        choices  = sort(unique(na.omit(si$Comparison))),
                        selected = character(0))
    })

    filtered_deg_assays <- reactive({
      si      <- study_info_r()
      si_uniq <- si[!duplicated(si$Assay), ]
      if (length(input$deg_filter_gene)  > 0)
        si_uniq <- si_uniq[si_uniq$Gene         %in% input$deg_filter_gene,  ]
      if (length(input$deg_filter_model) > 0)
        si_uniq <- si_uniq[si_uniq$Model_system %in% input$deg_filter_model, ]
      if (length(input$deg_filter_comp)  > 0)
        si_uniq <- si_uniq[si_uniq$Comparison   %in% input$deg_filter_comp, ]
      sort(si_uniq$Assay)
    })

    observeEvent(filtered_deg_assays(), {
      ch   <- filtered_deg_assays()
      cur  <- isolate(input$deg_assay)
      kept <- if (!is.null(cur) && nzchar(cur) && cur %in% ch) cur else ""
      updateSelectizeInput(session, "deg_assay",
                           choices = ch, selected = kept, server = TRUE)
    })

    output$deg_no_assays <- renderUI({
      if (length(filtered_deg_assays()) == 0)
        tags$p(class = "text-muted small fst-italic mb-0 mt-1",
               "No assays match the current filters.")
    })

    # ── Match count badge ───────────────────────────────────────────────────
    output$deg_match_count <- renderUI({
      n <- length(filtered_deg_assays())
      tags$span(
        class = "badge rounded-pill",
        style = "font-size: 0.7rem; background: #6c757d; color: white;",
        paste0(n, " assay", if (n != 1) "s")
      )
    })

    # ── Active filter count badge ─────────────────────────────────────────────
    output$deg_filter_badge <- renderUI({
      n <- sum(
        length(input$deg_filter_gene)  > 0,
        length(input$deg_filter_model) > 0,
        length(input$deg_filter_comp)  > 0
      )
      if (n > 0) {
        tags$span(class = "badge rounded-pill bg-primary",
                  style = "font-size:0.65rem;", n)
      }
    })

    # ── Clear all filters ────────────────────────────────────────────────────
    observeEvent(input$deg_clear_filters, {
      updateSelectInput(session, "deg_filter_gene",  selected = character(0))
      updateSelectInput(session, "deg_filter_model", selected = character(0))
      updateSelectInput(session, "deg_filter_comp",  selected = character(0))
    })

    # ── Clear assay selection ──────────────────────────────────────────────────
    observeEvent(input$deg_clear_assay, {
      ch <- filtered_deg_assays()
      updateSelectizeInput(session, "deg_assay",
                           choices = ch, selected = "", server = TRUE)
    })

    # ── Empty state when no assay selected ────────────────────────────────────
    output$deg_no_selection <- renderUI({
      assay <- input$deg_assay
      if (is.null(assay) || !nzchar(assay)) {
        .empty_state("diagram-3",
                     "Select an assay to view plots",
                     "Use the assay picker above to choose an assay.")
      }
    })

    # ── Assay metadata card ──────────────────────────────────────────────────
    output$deg_meta_above <- renderUI({
      assay <- input$deg_assay
      req(assay, nzchar(assay))
      si <- study_info_r()
      r  <- si[si$Assay == assay, ][1, ]

      bdg <- function(val, bg_col) {
        v <- as.character(val)
        if (is.null(val) || is.na(val) || !nzchar(v)) return(NULL)
        tags$span(class = "badge me-1", style = paste0("background:", bg_col,
                  "; color:white; font-size:0.72rem;"), v)
      }

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

      tagList(
        div(class = "d-flex flex-wrap gap-2 mb-2",
            div(class = "border rounded p-2",
                style = "background:#f8f9fa; flex:1; min-width:180px;",
                tags$div(class = "small fw-semibold text-break mb-1",
                         style = "font-size:0.68rem; word-break:break-all; color:#495057;", assay),
                tags$div(bdg(r$Gene, "#198754"), bdg(r$Model_system, "#0d6efd"),
                         bdg(r$KO_strat, "#0dcaf0"), bdg(r$DPC, "#e67e22"),
                         bdg(r$Cell_Line, "#6c757d"), extra))),
        tags$hr(class = "my-2")
      )
    })

    assay_r <- reactive({ input$deg_assay })

    # KO gene for the selected assay (from study_info)
    assay_gene_r <- reactive({
      assay <- assay_r(); req(assay, nzchar(assay))
      si    <- study_info_r()
      gene  <- na.omit(unique(si$Gene[si$Assay == assay]))
      if (length(gene) == 0) return(NULL)
      gene[1]
    })

    # Assay metadata (model system, KO strategy, conditions, diff time point)
    assay_meta_r <- reactive({
      assay <- assay_r(); req(assay, nzchar(assay))
      si    <- study_info_r()
      row   <- si[si$Assay == assay, ][1, ]
      .safe <- function(col) {
        if (col %in% names(row) && !is.na(row[[col]])) as.character(row[[col]]) else NULL
      }
      list(
        model_system = .safe("Model_system"),
        ko_strat     = .safe("KO_strat"),
        conditions   = .safe("Condition_levels"),
        diff_tp      = .safe("Differentation_time_point")
      )
    })

    # =========================================================================
    # Shared: DEG symbols for phenotype tabs
    # =========================================================================

    deg_symbols_r <- reactive({
      assay <- assay_r(); req(assay, nzchar(assay))
      con <- con_r()
      tryCatch(
        dbGetQuery(con, sprintf(
          "SELECT COALESCE(NULLIF(TRIM(hgnc_symbol), ''), gene_ID) AS symbol, DEG
             FROM main.\"%s\" WHERE DEG IN ('up', 'down')", assay)),
        error = function(e) NULL)
    })

    # =========================================================================
    # DEG overlap factory helpers (Tabs 1–3 + PanelApp in Tab 5)
    # =========================================================================

    # Build a reactive that computes DEG overlap on-the-fly via ref tables.
    .make_overlap_data <- function(src, dir_pfx, topn_pfx, minsize_pfx) {
      ref_tbl <- .HP_REF_TBL[[src]]
      reactive({
        assay   <- assay_r();            req(assay, nzchar(assay))
        dir     <- input[[dir_pfx]];     req(dir)
        topn    <- input[[topn_pfx]];    req(topn)
        minsize <- input[[minsize_pfx]]; req(minsize)
        con     <- con_r()
        compute_deg_overlap(con, assay, ref_tbl, dir, topn, minsize)
      })
    }

    data_ro   <- .make_overlap_data("reactome",      "ro_dir", "ro_topn", "ro_minsize")
    data_gobp <- .make_overlap_data("go_bp",         "gobp_dir", "gobp_topn", "gobp_minsize")
    data_gomf <- .make_overlap_data("go_mf",         "gomf_dir", "gomf_topn", "gomf_minsize")
    data_pc   <- .make_overlap_data("panther_class", "pc_dir", "pc_topn", "pc_minsize")

    # Bar chart renderer (shared across all DEG overlap tabs)
    .render_overlap_bar <- function(data_r, dir_pfx, src_label) {
      renderPlotly({
        df  <- data_r(); req(df, nrow(df) > 0)
        dir <- input[[dir_pfx]]; if (is.null(dir)) dir <- "updn"
        dir_clr <- switch(dir,
          updn = .HP_CLRS$all, up = .HP_CLRS$up, down = .HP_CLRS$down, .HP_CLRS$all)
        df$pct      <- round(df$pct, 1)
        df$category <- ifelse(nchar(df$category) > 55,
                              paste0(substr(df$category, 1, 52), "\u2026"),
                              df$category)
        df$category <- factor(df$category, levels = rev(df$category))
        plot_ly(df, x = ~pct, y = ~category, type = "bar", orientation = "h",
                marker = list(color = dir_clr),
                hovertemplate = "<b>%{y}</b><br>Overlap: %{x}%<extra></extra>") |>
          layout(xaxis  = list(title = "% DEGs in Category"),
                 yaxis  = list(title = "", automargin = TRUE),
                 margin = list(l = 5, r = 5, t = 35, b = 5),
                 title  = list(text = paste0("DEG Overlap \u2014 ", src_label),
                               font = list(size = 13))) |>
          config(displayModeBar = FALSE)
      })
    }

    # Table renderer (shared; allow_select enables row-click gene chips)
    .render_overlap_tbl <- function(data_r, allow_select = FALSE) {
      DT::renderDT({
        df <- data_r(); req(df, nrow(df) > 0)
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
        DT::datatable(
          tbl, rownames = FALSE, filter = "top",
          selection = if (allow_select) "single" else "none",
          options = list(pageLength = 15, scrollX = TRUE,
                         scrollY = "300px", scrollCollapse = TRUE,
                         dom = "frtip",
                         initComplete = .dt_header_js),
          class = "compact row-border hover"
        ) |>
          DT::formatStyle(
            "% Overlap",
            background = DT::styleColorBar(c(0, max(df$pct, na.rm = TRUE)), "#b3cde3"),
            backgroundSize = "100% 60%", backgroundRepeat = "no-repeat",
            backgroundPosition = "center")
      })
    }

    # =========================================================================
    # Gene chip helpers (used by all tabs with row-click gene panels)
    # =========================================================================

    .parse_genes <- function(x) {
      if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(character(0))
      v <- trimws(strsplit(x, ";")[[1]]); v[nzchar(v)]
    }

    .fetch_logfc <- function(con, assay, all_genes) {
      col_names <- tryCatch(dbListFields(con, assay), error = function(e) character(0))
      lfc_col   <- col_names[grepl("log2FoldChange$", col_names)]
      if (length(lfc_col) == 0) return(NULL)
      gene_sql <- paste0("'", paste(gsub("'", "''", all_genes), collapse = "','"), "'")
      lfc_df   <- tryCatch(
        dbGetQuery(con, sprintf(
          'SELECT hgnc_symbol, "%s" AS log2fc FROM main."%s"
           WHERE hgnc_symbol IN (%s)', lfc_col[1], assay, gene_sql)),
        error = function(e) NULL)
      if (is.null(lfc_df) || nrow(lfc_df) == 0) return(NULL)
      setNames(lfc_df$log2fc, lfc_df$hgnc_symbol)
    }

    .make_gene_info_r <- function(sel_row_rv) {
      reactive({
        row   <- sel_row_rv(); req(row)
        assay <- assay_r();    req(assay, nzchar(assay))
        up_genes  <- .parse_genes(row$genes_up[1])
        dn_genes  <- .parse_genes(row$genes_down[1])
        all_genes <- unique(c(up_genes, dn_genes))
        if (length(all_genes) == 0) return(NULL)
        lfc_map <- .fetch_logfc(con_r(), assay, all_genes)
        list(up = up_genes, down = dn_genes, lfc_map = lfc_map)
      })
    }

    .render_gene_chip_panel <- function(sel_row_rv, gene_info_r) {
      renderUI({
        row <- sel_row_rv()
        if (is.null(row)) return(NULL)
        info <- gene_info_r(); req(info)
        up_genes <- info$up; dn_genes <- info$down; lfc_map <- info$lfc_map

        make_chips <- function(genes, is_up) {
          if (length(genes) == 0) return(NULL)
          lfc_vals <- if (!is.null(lfc_map)) lfc_map[genes] else rep(NA_real_, length(genes))
          lfc_vals[is.na(lfc_vals)] <- 0
          ord   <- if (is_up) order(-lfc_vals) else order(lfc_vals)
          genes <- genes[ord]; lfc_vals <- lfc_vals[ord]
          lapply(seq_along(genes), function(i) {
            g <- genes[i]; lfc <- lfc_vals[i]; lfc_abs <- abs(lfc)
            bg    <- .hp_chip_colour(is_up, lfc_abs)
            arrow <- if (is_up) "\u2191" else "\u2193"
            lfc_txt <- if (!is.na(lfc) && lfc != 0)
              paste0(" (", ifelse(lfc >= 0, "+", ""), round(lfc, 2), ")") else ""
            tags$span(style = .hp_chip_style(bg), paste0(g, " ", arrow, lfc_txt))
          })
        }

        up_chips <- make_chips(up_genes, TRUE)
        dn_chips <- make_chips(dn_genes, FALSE)
        n_up     <- length(up_genes); n_dn <- length(dn_genes)
        n_total  <- as.integer(row$n_pathway[1])

        section <- function(label, n, chips, border_col) {
          if (n == 0) return(NULL)
          div(class = "mb-3 p-2 rounded",
              style = paste0("border-left: 4px solid ", border_col, "; background:#fafafa;"),
              tags$div(class = "fw-bold mb-1", style = paste0("color:", border_col, ";"),
                       label, tags$span(class = "fw-normal text-muted ms-1",
                                        paste0("(", n, " gene", if (n == 1) "" else "s",
                                               " \u2014 sorted by |logFC|)"))),
              div(class = "d-flex flex-wrap", chips))
        }

        legend <- div(
          class = "text-muted small mb-2", "Colour intensity = |logFC|:",
          tags$span(style = .hp_chip_style("#F1948A"), "light < 1"),
          tags$span(style = .hp_chip_style("#E74C3C"), "medium 1\u20132"),
          tags$span(style = .hp_chip_style("#C0392B"), "dark \u22652"), " | ",
          tags$span(style = .hp_chip_style("#85C1E9"), "light < 1"),
          tags$span(style = .hp_chip_style("#5DADE2"), "medium 1\u20132"),
          tags$span(style = .hp_chip_style("#2C7BB6"), "dark \u22652")
        )

        n_label <- if (n_total > 0 && n_total != (n_up + n_dn))
          paste0(" / ", n_total, " total") else ""

        tagList(
          tags$hr(),
          div(class = "px-2 pb-2",
              tags$h6(class = "fw-bold mb-1",
                      paste0(row$category[1], " \u2014 ",
                             n_up + n_dn, " DEGs", n_label)),
              legend,
              section("\u2191 Up-regulated",   n_up, up_chips, "#C0392B"),
              section("\u2193 Down-regulated", n_dn, dn_chips, "#2C7BB6"))
        )
      })
    }

    # ── Tab 1: Reactome ───────────────────────────────────────────────────────

    output$ro_bar <- .render_overlap_bar(data_ro, "ro_dir", "Reactome Pathways")
    output$ro_tbl <- .render_overlap_tbl(data_ro, allow_select = TRUE)

    sel_row_ro <- reactiveVal(NULL)
    observeEvent(input$ro_tbl_rows_selected, {
      idx <- input$ro_tbl_rows_selected; df <- data_ro()
      if (length(idx) == 1 && !is.null(df) && idx <= nrow(df))
        sel_row_ro(df[idx, ]) else sel_row_ro(NULL)
    })
    observeEvent(list(input$deg_assay, input$ro_dir, input$ro_topn, input$ro_minsize),
                 { sel_row_ro(NULL) }, ignoreInit = TRUE)

    ro_gene_info_r       <- .make_gene_info_r(sel_row_ro)
    output$ro_gene_panel <- .render_gene_chip_panel(sel_row_ro, ro_gene_info_r)

    # ── Tab 2: GO BP ──────────────────────────────────────────────────────────

    output$gobp_bar <- .render_overlap_bar(data_gobp, "gobp_dir", "GO Biological Process")
    output$gobp_tbl <- .render_overlap_tbl(data_gobp, allow_select = TRUE)

    sel_row_gobp <- reactiveVal(NULL)
    observeEvent(input$gobp_tbl_rows_selected, {
      idx <- input$gobp_tbl_rows_selected; df <- data_gobp()
      if (length(idx) == 1 && !is.null(df) && idx <= nrow(df))
        sel_row_gobp(df[idx, ]) else sel_row_gobp(NULL)
    })
    observeEvent(list(input$deg_assay, input$gobp_dir, input$gobp_topn, input$gobp_minsize),
                 { sel_row_gobp(NULL) }, ignoreInit = TRUE)

    gobp_gene_info_r       <- .make_gene_info_r(sel_row_gobp)
    output$gobp_gene_panel <- .render_gene_chip_panel(sel_row_gobp, gobp_gene_info_r)

    # ── Tab 2: GO MF ──────────────────────────────────────────────────────────

    output$gomf_bar <- .render_overlap_bar(data_gomf, "gomf_dir", "GO Molecular Function")
    output$gomf_tbl <- .render_overlap_tbl(data_gomf, allow_select = TRUE)

    sel_row_gomf <- reactiveVal(NULL)
    observeEvent(input$gomf_tbl_rows_selected, {
      idx <- input$gomf_tbl_rows_selected; df <- data_gomf()
      if (length(idx) == 1 && !is.null(df) && idx <= nrow(df))
        sel_row_gomf(df[idx, ]) else sel_row_gomf(NULL)
    })
    observeEvent(list(input$deg_assay, input$gomf_dir, input$gomf_topn, input$gomf_minsize),
                 { sel_row_gomf(NULL) }, ignoreInit = TRUE)

    gomf_gene_info_r       <- .make_gene_info_r(sel_row_gomf)
    output$gomf_gene_panel <- .render_gene_chip_panel(sel_row_gomf, gomf_gene_info_r)

    # ── Tab 3: PANTHER Class ──────────────────────────────────────────────────

    output$pc_bar <- .render_overlap_bar(data_pc, "pc_dir", "PANTHER Protein Class")
    output$pc_tbl <- .render_overlap_tbl(data_pc, allow_select = TRUE)

    sel_row_pc <- reactiveVal(NULL)
    observeEvent(input$pc_tbl_rows_selected, {
      idx <- input$pc_tbl_rows_selected; df <- data_pc()
      if (length(idx) == 1 && !is.null(df) && idx <= nrow(df))
        sel_row_pc(df[idx, ]) else sel_row_pc(NULL)
    })
    observeEvent(list(input$deg_assay, input$pc_dir, input$pc_topn, input$pc_minsize),
                 { sel_row_pc(NULL) }, ignoreInit = TRUE)

    pc_gene_info_r       <- .make_gene_info_r(sel_row_pc)
    output$pc_gene_panel <- .render_gene_chip_panel(sel_row_pc, pc_gene_info_r)

    # =========================================================================
    # Gene-level annotation summary banners (top of each tab)
    # =========================================================================

    .ko_gene_annots <- function(src) {
      ref_tbl <- .HP_REF_TBL[[src]]
      reactive({
        gene <- assay_gene_r(); req(gene)
        con  <- con_r()
        tryCatch({
          df <- dbGetQuery(con, sprintf(
            "SELECT DISTINCT category FROM %s WHERE symbol = '%s' LIMIT 8",
            ref_tbl, gsub("'", "''", gene)))
          df$category
        }, error = function(e) character(0))
      })
    }

    .pheno_gene_annots <- function(ds_r) {
      reactive({
        ds   <- ds_r(); req(ds)
        gene <- assay_gene_r(); req(gene)
        if (ds == "panelapp") return(character(0))
        ph_col   <- .PB_COL[[ds]]; tbl_name <- .PB_TABLE[[ds]]
        if (is.null(ph_col) || is.null(tbl_name)) return(character(0))
        con <- con_r()
        tryCatch({
          df <- dbGetQuery(con, sprintf(
            'SELECT DISTINCT "%s" AS name FROM "%s"
             WHERE symbol = \'%s\' AND "%s" IS NOT NULL',
            ph_col, tbl_name, gsub("'", "''", gene), ph_col))
          head(na.omit(df$name), 8)
        }, error = function(e) character(0))
      })
    }

    gs_ro_r   <- .ko_gene_annots("reactome")
    gs_gobp_r <- .ko_gene_annots("go_bp")
    gs_gomf_r <- .ko_gene_annots("go_mf")
    gs_pc_r   <- .ko_gene_annots("panther_class")
    gs_ph_r   <- .pheno_gene_annots(reactive("impc"))
    gs_hpo_r  <- .pheno_gene_annots(reactive("hpo"))
    gs_di_r   <- .pheno_gene_annots(reactive(input$di_ds))

    # ── Per-tab annotation chip strip ───────────────────────────────────────

    .render_gs <- function(annots_r, annot_label_suffix, badge_bg = "#495057") {
      renderUI({
        gene <- assay_gene_r()
        if (is.null(gene)) return(NULL)
        items <- annots_r()
        chips <- if (length(items) > 0) {
          lapply(items, function(nm) {
            tags$span(
              class = "badge rounded-pill fw-normal me-1 mb-1",
              style = sprintf(
                "background:%s; font-size:0.7rem; opacity:0.88; color:white;",
                badge_bg),
              nm)
          })
        } else {
          list(tags$span(class = "text-muted small fst-italic", "No annotations found"))
        }
        lbl <- paste0("Gene ", gene, " Associated ", annot_label_suffix)
        div(
          class = "d-flex align-items-center flex-wrap gap-1 px-3 py-2 mb-3",
          style = "background:#f0f4f8; border-left:3px solid #adb5bd; border-radius:0 4px 4px 0;",
          tags$span(class = "text-secondary fw-semibold me-2",
                    style = "font-size:0.75rem; white-space:nowrap;", lbl),
          tagList(chips)
        )
      })
    }

    output$gs_ro   <- .render_gs(gs_ro_r,   "Reactome pathways:",      "#4E79A7")
    output$gs_gobp <- .render_gs(gs_gobp_r, "GO Biological Process:",  "#59A14F")
    output$gs_gomf <- .render_gs(gs_gomf_r, "GO Molecular Function:",  "#2E86C1")
    output$gs_pc   <- .render_gs(gs_pc_r,   "PANTHER class:",          "#B07AA1")
    output$gs_ph   <- .render_gs(gs_ph_r,   "IMPC phenotypes:",        "#E63946")
    output$gs_hpo  <- .render_gs(gs_hpo_r,  "HPO phenotypes:",         "#AF7AC5")
    output$gs_di   <- .render_gs(gs_di_r,   "Disease annotations:",    "#457B9D")

    # =========================================================================
    # Phenotype browser logic (Tabs 4 + 5 non-PanelApp)
    # =========================================================================

    .load_pheno_tbl <- function(ds_r, zyg_r = NULL) {
      reactive({
        ds  <- ds_r(); req(ds)
        if (ds == "panelapp") return(NULL)
        con <- con_r()
        tbl_name <- .PB_TABLE[[ds]]
        if (is.null(tbl_name)) return(NULL)
        tryCatch({
          df <- dbReadTable(con, tbl_name)
          if (ds == "impc" && !is.null(zyg_r)) {
            zyg <- zyg_r(); req(length(zyg) > 0)
            df  <- df[!is.na(df$impc_zygosity) & df$impc_zygosity %in% zyg, ]
          }
          df
        }, error = function(e) NULL)
      })
    }

    ph_pheno_r  <- .load_pheno_tbl(reactive("impc"), reactive(input$ph_zyg))
    hpo_pheno_r <- .load_pheno_tbl(reactive("hpo"))
    di_pheno_r  <- .load_pheno_tbl(reactive(input$di_ds))

    .make_joined_r <- function(pheno_r, ds_r) {
      reactive({
        ds <- ds_r(); req(ds)
        if (ds == "panelapp") return(NULL)
        deg    <- deg_symbols_r(); req(!is.null(deg), nrow(deg) > 0)
        ph     <- pheno_r();       req(!is.null(ph))
        ph_col <- .PB_COL[[ds]];   req(!is.null(ph_col))
        keep   <- unique(c("symbol", ph_col, if (ds == "impc") "impc_zygosity"))
        keep   <- keep[keep %in% colnames(ph)]
        ph_sub <- ph[, keep, drop = FALSE]
        ph_sub <- ph_sub[!is.na(ph_sub[[ph_col]]) &
                           nzchar(trimws(as.character(ph_sub[[ph_col]]))), ]
        merged <- merge(deg, ph_sub, by = "symbol", all = FALSE)
        if (nrow(merged) == 0L) return(NULL)
        merged
      })
    }

    ph_joined_r  <- .make_joined_r(ph_pheno_r,  reactive("impc"))
    hpo_joined_r <- .make_joined_r(hpo_pheno_r, reactive("hpo"))
    di_joined_r  <- .make_joined_r(di_pheno_r,  reactive(input$di_ds))

    .make_coverage_r <- function(pheno_r, ds_r) {
      reactive({
        ds <- ds_r(); req(ds)
        if (ds == "panelapp") return(NULL)
        deg <- deg_symbols_r(); req(!is.null(deg), nrow(deg) > 0)
        ph  <- pheno_r();       req(!is.null(ph))
        ann_syms <- unique(ph$symbol[!is.na(ph$symbol) & nzchar(ph$symbol)])
        data.frame(symbol = deg$symbol, DEG = deg$DEG,
                   annotated = deg$symbol %in% ann_syms) |>
          dplyr::group_by(DEG, annotated) |>
          dplyr::summarise(n = dplyr::n_distinct(symbol), .groups = "drop") |>
          dplyr::mutate(
            dir_lbl = ifelse(DEG == "up", "\u2191 Up", "\u2193 Down"),
            ann_lbl = ifelse(annotated, "Annotated", "No annotation"))
      })
    }

    ph_coverage_r  <- .make_coverage_r(ph_pheno_r,  reactive("impc"))
    hpo_coverage_r <- .make_coverage_r(hpo_pheno_r, reactive("hpo"))
    di_coverage_r  <- .make_coverage_r(di_pheno_r,  reactive(input$di_ds))

    .make_top_pheno_r <- function(joined_r, ds_r, topn_r, dir_r) {
      reactive({
        ds <- ds_r(); req(ds)
        if (ds == "panelapp") return(NULL)
        jd     <- joined_r(); req(!is.null(jd), nrow(jd) > 0)
        ph_col <- .PB_COL[[ds]]
        n_top  <- topn_r(); dir <- dir_r()
        if (dir != "both") { jd <- jd[jd$DEG == dir, ]; req(nrow(jd) > 0) }
        grp_cols <- unique(c(ph_col, "DEG", if (ds == "impc") "impc_zygosity"))
        summ <- jd |>
          dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) |>
          dplyr::summarise(
            n_genes = dplyr::n_distinct(symbol),
            genes   = paste(sort(unique(symbol)), collapse = "; "),
            .groups = "drop")
        top_labels <- summ |>
          dplyr::group_by(.data[[ph_col]]) |>
          dplyr::summarise(total = sum(n_genes), .groups = "drop") |>
          dplyr::slice_max(total, n = n_top, with_ties = FALSE) |>
          dplyr::arrange(total) |>
          dplyr::pull(.data[[ph_col]])
        list(summ = summ, top_labels = top_labels, ph_col = ph_col)
      })
    }

    ph_top_r  <- .make_top_pheno_r(ph_joined_r,  reactive("impc"),
                                    reactive(input$ph_topn),  reactive(input$ph_dir))
    hpo_top_r <- .make_top_pheno_r(hpo_joined_r, reactive("hpo"),
                                    reactive(input$hpo_topn), reactive(input$hpo_dir))
    di_top_r  <- .make_top_pheno_r(di_joined_r,  reactive(input$di_ds),
                                    reactive(input$di_topn),  reactive(input$di_dir))

    # ── Wide-data reactives (shared between DT render and gene chip) ──────────

    # Produces a data frame with raw genes_up / genes_down columns
    .make_pheno_wide_r <- function(top_r, ds_r) {
      reactive({
        td <- top_r(); req(!is.null(td))
        ds <- ds_r()
        summ   <- td$summ; labels <- td$top_labels; ph_col <- td$ph_col
        top    <- summ[summ[[ph_col]] %in% labels, ]
        if (ds == "impc") {
          top |>
            dplyr::group_by(.data[[ph_col]], impc_zygosity) |>
            dplyr::summarise(
              n_up       = sum(n_genes[DEG == "up"],   na.rm = TRUE),
              n_down     = sum(n_genes[DEG == "down"], na.rm = TRUE),
              n_total    = n_up + n_down,
              genes_up   = paste(na.omit(genes[DEG == "up"]),   collapse = "; "),
              genes_down = paste(na.omit(genes[DEG == "down"]), collapse = "; "),
              .groups    = "drop") |>
            dplyr::arrange(dplyr::desc(n_total))
        } else {
          top |>
            dplyr::group_by(.data[[ph_col]]) |>
            dplyr::summarise(
              n_up       = sum(n_genes[DEG == "up"],   na.rm = TRUE),
              n_down     = sum(n_genes[DEG == "down"], na.rm = TRUE),
              n_total    = n_up + n_down,
              genes_up   = paste(na.omit(genes[DEG == "up"]),   collapse = "; "),
              genes_down = paste(na.omit(genes[DEG == "down"]), collapse = "; "),
              .groups    = "drop") |>
            dplyr::arrange(dplyr::desc(n_total))
        }
      })
    }

    ph_wide_r  <- .make_pheno_wide_r(ph_top_r,  reactive("impc"))
    hpo_wide_r <- .make_pheno_wide_r(hpo_top_r, reactive("hpo"))
    di_wide_r  <- .make_pheno_wide_r(di_top_r,  reactive(input$di_ds))

    # Render DT from the shared wide reactive (renames columns for display only)
    .render_pheno_tbl <- function(wide_r, ds_r) {
      DT::renderDT({
        wide <- wide_r(); req(!is.null(wide), nrow(wide) > 0)
        ds   <- ds_r()
        ph_col <- .PB_COL[[ds]]; ph_lbl <- .PB_LABEL[[ds]]
        display <- as.data.frame(wide)
        names(display)[names(display) == ph_col]          <- ph_lbl
        if (ds == "impc")
          names(display)[names(display) == "impc_zygosity"] <- "Zygosity"
        names(display) <- gsub("^n_up$",      "\u2191 N Up",       names(display))
        names(display) <- gsub("^n_down$",    "\u2193 N Down",     names(display))
        names(display) <- gsub("^n_total$",   "N Total",           names(display))
        names(display) <- gsub("^genes_up$",  "Genes \u2191 Up",   names(display))
        names(display) <- gsub("^genes_down$","Genes \u2193 Down", names(display))
        DT::datatable(
          display, class = "compact row-border hover", rownames = FALSE,
          selection  = "single",
          extensions = "Buttons",
          options    = list(pageLength = 25, scrollX = TRUE,
                            dom = "Bfrtip", buttons = c("csv", "excel"),
                            initComplete = .dt_header_js))
      })
    }

    output$ph_tbl  <- .render_pheno_tbl(ph_wide_r,  reactive("impc"))
    output$hpo_tbl <- .render_pheno_tbl(hpo_wide_r, reactive("hpo"))
    output$di_tbl  <- .render_pheno_tbl(di_wide_r,  reactive(input$di_ds))

    # ── Phenotype tab gene chip ───────────────────────────────────────────────

    sel_row_ph <- reactiveVal(NULL)
    observeEvent(input$ph_tbl_rows_selected, {
      idx <- input$ph_tbl_rows_selected
      df  <- ph_wide_r()
      if (length(idx) == 1 && !is.null(df) && idx <= nrow(df)) {
        row <- as.data.frame(df[idx, ])
        ph_col <- .PB_COL[["impc"]]
        row$category  <- as.character(row[[ph_col]])
        row$n_pathway <- as.integer(row$n_total)
        sel_row_ph(row)
      } else {
        sel_row_ph(NULL)
      }
    })
    observeEvent(list(input$deg_assay, input$ph_dir, input$ph_topn, input$ph_zyg),
                 { sel_row_ph(NULL) }, ignoreInit = TRUE)

    ph_gene_info_r       <- .make_gene_info_r(sel_row_ph)
    output$ph_gene_panel <- .render_gene_chip_panel(sel_row_ph, ph_gene_info_r)

    # ── HPO tab gene chip ─────────────────────────────────────────────────

    sel_row_hpo <- reactiveVal(NULL)
    observeEvent(input$hpo_tbl_rows_selected, {
      idx <- input$hpo_tbl_rows_selected
      df  <- hpo_wide_r()
      if (length(idx) == 1 && !is.null(df) && idx <= nrow(df)) {
        row <- as.data.frame(df[idx, ])
        ph_col <- .PB_COL[["hpo"]]
        row$category  <- as.character(row[[ph_col]])
        row$n_pathway <- as.integer(row$n_total)
        sel_row_hpo(row)
      } else {
        sel_row_hpo(NULL)
      }
    })
    observeEvent(list(input$deg_assay, input$hpo_dir, input$hpo_topn),
                 { sel_row_hpo(NULL) }, ignoreInit = TRUE)

    hpo_gene_info_r       <- .make_gene_info_r(sel_row_hpo)
    output$hpo_gene_panel <- .render_gene_chip_panel(sel_row_hpo, hpo_gene_info_r)

    # ── Disease tab (OMIM / Orphanet) gene chip ───────────────────────────────

    sel_row_di <- reactiveVal(NULL)
    observeEvent(input$di_tbl_rows_selected, {
      idx <- input$di_tbl_rows_selected
      df  <- di_wide_r()
      if (length(idx) == 1 && !is.null(df) && idx <= nrow(df)) {
        row <- as.data.frame(df[idx, ])
        ph_col <- .PB_COL[[isolate(input$di_ds)]]
        row$category  <- as.character(row[[ph_col]])
        row$n_pathway <- as.integer(row$n_total)
        sel_row_di(row)
      } else {
        sel_row_di(NULL)
      }
    })
    observeEvent(list(input$deg_assay, input$di_ds, input$di_dir, input$di_topn),
                 { sel_row_di(NULL) }, ignoreInit = TRUE)

    di_gene_info_r       <- .make_gene_info_r(sel_row_di)
    output$di_gene_panel <- .render_gene_chip_panel(sel_row_di, di_gene_info_r)

    # ── Coverage + top-N charts ───────────────────────────────────────────────

    .render_coverage <- function(cov_r, dir_r) {
      renderPlotly({
        cov <- cov_r(); req(!is.null(cov), nrow(cov) > 0)
        dir <- dir_r(); if (is.null(dir)) dir <- "both"
        if (dir != "both") cov <- cov[cov$DEG == dir, ]
        fig <- plot_ly()
        for (a in c(FALSE, TRUE)) {
          sub <- cov[cov$annotated == a, ]; if (nrow(sub) == 0L) next
          lbl <- if (a) "Annotated" else "No annotation"
          fig <- fig |>
            add_bars(data = sub, x = ~n, y = ~dir_lbl, name = lbl, orientation = "h",
                     marker = list(color = .PB_ANN_CLRS[[lbl]]),
                     hovertemplate = paste0("<b>%{y}</b><br>", lbl, ": %{x} DEGs<extra></extra>"))
        }
        fig |>
          layout(barmode = "stack",
                 xaxis   = list(title = "N DEGs"),
                 yaxis   = list(title = "", categoryorder = "array",
                                categoryarray = c("\u2193 Down", "\u2191 Up")),
                 legend  = list(orientation = "h", x = 0, y = -0.45),
                 margin  = list(l = 70, t = 5, b = 80, r = 10)) |>
          config(displayModeBar = FALSE)
      })
    }

    output$ph_coverage  <- .render_coverage(ph_coverage_r,  reactive(input$ph_dir))
    output$hpo_coverage <- .render_coverage(hpo_coverage_r, reactive(input$hpo_dir))
    output$di_coverage  <- .render_coverage(di_coverage_r,  reactive(input$di_dir))

    .render_top_pheno <- function(top_r, ds_r, dir_r) {
      renderPlotly({
        td <- top_r(); req(!is.null(td))
        ds     <- ds_r()
        summ   <- td$summ; labels <- td$top_labels; ph_col <- td$ph_col
        dir    <- dir_r(); if (is.null(dir)) dir <- "both"

        build_bars <- function(fig, df, ph_labels, show_legend = TRUE) {
          df$phenotype <- factor(df[[ph_col]], levels = ph_labels)
          df <- df[!is.na(df$phenotype), ]
          dirs <- if (dir == "both") c("down", "up") else dir
          for (d in dirs) {
            sub <- df[df$DEG == d, ]; if (nrow(sub) == 0L) next
            nm  <- if (d == "up") "\u2191 Up" else "\u2193 Down"
            fig <- fig |>
              add_bars(data = sub, x = ~n_genes, y = ~phenotype, name = nm,
                       orientation = "h", marker = list(color = .PB_CLRS[[d]]),
                       legendgroup = nm, showlegend = show_legend,
                       hovertemplate = paste0("<b>%{y}</b><br>", nm,
                                              ": %{x} gene(s)<extra></extra>"))
          }
          fig
        }

        if (ds == "impc") {
          zyg_lvls <- sort(unique(summ$impc_zygosity))
          subs <- lapply(seq_along(zyg_lvls), function(i) {
            z    <- zyg_lvls[[i]]
            df_z <- summ[summ$impc_zygosity == z & summ[[ph_col]] %in% labels, ]
            z_labels <- df_z |>
              dplyr::group_by(.data[[ph_col]]) |>
              dplyr::summarise(total = sum(n_genes), .groups = "drop") |>
              dplyr::arrange(total) |> dplyr::pull(.data[[ph_col]])
            build_bars(plot_ly(), df_z, ph_labels = z_labels, show_legend = (i == 1L)) |>
              layout(barmode = "stack",
                     xaxis = list(title = "N DEGs"),
                     yaxis = list(title = "", tickfont = list(size = 8)),
                     annotations = list(list(
                       text = tools::toTitleCase(z), showarrow = FALSE,
                       xref = "paper", yref = "paper", x = 0.5, y = 1.05,
                       xanchor = "center", font = list(size = 11, color = "#444"))))
          })
          subplot(subs, nrows = length(zyg_lvls), shareX = FALSE, shareY = FALSE,
                  titleX = TRUE, titleY = TRUE) |>
            layout(legend = list(orientation = "h", x = 0, y = -0.04),
                   margin = list(l = 220, t = 50, b = 40, r = 10))
        } else {
          df <- summ[summ[[ph_col]] %in% labels, ]
          build_bars(plot_ly(), df, ph_labels = labels) |>
            layout(barmode = "stack",
                   xaxis  = list(title = "N DEGs with annotation"),
                   yaxis  = list(title = "", tickfont = list(size = 9)),
                   legend = list(orientation = "h", x = 0, y = -0.08),
                   margin = list(l = 220, t = 20, b = 60, r = 10))
        }
      })
    }

    output$ph_top  <- .render_top_pheno(ph_top_r,  reactive("impc"),        reactive(input$ph_dir))
    output$hpo_top <- .render_top_pheno(hpo_top_r, reactive("hpo"),         reactive(input$hpo_dir))
    output$di_top  <- .render_top_pheno(di_top_r,  reactive(input$di_ds),   reactive(input$di_dir))

    # =========================================================================
    # Suspend all outputs when hidden
    # =========================================================================

    invisible(lapply(
      c("deg_no_assays", "deg_match_count", "deg_filter_badge", "deg_no_selection", "deg_meta_above",
        "gs_ro", "gs_gobp", "gs_gomf", "gs_pc", "gs_ph", "gs_hpo", "gs_di",
        "ro_bar", "ro_tbl", "ro_gene_panel",
        "gobp_bar", "gobp_tbl", "gobp_gene_panel",
        "gomf_bar", "gomf_tbl", "gomf_gene_panel",
        "pc_bar", "pc_tbl", "pc_gene_panel",
        "ph_coverage", "ph_top", "ph_tbl", "ph_gene_panel",
        "hpo_coverage", "hpo_top", "hpo_tbl", "hpo_gene_panel",
        "di_coverage", "di_top", "di_tbl", "di_gene_panel"),
      function(oid) outputOptions(output, oid, suspendWhenHidden = TRUE)
    ))

  })
}
