# =============================================================================
# Module: assay_browser
# Description: Interactive browser for all MorPhiC assay metadata.
#   • Filter row: KO Gene / Model System / KO Strategy (multi-select)
#   • Value boxes that update with current filters
#   • Rich DT grouped by KO Gene — shows all assay details
#   • Click a row → detail card showing full assay metadata
# =============================================================================

library(DT)
library(bslib)
library(bsicons)
library(dplyr)

# ── UI ────────────────────────────────────────────────────────────────────────

assay_browserUI <- function(id) {
  ns <- NS(id)

  tagList(

    # Lift selectize dropdowns above card boundaries
    tags$style(HTML(".ab-filters .selectize-dropdown { z-index: 1080 !important; }")),

    # ── Filter row ──────────────────────────────────────────────────────────
    div(
      class = "ab-filters",
      layout_columns(
        col_widths = c(4, 4, 4),

        card(
          style = "overflow: visible;",
          card_header(tagList(bs_icon("diagram-3"), " KO Gene")),
          selectInput(ns("f_gene"), NULL,
                      choices  = NULL,
                      multiple = TRUE, selected = NULL, width = "100%")
        ),
        card(
          style = "overflow: visible;",
          card_header(tagList(bs_icon("grid-3x3-gap-fill"), " Model System")),
          selectInput(ns("f_model"), NULL,
                      choices  = NULL,
                      multiple = TRUE, selected = NULL, width = "100%")
        ),
        card(
          style = "overflow: visible;",
          card_header(tagList(bs_icon("scissors"), " KO Strategy")),
          selectInput(ns("f_kostrat"), NULL,
                      choices  = NULL,
                      multiple = TRUE, selected = NULL, width = "100%")
        )
      )
    ),

    # ── Value boxes ──────────────────────────────────────────────────────────
    layout_columns(
      col_widths = c(3, 3, 3, 3),

      value_box(
        title    = "Assays",
        value    = textOutput(ns("vb_assays")),
        showcase = bs_icon("bar-chart-fill"),
        theme    = "primary"
      ),
      value_box(
        title    = "KO Genes",
        value    = textOutput(ns("vb_genes")),
        showcase = bs_icon("diagram-3"),
        theme    = "secondary"
      ),
      value_box(
        title    = "Model Systems",
        value    = textOutput(ns("vb_models")),
        showcase = bs_icon("grid-3x3-gap-fill"),
        theme    = "info"
      ),
      value_box(
        title    = "KO Strategies",
        value    = textOutput(ns("vb_strats")),
        showcase = bs_icon("scissors"),
        theme    = "dark"
      )
    ),

    # ── Main table + detail card ──────────────────────────────────────────
    layout_columns(
      col_widths = c(8, 4),

      card(
        full_screen = TRUE,
        card_header(
          class = "d-flex justify-content-between align-items-center",
          tagList(bs_icon("table"), " Assay Metadata"),
          actionLink(ns("clear_filters"), "Clear filters",
                     class = "small text-muted",
                     icon  = icon("xmark"))
        ),
        DTOutput(ns("assay_tbl"))
      ),

      card(
        card_header(tagList(bs_icon("card-list"), " Assay Details")),
        uiOutput(ns("detail_panel"))
      )
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

assay_browserServer <- function(id, study_info_r) {
  moduleServer(id, function(input, output, session) {

    # ── Populate filter choices on first load ──────────────────────────────
    observe({
      si <- study_info_r()
      updateSelectInput(session, "f_gene",
                        choices  = sort(unique(na.omit(si$Gene))),
                        selected = character(0))
      updateSelectInput(session, "f_model",
                        choices  = sort(unique(na.omit(si$Model_system))),
                        selected = character(0))
      updateSelectInput(session, "f_kostrat",
                        choices  = sort(unique(na.omit(si$KO_strat))),
                        selected = character(0))
    })

    observeEvent(input$clear_filters, {
      updateSelectInput(session, "f_gene",    selected = character(0))
      updateSelectInput(session, "f_model",   selected = character(0))
      updateSelectInput(session, "f_kostrat", selected = character(0))
    })

    # ── Filtered study_info ────────────────────────────────────────────────
    filtered_si <- reactive({
      si <- study_info_r()
      if (length(input$f_gene)    > 0) si <- si |> filter(Gene         %in% input$f_gene)
      if (length(input$f_model)   > 0) si <- si |> filter(Model_system %in% input$f_model)
      if (length(input$f_kostrat) > 0) si <- si |> filter(KO_strat     %in% input$f_kostrat)
      si
    })

    # ── Value boxes ────────────────────────────────────────────────────────
    output$vb_assays  <- renderText({ n_distinct(filtered_si()$Assay)        })
    output$vb_genes   <- renderText({ n_distinct(na.omit(filtered_si()$Gene)) })
    output$vb_models  <- renderText({ n_distinct(na.omit(filtered_si()$Model_system)) })
    output$vb_strats  <- renderText({ n_distinct(na.omit(filtered_si()$KO_strat))     })

    # ── Build display table ────────────────────────────────────────────────
    display_tbl <- reactive({
      keep_cols <- c("Gene", "Model_system", "KO_strat", "DPC",
                     "Cell_Line", "Replicate", "Assay")
      si <- filtered_si()
      cols <- keep_cols[keep_cols %in% colnames(si)]
      si <- si[!duplicated(si$Assay), cols, drop = FALSE]
      si <- si[order(si$Gene, si$Model_system, si$Assay), ]
      si
    })

    # ── Render DT ─────────────────────────────────────────────────────────
    output$assay_tbl <- DT::renderDT({
      tbl <- display_tbl()
      req(nrow(tbl) > 0)

      pretty_names <- gsub("_", " ", colnames(tbl))

      DT::datatable(
        tbl,
        colnames   = pretty_names,
        rownames   = FALSE,
        selection  = "single",
        filter     = "top",
        extensions = c("RowGroup", "Buttons"),
        options    = list(
          rowGroup       = list(dataSrc = 0),        # group rows by Gene col (index 0)
          columnDefs     = list(
            list(visible = FALSE, targets = 0),      # hide Gene column in body (shown as group header)
            list(className = "dt-left", targets = "_all")
          ),
          pageLength     = 25,
          lengthMenu     = list(c(25, 50, 100, -1), c("25", "50", "100", "All")),
          scrollX        = TRUE,
          scrollY        = "520px",
          scrollCollapse = TRUE,
          dom            = "Bft",
          buttons        = list("csv", "excel")
        ),
        class = "display compact"
      )
    })

    # ── Detail panel on row click ──────────────────────────────────────────
    output$detail_panel <- renderUI({
      row_idx <- input$assay_tbl_rows_selected
      if (is.null(row_idx) || length(row_idx) == 0) {
        return(
          div(
            class = "text-muted small p-3",
            bs_icon("hand-index"), " Click any row in the table to see full assay details."
          )
        )
      }

      tbl <- display_tbl()
      row <- tbl[row_idx, , drop = FALSE]

      # Full row from study_info (may have more fields)
      full <- study_info_r() |> filter(Assay == row$Assay[1]) |> head(1)

      # Helper to render one labelled detail item
      detail_item <- function(label, value, badge_class = "bg-secondary") {
        if (is.null(value) || is.na(value) || !nzchar(as.character(value))) return(NULL)
        div(class = "mb-2",
          tags$span(class = "text-muted small me-1", paste0(label, ":")),
          tags$span(class = paste("badge", badge_class), as.character(value))
        )
      }

      tagList(
        div(class = "p-3",
          tags$h6(class = "fw-bold text-truncate mb-3",
                  style = "font-size:0.8rem; word-break:break-all;",
                  full$Assay),
          tags$hr(class = "my-2"),
          detail_item("KO Gene",       full$Gene,          "bg-success"),
          detail_item("Model System",  full$Model_system,  "bg-primary"),
          detail_item("KO Strategy",   full$KO_strat,      "bg-info text-dark"),
          detail_item("Cell Line",     full$Cell_Line,     "bg-secondary"),
          detail_item("DPC",           full$DPC,           "bg-warning text-dark"),
          if ("DPC_days" %in% names(full))
            detail_item("DPC days",    full$DPC_days,      "bg-light text-dark border"),
          if ("Replicate" %in% names(full))
            detail_item("Replicate",   full$Replicate,     "bg-light text-dark border"),
          if ("Differentation_time_point" %in% names(full))
            detail_item("Diff. time",  full$Differentation_time_point, "bg-light text-dark border")
        )
      )
    })

    # Suspend all outputs until tab is visited
    invisible(lapply(
      c("vb_assays", "vb_genes", "vb_models", "vb_strats",
        "assay_tbl", "detail_panel"),
      function(oid) outputOptions(output, oid, suspendWhenHidden = TRUE)
    ))

  })
}
