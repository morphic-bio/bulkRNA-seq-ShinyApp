# =============================================================================
# Module: compare_studies
# Page 2 — Compare Studies
#
# Sections:
#   5. Assay Comparison — filter + picker, metadata header, tabs:
#                         Gene Matrix | Overlapping DEGs
#
# Requires helpers.R to be sourced first (.HP_CLRS, load_deg_genes,
#   .LFC_ROW_CB, .fmt_dir_cell).
# Data contract:
#   con_r        — reactive: open read-only DuckDB connection
#   study_info_r — reactive: study_info data.frame
# =============================================================================

library(plotly)
library(DT)
library(UpSetR)
library(ggplot2)
library(shinyWidgets)

# ── UI ────────────────────────────────────────────────────────────────────────

compare_studiesUI <- function(id) {
  ns <- NS(id)

  tagList(
    card(
      card_header(tagList(bsicons::bs_icon("grid-3x3-gap-fill"), " Assay Comparison",
                          .info_tip("Compare differentially expressed genes across multiple assays side-by-side."))),
      card_body(
        # ── Top bar: picker + filters (left) | selected assay cards (right) ──
        layout_columns(
          col_widths = c(6, 6),

          # Left: assay picker + controls row
          div(
            tags$span(class = "fw-semibold small mb-1 d-block", "Select assays"),
            selectizeInput(ns("cmp_assays"), NULL,
                           choices = NULL, selected = NULL, multiple = TRUE,
                           options = list(placeholder = "Choose assays\u2026",
                                          plugins     = list("remove_button"),
                                          maxItems    = NULL),
                           width = "100%"),
            uiOutput(ns("cmp_no_assays")),
            div(class = "d-flex align-items-center gap-2 mt-1",
                actionButton(ns("cmp_clear_assays"),
                             tagList(bsicons::bs_icon("x-lg", size = "0.6rem"), " Clear selection"),
                             class = "btn btn-outline-secondary btn-sm",
                             style = "font-size: 0.7rem;"),
                uiOutput(ns("cmp_match_count"), inline = TRUE),
                dropdownButton(
                  tags$h6(class = "fw-semibold mb-3", "Filter assays"),
                  selectInput(ns("cmp_filter_gene"),  "KO Gene",
                              choices = NULL, multiple = TRUE, width = "100%"),
                  selectInput(ns("cmp_filter_model"), "Model System",
                              choices = NULL, multiple = TRUE, width = "100%"),
                  selectInput(ns("cmp_filter_comp"),  "Comparison",
                              choices = NULL, multiple = TRUE, width = "100%"),
                  actionButton(ns("cmp_clear_filters"), "Clear all filters",
                               class = "btn btn-outline-secondary btn-sm w-100 mt-1"),
                  circle = FALSE, status = "default", size = "sm",
                  icon = tagList(bsicons::bs_icon("funnel", size = "0.75rem"),
                                 uiOutput(ns("cmp_filter_badge"), inline = TRUE)),
                  label = "Filters", width = "320px",
                  inline = TRUE,
                  inputId = ns("cmp_filter_dropdown")
                )
            )
          ),

          # Right: selected assay metadata cards
          uiOutput(ns("cmp_meta_above"))
        ),

        tags$hr(class = "my-2"),

        # ── Results ─────────────────────────────────────────────────────
        navset_tab(
          nav_panel(tagList("Gene Matrix",
                          .info_tip("Matrix showing DEG direction per gene across selected assays. Sorted by number of assays.")),
                    div(class = "pt-2",
                        uiOutput(ns("cmp_matrix_ui")))),
          nav_panel(tagList("Overlapping DEGs",
                          .info_tip("Intersection plot showing overlap and unique DEG sets across selected assays.")),
                    div(class = "pt-2",
                        uiOutput(ns("cmp_upset_ui")))),
          nav_spacer(),
          .dl_dropdown_nav(ns, "cmp"),
          nav_item(
            dropdownButton(
              radioButtons(ns("cmp_dir"), "Direction",
                           choices  = c("\u2191\u2193 All DEGs" = "all",
                                        "\u2191 Up" = "up",
                                        "\u2193 Down" = "down"),
                           selected = "all", inline = TRUE),
              numericInput(ns("cmp_lfc_min"), HTML("Min log<sub>2</sub>FC"),
                           value = NA, step = 0.5, width = "100%"),
              numericInput(ns("cmp_lfc_max"), HTML("Max log<sub>2</sub>FC"),
                           value = NA, step = 0.5, width = "100%"),
              actionButton(ns("cmp_clear_deg_filters"), "Clear filters",
                           class = "btn btn-outline-secondary btn-sm w-100 mt-2"),
              circle = FALSE, status = "outline-secondary", size = "sm",
              icon = bsicons::bs_icon("gear", size = "0.8rem"),
              label = "Options", width = "300px",
              inputId = ns("cmp_opts_dd")
            )
          )
        )
      ),
      full_screen = TRUE
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

compare_studiesServer <- function(id, con_r, study_info_r) {
  moduleServer(id, function(input, output, session) {

    # ── Populate filter choices ───────────────────────────────────────────────
    observe({
      si <- study_info_r()
      updateSelectInput(session, "cmp_filter_gene",
                        choices  = sort(unique(na.omit(si$Gene))),
                        selected = character(0))
      updateSelectInput(session, "cmp_filter_model",
                        choices  = sort(unique(na.omit(si$Model_system))),
                        selected = character(0))
      updateSelectInput(session, "cmp_filter_comp",
                        choices  = sort(unique(na.omit(si$Comparison))),
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
      if (length(input$cmp_filter_comp)  > 0)
        si_uniq <- si_uniq[si_uniq$Comparison   %in% input$cmp_filter_comp, ]
      sort(si_uniq$Assay)
    })

    observeEvent(filtered_cmp_assays(), {
      ch       <- filtered_cmp_assays()
      selected <- isolate(input$cmp_assays)
      # Keep currently selected assays even if they don't match filters
      all_choices <- sort(unique(c(selected, ch)))
      updateSelectizeInput(session, "cmp_assays",
                           choices = all_choices, selected = selected, server = TRUE)
    })

    output$cmp_no_assays <- renderUI({
      if (length(filtered_cmp_assays()) == 0)
        tags$p(class = "text-muted small fst-italic mb-0 mt-1",
               "No assays match the current filters.")
    })

    # ── Match count badge ───────────────────────────────────────────────────
    output$cmp_match_count <- renderUI({
      n     <- length(filtered_cmp_assays())
      total <- length(unique(study_info_r()$Assay))
      tags$span(
        class = "badge rounded-pill",
        style = "font-size: 0.7rem; background: #6c757d; color: white;",
        paste0(n, "/", total, " assays")
      )
    })

    # ── Active filter count badge ─────────────────────────────────────────────
    output$cmp_filter_badge <- renderUI({
      n <- sum(
        length(input$cmp_filter_gene)  > 0,
        length(input$cmp_filter_model) > 0,
        length(input$cmp_filter_comp)  > 0
      )
      if (n > 0) {
        tags$span(class = "badge rounded-pill bg-primary",
                  style = "font-size:0.65rem;", n)
      }
    })

    # ── Clear all filters ────────────────────────────────────────────────────
    observeEvent(input$cmp_clear_filters, {
      updateSelectInput(session, "cmp_filter_gene",  selected = character(0))
      updateSelectInput(session, "cmp_filter_model", selected = character(0))
      updateSelectInput(session, "cmp_filter_comp",  selected = character(0))
    })

    # ── Clear assay selection ──────────────────────────────────────────────────
    observeEvent(input$cmp_clear_assays, {
      ch <- filtered_cmp_assays()
      updateSelectizeInput(session, "cmp_assays",
                           choices = ch, selected = character(0), server = TRUE)
    })

    # ── Clear DEG filters (direction + LFC) ──────────────────────────────────
    observeEvent(input$cmp_clear_deg_filters, {
      updateRadioButtons(session, "cmp_dir", selected = "all")
      updateNumericInput(session, "cmp_lfc_min", value = NA)
      updateNumericInput(session, "cmp_lfc_max", value = NA)
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
                               "Condition"),
                             colnames(study_info_r()))
      meta <- study_info_r()[!duplicated(study_info_r()$Assay), meta_cols]
      merge(rows, meta, by.x = "assay", by.y = "Assay", all.x = TRUE)
    })

    # ── Shared filtered data (direction + LFC) ─────────────────────────────────
    cmp_sub_data <- reactive({
      genes_df <- cmp_genes_data(); req(genes_df)
      assays   <- input$cmp_assays; req(assays, length(assays) >= 2)
      dir      <- input$cmp_dir

      sub <- genes_df[genes_df$assay %in% assays, ]

      # Direction filter
      if (dir == "up")   sub <- sub[sub$DEG == "up", ]
      if (dir == "down") sub <- sub[sub$DEG == "down", ]

      # LFC filter (signed log2FC range)
      lfc_min <- input$cmp_lfc_min
      lfc_max <- input$cmp_lfc_max
      if (!is.null(lfc_min) && !is.na(lfc_min))
        sub <- sub[!is.na(sub$log2fc) & sub$log2fc >= lfc_min, ]
      if (!is.null(lfc_max) && !is.na(lfc_max))
        sub <- sub[!is.na(sub$log2fc) & sub$log2fc <= lfc_max, ]

      sub
    })

    # ── Metadata header ───────────────────────────────────────────────────────
    output$cmp_meta_above <- renderUI({
      assays <- input$cmp_assays; req(assays, length(assays) >= 1)
      si <- study_info_r()

      info_row <- function(label, val, pool) {
        if (is.null(val) || is.na(val) || !nzchar(as.character(val))) return(NULL)
        tags$div(class = "small d-flex align-items-center gap-1", style = "font-size:0.75rem;",
                 tags$span(class = "text-muted", paste0(label, ":")),
                 .meta_badge(val, pool))
      }

      cards <- lapply(assays, function(a) {
        r <- si[si$Assay == a, ][1, ]
        div(class = "border rounded p-2",
            style = "background:#f8f9fa; min-width:180px; flex:1;",
            tags$div(class = "small fw-semibold text-break mb-1",
                     style = "font-size:0.68rem; word-break:break-all; color:#495057;", a),
            info_row("Gene",         r$Gene,         si$Gene),
            info_row("Model System", r$Model_system, si$Model_system),
            info_row("Comparison",   r$Comparison,   si$Comparison),
            info_row("DPC",          r$DPC,          si$DPC),
            info_row("Cell Line",    r$Cell_Line,    si$Cell_Line),
            info_row("Condition",    r$Condition,    si$Condition),
            if ("Differentation_time_point" %in% names(r))
              info_row("Differentiation", r$Differentation_time_point,
                       si$Differentation_time_point),
            if ("Replicate" %in% names(r))
              info_row("Replicate", r$Replicate, si$Replicate))
      })
      cards_div <- div(class = "d-flex flex-wrap gap-2", cards)

      if (length(assays) > 3) {
        collapse_id <- session$ns("cmp_meta_collapse")
        tagList(
          div(class = "d-flex align-items-center gap-2 mb-1",
              tags$a(class = "fw-semibold small text-decoration-none d-flex align-items-center gap-1",
                     href = paste0("#", collapse_id),
                     `data-bs-toggle` = "collapse",
                     role = "button", `aria-expanded` = "false",
                     `aria-controls` = collapse_id,
                     bsicons::bs_icon("chevron-right", size = "0.65rem",
                                      class = "collapse-chevron"),
                     paste0("Selected Assays (", length(assays), ")"))),
          div(id = collapse_id, class = "collapse", cards_div),
          tags$style(HTML(sprintf(
            "#%s.collapse.show ~ * .collapse-chevron,
             #%s.collapsing ~ * .collapse-chevron { /* no-op */ }
             [aria-expanded='true'] .collapse-chevron {
               transform: rotate(90deg);
               transition: transform 0.2s ease;
             }
             [aria-expanded='false'] .collapse-chevron {
               transform: rotate(0deg);
               transition: transform 0.2s ease;
             }", collapse_id, collapse_id)))
        )
      } else {
        tagList(
          tags$span(class = "fw-semibold small d-block mb-1", "Selected Assays"),
          cards_div
        )
      }
    })

    # ── Gene matrix (placeholder / DT switcher) ─────────────────────────────
    output$cmp_matrix_ui <- renderUI({
      if (is.null(input$cmp_assays) || length(input$cmp_assays) < 2) {
        .empty_state("grid-3x3-gap-fill",
                     "Select 2 or more assays to view DEG comparison",
                     "Use the filters and assay picker above to get started.")
      } else {
        DTOutput(session$ns("cmp_matrix_tbl"))
      }
    })

    # ── Gene matrix ───────────────────────────────────────────────────────────
    output$cmp_matrix_tbl <- DT::renderDT({
      sub    <- cmp_sub_data(); req(sub, nrow(sub) > 0)
      assays <- input$cmp_assays; req(assays, length(assays) >= 2)
      all_syms <- sort(unique(na.omit(sub$symbol)))
      if (length(all_syms) == 0)
        return(DT::datatable(data.frame(Message = "No DEGs found for selected assays."),
                             rownames = FALSE, options = list(dom = "t")))
      mat <- data.frame(Gene = all_syms, stringsAsFactors = FALSE)
      for (a in assays) {
        a_sub <- sub[sub$assay == a, ]
        is_up <- all_syms %in% a_sub$symbol[a_sub$DEG == "up"]
        is_dn <- all_syms %in% a_sub$symbol[a_sub$DEG == "down"]
        lbl   <- if (nchar(a) > 40) paste0("\u2026", substr(a, nchar(a)-39, nchar(a))) else a
        lfc_v <- if ("log2fc" %in% names(a_sub))
          setNames(a_sub$log2fc, a_sub$symbol) else setNames(numeric(0), character(0))
        mat[[lbl]] <- mapply(function(g, up, dn) {
          lfc <- lfc_v[g]; .fmt_dir_cell(up, dn, if (length(lfc) == 0) NA_real_ else lfc)
        }, all_syms, is_up, is_dn, SIMPLIFY = TRUE, USE.NAMES = FALSE)
      }
      ko_cols <- names(mat)[-1]
      n_assays <- rowSums(mat[, ko_cols, drop = FALSE] != "")
      mat <- mat[order(-n_assays), ]
      DT::datatable(mat, rownames = FALSE, filter = "top",
                    options = list(pageLength = 25, scrollX = TRUE,
                                   scrollY = "460px", scrollCollapse = TRUE,
                                   dom = "rtip",
                                   fixedColumns = list(leftColumns = 1),
                                   columnDefs = list(
                                     list(className = "dt-center",
                                          targets = seq_along(ko_cols)),
                                     list(width = "90px",
                                          targets = seq_along(ko_cols))),
                                   rowCallback = .LFC_ROW_CB,
                                   initComplete = .dt_header_js),
                    class = "compact row-border hover")
    })

    # ── UpSet plot (placeholder / empty-state switcher) ─────────────────────
    output$cmp_upset_ui <- renderUI({
      if (is.null(input$cmp_assays) || length(input$cmp_assays) < 2) {
        .empty_state("grid-3x3-gap-fill",
                     "Select 2 or more assays to view overlapping DEGs",
                     "Use the filters and assay picker above to get started.")
      } else {
        plotOutput(session$ns("cmp_upset"), height = "460px")
      }
    })

    # ── UpSet plot ────────────────────────────────────────────────────────────
    output$cmp_upset <- renderPlot({
      sub    <- cmp_sub_data(); req(sub)
      assays <- input$cmp_assays; req(assays, length(assays) >= 2)
      dir    <- input$cmp_dir
      if (nrow(sub) == 0 || length(unique(na.omit(sub$symbol))) == 0) {
        dir_lbl <- switch(dir, up = "up-regulated", down = "down-regulated", all = "")
        plot.new()
        text(0.5, 0.5, paste0("No ", if (nzchar(dir_lbl)) paste0(dir_lbl, " ") else "",
                               "DEGs found for selected assays."),
             cex = 1.2, col = "#6c757d")
        return(invisible())
      }
      all_syms <- sort(unique(na.omit(sub$symbol)))
      mat <- as.data.frame(lapply(setNames(assays, assays), function(a) {
        as.integer(all_syms %in% sub$symbol[sub$assay == a])
      }))
      rownames(mat) <- all_syms
      safe <- gsub("[^A-Za-z0-9_]", "_", names(mat))
      safe <- sub("^([0-9])", "X\\1", safe)
      safe <- make.unique(safe, sep = "_")
      names(mat) <- safe
      clr <- switch(dir, up = .HP_CLRS$up, down = .HP_CLRS$down, all = .HP_CLRS$all)
      UpSetR::upset(mat, sets = rev(names(mat)), order.by = "freq",
                    main.bar.color = clr, sets.bar.color = clr,
                    text.scale = c(2.5, 1.8, 1.5, 2.0, 1.8, 2.0),
                    keep.order = TRUE, mb.ratio = c(0.60, 0.40))
    })

    # ── Download: Table CSV ─────────────────────────────────────────────────────
    output$cmp_dl_csv <- downloadHandler(
      filename = function() {
        paste0("gene_matrix_", input$cmp_dir, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        sub <- cmp_sub_data(); assays <- input$cmp_assays
        if (is.null(sub) || nrow(sub) == 0 ||
            is.null(assays) || length(assays) < 2) {
          write.csv(data.frame(Message = "No data available"),
                    file, row.names = FALSE)
          return()
        }
        all_syms <- sort(unique(na.omit(sub$symbol)))
        mat <- data.frame(Gene = all_syms, stringsAsFactors = FALSE)
        for (a in assays) {
          a_sub <- sub[sub$assay == a, ]
          is_up <- all_syms %in% a_sub$symbol[a_sub$DEG == "up"]
          is_dn <- all_syms %in% a_sub$symbol[a_sub$DEG == "down"]
          lfc_v <- if ("log2fc" %in% names(a_sub))
            setNames(a_sub$log2fc, a_sub$symbol) else setNames(numeric(0), character(0))
          mat[[a]] <- mapply(function(g, up, dn) {
            lfc <- lfc_v[g]
            .fmt_dir_cell(up, dn, if (length(lfc) == 0) NA_real_ else lfc)
          }, all_syms, is_up, is_dn, SIMPLIFY = TRUE, USE.NAMES = FALSE)
        }
        n_assays <- rowSums(mat[, -1, drop = FALSE] != "")
        mat <- mat[order(-n_assays), ]
        write.csv(mat, file, row.names = FALSE)
      }
    )

    # ── Download: Plot PNG ───────────────────────────────────────────────────────
    output$cmp_dl_png <- downloadHandler(
      filename = function() {
        paste0("upset_plot_", input$cmp_dir, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        sub <- cmp_sub_data(); assays <- input$cmp_assays
        dir <- input$cmp_dir
        if (is.null(sub) || nrow(sub) == 0 ||
            is.null(assays) || length(assays) < 2) {
          png(file, width = 800, height = 500)
          plot.new(); text(0.5, 0.5, "No data available",
                           cex = 1.2, col = "#6c757d")
          dev.off(); return()
        }
        all_syms <- sort(unique(na.omit(sub$symbol)))
        if (length(all_syms) == 0) {
          png(file, width = 800, height = 500)
          plot.new(); text(0.5, 0.5, "No DEGs found",
                           cex = 1.2, col = "#6c757d")
          dev.off(); return()
        }
        mat <- as.data.frame(lapply(setNames(assays, assays), function(a) {
          as.integer(all_syms %in% sub$symbol[sub$assay == a])
        }))
        rownames(mat) <- all_syms
        safe <- gsub("[^A-Za-z0-9_]", "_", names(mat))
        safe <- sub("^([0-9])", "X\\1", safe)
        safe <- make.unique(safe, sep = "_")
        names(mat) <- safe
        clr <- switch(dir,
          up = .HP_CLRS$up, down = .HP_CLRS$down, all = .HP_CLRS$all)
        png(file, width = 1200, height = 700, res = 120)
        print(UpSetR::upset(mat, sets = rev(names(mat)), order.by = "freq",
                            main.bar.color = clr, sets.bar.color = clr,
                            text.scale = c(2.5, 1.8, 1.5, 2.0, 1.8, 2.0),
                            keep.order = TRUE, mb.ratio = c(0.60, 0.40)))
        dev.off()
      }
    )

    # ── Suspend outputs until tab is active ───────────────────────────────────
    invisible(lapply(
      c("cmp_no_assays", "cmp_match_count", "cmp_filter_badge",
        "cmp_meta_above", "cmp_matrix_ui", "cmp_matrix_tbl",
        "cmp_upset_ui", "cmp_upset"),
      function(oid) outputOptions(output, oid, suspendWhenHidden = TRUE)
    ))

  })
}
