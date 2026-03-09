# =============================================================================
# Module: compare_studies
# Page 2 — Compare Studies
#
# Sections:
#   5. Assay Comparison — filter + picker, metadata header, tabs:
#                         Gene Matrix | UpSet | Common DEGs
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
                        div(class = "d-flex align-items-center gap-3 mb-2",
                            tags$label("Direction:", class = "mb-0 small fw-semibold"),
                            radioButtons(ns("cmp_matrix_dir"), NULL,
                                         choices  = c("\u2191 Up" = "up",
                                                      "\u2193 Down" = "down",
                                                      "\u2191\u2193 All DEGs" = "all"),
                                         selected = "all", inline = TRUE)),
                        uiOutput(ns("cmp_matrix_ui")))),
          nav_panel(tagList("UpSet Plot",
                          .info_tip("Intersection plot showing overlap and unique DEG sets across selected assays.")),
                    div(class = "pt-2",
                        div(class = "d-flex align-items-center gap-3 mb-2",
                            tags$label("Direction:", class = "mb-0 small fw-semibold"),
                            radioButtons(ns("cmp_upset_dir"), NULL,
                                         choices  = c("\u2191 Up" = "up",
                                                      "\u2193 Down" = "down",
                                                      "\u2191\u2193 All DEGs" = "all"),
                                         selected = "up", inline = TRUE)),
                        uiOutput(ns("cmp_upset_ui")))),
          nav_panel(tagList("Common DEGs",
                          .info_tip("Genes differentially expressed in the same direction across all selected assays.")),
                    div(class = "pt-2", uiOutput(ns("cmp_consensus_ui"))))
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

    # ── Metadata header ───────────────────────────────────────────────────────
    output$cmp_meta_above <- renderUI({
      assays <- input$cmp_assays; req(assays, length(assays) >= 1)
      si <- study_info_r()

      info_row <- function(label, val) {
        if (is.null(val) || is.na(val) || !nzchar(as.character(val))) return(NULL)
        tags$div(class = "small", style = "font-size:0.75rem;",
                 tags$span(class = "text-muted", paste0(label, ": ")),
                 tags$span(class = "fw-semibold", as.character(val)))
      }

      cards <- lapply(assays, function(a) {
        r <- si[si$Assay == a, ][1, ]
        div(class = "border rounded p-2",
            style = "background:#f8f9fa; min-width:180px; flex:1;",
            tags$div(class = "small fw-semibold text-break mb-1",
                     style = "font-size:0.68rem; word-break:break-all; color:#495057;", a),
            info_row("Gene", r$Gene),
            info_row("Model System", r$Model_system),
            info_row("Comparison", r$Comparison),
            info_row("DPC", r$DPC),
            info_row("Cell Line", r$Cell_Line),
            info_row("Condition", r$Condition),
            if ("Differentation_time_point" %in% names(r))
              info_row("Differentiation", r$Differentation_time_point),
            if ("Replicate" %in% names(r))
              info_row("Replicate", r$Replicate))
      })
      tagList(
        tags$span(class = "fw-semibold small d-block mb-1", "Selected Assays"),
        div(class = "d-flex flex-wrap gap-2", cards)
      )
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
      genes_df <- cmp_genes_data(); req(genes_df)
      assays   <- input$cmp_assays; req(assays, length(assays) >= 2)
      dir      <- input$cmp_matrix_dir
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
        lbl   <- if (nchar(a) > 40) paste0("\u2026", substr(a, nchar(a)-39, nchar(a))) else a
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
      DT::datatable(mat, rownames = FALSE, filter = "top",
                    options = list(pageLength = 25, scrollX = TRUE,
                                   scrollY = "460px", scrollCollapse = TRUE,
                                   dom = "frtip",
                                   fixedColumns = list(leftColumns = 2),
                                   columnDefs = list(
                                     list(className = "dt-center",
                                          targets = seq_along(ko_cols) + 1),
                                     list(width = "90px",
                                          targets = seq_along(ko_cols) + 1)),
                                   rowCallback = .LFC_ROW_CB,
                                   initComplete = .dt_header_js),
                    class = "compact row-border hover") |>
        DT::formatStyle("N_Assays",
                        background = DT::styleColorBar(c(0, length(assays)), "#d4e6f1"),
                        backgroundSize = "100% 80%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "center")
    })

    # ── UpSet plot (placeholder / empty-state switcher) ─────────────────────
    output$cmp_upset_ui <- renderUI({
      if (is.null(input$cmp_assays) || length(input$cmp_assays) < 2) {
        .empty_state("grid-3x3-gap-fill",
                     "Select 2 or more assays to view UpSet plot",
                     "Use the filters and assay picker above to get started.")
      } else {
        plotOutput(session$ns("cmp_upset"), height = "460px")
      }
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
                    text.scale = c(1.2, 1, 1, 1, 1, 0.9),
                    keep.order = TRUE, mb.ratio = c(0.60, 0.40))
    })

    # ── Consensus DEGs data ──────────────────────────────────────────────────
    cmp_consensus_data <- reactive({
      genes_df <- cmp_genes_data(); req(genes_df)
      assays   <- input$cmp_assays; req(assays, length(assays) >= 2)

      up_per <- lapply(assays, function(a)
        unique(na.omit(genes_df$symbol[genes_df$assay == a & genes_df$DEG == "up"])))
      dn_per <- lapply(assays, function(a)
        unique(na.omit(genes_df$symbol[genes_df$assay == a & genes_df$DEG == "down"])))
      up_all <- sort(Reduce(intersect, up_per))
      dn_all <- sort(Reduce(intersect, dn_per))

      build_tbl <- function(genes, direction) {
        if (length(genes) == 0) return(data.frame(Gene = character(0)))
        sub <- genes_df[genes_df$symbol %in% genes & genes_df$DEG == direction &
                        genes_df$assay %in% assays, ]
        df <- data.frame(Gene = genes, stringsAsFactors = FALSE)
        for (a in assays) {
          lbl <- if (nchar(a) > 35) paste0("\u2026", substr(a, nchar(a)-34, nchar(a))) else a
          a_sub <- sub[sub$assay == a, ]
          if ("log2fc" %in% names(a_sub)) {
            lfc_map <- setNames(a_sub$log2fc, a_sub$symbol)
            df[[lbl]] <- round(lfc_map[genes], 2)
          }
        }
        df
      }

      list(up_tbl = build_tbl(up_all, "up"),
           dn_tbl = build_tbl(dn_all, "down"),
           n_assays = length(assays))
    })

    # ── Consensus DEGs (placeholder / empty-state switcher) ─────────────────
    output$cmp_consensus_ui <- renderUI({
      if (is.null(input$cmp_assays) || length(input$cmp_assays) < 2) {
        .empty_state("grid-3x3-gap-fill",
                     "Select 2 or more assays to view common DEGs",
                     "Use the filters and assay picker above to get started.")
      } else {
        uiOutput(session$ns("cmp_consensus"))
      }
    })

    # ── Consensus DEGs layout ────────────────────────────────────────────────
    output$cmp_consensus <- renderUI({
      cd <- cmp_consensus_data(); req(cd)
      n    <- cd$n_assays
      n_up <- nrow(cd$up_tbl)
      n_dn <- nrow(cd$dn_tbl)

      up_count <- tags$span(class = "fw-normal text-muted ms-1",
                            paste0("(", n_up, " gene", if (n_up != 1) "s", ")"))
      dn_count <- tags$span(class = "fw-normal text-muted ms-1",
                            paste0("(", n_dn, " gene", if (n_dn != 1) "s", ")"))

      tagList(
        tags$p(class = "text-muted small mb-3",
               paste0("Genes regulated in the same direction across all ",
                      n, " selected assays.")),
        layout_columns(
          col_widths = c(6, 6),
          div(
            tags$div(class = "fw-bold mb-2", style = "color: #C0392B;",
                     "\u2191 Up-regulated in all assays", up_count),
            if (n_up == 0)
              tags$p(class = "text-muted small", "No genes shared across all assays.")
            else
              DTOutput(session$ns("cmp_consensus_up_tbl"))
          ),
          div(
            tags$div(class = "fw-bold mb-2", style = "color: #2C7BB6;",
                     "\u2193 Down-regulated in all assays", dn_count),
            if (n_dn == 0)
              tags$p(class = "text-muted small", "No genes shared across all assays.")
            else
              DTOutput(session$ns("cmp_consensus_dn_tbl"))
          )
        )
      )
    })

    # ── Consensus Up table ───────────────────────────────────────────────────
    output$cmp_consensus_up_tbl <- DT::renderDT({
      cd <- cmp_consensus_data(); req(cd, nrow(cd$up_tbl) > 0)
      DT::datatable(cd$up_tbl, rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE,
                                   dom = "ftip", initComplete = .dt_header_js),
                    class = "compact row-border hover")
    })

    # ── Consensus Down table ─────────────────────────────────────────────────
    output$cmp_consensus_dn_tbl <- DT::renderDT({
      cd <- cmp_consensus_data(); req(cd, nrow(cd$dn_tbl) > 0)
      DT::datatable(cd$dn_tbl, rownames = FALSE,
                    options = list(pageLength = 15, scrollX = TRUE,
                                   dom = "ftip", initComplete = .dt_header_js),
                    class = "compact row-border hover")
    })

    # ── Suspend outputs until tab is active ───────────────────────────────────
    invisible(lapply(
      c("cmp_no_assays", "cmp_match_count", "cmp_filter_badge", "cmp_meta_above", "cmp_matrix_ui", "cmp_matrix_tbl", "cmp_upset_ui", "cmp_upset", "cmp_consensus_ui", "cmp_consensus", "cmp_consensus_up_tbl", "cmp_consensus_dn_tbl"),
      function(oid) outputOptions(output, oid, suspendWhenHidden = TRUE)
    ))

  })
}
