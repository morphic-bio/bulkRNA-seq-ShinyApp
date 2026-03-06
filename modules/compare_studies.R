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

# ── UI ────────────────────────────────────────────────────────────────────────

compare_studiesUI <- function(id) {
  ns <- NS(id)

  tagList(
    card(
      card_header(tagList(bsicons::bs_icon("grid-3x3-gap-fill"), " Assay Comparison",
                          .info_tip("Compare differentially expressed genes across multiple assays side-by-side."))),
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
            nav_panel(tagList("Gene Matrix",
                            .info_tip("Matrix showing DEG direction per gene across selected assays. Sorted by number of assays.")),
                      div(class = "pt-2", uiOutput(ns("cmp_matrix_ui")))),
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
                          plotOutput(ns("cmp_upset"), height = "460px"))),
            nav_panel(tagList("Common DEGs",
                            .info_tip("Genes differentially expressed in the same direction across all selected assays.")),
                      div(class = "pt-2", uiOutput(ns("cmp_consensus"))))
          )
        ))
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

    # ── Metadata header ───────────────────────────────────────────────────────
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
        div(class = "border rounded p-2 mb-1",
            style = "background:#f8f9fa; flex:1; min-width:180px;",
            tags$div(class = "small fw-semibold text-break mb-1",
                     style = "font-size:0.68rem; word-break:break-all; color:#495057;", a),
            tags$div(bdg(r$Gene, "#198754"), bdg(r$Model_system, "#0d6efd"),
                     bdg(r$KO_strat, "#0dcaf0"), bdg(r$DPC, "#e67e22"),
                     bdg(r$Cell_Line, "#6c757d"), extra))
      })
      tagList(div(class = "d-flex flex-wrap gap-2 mb-2", cards), tags$hr(class = "my-2"))
    })

    # ── Gene matrix (placeholder / DT switcher) ─────────────────────────────
    output$cmp_matrix_ui <- renderUI({
      if (is.null(input$cmp_assays) || length(input$cmp_assays) < 2) {
        .empty_state("grid-3x3-gap-fill",
                     "Select 2 or more assays to view DEG comparison",
                     "Use the filters and assay picker on the left to get started.")
      } else {
        DTOutput(session$ns("cmp_matrix_tbl"))
      }
    })

    # ── Gene matrix ───────────────────────────────────────────────────────────
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
                                   rowCallback = .LFC_ROW_CB,
                                   initComplete = .dt_header_js),
                    class = "compact row-border hover") |>
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

    # ── Consensus DEGs ────────────────────────────────────────────────────────
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

      chip <- function(g, bg)
        tags$span(style = paste0("background-color:", bg, "; color:white; font-size:0.78rem;",
                                 " font-weight:600; padding:0.28em 0.6em; border-radius:4px;",
                                 " display:inline-block; margin:2px 3px 2px 0;"), g)

      make_section <- function(genes, label, bg, border_col) {
        n_g <- length(genes)
        hdr  <- tags$div(class = "fw-bold mb-1", style = paste0("color:", border_col, ";"),
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

    # ── Suspend outputs until tab is active ───────────────────────────────────
    invisible(lapply(
      c("cmp_meta_above", "cmp_matrix_ui", "cmp_matrix_tbl", "cmp_upset", "cmp_consensus"),
      function(oid) outputOptions(output, oid, suspendWhenHidden = TRUE)
    ))

  })
}
