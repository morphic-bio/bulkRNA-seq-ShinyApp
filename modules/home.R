# =============================================================================
# Module: home
# Landing page — 3 full-width cards, each with description + live preview.
# =============================================================================

library(DT)

homeUI <- function(id) {
  ns <- NS(id)

  div(
    class = "py-4",
    style = "max-width: 1140px; margin: 0 auto;",

    # ── Title ───────────────────────────────────────────────────────────────
    div(
      class = "text-center mb-4",
      h1("MorPhiC Bulk RNA-seq Explorer", class = "fw-bold mb-1"),
      p(class = "text-muted fs-5 mb-0",
        "Multi-study differential expression analysis across KO assays")
    ),

    # ── Card 1: Perturbed Genes Overview ─────────────────────────────────────
    card(
      class = "mb-4 shadow-sm",
      card_body(
        layout_columns(
          col_widths = c(4, 8),

          # Left — description
          div(
            class = "d-flex flex-column justify-content-between h-100 pe-3",
            div(
              div(class = "d-flex align-items-center gap-2 mb-2",
                  bsicons::bs_icon("bar-chart-fill", size = "1.4rem"),
                  tags$h4("Perturbed Genes Overview", class = "mb-0 fw-semibold")),
              tags$p(class = "text-muted mb-3",
                     "Browse all MorPhiC assays across KO genes, model systems, and cell lines. ",
                     "Explore study breakdown statistics, N assays per KO gene, ",
                     "and pathway overlap for KO genes across Reactome, GO, and PANTHER.")
            ),
            actionButton(
              ns("go_page1"),
              label = tagList(bsicons::bs_icon("arrow-right-circle"), " Open"),
              class = "btn btn-primary mt-2",
              style = "width: fit-content;"
            )
          ),

          # Right — actual assay overview table
          div(
            class = "border rounded p-2",
            style = "background:#fafafa; max-height:260px; overflow:hidden;",
            tags$p(class = "text-muted small mb-1 fw-semibold", "Assay overview"),
            DTOutput(ns("preview_p1"), height = "210px")
          )
        )
      )
    ),

    # ── Card 2: Assay Comparison ──────────────────────────────────────────────
    card(
      class = "mb-4 shadow-sm",
      card_body(
        layout_columns(
          col_widths = c(4, 8),

          # Left — description
          div(
            class = "d-flex flex-column justify-content-between h-100 pe-3",
            div(
              div(class = "d-flex align-items-center gap-2 mb-2",
                  bsicons::bs_icon("grid-3x3-gap-fill", size = "1.4rem"),
                  tags$h4("Assay Comparison", class = "mb-0 fw-semibold")),
              tags$p(class = "text-muted mb-3",
                     "Select assays to compare differentially expressed genes side-by-side. ",
                     "View a gene-by-assay presence matrix, UpSet plot of overlapping DEG sets, ",
                     "and a table of genes shared across multiple assays.")
            ),
            actionButton(
              ns("go_page2"),
              label = tagList(bsicons::bs_icon("arrow-right-circle"), " Open"),
              class = "btn btn-primary mt-2",
              style = "width: fit-content;"
            )
          ),

          # Right — gene × assay matrix (FEZF2 CE vs KO)
          div(
            class = "border rounded p-2",
            style = "background:#fafafa; max-height:260px; overflow:hidden;",
            tags$p(class = "text-muted small mb-1 fw-semibold",
                   "Example: FEZF2 CE vs KO — shared DEGs"),
            DTOutput(ns("preview_p2"), height = "210px")
          )
        )
      )
    ),

    # ── Card 3: Assay Analysis ────────────────────────────────────────────────
    card(
      class = "mb-4 shadow-sm",
      card_body(
        layout_columns(
          col_widths = c(4, 8),

          # Left — description
          div(
            class = "d-flex flex-column justify-content-between h-100 pe-3",
            div(
              div(class = "d-flex align-items-center gap-2 mb-2",
                  bsicons::bs_icon("diagram-3", size = "1.4rem"),
                  tags$h4("Assay Analysis", class = "mb-0 fw-semibold")),
              tags$p(class = "text-muted mb-3",
                     "For any assay, explore DEG overlap with Reactome pathways, GO terms, ",
                     "and PANTHER protein classes. Cross-reference with IMPC mouse phenotypes, ",
                     "HPO human phenotypes, OMIM/Orphanet disease associations, and PanelApp gene panels.")
            ),
            actionButton(
              ns("go_page3"),
              label = tagList(bsicons::bs_icon("arrow-right-circle"), " Open"),
              class = "btn btn-primary mt-2",
              style = "width: fit-content;"
            )
          ),

          # Right — Reactome pathway overlap table (ISL1 KO example)
          div(
            class = "border rounded p-2",
            style = "background:#fafafa; max-height:260px; overflow:hidden;",
            tags$p(class = "text-muted small mb-1 fw-semibold",
                   "Example: ISL1 KO \u2014 Reactome pathway overlap"),
            DTOutput(ns("preview_p3"), height = "210px")
          )
        )
      )
    )
  )
}

homeServer <- function(id, parent_session, con_r, study_info_r) {
  moduleServer(id, function(input, output, session) {

    # ── Navigation ────────────────────────────────────────────────────────────
    observeEvent(input$go_page1, {
      nav_select("main_nav", selected = "page1", session = parent_session)
    })
    observeEvent(input$go_page2, {
      nav_select("main_nav", selected = "page2", session = parent_session)
    })
    observeEvent(input$go_page3, {
      nav_select("main_nav", selected = "page3", session = parent_session)
    })

    # ── Preview 1: full assay overview table ──────────────────────────────────
    output$preview_p1 <- DT::renderDT({
      si   <- study_info_r()
      want <- c("Gene", "Assay", "Model_system", "KO_strat", "DPC")
      cols <- want[want %in% colnames(si)]
      tbl  <- si[!duplicated(si$Assay), cols, drop = FALSE]
      tbl  <- tbl[order(tbl$Gene), ]
      pretty <- c(Gene = "KO Gene", Assay = "Assay", Model_system = "Model System",
                  KO_strat = "KO Strategy", DPC = "DPC")
      colnames(tbl) <- pretty[colnames(tbl)]
      DT::datatable(
        tbl,
        rownames  = FALSE,
        selection = "none",
        options   = list(
          dom            = "t",
          pageLength     = nrow(tbl),
          scrollY        = "185px",
          scrollCollapse = TRUE,
          ordering       = FALSE,
          columnDefs     = list(list(className = "dt-left", targets = "_all")),
          initComplete   = .dt_header_js
        ),
        class = "compact row-border hover"
      )
    })

    # ── Preview 2: gene × assay matrix — FEZF2 CE vs KO ──────────────────────
    output$preview_p2 <- DT::renderDT({
      con    <- con_r()
      assays <- c(
        "T2024_06_JAX_RNAseq_Neuroectoderm_CBO_day_6_FEZF2_CE",
        "T2024_06_JAX_RNAseq_Neuroectoderm_CBO_day_6_FEZF2_KO"
      )
      raw <- tryCatch(load_deg_genes(con, assays), error = function(e) NULL)
      req(!is.null(raw), nrow(raw) > 0)

      ce_rows <- raw[raw$assay == assays[1], ]
      ko_rows <- raw[raw$assay == assays[2], ]
      all_genes <- unique(raw$symbol)

      mat <- data.frame(
        Gene = all_genes,
        CE   = vapply(all_genes, function(g) {
          r <- ce_rows[ce_rows$symbol == g, ]
          if (nrow(r) == 0) return("")
          .fmt_dir_cell(any(r$DEG == "up"), any(r$DEG == "down"))
        }, character(1)),
        KO   = vapply(all_genes, function(g) {
          r <- ko_rows[ko_rows$symbol == g, ]
          if (nrow(r) == 0) return("")
          .fmt_dir_cell(any(r$DEG == "up"), any(r$DEG == "down"))
        }, character(1)),
        stringsAsFactors = FALSE
      )

      # genes shared by both assays first
      in_both <- mat$CE != "" & mat$KO != ""
      mat     <- rbind(mat[in_both, ], mat[!in_both, ])
      colnames(mat) <- c("Gene", "FEZF2 CE", "FEZF2 KO")

      DT::datatable(
        mat,
        rownames  = FALSE,
        selection = "none",
        options   = list(
          dom            = "t",
          pageLength     = nrow(mat),
          scrollY        = "185px",
          scrollCollapse = TRUE,
          ordering       = FALSE,
          columnDefs     = list(
            list(className = "dt-left",   targets = 0),
            list(className = "dt-center", targets = c(1, 2))
          ),
          rowCallback    = .LFC_ROW_CB,
          initComplete   = .dt_header_js
        ),
        class = "compact row-border hover"
      )
    })

    # ── Preview 3: Reactome pathway overlap — ISL1 KO ─────────────────────────
    output$preview_p3 <- DT::renderDT({
      con   <- con_r()
      assay <- "T2023_12_JAX_RNAseq_ExtraEmbryonic_ExM_day_6_nor_ISL1_KO"
      df    <- compute_deg_overlap(con, assay, "ref_reactome", "updn", 15L, 1L)
      req(!is.null(df), nrow(df) > 0)

      # truncate long pathway names
      df$category <- ifelse(
        nchar(df$category) > 48,
        paste0(substr(df$category, 1, 45), "\u2026"),
        df$category
      )
      df$pct <- round(df$pct, 1)
      tbl    <- df[, c("category", "n_pathway", "n_overlap", "pct")]
      names(tbl) <- c("Pathway", "Path. Size", "N DEGs", "% Overlap")

      DT::datatable(
        tbl,
        rownames  = FALSE,
        selection = "none",
        options   = list(
          dom            = "t",
          pageLength     = 15,
          scrollY        = "185px",
          scrollCollapse = TRUE,
          ordering       = FALSE,
          columnDefs     = list(
            list(className = "dt-left",   targets = 0),
            list(className = "dt-center", targets = c(1, 2, 3))
          ),
          initComplete   = .dt_header_js
        ),
        class = "compact row-border hover"
      ) |>
        DT::formatStyle(
          "% Overlap",
          background         = DT::styleColorBar(c(0, 100), "#4E79A7"),
          backgroundSize     = "98% 70%",
          backgroundRepeat   = "no-repeat",
          backgroundPosition = "center"
        )
    })

    invisible(lapply(
      c("preview_p1", "preview_p2", "preview_p3"),
      function(oid) outputOptions(output, oid, suspendWhenHidden = TRUE)
    ))
  })
}
