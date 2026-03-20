# =============================================================================
# Module: home
# Landing page — 3 clickable cards side-by-side, each with a live preview table.
# =============================================================================

library(DT)

homeUI <- function(id) {
  ns <- NS(id)

  # Shared CSS for clickable cards
  card_css <- paste0(
    "cursor:pointer; transition: box-shadow 0.15s, transform 0.15s;",
    "height:100%;"
  )
  card_hover_css <- HTML(
    "<style>
      .home-card:hover { box-shadow: 0 .5rem 1rem rgba(0,0,0,.12) !important;
                         transform: translateY(-2px); }
      .home-card .card-arrow { transition: transform 0.15s; }
      .home-card:hover .card-arrow { transform: translateX(4px); }
    </style>"
  )

  div(
    class = "py-4",
    style = "max-width: 1400px; margin: 0 auto;",
    card_hover_css,

    # ── Title ───────────────────────────────────────────────────────────────
    div(
      class = "text-center mb-4",
      h1("MorPhiC Bulk RNA-seq Explorer", class = "fw-bold mb-1"),
      uiOutput(ns("home_chips"), inline = TRUE)
    ),

    # ── Three cards side-by-side ──────────────────────────────────────────
    layout_columns(
      col_widths = c(4, 4, 4),
      fill       = TRUE,

      # ── Card 1: Assays Overview ──────────────────────────────────────
      div(
        id    = ns("card_page1"),
        class = "home-card",
        style = card_css,
        card(
          class = "shadow-sm h-100",
          card_header(
            class = "d-flex align-items-center gap-2 pb-1",
            bsicons::bs_icon("bar-chart-fill", size = "1.2rem"),
            tags$span("Assays Overview", class = "fw-semibold fs-6"),
            div(style = "flex:1;"),
            bsicons::bs_icon("arrow-right", size = "1rem", class = "card-arrow text-muted")
          ),
          card_body(
            class = "pt-1 pb-2 px-3",
            tags$p(class = "text-muted small mb-2",
                   "Browse assays across KO genes, model systems, and cell lines."),
            div(
              style = "max-height:280px; overflow:hidden;",
              DTOutput(ns("preview_p1"), height = "260px")
            )
          )
        )
      ),

      # ── Card 2: Assay Comparison ──────────────────────────────────────
      div(
        id    = ns("card_page2"),
        class = "home-card",
        style = card_css,
        card(
          class = "shadow-sm h-100",
          card_header(
            class = "d-flex align-items-center gap-2 pb-1",
            bsicons::bs_icon("grid-3x3-gap-fill", size = "1.2rem"),
            tags$span("Assay Comparison", class = "fw-semibold fs-6"),
            div(style = "flex:1;"),
            bsicons::bs_icon("arrow-right", size = "1rem", class = "card-arrow text-muted")
          ),
          card_body(
            class = "pt-1 pb-2 px-3",
            tags$p(class = "text-muted small mb-2",
                   "Compare DEGs across selected assays."),
            div(
              style = "max-height:280px; overflow:hidden;",
              DTOutput(ns("preview_p2"), height = "260px")
            )
          )
        )
      ),

      # ── Card 3: Assay Analysis ────────────────────────────────────────
      div(
        id    = ns("card_page3"),
        class = "home-card",
        style = card_css,
        card(
          class = "shadow-sm h-100",
          card_header(
            class = "d-flex align-items-center gap-2 pb-1",
            bsicons::bs_icon("diagram-3", size = "1.2rem"),
            tags$span("Assay Analysis", class = "fw-semibold fs-6"),
            div(style = "flex:1;"),
            bsicons::bs_icon("arrow-right", size = "1rem", class = "card-arrow text-muted")
          ),
          card_body(
            class = "pt-1 pb-2 px-3",
            tags$p(class = "text-muted small mb-2",
                   "Explore DEG functional/disease/phenotype annotations."),
            div(
              style = "max-height:280px; overflow:hidden;",
              DTOutput(ns("preview_p3"), height = "260px")
            )
          )
        )
      )
    ) # end layout_columns
  )
}

homeServer <- function(id, parent_session, con_r, study_info_r) {
  moduleServer(id, function(input, output, session) {

    # ── Make entire cards clickable via JS ──────────────────────────────────
    observe({
      ns <- session$ns
      shinyjs_click <- function(card_id, page_val) {
        shiny::insertUI(
          selector  = "body",
          where     = "beforeEnd",
          immediate = TRUE,
          ui = tags$script(HTML(sprintf(
            "$(document).on('click', '#%s', function() {
               Shiny.setInputValue('%s', Math.random());
             });",
            card_id, ns(page_val)
          )))
        )
      }
      shinyjs_click(ns("card_page1"), "go_page1")
      shinyjs_click(ns("card_page2"), "go_page2")
      shinyjs_click(ns("card_page3"), "go_page3")
    })

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

    # ── Summary chips ──────────────────────────────────────────────────────────
    output$home_chips <- renderUI({
      si <- study_info_r()
      n_studies <- length(unique(na.omit(si$Study_title)))
      n_assays  <- length(unique(si$Assay))
      n_genes   <- length(unique(na.omit(si$Gene)))
      n_models  <- length(unique(na.omit(si$Model_system)))
      n_comps   <- length(unique(na.omit(si$Comparison)))
      chip <- function(n, label, bg) {
        tags$span(class = "badge rounded-pill",
                  style = sprintf("background:%s; font-size:0.85rem; font-weight:500; padding:6px 14px;", bg),
                  paste(n, label))
      }
      div(
        class = "d-flex justify-content-center gap-2 mt-2 flex-wrap",
        chip(n_genes,   "genes",       "#E15759"),
        chip(n_studies, "studies",     "#4E79A7"),
        chip(n_assays,  "assays",      "#59A14F"),
        chip(n_models,  "model systems", "#F28E2B"),
        chip(n_comps,   "comparisons", "#76B7B2")
      )
    })

    # ── Preview 1: assay overview (matches current Studies Overview columns) ──
    output$preview_p1 <- DT::renderDT({
      si   <- study_info_r()
      want <- c("Gene", "Assay", "Comparison", "Model_system", "DPC", "Cell_Line")
      cols <- want[want %in% colnames(si)]
      tbl  <- si[!duplicated(si$Assay), cols, drop = FALSE]
      tbl  <- tbl[order(tbl$Gene), ]

      # Count DEGs per assay
      con <- con_r()
      deg_counts <- do.call(rbind, lapply(tbl$Assay, function(a) {
        tryCatch({
          df <- dbGetQuery(con, sprintf(
            "SELECT DEG, COUNT(*) AS n FROM main.\"%s\" WHERE DEG IN ('up','down') GROUP BY DEG", a))
          data.frame(Assay = a,
                     n_up   = sum(df$n[df$DEG == "up"],   na.rm = TRUE),
                     n_down = sum(df$n[df$DEG == "down"], na.rm = TRUE),
                     stringsAsFactors = FALSE)
        }, error = function(e) {
          data.frame(Assay = a, n_up = 0L, n_down = 0L, stringsAsFactors = FALSE)
        })
      }))
      tbl <- merge(tbl, deg_counts, by = "Assay", all.x = TRUE)
      tbl$n_up[is.na(tbl$n_up)]     <- 0L
      tbl$n_down[is.na(tbl$n_down)] <- 0L

      # Build HTML bar for N DEGs
      max_deg <- max(tbl$n_up + tbl$n_down, 1L, na.rm = TRUE)
      tbl$N_DEGs <- mapply(function(up, dn) {
        total <- up + dn
        if (total == 0) return('<span style="color:#adb5bd; font-size:0.82em;">0</span>')
        w_up <- round(up / max_deg * 80)
        w_dn <- round(dn / max_deg * 80)
        sprintf(
          paste0(
            '<div style="display:flex; align-items:center; gap:3px;">',
              '<div style="display:flex; height:12px; border-radius:2px; overflow:hidden;">',
                '<div style="width:%dpx; background:#E63946;"></div>',
                '<div style="width:%dpx; background:#457B9D;"></div>',
              '</div>',
              '<span style="font-size:0.7em; color:#495057;">',
                '%s\u2191 %s\u2193',
              '</span>',
            '</div>'),
          w_up, w_dn,
          formatC(up, big.mark = ","), formatC(dn, big.mark = ","))
      }, tbl$n_up, tbl$n_down)

      # Select display columns — keep compact to fit card width without scrolling
      tbl <- tbl[, c("Gene", "Comparison", "Model_system", "N_DEGs"),
                 drop = FALSE]

      pretty <- c(Gene = "Gene", Comparison = "Comparison",
                  Model_system = "Model System", N_DEGs = "N DEGs")
      colnames(tbl) <- pretty[colnames(tbl)]

      DT::datatable(
        tbl,
        rownames  = FALSE,
        selection = "none",
        escape    = FALSE,
        width     = "100%",
        options   = list(
          dom            = "t",
          pageLength     = nrow(tbl),
          scrollY        = "230px",
          scrollCollapse = TRUE,
          ordering       = FALSE,
          autoWidth      = TRUE,
          columnDefs     = list(
            list(className = "dt-left", targets = "_all")
          ),
          initComplete = .dt_header_js
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
          scrollY        = "230px",
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

    # ── Preview 3: Reactome pathway overlap — ISL1 KO (N DEGs view) ───────────
    output$preview_p3 <- DT::renderDT({
      con   <- con_r()
      assay <- "T2023_12_JAX_RNAseq_ExtraEmbryonic_ExM_day_6_nor_ISL1_KO"
      df    <- compute_deg_overlap(con, assay, "ref_reactome", "all", 15L, 500L,
                                   order_col = "n_overlap")
      req(!is.null(df), nrow(df) > 0)

      # truncate long pathway names
      df$category <- ifelse(
        nchar(df$category) > 40,
        paste0(substr(df$category, 1, 37), "\u2026"),
        df$category
      )
      tbl    <- df[, c("category", "n_pathway", "n_overlap")]
      names(tbl) <- c("Pathway", "Size", "N DEGs")

      DT::datatable(
        tbl,
        rownames  = FALSE,
        selection = "none",
        options   = list(
          dom            = "t",
          pageLength     = 15,
          scrollY        = "230px",
          scrollCollapse = TRUE,
          ordering       = FALSE,
          columnDefs     = list(
            list(className = "dt-left",   targets = 0),
            list(className = "dt-center", targets = c(1, 2))
          ),
          initComplete = .dt_header_js
        ),
        class = "compact row-border hover"
      ) |>
        DT::formatStyle(
          "N DEGs",
          background         = DT::styleColorBar(range(tbl[["N DEGs"]], na.rm = TRUE), "#4E79A7"),
          backgroundSize     = "98% 70%",
          backgroundRepeat   = "no-repeat",
          backgroundPosition = "center"
        )
    })

    invisible(lapply(
      c("home_chips", "preview_p1", "preview_p2", "preview_p3"),
      function(oid) outputOptions(output, oid, suspendWhenHidden = TRUE)
    ))
  })
}
