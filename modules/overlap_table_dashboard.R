library(shiny)
library(bslib)
library(bsicons)
library(plotly)
library(DT)
library(dplyr)
library(DBI)

# ── Source → column name maps (morphic_overlap tables only) ───────────────────

.OT_SOURCES <- c(
  "Reactome Pathways"      = "reactome",
  "GO Biological Process"  = "go_bp",
  "GO Molecular Function"  = "go_mf",
  "PANTHER Protein Class"  = "panther_class",
  "PANTHER Protein Family" = "panther_family",
  "PanelApp Panels"        = "panelapp"
)

# Name of the category label column per source (morphic_overlap)
.OT_GRP <- c(
  reactome       = "path_name",
  go_bp          = "go_term_BP",
  go_mf          = "go_term_MF",
  panther_class  = "CLASS_TERM",
  panther_family = "FAMILY_TERM",
  panelapp       = "panel_name"
)

# Name of the "N genes in pathway" column per source (morphic_overlap)
.OT_N <- c(
  reactome       = "pathway_genes_n",
  go_bp          = "biological_process_genes_n",
  go_mf          = "molecular_function_genes_n",
  panther_class  = "class_genes_n",
  panther_family = "family_genes_n",
  panelapp       = "panel_genes_n"
)

# deg_overlap tables store source values with different casing, and PANTHER
# further splits by category_type.  These maps translate the UI key → DB value.
.OT_DEG_SOURCE <- c(
  reactome       = "Reactome",
  go_bp          = "GO_BP",
  go_mf          = "GO_MF",
  panther_class  = "PANTHER",
  panther_family = "PANTHER",
  panelapp       = "PanelApp"
)

# category_type filter needed for PANTHER (to split class vs family); NA = no filter
.OT_DEG_CATTYPE <- c(
  reactome       = NA_character_,
  go_bp          = NA_character_,
  go_mf          = NA_character_,
  panther_class  = "protein_class",
  panther_family = "protein_family",
  panelapp       = NA_character_
)

# Human-readable column labels per source (used in DT headers & scatter axis)
.OT_LABEL <- list(
  reactome       = list(cat = "Pathway",        size = "Pathway genes"),
  go_bp          = list(cat = "GO Term (BP)",   size = "GO term genes"),
  go_mf          = list(cat = "GO Term (MF)",   size = "GO term genes"),
  panther_class  = list(cat = "Protein Class",  size = "Class genes"),
  panther_family = list(cat = "Protein Family", size = "Family genes"),
  panelapp       = list(cat = "Panel",          size = "Panel genes")
)

# ── UI ────────────────────────────────────────────────────────────────────────

overlap_table_dashboardUI <- function(id) {
  ns <- NS(id)

  tagList(

    # Lift all selectize dropdowns above card boundaries
    tags$style(HTML(".selectize-dropdown { z-index: 1080 !important; }")),

    # ── Row 1: Mode / Source / Direction / Assay ───────────────────────────
    layout_columns(
      col_widths = c(3, 3, 3, 3),

      card(
        card_header(tagList(bs_icon("layers"), " Mode")),
        radioButtons(
          ns("mode"), NULL,
          choices  = c("KO Gene Overlap" = "ko", "DEG Overlap" = "deg"),
          selected = "ko"
        )
      ),

      card(
        style = "overflow: visible;",
        card_header(tagList(bs_icon("funnel"), " Source")),
        selectInput(ns("source"), NULL,
                    choices  = .OT_SOURCES,
                    selected = "reactome",
                    width    = "100%")
      ),

      # Direction — only meaningful in DEG mode
      card(
        card_header(tagList(bs_icon("arrow-down-up"), " DEG Direction")),
        uiOutput(ns("direction_ui"))
      ),

      # Assay picker — server-rendered; shows a note in KO mode
      card(
        style = "overflow: visible;",
        card_header(tagList(bs_icon("database"), " Assay")),
        uiOutput(ns("assay_ui"))
      )
    ),

    # ── Assay info row — big cards shown when a single assay is selected ────
    uiOutput(ns("assay_info_row")),

    # ── Row 2: Filter controls ─────────────────────────────────────────────
    layout_columns(
      col_widths = c(4, 4, 4),

      card(
        card_header(
          tagList(bs_icon("filter-square"), " Min pathway size"),
          tooltip(bs_icon("info-circle"), "Exclude pathways with fewer than this many genes")
        ),
        sliderInput(ns("min_size"), NULL, min = 1, max = 200, value = 5, step = 1)
      ),

      card(
        card_header(
          tagList(bs_icon("intersect"), " Min overlap genes"),
          tooltip(bs_icon("info-circle"), "Only show pathways where at least this many MorPhiC genes overlap")
        ),
        sliderInput(ns("min_overlap"), NULL, min = 1, max = 50, value = 1, step = 1)
      ),

      card(
        card_header(
          tagList(bs_icon("hash"), " Top N rows"),
          tooltip(bs_icon("info-circle"), "Show the top N categories sorted by overlap %")
        ),
        sliderInput(ns("top_n"), NULL, min = 10, max = 100, value = 30, step = 5)
      )
    ),

    # ── Row 2: KPI value boxes ─────────────────────────────────────────────
    layout_columns(
      col_widths = c(3, 3, 3, 3),

      value_box(
        title    = "Categories shown",
        value    = textOutput(ns("vb_cats")),
        showcase = bs_icon("list-ul"),
        theme    = "primary"
      ),
      value_box(
        title    = "Highest overlap",
        value    = textOutput(ns("vb_max_pct")),
        showcase = bs_icon("trophy"),
        theme    = "#2C7BB6"
      ),
      value_box(
        title    = "Median pathway size",
        value    = textOutput(ns("vb_med_size")),
        showcase = bs_icon("diagram-3"),
        theme    = "secondary"
      ),
      value_box(
        title    = "Unique overlapping genes",
        value    = textOutput(ns("vb_genes")),
        showcase = bs_icon("universal-access-circle"),
        theme    = "dark"
      )
    ),

    # ── Row 3: Scatter (left) + DT table (right) ──────────────────────────
    layout_columns(
      col_widths = c(4, 8),

      card(
        card_header(
          tagList(bs_icon("graph-up"), " Pathway size vs. Overlap %"),
          tooltip(
            bs_icon("info-circle"),
            paste(
              "Each bubble is one pathway/category.",
              "X-axis: total genes in that pathway (log scale).",
              "Y-axis: % overlap with MorPhiC KO genes (or DEGs).",
              "Bubble size: raw overlap count.",
              "Click any bubble to populate the gene list below."
            )
          )
        ),
        plotlyOutput(ns("scatter"), height = "500px")
      ),

      card(
        card_header(
          tagList(bs_icon("table"), " Overlap Table"),
          tooltip(
            bs_icon("info-circle"),
            paste(
              "Default sort: Overlap % (descending).",
              "Click any column header to re-sort.",
              "Inline colour bars show relative values.",
              "Click a row to expand overlapping genes below."
            )
          )
        ),
        DTOutput(ns("tbl"))
      )
    ),

    # ── Row 3b: Pathway visualisations (bar chart + cross-assay heatmap) ──────
    navset_card_tab(
      title       = tagList(bs_icon("bar-chart-fill"), " Pathway Visualisations"),
      full_screen = TRUE,
      nav_panel(
        tagList(bs_icon("bar-chart-fill"), " Top Pathways"),
        plotlyOutput(ns("bar_chart"), height = "420px")
      ),
      nav_panel(
        tagList(bs_icon("grid-3x3-gap-fill"), " Cross-assay Heatmap"),
        plotlyOutput(ns("heatmap"), height = "420px")
      )
    ),

    # ── Row 4: Gene detail expansion ───────────────────────────────────────
    card(
      card_header(
        tagList(bs_icon("intersect"), " Overlapping Genes"),
        uiOutput(ns("gene_caption"))
      ),
      uiOutput(ns("gene_panel"))
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

overlap_table_dashboardServer <- function(id, con_r, study_info_r) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── Assay picker (DEG mode) ──────────────────────────────────────────────

    output$direction_ui <- renderUI({
      is_ko  <- isTRUE(input$mode == "ko")
      is_all <- isTRUE(tryCatch(input$assay, error = function(e) "") == "__all__")
      # Dim the control in KO mode (no directionality) or All-assay mode
      # (always uses Up & Down for cross-assay comparability)
      disabled <- is_ko || is_all
      note_txt <- if (is_ko) "\u2014 N/A in KO mode" else "\u2014 fixed Up & Down in all-assay mode"
      div(
        style = if (disabled) "opacity:0.35; pointer-events:none;" else "",
        radioButtons(
          ns("direction"), NULL,
          choices  = c("Up & Down" = "updn", "Up only" = "up", "Down only" = "down"),
          selected = if (is.null(isolate(input$direction))) "updn" else isolate(input$direction)
        ),
        if (disabled)
          tags$p(class = "text-muted small mb-0 mt-1", note_txt)
      )
    })

    # Assay info row — big value_box cards shown below the control row.
    # Single assay: one box per metadata field.
    # __all__: compact count boxes showing how many assays / genes / etc. exist.
    # KO mode: nothing rendered.
    output$assay_info_row <- renderUI({
      req(input$mode == "deg")
      assay_sel <- tryCatch(input$assay, error = function(e) "")
      req(nzchar(assay_sel))

      si_all <- study_info_r()

      if (isTRUE(assay_sel == "__all__")) {
        # ── Summary counts for "all assays" ─────────────────────────────────
        si_f <- si_all |> distinct(Assay, Gene, Model_system, KO_strat)
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          value_box(
            title    = "Total assays",
            value    = n_distinct(si_f$Assay),
            showcase = bs_icon("bar-chart-fill"),
            theme    = "primary"
          ),
          value_box(
            title    = "KO genes",
            value    = n_distinct(si_f$Gene),
            showcase = bs_icon("diagram-3"),
            theme    = "secondary"
          ),
          value_box(
            title    = "Model systems",
            value    = n_distinct(si_f$Model_system),
            showcase = bs_icon("grid-3x3-gap-fill"),
            theme    = "info"
          ),
          value_box(
            title    = "KO strategies",
            value    = n_distinct(si_f$KO_strat),
            showcase = bs_icon("scissors"),
            theme    = "dark"
          )
        )

      } else {
        # ── Single assay: one value_box per metadata field ───────────────────
        row <- si_all |> filter(Assay == assay_sel) |> head(1)
        req(nrow(row) > 0)

        na_dash <- function(x) {
          if (is.null(x) || (length(x) == 1 && is.na(x)) || !nzchar(as.character(x))) "\u2014" else as.character(x)
        }

        # Subtitle line: DPC · Replicate · Differentiation time point
        extra <- c(
          if (!is.na(row$DPC) && nzchar(as.character(row$DPC)))
            paste0("DPC: ", row$DPC),
          if ("Replicate" %in% names(row) && !is.na(row$Replicate))
            paste0("Rep: ", row$Replicate),
          if ("Differentation_time_point" %in% names(row) &&
              !is.na(row$Differentation_time_point) &&
              nzchar(as.character(row$Differentation_time_point)))
            paste0("TP: ", row$Differentation_time_point)
        )
        subtitle_txt <- if (length(extra) > 0) paste(extra, collapse = "  \u00b7  ") else NULL

        tagList(
          layout_columns(
            col_widths = c(3, 3, 3, 3),
            value_box(
              title    = "KO Gene",
              value    = na_dash(row$Gene),
              showcase = bs_icon("diagram-3"),
              theme    = "success"
            ),
            value_box(
              title    = "Model System",
              value    = na_dash(row$Model_system),
              showcase = bs_icon("grid-3x3-gap-fill"),
              theme    = "primary"
            ),
            value_box(
              title    = "KO Strategy",
              value    = na_dash(row$KO_strat),
              showcase = bs_icon("scissors"),
              theme    = "info"
            ),
            value_box(
              title    = "Cell Line",
              value    = na_dash(row$Cell_Line),
              showcase = bs_icon("card-list"),
              theme    = "secondary"
            )
          ),
          if (!is.null(subtitle_txt))
            div(class = "text-muted small text-center mt-n1 mb-2", subtitle_txt)
        )
      }
    })

    # Assay picker — grouped by Model System
    output$assay_ui <- renderUI({
      if (isTRUE(input$mode == "deg")) {
        si_all <- study_info_r() |> select(Assay, Model_system) |> distinct()
        by_sys <- split(si_all$Assay, si_all$Model_system)
        choices <- c(list("\u2014 All DEG Assays \u2014" = "__all__"), by_sys)
        selectInput(ns("assay"), NULL,
                    choices  = choices,
                    selected = "__all__",
                    width    = "100%")
      } else {
        tags$p(class = "text-muted small mt-2 mb-0", "\u2014 N/A in KO mode")
      }
    })

    # ── Data loaders ─────────────────────────────────────────────────────────

    # KO overlap: morphic_overlap.{source}
    # Columns differ per source → rename to standard names here.
    ko_data <- reactive({
      req(input$source)
      src <- input$source

      .empty <- data.frame(
        category = character(), n_pathway = integer(),
        n_overlap = integer(), pct = numeric(), genes = character(),
        genes_up = character(), genes_down = character()
      )

      tryCatch({
        df <- dbGetQuery(
          con_r(),
          sprintf('SELECT * FROM morphic_overlap."%s"', src)
        )
        df <- df |>
          rename(
            category  = all_of(.OT_GRP[[src]]),
            n_pathway = all_of(.OT_N[[src]]),
            n_overlap = overlap_n,
            pct       = pct_morphic_overlap,
            genes     = overlap_genes
          ) |>
          select(category, n_pathway, n_overlap, pct, genes)
        # KO mode has no directionality — add NA placeholders so downstream
        # code can always rely on these columns existing
        df$genes_up   <- NA_character_
        df$genes_down <- NA_character_
        df
      }, error = function(e) {
        showNotification(
          paste0("Could not load morphic_overlap.", src, ": ", conditionMessage(e)),
          type = "error", duration = 12
        )
        .empty
      })
    })

    # DEG overlap: deg_overlap.{assay}  (source stored in a column)
    # Each assay table has columns: source, category_name, category_genes_n,
    #   n_up/down/updn, pct_up/down/updn, overlap_up/down/updn.
    deg_data <- reactive({
      req(input$mode == "deg", input$source, input$direction, input$assay)
      # Guard: __all__ is handled by all_assay_deg_data() — do not proceed here
      req(!isTRUE(input$assay == "__all__"))

      dir       <- input$direction          # "up" | "down" | "updn"
      src       <- input$source
      assay_tbl <- input$assay

      n_col    <- paste0("n_",       dir)
      pct_col  <- paste0("pct_",     dir)
      gene_col <- paste0("overlap_", dir)

      # Map UI source key → stored source value + optional category_type filter
      db_source <- .OT_DEG_SOURCE[[src]]
      cat_type  <- .OT_DEG_CATTYPE[[src]]

      where_sql <- sprintf("WHERE source = '%s'", db_source)
      if (!is.na(cat_type))
        where_sql <- paste0(where_sql, sprintf(" AND category_type = '%s'", cat_type))

      # Always fetch overlap_up and overlap_down alongside the selected direction
      # column so the gene chip panel can colour-code each gene's direction.
      sql <- sprintf(
        'SELECT category_name     AS category,
                category_genes_n  AS n_pathway,
                "%s"              AS n_overlap,
                "%s"              AS pct,
                "%s"              AS genes,
                "overlap_up"      AS genes_up,
                "overlap_down"    AS genes_down
         FROM   deg_overlap."%s"
         %s',
        n_col, pct_col, gene_col, assay_tbl, where_sql
      )

      tryCatch(
        dbGetQuery(con_r(), sql),
        error = function(e) {
          showNotification(
            paste0("Could not load DEG overlap for assay '", assay_tbl, "': ", conditionMessage(e)),
            type = "error", duration = 12
          )
          data.frame(
            category = character(), n_pathway = integer(),
            n_overlap = integer(), pct = numeric(), genes = character(),
            genes_up = character(), genes_down = character()
          )
        }
      )
    })

    # All-assay aggregate (DEG mode, input$assay == "__all__"):
    # UNIONs every deg_overlap table for the selected source, then aggregates
    # by pathway → max overlap %, which assay achieved it, and how many assays
    # had any overlap.  Uses the "updn" direction columns (all DEGs) for fairness.
    all_assay_deg_data <- reactive({
      req(input$mode == "deg", input$source)
      req(isTRUE(tryCatch(input$assay, error = function(e) "") == "__all__"))

      src       <- input$source
      db_source <- .OT_DEG_SOURCE[[src]]
      cat_type  <- .OT_DEG_CATTYPE[[src]]

      # Discover all assay table names in deg_overlap schema
      tbls <- tryCatch(
        dbGetQuery(con_r(),
          "SELECT table_name
           FROM   information_schema.tables
           WHERE  table_schema = 'deg_overlap'
           ORDER  BY table_name"
        )$table_name,
        error = function(e) character(0)
      )
      req(length(tbls) > 0)

      where_extra <- if (!is.na(cat_type))
        sprintf(" AND category_type = '%s'", cat_type)
      else ""

      # Build UNION ALL — one SELECT per assay table
      union_parts <- vapply(tbls, function(tbl) {
        sprintf(
          "SELECT '%s'           AS assay,
                  category_name,
                  category_genes_n,
                  n_updn           AS n_overlap,
                  pct_updn         AS pct,
                  overlap_updn     AS genes
           FROM   deg_overlap.\"%s\"
           WHERE  source = '%s'%s",
          gsub("'", "''", tbl), tbl, db_source, where_extra
        )
      }, character(1))

      union_sql <- paste(union_parts, collapse = "\nUNION ALL\n")

      # Aggregate: for each pathway keep the assay with the highest overlap %
      agg_sql <- sprintf(
        "SELECT category_name                         AS category,
                MAX(category_genes_n)                 AS n_pathway,
                MAX(n_overlap)                        AS n_overlap,
                MAX(pct)                              AS pct,
                arg_max(assay, pct)                   AS best_assay,
                arg_max(genes, pct)                   AS genes,
                COUNT(*) FILTER (WHERE n_overlap > 0) AS n_assays
         FROM   (%s) t
         GROUP  BY category_name",
        union_sql
      )

      tryCatch({
        df <- dbGetQuery(con_r(), agg_sql)
        df$genes_up   <- NA_character_
        df$genes_down <- NA_character_
        df
      }, error = function(e) {
        showNotification(
          paste0("Cross-assay query failed: ", conditionMessage(e)),
          type = "error", duration = 15
        )
        data.frame(
          category   = character(), n_pathway  = integer(),
          n_overlap  = integer(),   pct        = numeric(),
          best_assay = character(), genes      = character(),
          n_assays   = integer(),   genes_up   = character(),
          genes_down = character()
        )
      })
    })

    # Unified dataset — pick based on mode / assay selection
    raw_data <- reactive({
      if (input$mode == "ko") return(ko_data())
      assay_sel <- tryCatch(input$assay, error = function(e) "")
      if (isTRUE(assay_sel == "__all__")) return(all_assay_deg_data())
      deg_data()
    })

    # Filter by min pathway size + min overlap count, sort by pct desc, slice top N
    display_data <- reactive({
      df <- raw_data()
      req(nrow(df) > 0)

      min_ov <- tryCatch(input$min_overlap, error = function(e) 1L)
      if (is.null(min_ov) || is.na(min_ov)) min_ov <- 1L

      df |>
        filter(n_pathway >= input$min_size,
               n_overlap >= min_ov) |>
        arrange(desc(pct), desc(n_overlap)) |>
        head(input$top_n) |>
        mutate(
          pct      = round(pct, 2),
          # Truncate very long category names to keep layouts clean
          category = ifelse(
            nchar(category) > 80,
            paste0(substr(category, 1, 77), "\u2026"),
            category
          )
        )
    })

    # ── KPI value boxes ───────────────────────────────────────────────────────

    output$vb_cats <- renderText(nrow(display_data()))

    output$vb_max_pct <- renderText({
      paste0(round(max(display_data()$pct, na.rm = TRUE), 1), "%")
    })

    output$vb_med_size <- renderText({
      formatC(
        median(display_data()$n_pathway, na.rm = TRUE),
        format = "d", big.mark = ","
      )
    })

    output$vb_genes <- renderText({
      gs   <- display_data()$genes
      gs   <- gs[!is.na(gs) & nzchar(gs)]
      gvec <- unique(trimws(unlist(strsplit(gs, ";"))))
      length(gvec[nzchar(gvec)])
    })

    # ── Scatter plot ──────────────────────────────────────────────────────────

    output$scatter <- renderPlotly({
      df <- display_data()
      req(nrow(df) > 0)

      df <- df |>
        mutate(
          row_idx   = seq_len(n()),
          tip_label = ifelse(
            nchar(category) > 55,
            paste0(substr(category, 1, 52), "\u2026"),
            category
          )
        )

      plot_ly(
        df,
        x             = ~n_pathway,
        y             = ~pct,
        type          = "scatter",
        mode          = "markers",
        size          = ~pmax(n_overlap, 1),
        sizes         = c(6, 36),
        text          = ~tip_label,
        customdata    = ~row_idx,
        hovertemplate = paste0(
          "<b>%{text}</b><br>",
          .OT_LABEL[[input$source]]$size, ": %{x:,d}<br>",
          "Overlap: %{y:.2f}%<br>",
          "<extra></extra>"
        ),
        marker = list(
          color      = ~pct,
          colorscale = list(
            c(0,   "#deebf7"),
            c(0.5, "#6baed6"),
            c(1,   "#08519c")
          ),
          showscale  = TRUE,
          colorbar   = list(
            title     = list(text = "Overlap %", side = "right"),
            thickness = 12
          ),
          opacity = 0.85,
          line    = list(width = 0.5, color = "#ffffff")
        ),
        source = ns("scatter_click")
      ) |>
        layout(
          xaxis = list(
            title      = paste0(.OT_LABEL[[input$source]]$size, " (log scale)"),
            type       = "log",
            automargin = TRUE
          ),
          yaxis = list(
            title      = "Overlap %",
            automargin = TRUE
          ),
          plot_bgcolor  = "white",
          paper_bgcolor = "white",
          margin        = list(t = 5, r = 10)
        ) |>
        config(displayModeBar = FALSE)
    })

    # ── DT table ──────────────────────────────────────────────────────────────

    output$tbl <- renderDT({
      df  <- display_data()
      req(nrow(df) > 0)

      lbl      <- .OT_LABEL[[input$source]]
      cat_lbl  <- lbl$cat   # e.g. "Pathway", "GO Term (BP)", "Panel" …
      size_lbl <- lbl$size  # e.g. "Pathway genes", "Class genes" …

      is_all <- isTRUE(tryCatch(input$assay, error = function(e) "") == "__all__")

      disp <- df
      names(disp)[names(disp) == "category"]  <- cat_lbl
      names(disp)[names(disp) == "n_pathway"] <- size_lbl
      names(disp)[names(disp) == "n_overlap"] <- "Overlap genes"
      names(disp)[names(disp) == "pct"]       <- "Overlap %"

      if (is_all && "best_assay" %in% names(disp)) {
        names(disp)[names(disp) == "best_assay"] <- "Best assay"
        names(disp)[names(disp) == "n_assays"]   <- "Assays w/ overlap"
        col_order <- c(cat_lbl, size_lbl, "Overlap genes", "Overlap %",
                       "Best assay", "Assays w/ overlap")
      } else {
        col_order <- c(cat_lbl, size_lbl, "Overlap genes", "Overlap %")
      }
      disp <- disp[, col_order]

      pct_range  <- c(0, max(disp[["Overlap %"]], na.rm = TRUE))
      size_range <- range(disp[[size_lbl]],        na.rm = TRUE)
      n_cols     <- ncol(disp)

      datatable(
        disp,
        rownames   = FALSE,
        filter     = "top",
        selection  = "single",
        extensions = "Buttons",
        options    = list(
          dom       = "Bft",
          paging    = FALSE,          # all rows visible → row index is always reliable
          scrollY   = "430px",
          scrollX   = TRUE,
          buttons   = list(
            list(extend = "csv",   filename = "overlap_table", text = "\u2b07 CSV"),
            list(extend = "excel", filename = "overlap_table", text = "\u2b07 Excel")
          ),
          order       = list(list(3, "desc")),  # default: Overlap % descending
          columnDefs  = list(
            list(width = "38%",     targets = 0),
            list(className = "dt-right", targets = seq(1, n_cols - 1))
          )
        )
      ) |>
        formatStyle(
          "Overlap %",
          background         = styleColorBar(pct_range, "#6baed6"),
          backgroundSize     = "100% 55%",
          backgroundRepeat   = "no-repeat",
          backgroundPosition = "center"
        ) |>
        formatStyle(
          size_lbl,
          background         = styleColorBar(size_range, "#fdae6b"),
          backgroundSize     = "100% 55%",
          backgroundRepeat   = "no-repeat",
          backgroundPosition = "center"
        ) |>
        formatRound("Overlap %", digits = 2)
    })

    # ── Bar chart: top pathways ranked by overlap % ──────────────────────────
    output$bar_chart <- renderPlotly({
      df <- display_data()
      req(nrow(df) > 0)

      size_lbl <- .OT_LABEL[[input$source]]$size

      # Sort ascending so the highest bar ends up at the top in horizontal layout.
      # Use make.unique() before factor() — GO/pathway names often share long
      # prefixes that become identical after truncation, causing duplicate factor
      # levels which silently kills plotly rendering.
      df <- df |>
        arrange(pct) |>
        mutate(
          lbl = ifelse(nchar(category) > 55,
                       paste0(substr(category, 1, 52), "\u2026"),
                       category),
          lbl = make.unique(lbl, sep = " "),
          lbl = factor(lbl, levels = lbl)
        )

      hover_txt <- paste0(
        "<b>", df$category, "</b><br>",
        "Overlap %: ",    df$pct, "<br>",
        "Overlap genes: ", df$n_overlap, "<br>",
        size_lbl, ": ",   formatC(df$n_pathway, format = "d", big.mark = ","),
        "<extra></extra>"
      )

      plot_ly(
        df, type = "bar", orientation = "h",
        y             = ~lbl,
        x             = ~pct,
        text          = ~n_overlap,
        textposition  = "outside",
        hovertemplate = hover_txt,
        marker = list(
          color      = ~pct,
          colorscale = list(c(0, "#deebf7"), c(0.5, "#6baed6"), c(1, "#08519c")),
          showscale  = TRUE,
          colorbar   = list(title = list(text = "Overlap %", side = "right"),
                            thickness = 12)
        )
      ) |>
        layout(
          xaxis = list(title = "Overlap %", automargin = TRUE),
          yaxis = list(title = "", automargin = TRUE, tickfont = list(size = 9)),
          plot_bgcolor  = "white",
          paper_bgcolor = "white",
          margin = list(l = 5, r = 60, t = 10, b = 30),
          uniformtext  = list(minsize = 8, mode = "hide")
        ) |>
        config(displayModeBar = FALSE)
    })

    # ── Cross-assay heatmap ────────────────────────────────────────────────────
    # Fetches per-assay overlap % for the currently displayed top pathways.
    # Only fires in DEG "all assays" mode (requires input$assay == "__all__").
    multi_assay_detail <- reactive({
      req(input$mode == "deg")
      req(isTRUE(tryCatch(input$assay, error = function(e) "") == "__all__"))

      top_cats <- display_data()$category
      req(length(top_cats) > 0)

      src       <- input$source
      db_source <- .OT_DEG_SOURCE[[src]]
      cat_type  <- .OT_DEG_CATTYPE[[src]]

      # Same filtered table list as the aggregate query
      tbls <- tryCatch(
        dbGetQuery(con_r(),
          "SELECT table_name FROM information_schema.tables
           WHERE  table_schema = 'deg_overlap' ORDER BY table_name"
        )$table_name,
        error = function(e) character(0)
      )
      req(length(tbls) > 0)

      where_extra <- if (!is.na(cat_type))
        sprintf(" AND category_type = '%s'", cat_type) else ""

      safe_cats <- gsub("'", "''", top_cats)
      cats_sql  <- paste0("'", paste(safe_cats, collapse = "','"), "'")

      union_parts <- vapply(tbls, function(tbl) {
        sprintf(
          "SELECT '%s' AS assay, category_name AS category, pct_updn AS pct
           FROM   deg_overlap.\"%s\"
           WHERE  source = '%s'%s AND category_name IN (%s)",
          gsub("'", "''", tbl), tbl, db_source, where_extra, cats_sql
        )
      }, character(1))

      tryCatch(
        dbGetQuery(con_r(), paste(union_parts, collapse = "\nUNION ALL\n")),
        error = function(e) NULL
      )
    })

    output$heatmap <- renderPlotly({
      assay_sel <- tryCatch(input$assay, error = function(e) "")
      is_all    <- isTRUE(assay_sel == "__all__") &&
                   isTRUE(tryCatch(input$mode, error = function(e) "") == "deg")

      # Placeholder when not in all-assay DEG mode — message varies by context
      if (!is_all) {
        placeholder_msg <- if (isTRUE(tryCatch(input$mode, error = function(e) "") == "ko")) {
          "Cross-assay heatmap is only available in <b>DEG mode</b> \u2014 switch Mode above and select <b>\u2014 All DEG Assays \u2014</b>"
        } else {
          "Select <b>\u2014 All DEG Assays \u2014</b> from the Assay picker to see the cross-assay heatmap"
        }
        return(
          plotly_empty(type = "scatter", mode = "markers") |>
            layout(
              title = list(
                text = placeholder_msg,
                font = list(size = 13, color = "#6c757d"), x = 0.5, xanchor = "center"
              ),
              paper_bgcolor = "white", plot_bgcolor = "white"
            ) |>
            config(displayModeBar = FALSE)
        )
      }

      det <- multi_assay_detail()
      req(!is.null(det) && nrow(det) > 0)

      top_cats <- display_data()$category
      assays   <- sort(unique(det$assay))

      # Shorten assay labels: strip common prefix patterns for readability
      short_assay <- function(x) {
        if (nchar(x) > 32) paste0(substr(x, 1, 30), "\u2026") else x
      }
      assay_lbls <- vapply(assays, short_assay, character(1))

      # Build matrix: rows = pathways (top → bottom), cols = assays
      mat <- matrix(NA_real_, nrow = length(top_cats), ncol = length(assays),
                    dimnames = list(top_cats, assays))
      for (i in seq_len(nrow(det))) {
        r <- match(det$category[i], top_cats)
        c <- match(det$assay[i],    assays)
        if (!is.na(r) && !is.na(c)) mat[r, c] <- det$pct[i]
      }

      cat_lbls <- ifelse(nchar(top_cats) > 55,
                         paste0(substr(top_cats, 1, 52), "\u2026"),
                         top_cats)

      # NA → 0 for display (white cells = no overlap in that assay)
      mat_disp <- mat
      mat_disp[is.na(mat_disp)] <- 0

      plot_ly(
        x    = assay_lbls,
        y    = rev(cat_lbls),
        z    = mat_disp[rev(seq_len(nrow(mat_disp))), ],
        type = "heatmap",
        colorscale  = list(c(0, "#f7fbff"), c(0.3, "#9ecae1"),
                           c(0.7, "#3182bd"), c(1, "#08519c")),
        colorbar    = list(title = list(text = "Overlap %"), thickness = 12),
        hoverongaps = FALSE,
        hovertemplate = "<b>%{y}</b><br>Assay: %{x}<br>Overlap: %{z:.2f}%<extra></extra>"
      ) |>
        layout(
          xaxis = list(title = "", tickangle = -50,
                       tickfont = list(size = 8), automargin = TRUE),
          yaxis = list(title = "", tickfont = list(size = 8), automargin = TRUE),
          plot_bgcolor  = "white",
          paper_bgcolor = "white",
          margin        = list(l = 10, r = 30, b = 110, t = 10)
        ) |>
        config(displayModeBar = FALSE)
    })

    # ── Row / point selection → gene panel ────────────────────────────────────

    selected_row <- reactiveVal(NULL)

    # DT row click
    observeEvent(input$tbl_rows_selected, {
      sel <- input$tbl_rows_selected
      df  <- display_data()
      if (length(sel) == 1 && sel <= nrow(df)) {
        selected_row(df[sel, ])
      }
    })

    # Scatter point click
    observe({
      ev  <- event_data("plotly_click", source = ns("scatter_click"))
      req(ev)
      idx <- as.integer(ev$customdata)
      df  <- display_data()
      if (!is.null(idx) && idx >= 1 && idx <= nrow(df)) {
        selected_row(df[idx, ])
      }
    })

    # Clear gene panel whenever the underlying data changes
    observeEvent(display_data(), { selected_row(NULL) }, ignoreInit = TRUE)

    # ── logFC lookup for selected row (DEG mode only) ─────────────────────────
    # Queries main.{assay} for log2FoldChange of every gene in the selected row.
    # The fold-change column is named {GENE}_{STRATEGY}_log2FoldChange and
    # discovered dynamically via information_schema.
    selected_logfc <- reactive({
      req(input$mode == "deg")
      row <- selected_row()
      req(row)

      # In "all assays" mode use the best_assay stored in the selected row;
      # in single-assay mode use the picker value directly.
      assay_sel <- tryCatch(input$assay, error = function(e) "")
      assay_tbl <- if (isTRUE(assay_sel == "__all__")) {
        ba <- row$best_assay
        req(!is.null(ba) && !is.na(ba) && nzchar(ba))
        ba
      } else {
        assay_sel
      }
      req(assay_tbl, nzchar(assay_tbl))

      gs <- row$genes
      if (is.na(gs) || !nzchar(trimws(gs))) return(NULL)

      gene_vec <- unique(trimws(strsplit(gs, ";")[[1]]))
      gene_vec <- gene_vec[nzchar(gene_vec)]
      req(length(gene_vec) > 0)

      # Discover the log2FoldChange column name for this assay table
      col_df <- tryCatch(
        dbGetQuery(con_r(), sprintf(
          "SELECT column_name
           FROM   information_schema.columns
           WHERE  table_schema = 'main'
             AND  table_name   = '%s'
             AND  column_name  LIKE '%%log2FoldChange'",
          assay_tbl
        )),
        error = function(e) NULL
      )
      if (is.null(col_df) || nrow(col_df) == 0) return(NULL)
      logfc_col <- col_df$column_name[1]

      # Safely quote gene symbols (escape any embedded single-quotes)
      safe <- gsub("'", "''", gene_vec)
      gene_sql <- paste0("'", paste(safe, collapse = "','"), "'")

      tryCatch(
        dbGetQuery(con_r(), sprintf(
          'SELECT hgnc_symbol,
                  "%s" AS log2fc
           FROM   main."%s"
           WHERE  hgnc_symbol IN (%s)',
          logfc_col, assay_tbl, gene_sql
        )),
        error = function(e) NULL
      )
    })

    # Caption above gene chips
    output$gene_caption <- renderUI({
      row <- selected_row()
      if (is.null(row)) {
        tags$span(
          class = "text-muted small fst-italic",
          " \u2014 select a table row or scatter bubble to expand"
        )
      } else {
        # Extra badges for the all-assay mode
        best_badge <- if (!is.null(row$best_assay) && !is.na(row$best_assay) &&
                         nzchar(row$best_assay)) {
          lbl <- if (nchar(row$best_assay) > 40)
            paste0(substr(row$best_assay, 1, 38), "\u2026")
          else row$best_assay
          tags$span(class = "badge ms-1",
                    style = "background:#2d6a4f; color:white;",
                    paste0("Best: ", lbl))
        }
        n_assays_badge <- if (!is.null(row$n_assays) && !is.na(row$n_assays)) {
          tags$span(class = "badge bg-warning text-dark ms-1",
                    paste0(row$n_assays, " assay", if (row$n_assays != 1) "s"))
        }
        tagList(
          tags$span(class = "fw-bold", row$category),
          tags$span(class = "badge bg-primary   ms-2", paste0(row$n_overlap, " overlapping")),
          tags$span(class = "badge bg-secondary  ms-1", paste0(row$n_pathway, " in pathway")),
          tags$span(
            class = "badge ms-1",
            style = "background:#2C7BB6; color:white;",
            paste0(row$pct, "% overlap")
          ),
          best_badge,
          n_assays_badge
        )
      }
    })

    # Gene chip badges
    output$gene_panel <- renderUI({
      row <- selected_row()

      if (is.null(row)) {
        return(tags$p(
          class = "text-muted small p-3 mb-0",
          "Click a row in the table or a bubble in the scatter plot to see the overlapping genes here."
        ))
      }

      gs <- row$genes
      if (is.na(gs) || !nzchar(trimws(gs))) {
        return(tags$p(class = "text-muted small p-3 mb-0",
                      "No gene names recorded for this category."))
      }

      gene_vec <- sort(unique(trimws(strsplit(gs, ";")[[1]])))
      gene_vec <- gene_vec[nzchar(gene_vec)]

      # ── Direction sets (DEG mode) ────────────────────────────────────────────
      has_dir <- isTRUE(!is.na(row$genes_up)   && nzchar(trimws(row$genes_up)))  ||
                 isTRUE(!is.na(row$genes_down)  && nzchar(trimws(row$genes_down)))

      parse_genes <- function(x) {
        if (is.null(x) || is.na(x) || !nzchar(trimws(x))) return(character(0))
        v <- trimws(strsplit(x, ";")[[1]]); v[nzchar(v)]
      }
      up_set   <- parse_genes(row$genes_up)
      down_set <- parse_genes(row$genes_down)

      # ── logFC lookup ─────────────────────────────────────────────────────────
      lfc_df  <- tryCatch(selected_logfc(), error = function(e) NULL)
      lfc_map <- if (!is.null(lfc_df) && nrow(lfc_df) > 0)
        setNames(lfc_df$log2fc, lfc_df$hgnc_symbol)
      else
        setNames(numeric(0), character(0))

      has_lfc <- length(lfc_map) > 0

      # Biologically meaningful 3-level colour scale per direction:
      #   level 1 = |logFC| < 1  (< 2-fold)  → light shade
      #   level 2 = 1 ≤ |logFC| < 2 (2-4x)  → medium shade
      #   level 3 = |logFC| ≥ 2  (> 4-fold)  → full/dark shade
      dir_color <- function(is_up, is_down, lfc_abs) {
        lvl <- if (lfc_abs >= 2) 3L else if (lfc_abs >= 1) 2L else 1L
        if (is_up && is_down)
          c("#C39BD3", "#9B59B6", "#7B2D8B")[lvl]   # purple scale
        else if (is_up)
          c("#F1948A", "#E74C3C", "#C0392B")[lvl]   # red scale
        else
          c("#85C1E9", "#5DADE2", "#2C7BB6")[lvl]   # blue scale
      }

      chip_style <- function(bg)
        paste0("background-color:", bg, "; color:white;",
               "font-size:0.78rem; font-weight:500;",
               "padding:0.3em 0.65em; border-radius:0.3rem;",
               "white-space:nowrap;")

      # ── Build chips ──────────────────────────────────────────────────────────
      chips <- lapply(gene_vec, function(g) {
        is_up   <- g %in% up_set
        is_down <- g %in% down_set
        lfc     <- if (has_lfc) lfc_map[g] else NA_real_
        lfc_abs <- if (!is.na(lfc)) abs(lfc) else 0

        # Direction arrow + logFC text
        arrow <- if (has_dir) {
          if      (is_up && is_down) " \u2191\u2193"
          else if (is_up)            " \u2191"
          else                       " \u2193"
        } else ""

        lfc_txt <- if (!is.na(lfc))
          paste0(" (", ifelse(lfc >= 0, "+", ""), round(lfc, 2), ")")
        else ""

        bg <- if (has_dir) dir_color(is_up, is_down, lfc_abs) else "#2C7BB6"

        tags$span(class = "badge me-1 mb-1",
                  style = chip_style(bg),
                  paste0(g, arrow, lfc_txt))
      })

      # ── Legend ───────────────────────────────────────────────────────────────
      legend <- if (has_dir || has_lfc) {
        dir_items <- if (has_dir) list(
          tags$span(class = "badge", style = chip_style("#C0392B"), "\u2191 Up"),
          tags$span(class = "badge", style = chip_style("#2C7BB6"), "\u2193 Down"),
          tags$span(class = "badge", style = chip_style("#7B2D8B"), "\u2191\u2193 Both")
        ) else list()

        lfc_items <- if (has_lfc) list(
          tags$span(class = "text-muted small ms-3 me-1",
                    style = "font-weight:600;", "log\u2082FC:"),
          tags$span(class = "badge", style = chip_style("#F1948A"), "< 1  (light)"),
          tags$span(class = "badge", style = chip_style("#E74C3C"), "1\u20132 (medium)"),
          tags$span(class = "badge", style = chip_style("#C0392B"), "\u2265 2  (dark)")
        ) else list()

        do.call(div, c(
          list(class = "d-flex gap-2 align-items-center px-3 pt-2 pb-2 flex-wrap",
               style = "border-bottom:1px solid #dee2e6; font-size:0.78rem;",
               tags$span(style = "font-weight:600; color:#6c757d;", "Key:")),
          dir_items,
          lfc_items
        ))
      }

      tagList(
        legend,
        div(class = "p-3 d-flex flex-wrap", chips)
      )
    })

    # ── Suspend all outputs until the tab is actually visible ─────────────────
    # bslib renders every nav_panel's HTML at startup so Shiny subscribes to
    # ALL outputs immediately — even on inactive tabs — which would fire DB
    # queries before the user ever opens this tab.  suspendWhenHidden = TRUE
    # defers evaluation until the output enters the viewport.
    # Must be called AFTER all output$... assignments are registered above.
    invisible(lapply(
      c("direction_ui", "assay_ui", "assay_info_row",
        "vb_cats", "vb_max_pct", "vb_med_size", "vb_genes",
        "scatter", "tbl", "bar_chart", "heatmap",
        "gene_caption", "gene_panel"),
      function(id) outputOptions(output, id, suspendWhenHidden = TRUE)
    ))

  })
}
