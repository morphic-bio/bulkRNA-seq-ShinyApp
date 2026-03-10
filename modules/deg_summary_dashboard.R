# =============================================================================
# Module: deg_summary_dashboard
# Description: Interactive DEG summary across all bulk RNA-seq assays.
#
#   1. N DEGs per Assay — bar chart, assay on X (or group axis),
#      filterable by Up / Down / All, switchable X-axis grouping
#   2. Most Recurrent DEGs — horizontal bar of genes appearing across
#      most assays; filterable by Up / Down / All and by Model System
#
#   Gene labels: hgnc_symbol when available, gene_ID otherwise.
#
# Data contract
#   con_r        — reactive: open read-only DuckDB connection
#   study_info_r — reactive: study_info data.frame
# =============================================================================

library(plotly)
library(DT)
library(UpSetR)

# ── Colours ───────────────────────────────────────────────────────────────────
.DEG_CLRS <- list(
  up   = "#E63946",
  down = "#457B9D",
  all  = "#6A4C93"
)
.GRP_PAL <- c("#4E79A7","#F28E2B","#E15759","#76B7B2",
              "#59A14F","#EDC948","#B07AA1","#FF9DA7",
              "#9C755F","#BAB0AC")

# ── SQL helper: display label = hgnc_symbol, fallback to gene_ID ──────────────
.display_sym_sql <-
  "CASE WHEN hgnc_symbol IS NULL OR hgnc_symbol = ''
        THEN gene_ID ELSE hgnc_symbol END"

# ── logFC cell helpers ────────────────────────────────────────────────────────

# Format one direction-matrix cell: "↑ +2.31", "↓ -1.45", "↑↓", or "".
.fmt_dir_cell <- function(is_up, is_dn, log2fc = NA_real_) {
  if (!is_up && !is_dn) return("")
  if (is_up  && is_dn)  return("\u2191\u2193")
  dir <- if (is_up) "\u2191" else "\u2193"
  if (is.na(log2fc)) return(dir)
  paste0(dir, " ", ifelse(log2fc >= 0, "+", ""), round(log2fc, 2))
}

# JavaScript rowCallback that reads the first character of every cell and
# applies a 3-level direction+magnitude background colour.
#   "↑ +2.31"  →  3-level red   (|logFC| < 1 / 1–2 / ≥ 2)
#   "↓ -1.45"  →  3-level blue
#   "↑↓"       →  flat purple
.LFC_ROW_CB <- DT::JS(
  "function(row, data, index) {",
  "  var UP='\u2191', DN='\u2193', UPDN='\u2191\u2193';",
  "  $(row).find('td').each(function() {",
  "    var v=$(this).text().trim(); if(!v) return;",
  "    var both=v.substring(0,2)===UPDN, isUp=!both&&v.charAt(0)===UP, isDn=!both&&v.charAt(0)===DN;",
  "    if(!isUp&&!isDn&&!both) return;",
  "    var bg,col;",
  "    if(both){bg='#f3e8fd';col='#7B2D8B';}",
  "    else{",
  "      var lfc=v.indexOf(' ')>-1?Math.abs(parseFloat(v.split(' ')[1])):0;",
  "      var l=lfc>=2?2:(lfc>=1?1:0);",
  "      if(isUp){bg=['#fde8ea','#f8c1c4','#f4989b'][l];col=['#E74C3C','#C0392B','#922B21'][l];}",
  "      else    {bg=['#ddeef7','#b8d8ef','#92c4e8'][l];col=['#457B9D','#2C6E8A','#1A5276'][l];}",
  "    }",
  "    $(this).css({'background-color':bg,'color':col,'font-weight':'bold','text-align':'center'});",
  "  });",
  "}"
)

# ── Data loaders ──────────────────────────────────────────────────────────────

#' Per-assay DEG counts (n_up, n_down).
load_deg_counts <- function(con, assay_names) {
  rows <- lapply(assay_names, function(tbl) {
    res <- tryCatch(
      dbGetQuery(con, paste0(
        'SELECT DEG, COUNT(*) AS n
         FROM "', tbl, '"
         WHERE DEG IN (\'up\',\'down\')
         GROUP BY DEG'
      )),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res) == 0) return(NULL)
    up   <- res[res$DEG == "up",   , drop = FALSE]
    down <- res[res$DEG == "down", , drop = FALSE]
    data.frame(
      assay  = tbl,
      n_up   = if (nrow(up)   > 0) up$n[1]   else 0L,
      n_down = if (nrow(down) > 0) down$n[1] else 0L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' All DEG gene display-labels per assay (hgnc_symbol or gene_ID fallback).
#' Also fetches log2FoldChange — column name is {GENE}_{STRATEGY}_log2FoldChange
#' and is discovered dynamically via dbListFields().
load_deg_genes <- function(con, assay_names) {
  rows <- lapply(assay_names, function(tbl) {
    # Discover the log2FoldChange column for this assay table
    col_names <- tryCatch(dbListFields(con, tbl), error = function(e) character(0))
    lfc_col   <- col_names[grepl("log2FoldChange$", col_names)]
    lfc_clause <- if (length(lfc_col) > 0)
      sprintf(', "%s" AS log2fc', lfc_col[1])
    else
      ", CAST(NULL AS DOUBLE) AS log2fc"

    res <- tryCatch(
      dbGetQuery(con, paste0(
        'SELECT ', .display_sym_sql, ' AS symbol, DEG', lfc_clause,
        ' FROM "', tbl, '"',
        " WHERE DEG IN ('up','down')"
      )),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res$assay <- tbl
    res
  })
  do.call(rbind, rows)
}

# ── UI ────────────────────────────────────────────────────────────────────────

deg_summary_dashboardUI <- function(id) {
  ns <- NS(id)

  tagList(

    # ── Load button + status ──────────────────────────────────────────────────
    layout_columns(
      col_widths = c(3, 9),
      card(
        card_body(
          class = "py-2",
          actionButton(ns("load_data"), "Load / Refresh DEG Data",
                       class = "btn-primary w-100",
                       icon  = icon("rotate")),
          uiOutput(ns("load_status"))
        )
      ),
      uiOutput(ns("vb_row"))   # value boxes appear after load
    ),

    # ── Chart 1: N DEGs per Assay ─────────────────────────────────────────────
    card(
      card_header(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
        span("N DEGs per Assay"),
        div(
          class = "d-flex gap-2 align-items-center flex-wrap",
          # Direction filter
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("Direction:", class = "mb-0 small fw-semibold"),
            radioButtons(ns("deg_dir"), NULL,
                         choices  = c("All" = "all", "Up" = "up", "Down" = "down"),
                         selected = "all", inline = TRUE)
          ),
          # X-axis grouping
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("X-axis:", class = "mb-0 small fw-semibold"),
            selectInput(ns("x_axis"), NULL,
                        choices  = c("Assay" = "assay",
                                     "Model System" = "Model_system",
                                     "KO Strategy"  = "KO_strat",
                                     "DPC"           = "DPC"),
                        selected = "assay",
                        width    = "160px")
          )
        )
      ),
      plotlyOutput(ns("deg_bar"), height = "420px"),
      full_screen = TRUE
    ),

    # ── Chart 2: Most Recurrent DEGs ─────────────────────────────────────────
    card(
      card_header(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
        span("Most Recurrent DEGs Across Assays"),
        div(
          class = "d-flex gap-2 align-items-center flex-wrap",
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("Direction:", class = "mb-0 small fw-semibold"),
            selectInput(ns("overlap_dir"), NULL,
                        choices  = c("All" = "all", "Up only" = "up", "Down only" = "down"),
                        selected = "all", width = "120px")
          ),
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("Model System:", class = "mb-0 small fw-semibold"),
            uiOutput(ns("overlap_ms_ui"))
          ),
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("Top N genes:", class = "mb-0 small fw-semibold"),
            numericInput(ns("top_n"), NULL, value = 30,
                         min = 5, max = 150, step = 5, width = "80px")
          )
        )
      ),
      layout_columns(
        col_widths = c(8, 4),
        plotlyOutput(ns("overlap_bar"), height = "500px"),
        div(
          class = "p-2",
          h6("Genes in selection", class = "fw-bold"),
          DTOutput(ns("overlap_tbl"))
        )
      ),
      full_screen = TRUE
    ),

    # ── Card 2b: UpSet — DEG gene intersections grouped by study ─────────────
    card(
      card_header(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
        span("DEG Intersection UpSet Plot"),
        div(
          class = "d-flex gap-2 align-items-center flex-wrap",
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("Direction:", class = "mb-0 small fw-semibold"),
            selectInput(ns("upset_dir"), NULL,
                        choices  = c("All" = "all", "Up only" = "up",
                                     "Down only" = "down"),
                        selected = "all", width = "110px")
          ),
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("Group by:", class = "mb-0 small fw-semibold"),
            selectInput(ns("upset_grp"), NULL,
                        choices  = c("Model System" = "Model_system",
                                     "KO Strategy"  = "KO_strat",
                                     "KO Gene"      = "Gene",
                                     "DPC"          = "DPC"),
                        selected = "Model_system", width = "150px")
          ),
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("Top N sets:", class = "mb-0 small fw-semibold"),
            numericInput(ns("upset_n"), NULL, value = 10,
                         min = 2, max = 30, step = 1, width = "75px")
          )
        )
      ),
      tags$p(class = "text-muted small px-3 pt-2 mb-0",
             "Each bar = intersection size (genes that are DEGs in exactly those groups). ",
             "Left bars = total DEGs per group. Grouped by the chosen study metadata."),
      uiOutput(ns("upset_ui")),
      full_screen = TRUE
    ),

    # ── Card 3: Cross-assay DEG overlap table ─────────────────────────────────
    # Rows = DEG gene symbols. Columns = KO-target genes (57 targets).
    # Cell value = "↑" / "↓" / "↑↓" / "" showing direction(s) that DEG
    # appears in relative to each KO-target, across all its assays.
    card(
      card_header(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
        span("DEG Overlap Across KO Targets"),
        div(
          class = "d-flex gap-2 align-items-center flex-wrap",
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("Direction:", class = "mb-0 small fw-semibold"),
            selectInput(ns("ov_dir"), NULL,
                        choices  = c("All" = "all", "Up only" = "up",
                                     "Down only" = "down"),
                        selected = "all", width = "120px")
          ),
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("Min KO targets (≥):", class = "mb-0 small fw-semibold"),
            numericInput(ns("ov_min_ko"), NULL, value = 2,
                         min = 1, max = 57, step = 1, width = "70px")
          ),
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("Model System:", class = "mb-0 small fw-semibold"),
            uiOutput(ns("ov_ms_ui"))
          )
        )
      ),
      div(
        class = "p-2",
        tags$p(class = "text-muted small mb-2",
               "Each row is a DEG gene. Each column is a KO-target gene. ",
               tags$b("↑"), " = up-regulated, ",
               tags$b("↓"), " = down-regulated, ",
               tags$b("↑↓"), " = both (across different assays of that KO-target). ",
               "N = number of KO-targets the DEG appears in. ",
               "Use the column search boxes to filter to specific KO-targets."),
        DTOutput(ns("ov_table"))
      ),
      full_screen = TRUE
    ),

    # ── Card 4: Assay-level comparison ────────────────────────────────────────
    # Pick 2–N assays; see set-size bars + shared/unique breakdown, then drill
    # into the actual genes shared across all selected assays.
    card(
      card_header(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
        span("Assay-level DEG Comparison"),
        div(
          class = "d-flex gap-2 align-items-center flex-wrap",
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("Direction:", class = "mb-0 small fw-semibold"),
            selectInput(ns("cmp_dir"), NULL,
                        choices  = c("All" = "all", "Up only" = "up",
                                     "Down only" = "down"),
                        selected = "all", width = "120px")
          )
        )
      ),
      card_body(
        # Assay picker (grouped by Model System, rendered after load)
        uiOutput(ns("cmp_assay_ui")),
        # Summary bar chart: total DEGs per assay, colour = unique vs shared
        uiOutput(ns("cmp_summary_ui")),
        # Breakdown tabs: shared-all / per-assay-unique / full gene matrix
        uiOutput(ns("cmp_tabs_ui"))
      ),
      full_screen = TRUE
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

deg_summary_dashboardServer <- function(id, con_r, study_info_r) {
  moduleServer(id, function(input, output, session) {

    # ── Stored data after load ────────────────────────────────────────────────
    deg_counts_rv <- reactiveVal(NULL)
    deg_genes_rv  <- reactiveVal(NULL)

    # ── Load on button press ──────────────────────────────────────────────────
    observeEvent(input$load_data, {
      con    <- con_r()
      si     <- study_info_r()
      assays <- unique(si$Assay)

      withProgress(message = "Loading DEG data…", value = 0, {
        setProgress(0.1, detail = "Counting DEGs")
        counts <- load_deg_counts(con, assays)

        setProgress(0.55, detail = "Loading gene lists")
        genes <- load_deg_genes(con, assays)

        setProgress(0.9, detail = "Joining metadata")
        meta <- si[!duplicated(si$Assay),
                   c("Assay","Gene","Model_system","KO_strat","DPC",
                     "Differentation_time_point","Cell_Line")]

        counts <- merge(counts, meta, by.x = "assay", by.y = "Assay", all.x = TRUE)
        counts$n_all <- counts$n_up + counts$n_down

        genes  <- merge(genes, meta, by.x = "assay", by.y = "Assay", all.x = TRUE)

        deg_counts_rv(counts)
        deg_genes_rv(genes)
        setProgress(1)
      })
    })

    # ── Load status line ──────────────────────────────────────────────────────
    output$load_status <- renderUI({
      counts <- deg_counts_rv()
      if (is.null(counts))
        tags$p(class = "text-muted small mt-1 mb-0",
               "Click to load DEG data from all 167 assay tables.")
      else
        tags$p(class = "text-success small mt-1 mb-0",
               icon("check-circle"),
               sprintf(" %d assays loaded.", nrow(counts)))
    })

    # ── Value boxes ───────────────────────────────────────────────────────────
    output$vb_row <- renderUI({
      counts <- deg_counts_rv()
      genes  <- deg_genes_rv()
      if (is.null(counts)) return(div())
      layout_columns(
        col_widths = c(3, 3, 3, 3),
        value_box("Assays Loaded", nrow(counts),
                  showcase = bsicons::bs_icon("bar-chart-fill"), theme = "primary"),
        value_box("Total Up DEGs", format(sum(counts$n_up), big.mark = ","),
                  showcase = bsicons::bs_icon("arrow-up-circle-fill"),
                  theme = value_box_theme(bg = .DEG_CLRS$up,   fg = "#fff")),
        value_box("Total Down DEGs", format(sum(counts$n_down), big.mark = ","),
                  showcase = bsicons::bs_icon("arrow-down-circle-fill"),
                  theme = value_box_theme(bg = .DEG_CLRS$down, fg = "#fff")),
        value_box("Unique DEG Genes",
                  format(length(unique(na.omit(genes$symbol))), big.mark = ","),
                  showcase = bsicons::bs_icon("diagram-3"), theme = "secondary")
      )
    })

    # ── Chart 1: N DEGs per Assay / Group ────────────────────────────────────
    output$deg_bar <- renderPlotly({
      counts <- deg_counts_rv()
      req(counts)

      dir   <- input$deg_dir
      x_var <- input$x_axis

      y_col  <- switch(dir, all = "n_all", up = "n_up", down = "n_down")
      y_lab  <- switch(dir,
                       all  = "N DEGs (up + down)",
                       up   = "N Up DEGs",
                       down = "N Down DEGs")
      bar_colour <- switch(dir,
                           all  = .DEG_CLRS$all,
                           up   = .DEG_CLRS$up,
                           down = .DEG_CLRS$down)

      if (x_var == "assay") {
        # ── Individual assay bars ─────────────────────────────────────────
        # Colour by Model System; sort descending
        counts <- counts[order(-counts[[y_col]]), ]
        cats   <- sort(unique(na.omit(counts$Model_system)))
        pal    <- setNames(.GRP_PAL[seq_along(cats)], cats)

        p <- plot_ly()
        for (ms in cats) {
          sub <- counts[!is.na(counts$Model_system) & counts$Model_system == ms, ]
          if (nrow(sub) == 0) next
          p <- p |> add_bars(
            x = sub$assay, y = sub[[y_col]],
            name          = ms,
            marker        = list(color = pal[[ms]]),
            hovertemplate = "<b>%{x}</b><br>%{y} DEGs<extra></extra>"
          )
        }
        p |> layout(
          barmode = "stack",
          xaxis   = list(title = "", tickangle = -55, automargin = TRUE,
                         showticklabels = nrow(counts) <= 60),
          yaxis   = list(title = y_lab),
          legend  = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.25),
          margin  = list(l = 5, r = 5, t = 5, b = 5)
        ) |> config(displayModeBar = FALSE)

      } else {
        # ── Grouped / aggregated bars ────────────────────────────────────
        # When direction is "all": show stacked up+down per group
        # Otherwise: one bar per group
        grp_col <- x_var

        if (dir == "all") {
          # Aggregate up and down separately then stack
          up_agg <- tapply(counts$n_up,   counts[[grp_col]], sum, na.rm = TRUE)
          dn_agg <- tapply(counts$n_down, counts[[grp_col]], sum, na.rm = TRUE)
          grps   <- sort(unique(na.omit(counts[[grp_col]])))
          agg    <- data.frame(
            group  = grps,
            n_up   = as.integer(up_agg[grps]),
            n_down = as.integer(dn_agg[grps]),
            stringsAsFactors = FALSE
          )
          agg$n_all <- agg$n_up + agg$n_down
          # Sort by total
          agg <- agg[order(-agg$n_all), ]

          plot_ly(agg) |>
            add_bars(x = ~group, y = ~n_up,   name = "Up",
                     marker = list(color = .DEG_CLRS$up),
                     hovertemplate = "<b>%{x}</b><br>Up DEGs: %{y}<extra></extra>") |>
            add_bars(x = ~group, y = ~n_down, name = "Down",
                     marker = list(color = .DEG_CLRS$down),
                     hovertemplate = "<b>%{x}</b><br>Down DEGs: %{y}<extra></extra>") |>
            layout(
              barmode = "stack",
              xaxis   = list(title = "", automargin = TRUE, tickangle = -30),
              yaxis   = list(title = "Total DEGs"),
              legend  = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.25),
              margin  = list(l = 5, r = 5, t = 5, b = 5)
            ) |> config(displayModeBar = FALSE)

        } else {
          agg_vals <- tapply(counts[[y_col]], counts[[grp_col]], sum, na.rm = TRUE)
          grps     <- sort(unique(na.omit(counts[[grp_col]])))
          agg      <- data.frame(
            group = grps,
            n     = as.integer(agg_vals[grps]),
            stringsAsFactors = FALSE
          )
          agg <- agg[order(-agg$n), ]

          plot_ly(agg, x = ~group, y = ~n, type = "bar",
                  marker = list(color = bar_colour),
                  hovertemplate = "<b>%{x}</b><br>%{y} DEGs<extra></extra>") |>
            layout(
              xaxis  = list(title = "", automargin = TRUE, tickangle = -30),
              yaxis  = list(title = y_lab),
              margin = list(l = 5, r = 5, t = 5, b = 5)
            ) |> config(displayModeBar = FALSE)
        }
      }
    })

    # ── Model System filter for overlap chart ─────────────────────────────────
    output$overlap_ms_ui <- renderUI({
      counts <- deg_counts_rv()
      req(counts)
      ms_opts <- c("All" = "all",
                   setNames(sort(unique(na.omit(counts$Model_system))),
                            sort(unique(na.omit(counts$Model_system)))))
      selectInput(session$ns("overlap_ms"), NULL,
                  choices = ms_opts, selected = "all", width = "200px")
    })

    # ── Filtered gene set for overlap chart ───────────────────────────────────
    overlap_genes_filtered <- reactive({
      genes <- deg_genes_rv()
      req(genes)

      # Direction filter
      dir <- input$overlap_dir
      if (dir == "up")   genes <- genes[genes$DEG == "up",   ]
      if (dir == "down") genes <- genes[genes$DEG == "down", ]

      # Model system filter
      ms <- input$overlap_ms
      if (!is.null(ms) && ms != "all")
        genes <- genes[!is.na(genes$Model_system) & genes$Model_system == ms, ]

      genes
    })

    # ── Chart 2: Most Recurrent DEGs ─────────────────────────────────────────
    output$overlap_bar <- renderPlotly({
      genes <- overlap_genes_filtered()
      req(genes, nrow(genes) > 0)

      top_n <- max(5L, as.integer(input$top_n))
      dir   <- input$overlap_dir

      if (dir == "all") {
        # Stacked up+down
        up_cnt   <- table(genes$symbol[genes$DEG == "up"])
        dn_cnt   <- table(genes$symbol[genes$DEG == "down"])
        all_syms <- union(names(up_cnt), names(dn_cnt))
        freq <- data.frame(
          gene   = all_syms,
          n_up   = as.integer(up_cnt[all_syms]),
          n_down = as.integer(dn_cnt[all_syms]),
          stringsAsFactors = FALSE
        )
        freq[is.na(freq$n_up),   "n_up"]   <- 0L
        freq[is.na(freq$n_down), "n_down"] <- 0L
        freq$n_total <- freq$n_up + freq$n_down
        freq <- head(freq[order(-freq$n_total), ], top_n)
        freq$gene <- factor(freq$gene, levels = rev(freq$gene))

        plot_ly(freq) |>
          add_bars(x = ~n_up,   y = ~gene, name = "Up",
                   orientation = "h",
                   marker = list(color = .DEG_CLRS$up),
                   hovertemplate = "<b>%{y}</b><br>Up in %{x} assays<extra></extra>") |>
          add_bars(x = ~n_down, y = ~gene, name = "Down",
                   orientation = "h",
                   marker = list(color = .DEG_CLRS$down),
                   hovertemplate = "<b>%{y}</b><br>Down in %{x} assays<extra></extra>") |>
          layout(
            barmode = "stack",
            xaxis   = list(title = "Number of assays"),
            yaxis   = list(title = "", automargin = TRUE),
            legend  = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.1),
            margin  = list(l = 0, r = 5, t = 5, b = 5)
          ) |> config(displayModeBar = FALSE)

      } else {
        freq <- as.data.frame(sort(table(genes$symbol), decreasing = TRUE))
        colnames(freq) <- c("gene", "n")
        freq <- head(freq, top_n)
        freq$gene <- factor(freq$gene, levels = rev(freq$gene))
        bar_col <- if (dir == "up") .DEG_CLRS$up else .DEG_CLRS$down

        plot_ly(freq, x = ~n, y = ~gene, type = "bar",
                orientation = "h",
                marker = list(color = bar_col),
                hovertemplate = "<b>%{y}</b><br>DE in %{x} assays<extra></extra>") |>
          layout(
            xaxis  = list(title = "Number of assays"),
            yaxis  = list(title = "", automargin = TRUE),
            margin = list(l = 0, r = 5, t = 5, b = 5)
          ) |> config(displayModeBar = FALSE)
      }
    })

    # ── Card 2b: UpSet plot ────────────────────────────────────────────────────

    # Render the plotOutput with a dynamic height that scales with set count
    output$upset_ui <- renderUI({
      req(deg_genes_rv())
      grp_var <- input$upset_grp
      top_n   <- max(2L, as.integer(input$upset_n))
      n_sets  <- min(top_n, length(unique(na.omit(deg_genes_rv()[[grp_var]]))))
      # ~40px per set plus padding for the bar chart panel above the dot matrix
      height_px <- max(380L, 180L + n_sets * 42L)
      plotOutput(session$ns("upset_plot"), height = paste0(height_px, "px"))
    })

    output$upset_plot <- renderPlot({
      genes <- deg_genes_rv()
      req(genes, nrow(genes) > 0)

      dir     <- input$upset_dir
      grp_var <- input$upset_grp
      top_n   <- max(2L, as.integer(input$upset_n))

      # Direction filter
      if (dir == "up")   genes <- genes[genes$DEG == "up",   ]
      if (dir == "down") genes <- genes[genes$DEG == "down", ]
      req(nrow(genes) > 0)

      # Rank groups by how many unique DEGs they contain; keep top_n
      valid <- genes[!is.na(genes[[grp_var]]) & !is.na(genes$symbol), ]
      req(nrow(valid) > 0)

      grp_sizes <- sort(
        tapply(valid$symbol, valid[[grp_var]], function(x) length(unique(x))),
        decreasing = TRUE
      )
      groups <- names(grp_sizes)[seq_len(min(top_n, length(grp_sizes)))]
      req(length(groups) >= 2)

      sub      <- valid[valid[[grp_var]] %in% groups, ]
      all_syms <- unique(sub$symbol)
      req(length(all_syms) > 0)

      # Build binary membership matrix (gene × group)
      # Column names: replace leading digits and non-alphanumeric chars so
      # UpSetR is happy, but keep labels as human-readable as possible.
      raw_names   <- groups
      clean_names <- gsub("[^A-Za-z0-9_]", "_", raw_names)
      clean_names <- ifelse(
        grepl("^[0-9]", clean_names),
        paste0("G", clean_names),
        clean_names
      )
      clean_names <- make.unique(clean_names)

      mat <- as.data.frame(
        matrix(0L, nrow = length(all_syms), ncol = length(groups)),
        stringsAsFactors = FALSE
      )
      colnames(mat) <- clean_names

      for (i in seq_along(groups)) {
        g_syms <- unique(sub$symbol[sub[[grp_var]] == groups[i]])
        mat[all_syms %in% g_syms, clean_names[i]] <- 1L
      }

      # Drop genes not in any group (should not happen but be safe)
      mat <- mat[rowSums(mat) > 0, , drop = FALSE]
      req(nrow(mat) > 0)

      # Colour palette consistent with the rest of the module
      bar_col <- switch(dir,
        all  = .DEG_CLRS$all,
        up   = .DEG_CLRS$up,
        down = .DEG_CLRS$down
      )

      UpSetR::upset(
        mat,
        sets             = rev(clean_names),   # bottom → top ordering
        order.by         = "freq",
        decreasing       = TRUE,
        mb.ratio         = c(0.60, 0.40),
        text.scale       = c(1.4, 1.2, 1.1, 1.0, 1.3, 0.9),
        # c(intersection size title, intersection size ticks,
        #   set size title, set size ticks, set names, numbers above bars)
        point.size       = 3.2,
        line.size        = 0.9,
        sets.bar.color   = bar_col,
        main.bar.color   = bar_col,
        sets.x.label     = "DEGs in group",
        mainbar.y.label  = "Intersection size (genes)",
        show.numbers     = "yes",
        number.angles    = 0,
        set_size.show    = TRUE,
        set_size.numbers_size = 3.5
      )
    }, bg = "white")

    # ── Card 3: Model System filter UI ───────────────────────────────────────
    output$ov_ms_ui <- renderUI({
      counts <- deg_counts_rv()
      req(counts)
      ms_opts <- c("All" = "all",
                   setNames(sort(unique(na.omit(counts$Model_system))),
                            sort(unique(na.omit(counts$Model_system)))))
      selectInput(session$ns("ov_ms"), NULL,
                  choices = ms_opts, selected = "all", width = "200px")
    })

    # ── Card 3: Build the DEG × KO-target overlap matrix ─────────────────────
    # Computed once per load (cached in a reactive), then filtered by UI inputs.
    # Structure: each row = a DEG symbol; each KO-target gene gets one column
    # containing "" / "↑" / "↓" / "↑↓" depending on direction in that target's assays.
    ov_matrix_base <- reactive({
      genes <- deg_genes_rv()
      req(genes)

      # Collapse to: symbol × KO_target × direction
      # For each (symbol, KO_target) pair, collect all directions across assays
      genes_valid <- genes[!is.na(genes$symbol) & genes$symbol != "" &
                             !is.na(genes$Gene), ]

      ko_targets <- sort(unique(genes_valid$Gene))
      all_syms   <- sort(unique(genes_valid$symbol))

      # Build a named list: ko_target -> list(up_syms, down_syms, lfc_vec)
      has_lfc <- "log2fc" %in% names(genes_valid)

      ko_sets <- lapply(ko_targets, function(ko) {
        sub <- genes_valid[genes_valid$Gene == ko, ]
        lfc_vec <- if (has_lfc && any(!is.na(sub$log2fc))) {
          tapply(sub$log2fc, sub$symbol, mean, na.rm = TRUE)
        } else setNames(numeric(0), character(0))
        list(
          up   = unique(sub$symbol[sub$DEG == "up"]),
          down = unique(sub$symbol[sub$DEG == "down"]),
          lfc  = lfc_vec
        )
      })
      names(ko_sets) <- ko_targets

      # For each symbol build a row: one cell per KO-target.
      # Cell format: "↑ +2.31", "↓ -1.45", "↑↓", or "".
      mat <- data.frame(
        DEG_Gene = all_syms,
        stringsAsFactors = FALSE
      )
      for (ko in ko_targets) {
        s     <- ko_sets[[ko]]
        is_up <- all_syms %in% s$up
        is_dn <- all_syms %in% s$down
        lfc_v <- s$lfc
        mat[[ko]] <- mapply(function(sym, up, dn) {
          lfc <- if (length(lfc_v) > 0) lfc_v[sym] else NA_real_
          .fmt_dir_cell(up, dn, if (length(lfc) == 0) NA_real_ else lfc)
        }, all_syms, is_up, is_dn, SIMPLIFY = TRUE, USE.NAMES = FALSE)
      }

      # N = how many KO-targets this DEG appears in (any direction)
      mat$N_KO_Targets <- rowSums(mat[, ko_targets, drop = FALSE] != "")
      # Re-order: DEG_Gene, N_KO_Targets, then KO columns
      mat <- mat[, c("DEG_Gene", "N_KO_Targets", ko_targets)]
      mat
    })

    # ── Card 3: Filtered view ─────────────────────────────────────────────────
    ov_matrix_filtered <- reactive({
      mat    <- ov_matrix_base()
      req(mat)
      genes  <- deg_genes_rv()
      dir    <- input$ov_dir
      min_ko <- max(1L, as.integer(input$ov_min_ko))
      ms     <- input$ov_ms

      ko_targets <- names(mat)[-(1:2)]   # all KO-target columns

      # Model System filter: restrict to KO-targets belonging to that model system
      if (!is.null(ms) && ms != "all") {
        counts <- deg_counts_rv()
        req(counts)
        allowed_ko <- unique(counts$Gene[!is.na(counts$Model_system) &
                                           counts$Model_system == ms])
        ko_keep    <- ko_targets[ko_targets %in% allowed_ko]
        if (length(ko_keep) == 0) return(data.frame())
        mat <- mat[, c("DEG_Gene", "N_KO_Targets", ko_keep), drop = FALSE]
        ko_targets <- ko_keep
        # Recompute N for remaining columns
        mat$N_KO_Targets <- rowSums(mat[, ko_targets, drop = FALSE] != "")
      }

      # Direction filter: blank out cells that don't match.
      # Cells now have format "↑ +2.31" / "↓ -1.45" / "↑↓" so we match by
      # first character rather than exact equality.
      if (dir == "up") {
        for (ko in ko_targets) {
          v    <- mat[[ko]]
          both <- v == "\u2191\u2193"
          mat[[ko]] <- ifelse(both, "\u2191",                           # ↑↓ → ↑ (no single logFC)
                       ifelse(startsWith(v, "\u2191") & !both, v, "")) # keep "↑ +2.31" as-is
        }
        mat$N_KO_Targets <- rowSums(mat[, ko_targets, drop = FALSE] != "")
      } else if (dir == "down") {
        for (ko in ko_targets) {
          v    <- mat[[ko]]
          both <- v == "\u2191\u2193"
          mat[[ko]] <- ifelse(both, "\u2193",
                       ifelse(startsWith(v, "\u2193") & !both, v, ""))
        }
        mat$N_KO_Targets <- rowSums(mat[, ko_targets, drop = FALSE] != "")
      }

      # Minimum KO-target threshold
      mat <- mat[mat$N_KO_Targets >= min_ko, ]
      # Sort by most shared
      mat[order(-mat$N_KO_Targets), ]
    })

    # ── Card 3: Render table ──────────────────────────────────────────────────
    output$ov_table <- DT::renderDT({
      mat <- ov_matrix_filtered()
      req(mat, nrow(mat) > 0)

      ko_cols <- names(mat)[-(1:2)]

      DT::datatable(
        mat,
        rownames   = FALSE,
        filter     = "top",
        extensions = "Buttons",
        options    = list(
          pageLength     = 25,
          scrollX        = TRUE,
          scrollY        = "520px",
          scrollCollapse = TRUE,
          dom            = "Bfrtip",
          buttons        = list("csv", "excel"),
          fixedColumns   = list(leftColumns = 2),
          columnDefs     = list(
            list(className = "dt-center", targets = seq_along(ko_cols)),
            list(width = "80px", targets = seq_along(ko_cols))
          ),
          rowCallback = .LFC_ROW_CB
        ),
        class = "display compact cell-border"
      ) |>
        DT::formatStyle(
          "N_KO_Targets",
          background         = DT::styleColorBar(range(mat$N_KO_Targets), "#d4e6f1"),
          backgroundSize     = "100% 80%",
          backgroundRepeat   = "no-repeat",
          backgroundPosition = "center"
        )
    })

    # ── Table: genes in current overlap selection ─────────────────────────────
    output$overlap_tbl <- DT::renderDT({
      genes <- overlap_genes_filtered()
      req(genes, nrow(genes) > 0)

      top_n <- max(5L, as.integer(input$top_n))
      dir   <- input$overlap_dir

      has_lfc <- "log2fc" %in% names(genes)

      if (dir == "all") {
        up_cnt   <- table(genes$symbol[genes$DEG == "up"])
        dn_cnt   <- table(genes$symbol[genes$DEG == "down"])
        all_syms <- union(names(up_cnt), names(dn_cnt))
        tbl <- data.frame(
          Gene        = all_syms,
          Up_assays   = as.integer(up_cnt[all_syms]),
          Down_assays = as.integer(dn_cnt[all_syms]),
          stringsAsFactors = FALSE
        )
        tbl[is.na(tbl)] <- 0L
        tbl$Total <- tbl$Up_assays + tbl$Down_assays
        tbl <- head(tbl[order(-tbl$Total), ], top_n)
        if (has_lfc) {
          up_lfc <- tapply(genes$log2fc[genes$DEG == "up"],
                           genes$symbol[genes$DEG == "up"],
                           function(x) round(median(abs(x), na.rm=TRUE), 2))
          dn_lfc <- tapply(genes$log2fc[genes$DEG == "down"],
                           genes$symbol[genes$DEG == "down"],
                           function(x) round(median(abs(x), na.rm=TRUE), 2))
          tbl[["Median Up |logFC|"]]   <- as.numeric(up_lfc[tbl$Gene])
          tbl[["Median Down |logFC|"]] <- as.numeric(dn_lfc[tbl$Gene])
        }
      } else {
        freq <- sort(table(genes$symbol), decreasing = TRUE)
        tbl  <- data.frame(Gene = names(freq), N_assays = as.integer(freq),
                           stringsAsFactors = FALSE)
        tbl  <- head(tbl, top_n)
        if (has_lfc) {
          lfc_med <- tapply(genes$log2fc, genes$symbol,
                            function(x) round(median(abs(x), na.rm=TRUE), 2))
          tbl[["Median |logFC|"]] <- as.numeric(lfc_med[tbl$Gene])
        }
      }

      lfc_cols <- grep("logFC", names(tbl), value = TRUE)
      bar_col  <- if (dir == "down") "#b8d8ef" else "#f8c1c4"

      dt <- DT::datatable(
        tbl, rownames = FALSE,
        options = list(pageLength = 15, dom = "tp", scrollY = "380px",
                       scrollCollapse = TRUE),
        class = "display compact"
      )
      for (col in lfc_cols) {
        rng <- range(tbl[[col]], na.rm = TRUE)
        dt  <- dt |> DT::formatStyle(
          col,
          background         = DT::styleColorBar(rng, bar_col),
          backgroundSize     = "100% 70%",
          backgroundRepeat   = "no-repeat",
          backgroundPosition = "center"
        )
      }
      dt
    })

    # =========================================================================
    # Card 4: Assay-level DEG comparison
    # =========================================================================

    # ── Assay picker UI (grouped by KO Gene, labelled with metadata) ─────────
    output$cmp_assay_ui <- renderUI({
      counts <- deg_counts_rv()
      req(counts)

      # Build a human-readable label: "ModelSystem · Strategy [DPC]"
      # Value = raw assay name (used internally); label = readable metadata.
      # Append a numeric suffix via make.unique() for the rare case where two
      # assays share the same Gene + Model + Strategy + DPC combination.
      make_lbl <- function(ms, ko, dpc) {
        base <- paste0(
          ifelse(!is.na(ms),  ms,  "?"), " \u00b7 ",
          ifelse(!is.na(ko),  ko,  "?")
        )
        if (!is.na(dpc) && nzchar(as.character(dpc)))
          paste0(base, " [", dpc, "]")
        else
          base
      }

      gene_order <- sort(unique(na.omit(counts$Gene)))
      choices_grouped <- lapply(gene_order, function(g) {
        sub  <- counts[!is.na(counts$Gene) & counts$Gene == g, ]
        sub  <- sub[order(sub$Model_system, sub$KO_strat, sub$assay), ]
        lbls <- mapply(make_lbl, sub$Model_system, sub$KO_strat, sub$DPC,
                       SIMPLIFY = TRUE)
        lbls <- make.unique(lbls, sep = " #")   # ensure labels unique within group
        setNames(sub$assay, lbls)
      })
      names(choices_grouped) <- gene_order

      tagList(
        selectizeInput(
          session$ns("cmp_assays"),
          label    = "Select assays to compare (2 or more):",
          choices  = choices_grouped,
          selected = NULL,
          multiple = TRUE,
          options  = list(
            placeholder = "Search by gene, model system, or KO strategy\u2026",
            plugins     = list("remove_button"),
            maxItems    = NULL
          )
        ),
        tags$p(class = "text-muted small mt-1",
               icon("circle-info"),
               " Groups = KO Gene. Labels show Model System \u00b7 Strategy [DPC].",
               " Select at least 2 assays.")
      )
    })

    # ── Helper: extract DEG gene sets for selected assays + direction ─────────
    cmp_sets <- reactive({
      genes   <- deg_genes_rv()
      req(genes)
      assays  <- input$cmp_assays
      req(assays, length(assays) >= 2)
      dir     <- input$cmp_dir

      sub <- genes[genes$assay %in% assays, ]
      if (dir == "up")   sub <- sub[sub$DEG == "up",   ]
      if (dir == "down") sub <- sub[sub$DEG == "down", ]

      # Per-assay gene sets (unique symbols)
      lapply(setNames(assays, assays), function(a) {
        unique(na.omit(sub$symbol[sub$assay == a]))
      })
    })

    # ── Summary bar chart: total + shared-with-all breakdown ─────────────────
    output$cmp_summary_ui <- renderUI({
      req(cmp_sets())
      plotlyOutput(session$ns("cmp_bar"), height = "260px")
    })

    output$cmp_bar <- renderPlotly({
      sets   <- cmp_sets()
      req(sets)
      assays <- names(sets)
      dir    <- input$cmp_dir

      # Core shared across ALL selected assays
      core      <- Reduce(intersect, sets)
      n_core    <- length(core)

      # Per assay: unique to that assay only (not in any other)
      bar_colour <- switch(dir, all = .DEG_CLRS$all,
                           up = .DEG_CLRS$up, down = .DEG_CLRS$down)

      rows <- lapply(assays, function(a) {
        others  <- sets[assays != a]
        in_any  <- length(intersect(sets[[a]], unique(unlist(others))))
        unique_n <- length(sets[[a]]) - in_any
        data.frame(
          assay    = a,
          total    = length(sets[[a]]),
          shared   = in_any,
          unique_n = unique_n,
          core     = n_core,
          stringsAsFactors = FALSE
        )
      })
      df <- do.call(rbind, rows)
      df <- df[order(-df$total), ]
      # Shorten assay labels for display (last 40 chars)
      df$label <- ifelse(nchar(df$assay) > 45,
                         paste0("…", substr(df$assay, nchar(df$assay)-44, nchar(df$assay))),
                         df$assay)
      df$label <- factor(df$label, levels = rev(df$label))

      plot_ly(df, y = ~label) |>
        add_bars(x = ~unique_n, name = "Unique to assay",
                 orientation = "h",
                 marker = list(color = bar_colour, opacity = 0.45),
                 hovertemplate = "<b>%{y}</b><br>Unique: %{x}<extra></extra>") |>
        add_bars(x = ~shared, name = "Shared with ≥1 other",
                 orientation = "h",
                 marker = list(color = bar_colour, opacity = 0.85),
                 hovertemplate = "<b>%{y}</b><br>Shared w/ ≥1: %{x}<extra></extra>") |>
        layout(
          barmode = "stack",
          xaxis   = list(title = "N DEGs"),
          yaxis   = list(title = "", automargin = TRUE),
          legend  = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.18),
          margin  = list(l = 5, r = 5, t = 30, b = 5),
          title   = list(
            text = sprintf("<b>Core (shared across ALL %d assays): %d gene%s</b>",
                           length(assays), n_core, if (n_core == 1) "" else "s"),
            x = 0, font = list(size = 13)
          )
        ) |>
        config(displayModeBar = FALSE)
    })

    # ── Tabs: shared-all genes / per-assay unique genes / full matrix ─────────
    output$cmp_tabs_ui <- renderUI({
      req(cmp_sets(), length(cmp_sets()) >= 2)
      navset_tab(
        nav_panel("Shared across ALL assays",
                  div(class = "p-2", DTOutput(session$ns("cmp_core_tbl")))),
        nav_panel("Unique per assay",
                  div(class = "p-2", DTOutput(session$ns("cmp_unique_tbl")))),
        nav_panel("Full gene matrix",
                  div(class = "p-2",
                      tags$p(class = "text-muted small",
                             "Rows = DEG genes present in ≥1 selected assay. ",
                             "↑ up, ↓ down, ↑↓ both (if gene appears multiple times). ",
                             "Sorted by number of assays gene appears in."),
                      DTOutput(session$ns("cmp_matrix_tbl"))))
      )
    })

    # ── Core (shared-ALL) gene table ──────────────────────────────────────────
    output$cmp_core_tbl <- DT::renderDT({
      sets  <- cmp_sets()
      req(sets, length(sets) >= 2)
      genes <- deg_genes_rv()
      core  <- Reduce(intersect, sets)

      if (length(core) == 0)
        return(DT::datatable(data.frame(Message = "No genes shared across ALL selected assays."),
                             rownames = FALSE, options = list(dom = "t")))

      assays <- names(sets)
      dir    <- input$cmp_dir
      sub    <- genes[genes$assay %in% assays &
                        genes$symbol %in% core, ]
      if (dir == "up")   sub <- sub[sub$DEG == "up",   ]
      if (dir == "down") sub <- sub[sub$DEG == "down", ]

      # One row per gene: show direction + logFC per assay
      tbl <- data.frame(Gene = sort(unique(core)), stringsAsFactors = FALSE)
      for (a in assays) {
        a_sub  <- sub[sub$assay == a, ]
        is_up  <- tbl$Gene %in% a_sub$symbol[a_sub$DEG == "up"]
        is_dn  <- tbl$Gene %in% a_sub$symbol[a_sub$DEG == "down"]
        lbl    <- ifelse(nchar(a) > 40, paste0("\u2026", substr(a, nchar(a)-39, nchar(a))), a)
        lfc_v  <- if ("log2fc" %in% names(a_sub))
          setNames(a_sub$log2fc, a_sub$symbol) else setNames(numeric(0), character(0))
        tbl[[lbl]] <- mapply(function(g, up, dn) {
          lfc <- lfc_v[g]; .fmt_dir_cell(up, dn, if (length(lfc)==0) NA_real_ else lfc)
        }, tbl$Gene, is_up, is_dn, SIMPLIFY = TRUE, USE.NAMES = FALSE)
      }

      ko_cols <- names(tbl)[-1]
      DT::datatable(tbl, rownames = FALSE, filter = "top",
                    extensions = "Buttons",
                    options = list(pageLength = 25, scrollX = TRUE,
                                   dom = "Bfrtip",
                                   buttons = list("csv", "excel"),
                                   columnDefs = list(list(width = "90px",
                                                          targets = seq_along(ko_cols))),
                                   rowCallback = .LFC_ROW_CB),
                    class = "display compact cell-border")
    })

    # ── Unique-per-assay gene table ───────────────────────────────────────────
    output$cmp_unique_tbl <- DT::renderDT({
      sets  <- cmp_sets()
      req(sets, length(sets) >= 2)
      assays <- names(sets)

      genes_df <- deg_genes_rv()
      rows <- lapply(assays, function(a) {
        others     <- unique(unlist(sets[assays != a]))
        uniq_genes <- setdiff(sets[[a]], others)
        if (length(uniq_genes) == 0) return(NULL)
        lbl    <- ifelse(nchar(a) > 55, paste0("\u2026", substr(a, nchar(a)-54, nchar(a))), a)
        sorted <- sort(uniq_genes)
        # Attach direction and logFC from the loaded gene table
        a_sub  <- genes_df[genes_df$assay == a & genes_df$symbol %in% sorted, ]
        dir_map <- setNames(a_sub$DEG,     a_sub$symbol)
        lfc_map <- if ("log2fc" %in% names(a_sub))
          setNames(a_sub$log2fc, a_sub$symbol) else setNames(numeric(0), character(0))
        dir_lbl <- ifelse(dir_map[sorted] == "up", "\u2191 up", "\u2193 down")
        lfc_val <- round(as.numeric(lfc_map[sorted]), 3)
        data.frame(Assay = lbl, Gene = sorted,
                   Direction = dir_lbl, logFC = lfc_val,
                   stringsAsFactors = FALSE)
      })
      tbl <- do.call(rbind, rows)

      if (is.null(tbl) || nrow(tbl) == 0)
        return(DT::datatable(data.frame(Message = "No genes are unique to a single assay."),
                             rownames = FALSE, options = list(dom = "t")))

      has_lfc_col <- "logFC" %in% names(tbl) && any(!is.na(tbl$logFC))
      lfc_rng     <- if (has_lfc_col) range(abs(tbl$logFC), na.rm = TRUE) else c(0,1)

      dt <- DT::datatable(tbl, rownames = FALSE, filter = "top",
                    extensions = "Buttons",
                    options = list(pageLength = 25, scrollX = TRUE,
                                   dom = "Bfrtip",
                                   buttons = list("csv", "excel"),
                                   rowCallback = .LFC_ROW_CB),
                    class = "display compact")
      if (has_lfc_col)
        dt <- dt |> DT::formatStyle(
          "logFC",
          background         = DT::styleColorBar(lfc_rng, "#e0e0e0"),
          backgroundSize     = "100% 70%",
          backgroundRepeat   = "no-repeat",
          backgroundPosition = "center"
        )
      dt
    })

    # ── Full gene matrix for selected assays ──────────────────────────────────
    output$cmp_matrix_tbl <- DT::renderDT({
      sets  <- cmp_sets()
      req(sets, length(sets) >= 2)
      genes_df <- deg_genes_rv()
      assays   <- names(sets)
      dir      <- input$cmp_dir

      sub <- genes_df[genes_df$assay %in% assays, ]
      if (dir == "up")   sub <- sub[sub$DEG == "up",   ]
      if (dir == "down") sub <- sub[sub$DEG == "down", ]

      all_syms <- sort(unique(na.omit(sub$symbol)))
      if (length(all_syms) == 0)
        return(DT::datatable(data.frame(Message = "No DEGs found."),
                             rownames = FALSE, options = list(dom = "t")))

      mat <- data.frame(Gene = all_syms, stringsAsFactors = FALSE)
      for (a in assays) {
        a_sub  <- sub[sub$assay == a, ]
        is_up  <- all_syms %in% a_sub$symbol[a_sub$DEG == "up"]
        is_dn  <- all_syms %in% a_sub$symbol[a_sub$DEG == "down"]
        lbl    <- ifelse(nchar(a) > 40, paste0("\u2026", substr(a, nchar(a)-39, nchar(a))), a)
        lfc_v  <- if ("log2fc" %in% names(a_sub))
          setNames(a_sub$log2fc, a_sub$symbol) else setNames(numeric(0), character(0))
        mat[[lbl]] <- mapply(function(g, up, dn) {
          lfc <- lfc_v[g]; .fmt_dir_cell(up, dn, if (length(lfc)==0) NA_real_ else lfc)
        }, all_syms, is_up, is_dn, SIMPLIFY = TRUE, USE.NAMES = FALSE)
      }

      ko_cols <- names(mat)[-1]
      mat$N_Assays <- rowSums(mat[, ko_cols, drop = FALSE] != "")
      mat <- mat[, c("Gene", "N_Assays", ko_cols)]
      mat <- mat[order(-mat$N_Assays), ]

      DT::datatable(mat, rownames = FALSE, filter = "top",
                    extensions = "Buttons",
                    options = list(pageLength = 25, scrollX = TRUE,
                                   scrollY = "480px", scrollCollapse = TRUE,
                                   dom = "Bfrtip",
                                   buttons = list("csv", "excel"),
                                   fixedColumns = list(leftColumns = 2),
                                   columnDefs = list(
                                     list(className = "dt-center",
                                          targets = seq_along(ko_cols) + 1),
                                     list(width = "90px",
                                          targets = seq_along(ko_cols) + 1)),
                                   rowCallback = .LFC_ROW_CB),
                    class = "display compact cell-border") |>
        DT::formatStyle("N_Assays",
                        background         = DT::styleColorBar(c(0, length(assays)), "#d4e6f1"),
                        backgroundSize     = "100% 80%",
                        backgroundRepeat   = "no-repeat",
                        backgroundPosition = "center")
    })

  })
}
