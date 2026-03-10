# =============================================================================
# Module: study_info_dashboard
# Description: Interactive overview of bulk RNA-seq study coverage.
#   Shows N unique assays & genes broken down by Model System, KO Strategy,
#   and DPC. Charts are cross-linked — clicking a bar in any chart filters
#   the counts shown in the other two.
#   Each card also contains a gene-list tab showing which genes belong to
#   each category visible in the bar chart.
# =============================================================================

library(plotly)
library(DT)
library(UpSetR)

# Colour palette shared across all charts
.CLRS <- list(assays = "#4E79A7", genes = "#F28E2B")

# ── UI ────────────────────────────────────────────────────────────────────────

study_info_dashboardUI <- function(id) {
  ns <- NS(id)

  tagList(

    # ── Value boxes ──────────────────────────────────────────────────────────
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      value_box(
        title    = "Total Assays",
        value    = textOutput(ns("vb_assays"), inline = TRUE),
        showcase = bsicons::bs_icon("bar-chart-fill"),
        theme    = "primary"
      ),
      value_box(
        title    = "Unique Genes",
        value    = textOutput(ns("vb_genes"), inline = TRUE),
        showcase = bsicons::bs_icon("diagram-3"),
        theme    = value_box_theme(bg = .CLRS$genes, fg = "#fff")
      ),
      value_box(
        title    = "Model Systems",
        value    = textOutput(ns("vb_models"), inline = TRUE),
        showcase = bsicons::bs_icon("grid-3x3-gap-fill"),
        theme    = "secondary"
      ),
      value_box(
        title    = "KO Strategies",
        value    = textOutput(ns("vb_ko"), inline = TRUE),
        showcase = bsicons::bs_icon("scissors"),
        theme    = "dark"
      )
    ),

    # ── Chart row 1: Model System (wide) + KO Strategy ───────────────────────
    layout_columns(
      col_widths = c(8, 4),

      card(
        card_header(
          class = "d-flex justify-content-between align-items-center",
          "Assays & Genes by Model System",
          actionLink(ns("clear_model"), "Clear", class = "small text-muted")
        ),
        navset_tab(
          nav_panel("Chart",
                    plotlyOutput(ns("model_plot"), height = "300px")),
          nav_panel("Genes",
                    div(class = "p-2",
                        uiOutput(ns("model_tbl_header")),
                        DTOutput(ns("model_gene_tbl"))))
        ),
        full_screen = TRUE
      ),

      card(
        card_header(
          class = "d-flex justify-content-between align-items-center",
          "Assays & Genes by KO Strategy",
          actionLink(ns("clear_ko"), "Clear", class = "small text-muted")
        ),
        navset_tab(
          nav_panel("Chart",
                    plotlyOutput(ns("ko_plot"), height = "300px")),
          nav_panel("Genes",
                    div(class = "p-2",
                        uiOutput(ns("ko_tbl_header")),
                        DTOutput(ns("ko_gene_tbl"))))
        ),
        full_screen = TRUE
      )
    ),

    # ── Chart row 2: DPC + selection summary ─────────────────────────────────
    layout_columns(
      col_widths = c(6, 6),

      card(
        card_header(
          class = "d-flex justify-content-between align-items-center",
          "Assays & Genes by Data Production Centre (DPC)",
          actionLink(ns("clear_dpc"), "Clear", class = "small text-muted")
        ),
        navset_tab(
          nav_panel("Chart",
                    plotlyOutput(ns("dpc_plot"), height = "230px")),
          nav_panel("Genes",
                    div(class = "p-2",
                        uiOutput(ns("dpc_tbl_header")),
                        DTOutput(ns("dpc_gene_tbl"))))
        ),
        full_screen = TRUE
      ),

      card(
        card_header("Active Cross-filter Selection"),
        uiOutput(ns("selection_summary")),
        actionButton(
          ns("clear_all"), "Clear All Selections",
          class = "btn-sm btn-outline-secondary mt-2"
        )
      )
    ),

    # ── UpSet: KO gene coverage across study groupings ────────────────────────
    card(
      card_header(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
        tagList(bsicons::bs_icon("intersect"), " KO Gene Coverage — UpSet Plot"),
        div(
          class = "d-flex gap-2 align-items-center flex-wrap",
          div(
            class = "d-flex align-items-center gap-1",
            tags$label("Group by:", class = "mb-0 small fw-semibold"),
            selectInput(
              ns("si_upset_grp"), NULL,
              choices  = c("Model System" = "Model_system",
                           "KO Strategy"  = "KO_strat",
                           "DPC"          = "DPC"),
              selected = "Model_system",
              width    = "160px"
            )
          )
        )
      ),
      tags$p(
        class = "text-muted small px-3 pt-2 mb-0",
        "Each vertical bar = KO genes studied in exactly that combination of groups. ",
        "Horizontal bars = total KO genes per group. ",
        tags$b("Responds to cross-filter selections above.")
      ),
      uiOutput(ns("si_upset_ui")),
      full_screen = TRUE
    ),

    # ── KO Gene Study Coverage panel ──────────────────────────────────────────
    card(
      card_header(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
        tagList(bsicons::bs_icon("card-list"), " KO Gene Study Coverage"),
        tags$span(
          class = "text-muted small fw-normal",
          "All studies, conditions, model systems and KO strategies per KO gene ",
          tags$b("— responds to cross-filter")
        )
      ),
      tags$p(
        class = "text-muted small px-3 pt-2 mb-1",
        "Rows are grouped by KO gene. Use the column filters or search box to narrow down.  ",
        "Click a column header to sort. Use ",
        tags$kbd("Shift+click"), " to multi-sort."
      ),
      DTOutput(ns("gene_panel_tbl")),
      full_screen = TRUE
    ),

    # ── Assay Coverage Matrix ──────────────────────────────────────────────────
    card(
      card_header(
        class = "d-flex justify-content-between align-items-center flex-wrap gap-2",
        tagList(bsicons::bs_icon("grid-3x3-gap-fill"), " Assay Coverage Matrix"),
        div(
          class = "d-flex align-items-center gap-2",
          tags$label("X axis:", class = "mb-0 small fw-semibold"),
          selectInput(
            ns("matrix_x"), NULL,
            choices  = c("Model System" = "Model_system",
                         "KO Strategy"  = "KO_strat",
                         "DPC"          = "DPC"),
            selected = "Model_system",
            width    = "160px"
          )
        )
      ),
      tags$p(
        class = "text-muted small px-3 pt-2 mb-1",
        "Colour intensity = number of assays. Hover a cell to see the assay list."
      ),
      uiOutput(ns("matrix_ui")),
      full_screen = TRUE
    )
  )
}


# ── Server ────────────────────────────────────────────────────────────────────

study_info_dashboardServer <- function(id, study_info_r) {
  moduleServer(id, function(input, output, session) {

    # ── Cross-filter state ────────────────────────────────────────────────────
    sel <- reactiveValues(
      model_system = NULL,
      ko_strat     = NULL,
      dpc          = NULL
    )

    # ── Helpers ───────────────────────────────────────────────────────────────

    # Apply cross-filter selections; `exclude` skips its own slot
    filtered_data <- function(exclude = character(0)) {
      df <- study_info_r()
      if (!"model_system" %in% exclude && !is.null(sel$model_system))
        df <- df[df$Model_system %in% sel$model_system, ]
      if (!"ko_strat" %in% exclude && !is.null(sel$ko_strat))
        df <- df[df$KO_strat %in% sel$ko_strat, ]
      if (!"dpc" %in% exclude && !is.null(sel$dpc))
        df <- df[df$DPC %in% sel$dpc, ]
      df
    }

    # Aggregate counts per group
    agg_counts <- function(df, group_col) {
      if (nrow(df) == 0)
        return(data.frame(group = character(0),
                          n_assays = integer(0), n_genes = integer(0)))
      groups   <- df[[group_col]]
      all_grp  <- sort(unique(na.omit(groups)))
      assays_n <- tapply(df$Assay, groups, function(x) length(unique(na.omit(x))))
      genes_n  <- tapply(df$Gene,  groups, function(x) length(unique(na.omit(x))))
      data.frame(
        group    = all_grp,
        n_assays = as.integer(assays_n[all_grp]),
        n_genes  = as.integer(genes_n[all_grp]),
        stringsAsFactors = FALSE
      )
    }

    # Build gene table per group category for a given grouping column.
    # N_Assays is the count of assays that gene appears in *within* that
    # category — not the total assays for the whole category.
    # Returns data.frame: Category, Gene, N_Assays (per gene within category)
    gene_table_for <- function(df, group_col) {
      if (nrow(df) == 0) return(data.frame())
      groups <- sort(unique(na.omit(df[[group_col]])))
      do.call(rbind, lapply(groups, function(g) {
        sub   <- df[!is.na(df[[group_col]]) & df[[group_col]] == g, ]
        genes <- sort(unique(na.omit(sub$Gene)))
        if (length(genes) == 0) return(NULL)
        # Per-gene assay count within this group
        per_gene_n <- tapply(sub$Assay, sub$Gene,
                             function(x) length(unique(na.omit(x))))
        data.frame(
          Category = g,
          Gene     = genes,
          N_Assays = as.integer(per_gene_n[genes]),
          stringsAsFactors = FALSE
        )
      }))
    }

    # Build a plotly grouped-bar chart
    make_bar_chart <- function(agg, selected_val = NULL,
                               orientation = "v", source_id = NULL) {
      if (nrow(agg) == 0) return(plotly_empty())

      opacity_v <- if (is.null(selected_val)) rep(0.9, nrow(agg)) else
        ifelse(agg$group == selected_val, 1.0, 0.3)

      if (orientation == "h") {
        p <- plot_ly(source = source_id) |>
          add_bars(data = agg, x = ~n_assays, y = ~group,
                   name = "Assays", orientation = "h",
                   marker = list(color = .CLRS$assays, opacity = opacity_v),
                   hovertemplate = "<b>%{y}</b><br>Assays: %{x}<extra></extra>",
                   customdata = ~group) |>
          add_bars(data = agg, x = ~n_genes, y = ~group,
                   name = "Genes", orientation = "h",
                   marker = list(color = .CLRS$genes, opacity = opacity_v),
                   hovertemplate = "<b>%{y}</b><br>Genes: %{x}<extra></extra>",
                   customdata = ~group) |>
          layout(barmode = "group",
                 xaxis = list(title = "Count", fixedrange = TRUE),
                 yaxis = list(title = "", automargin = TRUE,
                              categoryorder = "total ascending"),
                 legend = list(orientation = "h", x = 0.5,
                               xanchor = "center", y = -0.15),
                 margin = list(l = 0, r = 10, t = 10, b = 10))
      } else {
        p <- plot_ly(source = source_id) |>
          add_bars(data = agg, x = ~group, y = ~n_assays,
                   name = "Assays",
                   marker = list(color = .CLRS$assays, opacity = opacity_v),
                   hovertemplate = "<b>%{x}</b><br>Assays: %{y}<extra></extra>",
                   customdata = ~group) |>
          add_bars(data = agg, x = ~group, y = ~n_genes,
                   name = "Genes",
                   marker = list(color = .CLRS$genes, opacity = opacity_v),
                   hovertemplate = "<b>%{x}</b><br>Genes: %{y}<extra></extra>",
                   customdata = ~group) |>
          layout(barmode = "group",
                 xaxis = list(title = "", automargin = TRUE),
                 yaxis = list(title = "Count", fixedrange = TRUE),
                 legend = list(orientation = "h", x = 0.5,
                               xanchor = "center", y = -0.2),
                 margin = list(l = 0, r = 10, t = 10, b = 10))
      }
      p |> config(displayModeBar = FALSE)
    }

    # Render a gene DT nicely
    render_gene_dt <- function(tbl) {
      if (is.null(tbl) || nrow(tbl) == 0)
        return(DT::datatable(data.frame(Message = "No data"),
                             rownames = FALSE, options = list(dom = "t")))
      DT::datatable(
        tbl, rownames = FALSE, filter = "top",
        options = list(pageLength = 15, scrollY = "320px",
                       scrollCollapse = TRUE, dom = "ftip"),
        class = "display compact"
      )
    }

    # ── Aggregated reactives (each excludes its own filter) ───────────────────

    model_agg <- reactive({ agg_counts(filtered_data("model_system"), "Model_system") })
    ko_agg    <- reactive({ agg_counts(filtered_data("ko_strat"),     "KO_strat")     })
    dpc_agg   <- reactive({ agg_counts(filtered_data("dpc"),          "DPC")          })
    full_filt <- reactive({ filtered_data() })

    # Gene tables (use full filter so tables reflect all active selections)
    model_genes <- reactive({ gene_table_for(full_filt(), "Model_system") })
    ko_genes    <- reactive({ gene_table_for(full_filt(), "KO_strat")     })
    dpc_genes   <- reactive({ gene_table_for(full_filt(), "DPC")          })

    # ── Value boxes ───────────────────────────────────────────────────────────

    output$vb_assays <- renderText({ length(unique(na.omit(full_filt()$Assay)))       })
    output$vb_genes  <- renderText({ length(unique(na.omit(full_filt()$Gene)))        })
    output$vb_models <- renderText({ length(unique(na.omit(full_filt()$Model_system)))})
    output$vb_ko     <- renderText({ length(unique(na.omit(full_filt()$KO_strat)))    })

    # ── Chart outputs ─────────────────────────────────────────────────────────

    output$model_plot <- renderPlotly({
      make_bar_chart(model_agg(), sel$model_system, "h", session$ns("model_src"))
    })
    output$ko_plot <- renderPlotly({
      make_bar_chart(ko_agg(), sel$ko_strat, "v", session$ns("ko_src"))
    })
    output$dpc_plot <- renderPlotly({
      make_bar_chart(dpc_agg(), sel$dpc, "v", session$ns("dpc_src"))
    })

    # ── Gene table headers (show active filter context) ───────────────────────

    make_tbl_header <- function(filter_name, filter_val) {
      renderUI({
        if (!is.null(filter_val))
          tags$p(class = "text-muted small mb-1",
                 "Filtered to: ",
                 tags$span(class = "badge bg-primary", filter_val))
        else
          tags$p(class = "text-muted small mb-1",
                 "Showing all (no cross-filter active for this dimension)")
      })
    }

    output$model_tbl_header <- make_tbl_header("Model System", sel$model_system)
    output$ko_tbl_header    <- make_tbl_header("KO Strategy",  sel$ko_strat)
    output$dpc_tbl_header   <- make_tbl_header("DPC",          sel$dpc)

    # ── Gene table outputs ────────────────────────────────────────────────────

    output$model_gene_tbl <- DT::renderDT({ render_gene_dt(model_genes()) })
    output$ko_gene_tbl    <- DT::renderDT({ render_gene_dt(ko_genes())    })
    output$dpc_gene_tbl   <- DT::renderDT({ render_gene_dt(dpc_genes())   })

    # ── Click handlers (cross-filter) ─────────────────────────────────────────

    observeEvent(event_data("plotly_click", source = session$ns("model_src")), {
      val <- event_data("plotly_click", source = session$ns("model_src"))$customdata[[1]]
      sel$model_system <- if (identical(sel$model_system, val)) NULL else val
    })
    observeEvent(event_data("plotly_click", source = session$ns("ko_src")), {
      val <- event_data("plotly_click", source = session$ns("ko_src"))$customdata[[1]]
      sel$ko_strat <- if (identical(sel$ko_strat, val)) NULL else val
    })
    observeEvent(event_data("plotly_click", source = session$ns("dpc_src")), {
      val <- event_data("plotly_click", source = session$ns("dpc_src"))$customdata[[1]]
      sel$dpc <- if (identical(sel$dpc, val)) NULL else val
    })

    # ── Clear buttons ─────────────────────────────────────────────────────────

    observeEvent(input$clear_model, { sel$model_system <- NULL })
    observeEvent(input$clear_ko,    { sel$ko_strat     <- NULL })
    observeEvent(input$clear_dpc,   { sel$dpc          <- NULL })
    observeEvent(input$clear_all,   {
      sel$model_system <- NULL
      sel$ko_strat     <- NULL
      sel$dpc          <- NULL
    })

    # ── UpSet: KO gene coverage ───────────────────────────────────────────────

    # Dynamic plot height: scale with number of groups so the dot-matrix breathes
    output$si_upset_ui <- renderUI({
      df      <- full_filt()
      grp_var <- input$si_upset_grp
      req(df, nrow(df) > 0, !is.null(grp_var))
      n_sets    <- length(unique(na.omit(df[[grp_var]])))
      height_px <- max(360L, 180L + n_sets * 45L)
      plotOutput(session$ns("si_upset_plot"), height = paste0(height_px, "px"))
    })

    output$si_upset_plot <- renderPlot({
      df      <- full_filt()
      grp_var <- input$si_upset_grp
      req(df, nrow(df) > 0, !is.null(grp_var))

      # Elements = unique KO target genes (Gene column)
      # Membership: gene ∈ group G  iff  ≥1 assay has that Gene & that group value
      groups    <- sort(unique(na.omit(df[[grp_var]])))
      all_genes <- sort(unique(na.omit(df$Gene)))
      req(length(groups) >= 2, length(all_genes) >= 1)

      # Sanitise column names for UpSetR (no leading digits, no special chars)
      clean_names <- gsub("[^A-Za-z0-9_]", "_", as.character(groups))
      clean_names <- ifelse(grepl("^[0-9]", clean_names),
                            paste0("G", clean_names),
                            clean_names)
      clean_names <- make.unique(clean_names)

      # Build binary membership matrix: gene × group
      mat <- as.data.frame(
        matrix(0L, nrow = length(all_genes), ncol = length(groups)),
        stringsAsFactors = FALSE
      )
      colnames(mat) <- clean_names

      for (i in seq_along(groups)) {
        in_grp <- unique(na.omit(df$Gene[df[[grp_var]] == groups[i]]))
        mat[all_genes %in% in_grp, clean_names[i]] <- 1L
      }

      mat <- mat[rowSums(mat) > 0, , drop = FALSE]
      req(nrow(mat) > 0)

      UpSetR::upset(
        mat,
        sets             = rev(clean_names),  # bottom → top in dot matrix
        order.by         = "freq",
        decreasing       = TRUE,
        mb.ratio         = c(0.60, 0.40),
        text.scale       = c(1.4, 1.2, 1.1, 1.0, 1.3, 0.9),
        point.size       = 3.2,
        line.size        = 0.9,
        sets.bar.color   = .CLRS$genes,   # orange — consistent with gene bars above
        main.bar.color   = .CLRS$assays,  # blue  — consistent with assay bars above
        sets.x.label     = "KO genes in group",
        mainbar.y.label  = "Intersection size (KO genes)",
        show.numbers     = "yes",
        number.angles    = 0,
        set_size.show    = TRUE,
        set_size.numbers_size = 3.5
      )
    }, bg = "white")

    # ── Selection summary card ────────────────────────────────────────────────

    output$selection_summary <- renderUI({
      items <- list()
      if (!is.null(sel$model_system))
        items <- c(items, list(tags$p(tags$strong("Model System: "),
                                      tags$span(class = "badge bg-primary",
                                                sel$model_system))))
      if (!is.null(sel$ko_strat))
        items <- c(items, list(tags$p(tags$strong("KO Strategy: "),
                                      tags$span(class = "badge bg-primary",
                                                sel$ko_strat))))
      if (!is.null(sel$dpc))
        items <- c(items, list(tags$p(tags$strong("DPC: "),
                                      tags$span(class = "badge bg-primary",
                                                sel$dpc))))
      if (length(items) == 0)
        tags$p(class = "text-muted small",
               "No filters active. Click any bar to cross-filter the other charts.")
      else
        tagList(items)
    })

    # ── KO Gene Study Coverage table ──────────────────────────────────────────

    output$gene_panel_tbl <- DT::renderDT({
      df <- full_filt()
      req(nrow(df) > 0)

      # Keep only short categorical columns — drop Assay and any free-text fields.
      keep_cols <- c("Gene", "Model_system", "KO_strat", "DPC",
                     "Cell_Line", "DPC_days", "Replicate")
      cols      <- keep_cols[keep_cols %in% colnames(df)]

      df <- df[order(df$Gene, df$Model_system), cols, drop = FALSE]

      # Pretty column names for display
      pretty <- gsub("_", " ", cols)

      # Gene column is index 0 (JS) — used for RowGroup, then hidden in body.
      DT::datatable(
        df,
        colnames   = pretty,
        extensions = "RowGroup",
        rownames   = FALSE,
        filter     = "top",
        options    = list(
          # Group rows by Gene (column index 0 in JS)
          rowGroup       = list(dataSrc = 0),
          # Hide the Gene column in the row body (it appears as the group header)
          columnDefs     = list(
            list(visible = FALSE, targets = 0),
            list(className = "dt-left", targets = "_all")
          ),
          pageLength     = 25,
          lengthMenu     = list(c(25, 50, 100, -1), c("25", "50", "100", "All")),
          scrollX        = TRUE,
          scrollY        = "520px",
          scrollCollapse = TRUE,
          dom            = "lftip",   # length + filter + table + info + paging
          language       = list(
            info = "Showing _START_ to _END_ of _TOTAL_ assay entries"
          )
        ),
        class = "display compact stripe hover"
      )
    })

    # ── Assay Coverage Matrix ─────────────────────────────────────────────────

    # Dynamic height so each gene row has breathing room
    output$matrix_ui <- renderUI({
      df    <- full_filt()
      req(nrow(df) > 0)
      n_genes   <- length(unique(na.omit(df$Gene)))
      height_px <- max(360L, 120L + n_genes * 28L)
      plotlyOutput(session$ns("coverage_matrix"),
                   height = paste0(height_px, "px"))
    })

    output$coverage_matrix <- renderPlotly({
      df    <- full_filt()
      x_var <- req(input$matrix_x)
      req(nrow(df) > 0)

      genes  <- sort(unique(na.omit(df$Gene)))
      x_vals <- sort(unique(na.omit(df[[x_var]])))

      # Build n-assay count matrix and hover-text matrix
      n_mat   <- matrix(0L, nrow = length(genes), ncol = length(x_vals),
                        dimnames = list(genes, x_vals))
      txt_mat <- matrix("", nrow = length(genes), ncol = length(x_vals),
                        dimnames = list(genes, x_vals))

      for (g in genes) {
        for (x in x_vals) {
          sub    <- df[!is.na(df$Gene)      & df$Gene      == g &
                       !is.na(df[[x_var]]) & df[[x_var]] == x, ]
          assays <- sort(unique(na.omit(sub$Assay)))
          n_mat[g, x]   <- length(assays)
          txt_mat[g, x] <- if (length(assays) == 0)
            paste0("<b>", g, " \u00d7 ", x, "</b><br>No assays")
          else
            paste0("<b>", g, " \u00d7 ", x, "</b><br>",
                   length(assays), " assay(s):<br>",
                   paste(assays, collapse = "<br>"))
        }
      }

      plot_ly(
        x            = x_vals,
        y            = genes,
        z            = n_mat,
        text         = txt_mat,
        type         = "heatmap",
        colorscale   = list(c(0, "#eef2ff"), c(1, .CLRS$assays)),
        hovertemplate = "%{text}<extra></extra>",
        showscale    = TRUE,
        colorbar     = list(title = "# Assays", thickness = 14, len = 0.6)
      ) |>
        layout(
          xaxis  = list(title = gsub("_", " ", x_var), automargin = TRUE,
                        tickangle = -30),
          yaxis  = list(title = "KO Gene", automargin = TRUE,
                        autorange = "reversed"),
          margin = list(l = 10, r = 10, t = 10, b = 60)
        ) |>
        config(displayModeBar = FALSE)
    })

  })
}
