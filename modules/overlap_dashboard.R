library(shiny)
library(bslib)
library(bsicons)
library(plotly)
library(DT)
library(dplyr)

# ── Constants ──────────────────────────────────────────────────────────────────

.OVERLAP_SOURCES <- c(
  "Reactome Pathways"    = "reactome",
  "GO Biological Process" = "go_bp",
  "GO Molecular Function" = "go_mf",
  "PANTHER Protein Class" = "panther_class",
  "PANTHER Protein Family"= "panther_family",
  "PanelApp Panels"      = "panelapp"
)

# group-by column per source (for morphic_overlap and deg_overlap tables)
.GROUP_COL <- c(
  reactome      = "path_name",
  go_bp         = "go_term_BP",
  go_mf         = "go_term_MF",
  panther_class = "CLASS_TERM",
  panther_family= "FAMILY_TERM",
  panelapp      = "panel_name"
)

.N_COL <- c(
  reactome      = "pathway_genes_n",
  go_bp         = "biological_process_genes_n",
  go_mf         = "molecular_function_genes_n",
  panther_class = "class_genes_n",
  panther_family= "family_genes_n",
  panelapp      = "panel_genes_n"
)

# ── UI ────────────────────────────────────────────────────────────────────────

overlap_dashboardUI <- function(id) {
  ns <- NS(id)

  tagList(

    # ── Row 1: controls ──────────────────────────────────────────────────────
    layout_columns(
      col_widths = c(3, 3, 3, 3),

      card(
        card_header(tagList(bs_icon("funnel"), " Source")),
        selectInput(ns("source"), NULL,
                    choices  = .OVERLAP_SOURCES,
                    selected = "reactome",
                    width    = "100%")
      ),

      card(
        card_header(tagList(bs_icon("arrow-down-up"), " DEG Direction")),
        radioButtons(ns("direction"), NULL,
                     choices  = c("Up & Down" = "updn", "Up only" = "up", "Down only" = "down"),
                     selected = "updn", inline = TRUE)
      ),

      card(
        card_header(tagList(bs_icon("sliders"), " Min DEG hits")),
        sliderInput(ns("min_n"), NULL, min = 1, max = 20, value = 2, step = 1)
      ),

      card(
        card_header(tagList(bs_icon("hash"), " Top N categories")),
        sliderInput(ns("top_n"), NULL, min = 5, max = 50, value = 20, step = 5)
      )
    ),

    # ── Row 2: Morphic gene overlap (static — all 57 KO genes) ───────────────
    card(
      card_header(
        tagList(bs_icon("diagram-3"), " Morphic KO Gene Overlap"),
        tooltip(
          bs_icon("info-circle"),
          "Percentage of each category's genes that are MorPhiC KO targets"
        )
      ),
      navset_tab(
        nav_panel("Chart",
          plotlyOutput(ns("morphic_bar"), height = "420px")
        ),
        nav_panel("Table",
          DTOutput(ns("morphic_tbl"))
        )
      )
    ),

    # ── Row 3: DEG overlap across assays ─────────────────────────────────────
    layout_columns(
      col_widths = c(8, 4),

      card(
        card_header(
          tagList(bs_icon("activity"), " DEG Overlap — Category × Assay Heatmap"),
          tooltip(
            bs_icon("info-circle"),
            "% of category genes that are DEGs in each assay. Filtered by min DEG hits and top N categories."
          )
        ),
        # assay filter — populated server-side
        uiOutput(ns("assay_filter_ui")),
        plotlyOutput(ns("deg_heatmap"), height = "500px")
      ),

      card(
        card_header(tagList(bs_icon("bar-chart-steps"), " Top Categories (summed hits)")),
        plotlyOutput(ns("deg_bar"), height = "500px")
      )
    ),

    # ── Row 4: gene-level detail table for selected category ─────────────────
    card(
      card_header(
        tagList(bs_icon("table"), " Gene Detail — click a cell in the heatmap"),
      ),
      DTOutput(ns("gene_detail_tbl"))
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

overlap_dashboardServer <- function(id, con_r, study_info_r) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── helpers ──────────────────────────────────────────────────────────────

    # load morphic_overlap table for chosen source
    morphic_data <- reactive({
      src <- input$source
      req(src)
      dbGetQuery(con_r(),
        sprintf("SELECT * FROM morphic_overlap.%s", src)
      )
    })

    # load deg_overlap table for chosen source (sym version)
    deg_data_raw <- reactive({
      src   <- input$source
      req(src)
      tbl   <- paste0(src, "_sym")
      dbGetQuery(con_r(),
        sprintf("SELECT * FROM deg_overlap.%s", tbl)
      )
    })

    # join deg_overlap with study_info to get metadata columns
    deg_data <- reactive({
      dd  <- deg_data_raw()
      si  <- study_info_r() |>
               select(Assay, Gene, Model_system, KO_strat, DPC) |>
               distinct(Assay, .keep_all = TRUE)
      left_join(dd, si, by = c("assay" = "Assay"))
    })

    # dynamic assay filter choices — grouped by Model_system
    output$assay_filter_ui <- renderUI({
      si <- study_info_r() |>
              select(Assay, Model_system) |>
              distinct()
      choices <- split(si$Assay, si$Model_system)
      selectizeInput(
        ns("assay_filter"),
        label    = NULL,
        choices  = choices,
        selected = si$Assay,
        multiple = TRUE,
        options  = list(placeholder = "Filter assays…", plugins = list("remove_button"))
      )
    })

    # n_col and group_col for current source
    grp_col <- reactive(.GROUP_COL[[input$source]])
    n_col   <- reactive(.N_COL[[input$source]])

    # direction-aware column names
    dir_cols <- reactive({
      switch(input$direction,
        updn = list(n = "n_updn", pct = "pct_updn", genes = "overlap_updn"),
        up   = list(n = "n_up",   pct = "pct_up",   genes = "overlap_up"),
        down = list(n = "n_down", pct = "pct_down",  genes = "overlap_down")
      )
    })

    # ── Morphic bar chart ────────────────────────────────────────────────────

    output$morphic_bar <- renderPlotly({
      df  <- morphic_data()
      gc  <- grp_col()
      req(nrow(df) > 0, gc %in% colnames(df))

      df <- df |>
        arrange(desc(pct_morphic_overlap), desc(overlap_n)) |>
        head(input$top_n) |>
        mutate(label = ifelse(nchar(.data[[gc]]) > 55,
                              paste0(substr(.data[[gc]], 1, 52), "…"),
                              .data[[gc]]))

      plot_ly(df,
        x           = ~pct_morphic_overlap,
        y           = ~reorder(label, pct_morphic_overlap),
        type        = "bar",
        orientation = "h",
        text        = ~paste0(overlap_n, " / ", .data[[n_col()]], " genes"),
        hovertemplate = paste0(
          "<b>%{y}</b><br>",
          "Overlap: %{text}<br>",
          "%{x:.1f}%<extra></extra>"
        ),
        marker = list(color = "#2C7BB6")
      ) |>
        layout(
          xaxis  = list(title = "% MorPhiC KO genes in category", range = c(0, 105)),
          yaxis  = list(title = "", automargin = TRUE),
          margin = list(l = 10),
          plot_bgcolor  = "white",
          paper_bgcolor = "white"
        )
    })

    output$morphic_tbl <- renderDT({
      df  <- morphic_data()
      gc  <- grp_col()
      req(nrow(df) > 0, gc %in% colnames(df))
      df |>
        select(Category = all_of(gc),
               N_category_genes = all_of(n_col()),
               overlap_n,
               pct_morphic_overlap,
               overlap_genes) |>
        mutate(pct_morphic_overlap = round(pct_morphic_overlap, 2)) |>
        datatable(
          rownames  = FALSE,
          filter    = "top",
          options   = list(pageLength = 15, scrollX = TRUE)
        )
    })

    # ── DEG heatmap ──────────────────────────────────────────────────────────

    # filtered + aggregated for heatmap
    heatmap_data <- reactive({
      dd   <- deg_data()
      dc   <- dir_cols()
      gc   <- grp_col()
      req(nrow(dd) > 0, gc %in% colnames(dd))

      assays_sel <- input$assay_filter
      req(length(assays_sel) > 0)

      dd <- dd |> filter(assay %in% assays_sel)

      # pick top N categories by total hits across assays
      top_cats <- dd |>
        group_by(across(all_of(gc))) |>
        summarise(total = sum(.data[[dc$n]], na.rm = TRUE), .groups = "drop") |>
        filter(total >= input$min_n) |>
        arrange(desc(total)) |>
        head(input$top_n) |>
        pull(all_of(gc))

      req(length(top_cats) > 0)

      dd |>
        filter(.data[[gc]] %in% top_cats) |>
        select(assay, category = all_of(gc),
               n      = all_of(dc$n),
               pct    = all_of(dc$pct),
               genes  = all_of(dc$genes),
               Gene, Model_system, KO_strat)
    })

    output$deg_heatmap <- renderPlotly({
      hd   <- heatmap_data()
      req(nrow(hd) > 0)

      # short assay labels
      hd <- hd |>
        mutate(assay_short = ifelse(nchar(assay) > 45,
                                    paste0("…", substr(assay, nchar(assay)-42, nchar(assay))),
                                    assay))

      plot_ly(
        data        = hd,
        x           = ~assay_short,
        y           = ~category,
        z           = ~pct,
        customdata  = ~paste0(assay, "||", category),
        type        = "heatmap",
        colorscale  = list(c(0, "#f7fbff"), c(0.5, "#6baed6"), c(1, "#08306b")),
        hovertemplate = paste0(
          "<b>%{y}</b><br>",
          "Assay: %{x}<br>",
          "%{z:.1f}% DEGs overlap<br>",
          "<extra></extra>"
        ),
        source = ns("heatmap_click")
      ) |>
        layout(
          xaxis  = list(title = "", tickangle = -45, automargin = TRUE),
          yaxis  = list(title = "", automargin = TRUE),
          margin = list(b = 120),
          plot_bgcolor  = "white",
          paper_bgcolor = "white"
        )
    })

    # ── DEG bar chart (top categories by total hits) ──────────────────────────

    output$deg_bar <- renderPlotly({
      hd   <- heatmap_data()
      gc   <- "category"
      req(nrow(hd) > 0)

      agg <- hd |>
        group_by(category) |>
        summarise(total_hits = sum(n, na.rm = TRUE), .groups = "drop") |>
        arrange(desc(total_hits)) |>
        mutate(label = ifelse(nchar(category) > 45,
                              paste0(substr(category, 1, 42), "…"),
                              category))

      plot_ly(agg,
        x           = ~total_hits,
        y           = ~reorder(label, total_hits),
        type        = "bar",
        orientation = "h",
        marker      = list(color = "#31a354"),
        hovertemplate = "<b>%{y}</b><br>Total DEG hits across assays: %{x}<extra></extra>"
      ) |>
        layout(
          xaxis  = list(title = "Total DEG hits across assays"),
          yaxis  = list(title = "", automargin = TRUE),
          margin = list(l = 10),
          plot_bgcolor  = "white",
          paper_bgcolor = "white"
        )
    })

    # ── Gene detail table — populated on heatmap click ────────────────────────

    clicked_cell <- reactiveVal(NULL)

    observe({
      ev <- event_data("plotly_click", source = ns("heatmap_click"))
      req(ev)
      # customdata encodes "assay||category"
      parts <- strsplit(ev$customdata, "\\|\\|")[[1]]
      if (length(parts) == 2) clicked_cell(list(assay = parts[1], category = parts[2]))
    })

    output$gene_detail_tbl <- renderDT({
      cell <- clicked_cell()
      if (is.null(cell)) {
        return(datatable(
          data.frame(Info = "Click a cell in the heatmap to see gene details"),
          rownames = FALSE, options = list(dom = "t")
        ))
      }

      hd  <- heatmap_data()
      dc  <- dir_cols()

      row <- hd |>
        filter(assay == cell$assay, category == cell$category)

      req(nrow(row) > 0)

      genes_str <- row$genes[1]
      if (is.na(genes_str) || genes_str == "") {
        return(datatable(
          data.frame(Info = "No overlapping genes for this cell"),
          rownames = FALSE, options = list(dom = "t")
        ))
      }

      gene_vec <- trimws(strsplit(genes_str, ";")[[1]])
      df <- data.frame(
        Assay    = cell$assay,
        Category = cell$category,
        Gene     = gene_vec,
        stringsAsFactors = FALSE
      )

      datatable(
        df,
        rownames = FALSE,
        caption  = htmltools::tags$caption(
          style = "caption-side: top; font-weight: bold;",
          sprintf("Overlapping genes — %s × %s", cell$category, cell$assay)
        ),
        options  = list(pageLength = 20, scrollX = TRUE)
      )
    })
  })
}
