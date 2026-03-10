# ================================================================
############### PHENOTYPE BROWSER MODULE #########################
# ================================================================
# Allows selection of an assay + phenotype database (IMPC, HPO,
# OMIM, Orphanet), then shows:
#   • Value boxes  — total DEG counts + annotated fractions
#   • Coverage bar — annotated vs unannotated DEGs per direction
#   • Top-N phenotypes stacked bar (plotly), per-zygosity subplots
#     for IMPC (each subplot only shows its own phenotypes)
#   • Summary DT   — one row per phenotype (+ zygosity for IMPC)
#   • Assay Overlap — ranked bar of all assays by N annotated DEGs
# ================================================================

library(dplyr)
library(DT)

# ── Dataset metadata ──────────────────────────────────────────────────────────

.PB_DATASETS <- c(
  "IMPC Mouse Phenotypes" = "impc",
  "HPO Human Phenotypes"  = "hpo",
  "OMIM Phenotypes"       = "omim",
  "Orphanet Disorders"    = "orphanet"
)

.PB_TABLE <- c(
  impc     = "impc_phenotypes",
  hpo      = "hpo_phenotypes",
  omim     = "omim_phenotypes",
  orphanet = "orphanet_data"
)

# Column holding the phenotype / term name in each table
.PB_COL <- c(
  impc     = "impc_phenotypes",
  hpo      = "hpo_name",
  omim     = "omim_phenotype",
  orphanet = "disorder_name"
)

# Human-readable column label for the DT
.PB_LABEL <- c(
  impc     = "Phenotype",
  hpo      = "HPO Term",
  omim     = "OMIM Phenotype",
  orphanet = "Disorder"
)

.PB_CLRS     <- c(up = "#E63946", down = "#457B9D")
.PB_ANN_CLRS <- c("Annotated" = "#6AB187", "No annotation" = "#CBD5D9")

# ── UI ────────────────────────────────────────────────────────────────────────

phenotype_browserUI <- function(id) {
  ns <- NS(id)
  tagList(

    # ── Controls card ────────────────────────────────────────────────────────
    card(
      class = "mb-3",
      card_header(
        class = "bg-dark text-white d-flex align-items-center gap-2",
        bsicons::bs_icon("bar-chart-fill"), " Phenotype Browser"
      ),
      card_body(
        layout_columns(
          col_widths = c(3, 3, 4, 2),
          selectInput(
            ns("pb_gene"), "KO Gene",
            choices = NULL, multiple = TRUE, selectize = TRUE
          ),
          selectInput(
            ns("pb_model"), "Model System",
            choices = NULL, multiple = TRUE, selectize = TRUE
          ),
          selectizeInput(
            ns("pb_assay"), "Assay",
            choices = NULL,
            options = list(placeholder = "Select an assay \u2026")
          ),
          selectInput(
            ns("pb_dataset"), "Phenotype Database",
            choices = .PB_DATASETS, selected = "impc"
          )
        ),
        layout_columns(
          col_widths = c(4, 4, 4),
          # IMPC-only: zygosity selector
          conditionalPanel(
            condition = sprintf("input['%s'] === 'impc'", ns("pb_dataset")),
            checkboxGroupInput(
              ns("pb_zyg"), "IMPC Zygosity",
              choices  = c("Homozygote" = "homozygote",
                            "Heterozygote" = "heterozygote"),
              selected = c("homozygote", "heterozygote"),
              inline   = TRUE
            )
          ),
          sliderInput(
            ns("pb_n_top"), "Top N phenotypes shown",
            min = 5, max = 50, value = 20, step = 5
          ),
          radioButtons(
            ns("pb_dir"), "Direction",
            choices  = c("Both" = "both",
                          "\u2191 Up only" = "up",
                          "\u2193 Down only" = "down"),
            selected = "both",
            inline   = TRUE
          )
        )
      )
    ),

    # ── Value boxes ───────────────────────────────────────────────────────────
    uiOutput(ns("pb_value_boxes")),

    # ── Charts row ───────────────────────────────────────────────────────────
    layout_columns(
      col_widths = c(4, 8),
      card(
        full_screen = FALSE,
        card_header("Annotation Coverage"),
        card_body(plotlyOutput(ns("pb_coverage"), height = "180px"))
      ),
      card(
        full_screen = TRUE,
        card_header("Top Phenotypes"),
        card_body(plotlyOutput(ns("pb_top"), height = "500px"))
      )
    ),

    # ── Phenotype summary table ───────────────────────────────────────────────
    card(
      class = "mt-3",
      card_header(
        class = "d-flex align-items-center gap-2",
        bsicons::bs_icon("table"), " Phenotype Summary Table"
      ),
      card_body(DT::DTOutput(ns("pb_tbl")))
    ),

    # ── Assay annotation overlap ──────────────────────────────────────────────
    card(
      class = "mt-3",
      card_header(
        class = "bg-dark text-white d-flex align-items-center gap-2",
        bsicons::bs_icon("bar-chart-fill"),
        " Assay Annotation Overlap"
      ),
      card_body(
        div(
          class = "d-flex align-items-center gap-3 mb-3",
          actionButton(
            ns("pb_run_overlap"),
            "Compute Overlap Across Assays",
            class = "btn btn-outline-dark"
          ),
          span(class = "text-muted small",
               "Ranks all assays (within current Gene/Model filters) by N annotated DEGs")
        ),
        uiOutput(ns("pb_overlap_info")),
        plotlyOutput(ns("pb_assay_overlap"), height = "700px"),
        div(class = "mt-3", DT::DTOutput(ns("pb_assay_tbl")))
      )
    )
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

phenotype_browserServer <- function(id, con_r, study_info_r) {
  moduleServer(id, function(input, output, session) {

    # ── Populate filter dropdowns ─────────────────────────────────────────────
    observe({
      si <- study_info_r()
      updateSelectInput(session, "pb_gene",
                        choices = c("(All genes)"         = "",
                                    sort(unique(si$Gene))))
      updateSelectInput(session, "pb_model",
                        choices = c("(All model systems)" = "",
                                    sort(unique(si$Model_system))))
    })

    # ── Filtered assay list ───────────────────────────────────────────────────
    filtered_assays <- reactive({
      si     <- study_info_r()
      genes  <- Filter(nzchar, if (is.null(input$pb_gene))  character(0) else input$pb_gene)
      models <- Filter(nzchar, if (is.null(input$pb_model)) character(0) else input$pb_model)
      if (length(genes)  > 0) si <- si[si$Gene         %in% genes,  ]
      if (length(models) > 0) si <- si[si$Model_system %in% models, ]
      sort(unique(si$Assay))
    })

    observeEvent(filtered_assays(), {
      updateSelectizeInput(session, "pb_assay",
                           choices = filtered_assays(), server = TRUE)
    }, ignoreNULL = FALSE)

    # ── Load DEG data for the selected assay ──────────────────────────────────
    deg_data <- reactive({
      assay <- input$pb_assay
      req(assay, nzchar(assay))
      con <- con_r()
      tryCatch(
        dbGetQuery(con, sprintf(
          "SELECT COALESCE(NULLIF(TRIM(hgnc_symbol), ''), gene_ID) AS symbol, DEG
             FROM main.\"%s\"
            WHERE DEG IN ('up', 'down')",
          assay
        )),
        error = function(e) NULL
      )
    })

    # ── Load phenotype annotation table ───────────────────────────────────────
    pheno_data <- reactive({
      ds  <- input$pb_dataset; req(ds)
      con <- con_r()
      tryCatch({
        df <- dbReadTable(con, .PB_TABLE[[ds]])
        if (ds == "impc") {
          zyg <- input$pb_zyg; req(length(zyg) > 0)
          df  <- df[!is.na(df$impc_zygosity) & df$impc_zygosity %in% zyg, ]
        }
        df
      }, error = function(e) NULL)
    })

    # ── Inner join DEG × phenotype ────────────────────────────────────────────
    joined_data <- reactive({
      deg    <- deg_data();  req(!is.null(deg),  nrow(deg) > 0)
      ph     <- pheno_data(); req(!is.null(ph))
      ds     <- input$pb_dataset
      ph_col <- .PB_COL[[ds]]

      keep   <- unique(c("symbol", ph_col,
                          if (ds == "impc") "impc_zygosity"))
      ph_sub <- ph[, keep, drop = FALSE]
      ph_sub <- ph_sub[!is.na(ph_sub[[ph_col]]) &
                         nzchar(trimws(as.character(ph_sub[[ph_col]]))), ]

      merged <- merge(deg, ph_sub, by = "symbol", all = FALSE)
      if (nrow(merged) == 0L) return(NULL)
      merged
    })

    # ── Coverage (annotated vs unannotated counts per direction) ──────────────
    coverage_data <- reactive({
      deg <- deg_data();  req(!is.null(deg), nrow(deg) > 0)
      ph  <- pheno_data(); req(!is.null(ph))

      ann_syms <- unique(ph$symbol[!is.na(ph$symbol) & nzchar(ph$symbol)])

      data.frame(
        symbol    = deg$symbol,
        DEG       = deg$DEG,
        annotated = deg$symbol %in% ann_syms
      ) |>
        dplyr::group_by(DEG, annotated) |>
        dplyr::summarise(n = dplyr::n_distinct(symbol), .groups = "drop") |>
        dplyr::mutate(
          dir_lbl = ifelse(DEG == "up", "\u2191 Up", "\u2193 Down"),
          ann_lbl = ifelse(annotated, "Annotated", "No annotation")
        )
    })

    # ── Top-N phenotype summary — respects direction filter ───────────────────
    top_pheno_data <- reactive({
      jd     <- joined_data(); req(!is.null(jd), nrow(jd) > 0)
      ds     <- input$pb_dataset
      ph_col <- .PB_COL[[ds]]
      n_top  <- input$pb_n_top; if (is.null(n_top)) n_top <- 20L
      dir    <- input$pb_dir;   if (is.null(dir))   dir   <- "both"

      # Apply direction filter
      if (dir != "both") {
        jd <- jd[jd$DEG == dir, ]
        req(nrow(jd) > 0)
      }

      grp_cols <- unique(c(ph_col, "DEG",
                            if (ds == "impc") "impc_zygosity"))

      summ <- jd |>
        dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) |>
        dplyr::summarise(
          n_genes = dplyr::n_distinct(symbol),
          genes   = paste(sort(unique(symbol)), collapse = "; "),
          .groups = "drop"
        )

      # Global top-N labels (used for non-IMPC and to cap IMPC per-zygosity labels)
      top_labels <- summ |>
        dplyr::group_by(.data[[ph_col]]) |>
        dplyr::summarise(total = sum(n_genes), .groups = "drop") |>
        dplyr::slice_max(total, n = n_top, with_ties = FALSE) |>
        dplyr::arrange(total) |>
        dplyr::pull(.data[[ph_col]])

      list(summ = summ, top_labels = top_labels, ph_col = ph_col)
    })

    # ── Value boxes ───────────────────────────────────────────────────────────
    output$pb_value_boxes <- renderUI({
      deg <- deg_data()
      if (is.null(deg) || nrow(deg) == 0L)
        return(div(class = "text-muted text-center py-3",
                   "Select an assay above to load DEG data."))

      n_up   <- sum(deg$DEG == "up")
      n_dn   <- sum(deg$DEG == "down")
      jd     <- joined_data()
      ann_up <- if (!is.null(jd)) length(unique(jd$symbol[jd$DEG == "up"]))   else 0L
      ann_dn <- if (!is.null(jd)) length(unique(jd$symbol[jd$DEG == "down"])) else 0L
      pct_up <- if (n_up > 0) sprintf("%.0f%%", 100 * ann_up / n_up) else "\u2014"
      pct_dn <- if (n_dn > 0) sprintf("%.0f%%", 100 * ann_dn / n_dn) else "\u2014"

      ds_short <- switch(input$pb_dataset,
        impc = "IMPC", hpo = "HPO", omim = "OMIM", orphanet = "Orphanet",
        input$pb_dataset
      )

      layout_columns(
        col_widths = c(3, 3, 3, 3),
        value_box(title = "\u2191 Up DEGs",
                  value = format(n_up, big.mark = ","),
                  showcase = bsicons::bs_icon("bar-chart-fill"), theme = "danger"),
        value_box(title = "\u2193 Down DEGs",
                  value = format(n_dn, big.mark = ","),
                  showcase = bsicons::bs_icon("bar-chart-fill"), theme = "primary"),
        value_box(title = paste0("\u2191 Up w/ ", ds_short),
                  value = sprintf("%d / %d (%s)", ann_up, n_up, pct_up),
                  showcase = bsicons::bs_icon("table"), theme = "danger"),
        value_box(title = paste0("\u2193 Down w/ ", ds_short),
                  value = sprintf("%d / %d (%s)", ann_dn, n_dn, pct_dn),
                  showcase = bsicons::bs_icon("table"), theme = "primary")
      )
    })

    # ── Coverage plot ─────────────────────────────────────────────────────────
    output$pb_coverage <- renderPlotly({
      cov <- coverage_data(); req(!is.null(cov), nrow(cov) > 0)

      # Apply direction filter to coverage display
      dir <- input$pb_dir; if (is.null(dir)) dir <- "both"
      if (dir != "both") cov <- cov[cov$DEG == dir, ]

      fig <- plot_ly()
      for (a in c(FALSE, TRUE)) {
        sub <- cov[cov$annotated == a, ]
        if (nrow(sub) == 0L) next
        lbl <- if (a) "Annotated" else "No annotation"
        fig <- fig |>
          add_bars(
            data          = sub,
            x             = ~n,
            y             = ~dir_lbl,
            name          = lbl,
            orientation   = "h",
            marker        = list(color = .PB_ANN_CLRS[[lbl]]),
            hovertemplate = paste0("<b>%{y}</b><br>",
                                   lbl, ": %{x} DEGs<extra></extra>")
          )
      }
      fig |> layout(
        barmode = "stack",
        xaxis   = list(title = "N DEGs"),
        yaxis   = list(title = "", categoryorder = "array",
                       categoryarray = c("\u2193 Down", "\u2191 Up")),
        legend  = list(orientation = "h", x = 0, y = -0.45),
        margin  = list(l = 70, t = 5, b = 80, r = 10)
      )
    })

    # ── Top phenotypes plot ───────────────────────────────────────────────────
    output$pb_top <- renderPlotly({
      td <- top_pheno_data(); req(!is.null(td))
      ds     <- input$pb_dataset
      summ   <- td$summ
      labels <- td$top_labels    # global top-N (used for non-IMPC)
      ph_col <- td$ph_col
      dir    <- input$pb_dir; if (is.null(dir)) dir <- "both"

      # Helper: build stacked bars for a df slice using explicit y-axis levels
      build_bars <- function(fig, df, ph_labels, show_legend = TRUE) {
        df$phenotype <- factor(df[[ph_col]], levels = ph_labels)
        df <- df[!is.na(df$phenotype), ]          # drop any labels not in this slice
        dirs <- if (dir == "both") c("down", "up") else dir
        for (d in dirs) {
          sub <- df[df$DEG == d, ]
          if (nrow(sub) == 0L) next
          nm  <- if (d == "up") "\u2191 Up" else "\u2193 Down"
          fig <- fig |>
            add_bars(
              data          = sub,
              x             = ~n_genes,
              y             = ~phenotype,
              name          = nm,
              orientation   = "h",
              marker        = list(color = .PB_CLRS[[d]]),
              legendgroup   = nm,
              showlegend    = show_legend,
              hovertemplate = paste0("<b>%{y}</b><br>", nm,
                                     ": %{x} gene(s)<extra></extra>")
            )
        }
        fig
      }

      if (ds == "impc") {
        zyg_lvls <- sort(unique(summ$impc_zygosity))
        subs <- lapply(seq_along(zyg_lvls), function(i) {
          z    <- zyg_lvls[[i]]
          df_z <- summ[summ$impc_zygosity == z &
                         summ[[ph_col]] %in% labels, ]

          # ── FIX: use only phenotypes present in THIS zygosity, ordered by
          #         their own totals so each subplot has a clean, dense y-axis
          z_labels <- df_z |>
            dplyr::group_by(.data[[ph_col]]) |>
            dplyr::summarise(total = sum(n_genes), .groups = "drop") |>
            dplyr::arrange(total) |>
            dplyr::pull(.data[[ph_col]])

          build_bars(plot_ly(), df_z, ph_labels = z_labels,
                     show_legend = (i == 1L)) |>
            layout(
              barmode = "stack",
              xaxis   = list(title = "N DEGs"),
              yaxis   = list(title = "", tickfont = list(size = 8)),
              annotations = list(list(
                text = tools::toTitleCase(z), showarrow = FALSE,
                xref = "paper", yref = "paper",
                x = 0.5, y = 1.05, xanchor = "center",
                font = list(size = 11, color = "#444")
              ))
            )
        })

        # Dynamic height based on max n phenotypes across any zygosity
        subplot(subs,
                nrows  = length(zyg_lvls),
                shareX = FALSE, shareY = FALSE,
                titleX = TRUE, titleY = TRUE) |>
          layout(
            legend = list(orientation = "h", x = 0, y = -0.04),
            margin = list(l = 220, t = 50, b = 40, r = 10)
          )

      } else {
        df <- summ[summ[[ph_col]] %in% labels, ]
        build_bars(plot_ly(), df, ph_labels = labels) |>
          layout(
            barmode = "stack",
            xaxis   = list(title = "N DEGs with phenotype"),
            yaxis   = list(title = "", tickfont = list(size = 9)),
            legend  = list(orientation = "h", x = 0, y = -0.08),
            margin  = list(l = 220, t = 20, b = 60, r = 10)
          )
      }
    })

    # ── Phenotype summary DT ──────────────────────────────────────────────────
    output$pb_tbl <- DT::renderDT({
      td <- top_pheno_data(); req(!is.null(td))
      ds     <- input$pb_dataset
      summ   <- td$summ
      labels <- td$top_labels
      ph_col <- td$ph_col
      ph_lbl <- .PB_LABEL[[ds]]

      top <- summ[summ[[ph_col]] %in% labels, ]

      if (ds == "impc") {
        wide <- top |>
          dplyr::group_by(.data[[ph_col]], impc_zygosity) |>
          dplyr::summarise(
            n_up       = sum(n_genes[DEG == "up"],   na.rm = TRUE),
            n_down     = sum(n_genes[DEG == "down"], na.rm = TRUE),
            n_total    = n_up + n_down,
            genes_up   = paste(na.omit(genes[DEG == "up"]),   collapse = "; "),
            genes_down = paste(na.omit(genes[DEG == "down"]), collapse = "; "),
            .groups    = "drop"
          ) |> dplyr::arrange(dplyr::desc(n_total))
        names(wide)[names(wide) == ph_col]          <- ph_lbl
        names(wide)[names(wide) == "impc_zygosity"] <- "Zygosity"
      } else {
        wide <- top |>
          dplyr::group_by(.data[[ph_col]]) |>
          dplyr::summarise(
            n_up       = sum(n_genes[DEG == "up"],   na.rm = TRUE),
            n_down     = sum(n_genes[DEG == "down"], na.rm = TRUE),
            n_total    = n_up + n_down,
            genes_up   = paste(na.omit(genes[DEG == "up"]),   collapse = "; "),
            genes_down = paste(na.omit(genes[DEG == "down"]), collapse = "; "),
            .groups    = "drop"
          ) |> dplyr::arrange(dplyr::desc(n_total))
        names(wide)[names(wide) == ph_col] <- ph_lbl
      }

      names(wide) <- gsub("^n_up$",     "\u2191 N Up",        names(wide))
      names(wide) <- gsub("^n_down$",   "\u2193 N Down",       names(wide))
      names(wide) <- gsub("^n_total$",  "N Total",            names(wide))
      names(wide) <- gsub("^genes_up$",   "Genes \u2191 Up",  names(wide))
      names(wide) <- gsub("^genes_down$", "Genes \u2193 Down",names(wide))

      DT::datatable(wide,
        class = "compact row-border hover", rownames = FALSE,
        extensions = "Buttons",
        options = list(pageLength = 25, scrollX = TRUE,
                       dom = "Bfrtip", buttons = c("csv", "excel")))
    })

    # ── Assay overlap: compute across all filtered assays ─────────────────────
    assay_overlap_data <- eventReactive(input$pb_run_overlap, {
      assays <- filtered_assays()
      req(length(assays) > 0)
      ph  <- pheno_data(); req(!is.null(ph))
      ds  <- input$pb_dataset
      con <- con_r()

      # Annotated symbol set for this database + zygosity selection
      ph_col   <- .PB_COL[[ds]]
      ph_clean <- ph[!is.na(ph[[ph_col]]) & nzchar(trimws(as.character(ph[[ph_col]]))), ]
      ann_syms <- unique(ph_clean$symbol[!is.na(ph_clean$symbol) & nzchar(ph_clean$symbol)])

      # Build one UNION ALL query across all assays
      parts <- sprintf(
        "SELECT '%s' AS assay,
                COALESCE(NULLIF(TRIM(hgnc_symbol), ''), gene_ID) AS symbol,
                DEG
           FROM main.\"%s\"
          WHERE DEG IN ('up', 'down')",
        gsub("'", "''", assays), assays   # escape single quotes in names
      )
      all_degs <- tryCatch(
        dbGetQuery(con, paste(parts, collapse = "\nUNION ALL\n")),
        error = function(e) NULL
      )
      req(!is.null(all_degs), nrow(all_degs) > 0)

      all_degs$annotated <- all_degs$symbol %in% ann_syms

      all_degs |>
        dplyr::group_by(assay) |>
        dplyr::summarise(
          n_degs      = dplyr::n_distinct(symbol),
          n_up        = dplyr::n_distinct(symbol[DEG == "up"]),
          n_down      = dplyr::n_distinct(symbol[DEG == "down"]),
          n_ann       = dplyr::n_distinct(symbol[annotated]),
          n_ann_up    = dplyr::n_distinct(symbol[DEG == "up"   & annotated]),
          n_ann_down  = dplyr::n_distinct(symbol[DEG == "down" & annotated]),
          pct_ann     = round(100 * n_ann / pmax(n_degs, 1), 1),
          .groups     = "drop"
        ) |>
        dplyr::arrange(dplyr::desc(n_ann))
    })

    output$pb_overlap_info <- renderUI({
      aod <- assay_overlap_data()
      if (is.null(aod)) return(NULL)
      ds_short <- switch(input$pb_dataset,
        impc = "IMPC", hpo = "HPO", omim = "OMIM", orphanet = "Orphanet",
        input$pb_dataset
      )
      div(class = "alert alert-light border mb-2 py-2 small",
          sprintf("%d assays computed. Database: %s. Showing top 40 by N annotated DEGs in chart; full table below.",
                  nrow(aod), ds_short))
    })

    output$pb_assay_overlap <- renderPlotly({
      aod <- assay_overlap_data(); req(!is.null(aod), nrow(aod) > 0)
      dir <- input$pb_dir; if (is.null(dir)) dir <- "both"

      # Top 40 by total annotated
      top40 <- head(aod, 40)
      top40$assay_f <- factor(top40$assay, levels = rev(top40$assay))

      fig <- plot_ly()

      if (dir %in% c("both", "up")) {
        fig <- fig |>
          add_bars(
            data          = top40,
            x             = ~n_ann_up,
            y             = ~assay_f,
            name          = "\u2191 Up annotated",
            orientation   = "h",
            marker        = list(color = .PB_CLRS[["up"]]),
            hovertemplate = paste0("<b>%{y}</b><br>Up annotated: %{x}",
                                   " / ", top40$n_up, "<extra></extra>")
          )
      }
      if (dir %in% c("both", "down")) {
        fig <- fig |>
          add_bars(
            data          = top40,
            x             = ~n_ann_down,
            y             = ~assay_f,
            name          = "\u2193 Down annotated",
            orientation   = "h",
            marker        = list(color = .PB_CLRS[["down"]]),
            hovertemplate = paste0("<b>%{y}</b><br>Down annotated: %{x}",
                                   " / ", top40$n_down, "<extra></extra>")
          )
      }

      # % annotation as scatter on secondary x-axis
      pct_col <- switch(dir,
        up   = top40$n_ann_up   / pmax(top40$n_up,   1) * 100,
        down = top40$n_ann_down / pmax(top40$n_down, 1) * 100,
        top40$pct_ann
      )
      fig <- fig |>
        add_trace(
          data   = top40,
          x      = pct_col,
          y      = ~assay_f,
          type   = "scatter",
          mode   = "markers",
          name   = "% annotated",
          xaxis  = "x2",
          marker = list(color = "#F28E2B", symbol = "diamond", size = 7),
          hovertemplate = "<b>%{y}</b><br>% annotated: %{x:.1f}%<extra></extra>"
        )

      fig |> layout(
        barmode = "stack",
        xaxis   = list(title = "N annotated DEGs", side = "bottom"),
        xaxis2  = list(title = "% annotated", overlaying = "x",
                       side = "top", showgrid = FALSE,
                       range = c(0, 100)),
        yaxis   = list(title = "", tickfont = list(size = 8)),
        legend  = list(orientation = "h", x = 0, y = -0.04),
        margin  = list(l = 300, t = 50, b = 50, r = 10)
      )
    })

    output$pb_assay_tbl <- DT::renderDT({
      aod <- assay_overlap_data(); req(!is.null(aod), nrow(aod) > 0)
      dir <- input$pb_dir; if (is.null(dir)) dir <- "both"

      tbl <- aod
      names(tbl) <- c("Assay", "N DEGs", "\u2191 N Up", "\u2193 N Down",
                       "N Annotated", "\u2191 Annotated", "\u2193 Annotated",
                       "% Annotated")

      # Hide direction columns not relevant to filter
      hide_cols <- if (dir == "up")   c(3, 6)    # down cols (0-indexed later)
                   else if (dir == "down") c(2, 5)
                   else integer(0)

      DT::datatable(tbl,
        class = "compact row-border hover", rownames = FALSE,
        extensions = "Buttons",
        options = list(
          pageLength = 30, scrollX = TRUE,
          dom = "Bfrtip", buttons = c("csv", "excel"),
          columnDefs = if (length(hide_cols) > 0)
            list(list(visible = FALSE, targets = hide_cols))
          else list()
        )
      ) |>
        DT::formatStyle("\u2191 Annotated",
          background = DT::styleColorBar(aod$n_ann_up,   "#fde8ea"),
          backgroundSize = "98% 80%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center") |>
        DT::formatStyle("\u2193 Annotated",
          background = DT::styleColorBar(aod$n_ann_down, "#ddeef7"),
          backgroundSize = "98% 80%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center") |>
        DT::formatStyle("% Annotated",
          background = DT::styleColorBar(aod$pct_ann, "#e8f5e9"),
          backgroundSize = "98% 80%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center")
    })

    # ── outputOptions — MUST follow all output$ assignments ──────────────────
    for (o in c("pb_value_boxes", "pb_coverage", "pb_top", "pb_tbl",
                "pb_overlap_info", "pb_assay_overlap", "pb_assay_tbl"))
      outputOptions(output, o, suspendWhenHidden = TRUE)
  })
}
