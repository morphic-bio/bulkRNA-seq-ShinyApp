library(shiny)
library(duckdb)
library(DBI)
library(bslib)
library(bsicons)
library(plotly)

source("modules/helpers.R")           # shared constants + helper functions
source("modules/home.R")              # landing page
source("modules/studies_overview.R")  # page 1
source("modules/compare_studies.R")   # page 2
source("modules/deg_annotations.R")   # page 3 (all annotation tabs incl. phenotype)

DB_PATH <- "morphic_bulkRNA3.duckdb"

# Load study_info once at startup
study_info_data <- local({
  con <- dbConnect(duckdb(), DB_PATH, read_only = TRUE)
  on.exit(dbDisconnect(con))
  dbGetQuery(con, "SELECT * FROM study_info")
})

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- page_navbar(
  id       = "main_nav",
  title    = "MorPhiC Bulk RNA-seq Explorer",
  theme    = bs_theme(bootswatch = "litera", base_font = font_google("Inter")),
  fillable = FALSE,

  # Home
  nav_panel(
    title = tagList(bsicons::bs_icon("house-fill"), " Home"),
    value = "home",
    div(class = "container-fluid py-3",
        homeUI("home"))
  ),

  # Page 1 — Assays Overview
  nav_panel(
    title = tagList(bsicons::bs_icon("bar-chart-fill"), " Assays Overview"),
    value = "page1",
    div(class = "container-fluid py-3",
        studies_overviewUI("studies_overview"))
  ),

  # Page 2 — Assay Comparison
  nav_panel(
    title = tagList(bsicons::bs_icon("grid-3x3-gap-fill"), " Assay Comparison"),
    value = "page2",
    div(class = "container-fluid py-3",
        compare_studiesUI("compare_studies"))
  ),

  # Page 3 — Assay Analysis
  nav_panel(
    title = tagList(bsicons::bs_icon("diagram-3"), " Assay Analysis"),
    value = "page3",
    div(class = "container-fluid py-3",
        deg_annotationsUI("deg_annotations"))
  ),

  nav_spacer(),
  nav_item(
    tags$span(
      class = "badge rounded-pill bg-secondary",
      style = "font-size: 0.65rem; opacity: 0.7; vertical-align: middle;",
      "v1.0"
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  db_con <- dbConnect(duckdb(), DB_PATH, read_only = TRUE)
  onStop(function() {
    if (dbIsValid(db_con)) dbDisconnect(db_con)
  })

  study_info_r <- reactive(study_info_data)

  homeServer(
    id             = "home",
    parent_session = session,
    con_r          = reactive(db_con),
    study_info_r   = study_info_r
  )

  studies_overviewServer(
    id           = "studies_overview",
    con_r        = reactive(db_con),
    study_info_r = study_info_r
  )

  compare_studiesServer(
    id           = "compare_studies",
    con_r        = reactive(db_con),
    study_info_r = study_info_r
  )

  deg_annotationsServer(
    id           = "deg_annotations",
    con_r        = reactive(db_con),
    study_info_r = study_info_r
  )
}

shinyApp(ui, server)
