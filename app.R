# SHINY APP FOR THE MAX MULTI-OMICS STUDY: --------------------------------
#
# page_navbar + bslib theme + custom CSS. Data is loaded once here with vroom;
# the plotting logic lives in functions/functions.R (sourced below).
#
# Required packages: shiny, bslib, vroom, here, dplyr, ggplot2, plotly
#   (xlsx download also needs openxlsx). Install with, e.g.:
#   renv::install(c("vroom", "here", "dplyr", "ggplot2", "plotly"))


# Libraries: --------------------------------------------------------------

library(shiny)
library(bslib)
# Other packages are called with :: (vroom::, dplyr::, ggplot2::, plotly::)


# Data: -------------------------------------------------------------------
# Loaded once at startup with vroom (fast). Columns used downstream:
#   data_proteins      : Gene_name, omic_layer, PTM_collapse_key, Timepoint, Group, Value
#   limma_proteins_T2D : Gene_name, omics_layer, PTM_collapse_key, P.Value_*

data_proteins      <- vroom::vroom(here::here("data/data_files/data_proteins.txt"))
limma_proteins_T2D <- vroom::vroom(here::here("data/limma_outputs/limma_proteins_T2D.txt"))


# Functions: --------------------------------------------------------------
# Defines violin_T2D(). Data is already loaded above, so functions.R skips its
# own (guarded) load.

source(here::here("functions/functions.R"))


# Choices derived from the data: ------------------------------------------

genes        <- sort(unique(data_proteins$Gene_name))
layers       <- sort(unique(data_proteins$omic_layer))  # phospho / proteome / transcriptome
default_gene <- if ("IQGAP1" %in% genes) "IQGAP1" else genes[1]


# Helper - write a data.frame in the chosen download format: ---------------

write_table_by_type <- function(df, file, type) {
  switch(
    type,
    csv  = utils::write.csv(df, file, row.names = FALSE),
    tsv  = utils::write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE),
    xlsx = if (requireNamespace("openxlsx", quietly = TRUE)) {
      openxlsx::write.xlsx(df, file)
    } else {
      utils::write.csv(df, file, row.names = FALSE)  # fallback if openxlsx missing
    }
  )
}


# Theme + custom CSS: -----------------------------------------------------

# TODO: rebrand to your study's colours.
color_palette <- c(
  bg      = "#F4F7FB",  # app background
  primary = "#1B6CA8",  # navbar / accents
  panel   = "#D6E4F0",  # sidebar background
  border  = "#1B6CA8"   # card / sidebar borders
)

app_theme <- bs_theme(
  version = 5,
  bg = color_palette[["bg"]],
  fg = "black",
  primary = color_palette[["primary"]]
) |>
  bs_add_rules(sprintf("
    .navbar-brand { color: white !important; font-weight: bold; font-size: 24px; }
    .nav-link { font-size: 18px; color: white; }

    /* Sidebar */
    .custom-sidebar {
      background-color: %s;
      border: 3px solid %s;
      padding: 15px;
      border-radius: 5px;
    }
    .sidebar-section {
      margin-bottom: 20px;
      padding-bottom: 20px;
      border-bottom: 2px solid #ddd;
    }
    .section-header { font-weight: bold; margin-bottom: 10px; }

    /* Cards holding the plots / tables */
    .card .card-header {
      background-color: %s;
      color: white;
      font-size: 20px;
      font-weight: bold;
      padding: 10px;
    }
  ", color_palette[["panel"]], color_palette[["border"]], color_palette[["primary"]]))

# Shared style string for the plot/table panels
panel_style <- "background-color: #FFFFFF; border: 3px solid #1B6CA8; color: black;"


# Graphical user interface: -----------------------------------------------

ui <- page_navbar(
  title = "MAX app",          # TODO: app title
  bg = color_palette[["primary"]],
  theme = app_theme,

  # --- Home ----------------------------------------------------------------
  nav_panel(
    "Home",
    fluidPage(
      # Citation
      fluidRow(column(12, wellPanel(
        p("If you use our app, please cite our publication:"),
        tags$span(strong("Molecular Pathways of Exercise in Type 2 Diabetes Revealed by Multi-Omics"), style = "font-size: 24px;"),
        p("Ben Stocks*, Stephen P Ashcroft*, Signe Schmidt Kjølner Hansen, Jeppe Kjærgaard,
        Kirstin A MacGregor, Dimitrius Santiago Passos Simões Fróes Guimaraes, David Rizo-Roca,
        Marc Pielies Avelli, Simon Wengert, Konstantinos Makris, Amy M Ehrlich, Scott Frendo-Cumbo,
        Simone Jensen, Mladen Savikj, Roger Moreno-Justicia, Torkil Rogneflåten, Håvard Hamarsland,
        Daniel Hammarström, Dominik Lutter, Julia Otten, Tommy Olsson, Simon Rasmussen, Kenneth Caidahl,
        Harriet Wallberg-Henriksson, Anna Krook, Atul S Deshmukh#, and Juleen R Zierath#"),
        p("*Joint first authors #Co-corresponding authors"),
        p("<Journal / DOI>"),                                 # TODO
        p(strong("Summary")),
        p("<One-paragraph summary of the study and what this app lets users explore.>")  # TODO
      ))),

      # Graphical abstract (drop an image in www/ and point to it here)
      fluidRow(column(12, wellPanel(
        tags$div(
          style = "text-align: center;",
          # TODO: add www/graphical_abstract.png then uncomment:
          # tags$img(src = "graphical_abstract.png",
          #          style = "max-width: 100%; height: auto; max-height: 600px;"),
          tags$em("Add a graphical abstract image to www/ and reference it here.")
        )
      ))),

      # How to use
      fluidRow(column(12, wellPanel(
        strong("Introduction to the data:"),
        p("This app visualises the MAX multi-omics exercise study in people with normal
           glucose tolerance (NGT) and type 2 diabetes (T2D)."),
        p("The 'Expression' tab shows the abundance of a chosen gene / protein across the
           exercise time-course (Base, Post, Rec) in each group, for the transcriptome,
           proteome or phosphoproteome. Significant timepoint comparisons (p ≤ 0.05)
           are annotated with brackets."),
        p("The 'Differential abundance' tab (coming soon) will show volcano plots from the
           limma analysis."),
        p("Disclaimer: <add any preprint / preliminary-results disclaimer here>.")  # TODO
      ))),

      # Contact
      fluidRow(column(12, wellPanel(
        p("Questions or trouble using the tool? Contact us:"),
        p(tags$a(href = "mailto:roger.moreno.justicia@sund.ku.dk",
                 "roger.moreno.justicia@sund.ku.dk"))         # TODO: confirm contact(s)
      ))),

      # Institutional / funder logos
      fluidRow(column(12,
        # TODO: add logo PNGs to www/ and render via plotOutput() (see server).
        tags$em("Add institutional / funder logos to www/ and wire them up in the server.")
      ))
    )
  ),

  # --- Expression (violin_T2D) ---------------------------------------------
  nav_panel(
    "Expression",
    sidebarLayout(
      sidebarPanel(
        class = "custom-sidebar",
        div(
          class = "sidebar-section",
          p("Select a gene / protein and an omics layer to view its abundance across the
             exercise time-course."),
          selectizeInput("feature", "Gene / protein:", choices = NULL),
          radioButtons(
            "omics_layer", "Omics layer:",
            choices = layers, selected = if ("proteome" %in% layers) "proteome" else layers[1]
          ),
          # Phosphosite picker - only shown for the phosphoproteome layer.
          conditionalPanel(
            condition = "input.omics_layer == 'phosphoproteome'",
            selectizeInput("phosphosite", "Phosphosite:", choices = NULL)
          )
        ),
        div(
          class = "sidebar-section",
          radioButtons(
            "filetype", "Download the underlying data as:",
            choices = c("csv", "tsv", "xlsx"), selected = "csv"
          ),
          downloadButton("download_data", "Download data", class = "btn-primary")
        )
      ),
      mainPanel(
        navset_card_underline(
          full_screen = TRUE, title = "Expression",
          nav_panel(
            "Ins.Sensitivity",
            plotly::plotlyOutput("violin", height = "520px"),
            style = panel_style
          )
        )
      )
    )
  ),

  # --- Differential abundance (placeholder) --------------------------------
  nav_panel(
    "Differential abundance",
    fluidPage(fluidRow(column(12, wellPanel(
      h4("Differential abundance"),
      p("Volcano plots from the limma results (limma_proteins_T2D) will live here."),
      tags$em("Coming soon.")
    ))))
  )
)


# Server: -----------------------------------------------------------------

server <- function(input, output, session) {

  # --- Institutional logos (Home) -----------------------------------------
  # TODO: add PNGs to www/ then uncomment and add plotOutput("logo_*") in the UI.
  # output$logo_CBMR <- renderPlot({
  #   grid::grid.raster(png::readPNG("www/logo_CBMR.png"))
  # }, bg = "transparent")

  # --- Expression tab ------------------------------------------------------

  # 21k+ genes -> populate the selectize on the server for performance.
  updateSelectizeInput(session, "feature", choices = genes,
                       selected = default_gene, server = TRUE)

  # Phosphosites depend on the chosen gene; refresh them when gene / layer change.
  observeEvent(list(input$feature, input$omics_layer), {
    req(input$feature)
    if (identical(input$omics_layer, "phosphoproteome")) {
      sites <- data_proteins |>
        dplyr::filter(Gene_name == input$feature, omic_layer == "phosphoproteome") |>
        dplyr::pull(PTM_collapse_key) |>
        unique() |>
        sort()
      updateSelectizeInput(session, "phosphosite", choices = sites,
                           selected = if (length(sites)) sites[1] else NULL)
    }
  })

  # Data underlying the current selection (used by the download handler).
  selection <- reactive({
    req(input$feature, input$omics_layer)
    d <- data_proteins |>
      dplyr::filter(Gene_name == input$feature, omic_layer == input$omics_layer)
    if (identical(input$omics_layer, "phosphoproteome") &&
        isTRUE(nzchar(input$phosphosite))) {
      d <- d |> dplyr::filter(PTM_collapse_key == input$phosphosite)
    }
    d
  })

  # Main output: the violin from violin_T2D(), made interactive with ggplotly().
  output$violin <- plotly::renderPlotly({
    req(input$feature, input$omics_layer)

    # Guard the combinations violin_T2D() cannot handle, with a friendly message.
    validate(need(
      nrow(dplyr::filter(data_proteins,
                         Gene_name == input$feature,
                         omic_layer == input$omics_layer)) > 0,
      "No data for this gene in the selected omics layer."
    ))

    site <- NULL
    if (identical(input$omics_layer, "phosphoproteome")) {
      validate(need(isTRUE(nzchar(input$phosphosite)), "Select a phosphosite."))
      site <- input$phosphosite
    }

    p <- violin_T2D(data_proteins, limma_proteins_T2D,
                    feature = input$feature, omics_layer = input$omics_layer,
                    phosphosite = site)
    plotly::ggplotly(p, tooltip = "text")
  })

  # Download the underlying data for the current selection.
  output$download_data <- downloadHandler(
    filename = function() {
      base <- input$feature
      if (identical(input$omics_layer, "phosphoproteome") && isTRUE(nzchar(input$phosphosite))) {
        base <- input$phosphosite
      }
      paste0("MAX_", input$omics_layer, "_", base, ".", input$filetype)
    },
    content = function(file) write_table_by_type(selection(), file, input$filetype)
  )
}


# Run the application: ----------------------------------------------------

shinyApp(ui = ui, server = server)
