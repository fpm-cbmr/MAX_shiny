# SHINY APP FOR THE MAX MULTI-OMICS STUDY: --------------------------------
#
# page_navbar + bslib theme + custom CSS. Data is loaded once here with vroom;
# the plotting logic lives in functions/functions.R (sourced below).
#
# Required packages: shiny, bslib, data.table, vroom, here, dplyr, ggplot2, plotly, DT
#   (xlsx download also needs openxlsx). Install with, e.g.:
#   renv::install(c("data.table", "vroom", "here", "dplyr", "ggplot2", "plotly", "DT"))


# Libraries: --------------------------------------------------------------

library(shiny)
library(bslib)
library(data.table)   # fast load (fread) + keyed lookups for the large protein table
# Other packages are called with :: (vroom::, dplyr::, ggplot2::, plotly::)


# Data: -------------------------------------------------------------------
# Loaded once at startup with vroom (fast). Columns used downstream:
#   data_proteins      : Gene_name, omic_layer, PTM_collapse_key, Timepoint, Group, Value
#   limma_proteins_T2D : Gene_name, omics_layer, PTM_collapse_key, P.Value_*

# fread (~2.8 s vs ~9.6 s for vroom) + key on Gene_name + omic_layer, so the
# per-selection lookup in protein_subset() is a ~0 ms binary search.
data_proteins      <- data.table::fread(here::here("data/data_files/data_proteins.txt"))
data.table::setkey(data_proteins, Gene_name, omic_layer)
limma_proteins_T2D <- vroom::vroom(here::here("data/limma_outputs/limma_proteins_T2D.txt"))
limma_proteins_sex <- vroom::vroom(here::here("data/limma_outputs/limma_proteins_sex.txt"))

# Metabolites (single feature type: CHEMICAL_NAME; no omics layer / phosphosites).
data_metabolites      <- vroom::vroom(here::here("data/data_files/data_metabolites.txt"))
limma_metabolites_TD  <- vroom::vroom(here::here("data/limma_outputs/limma_metabolites_TD.txt"))
limma_metabolites_sex <- vroom::vroom(here::here("data/limma_outputs/limma_metabolites_sex.txt"))


# Functions: --------------------------------------------------------------
# Defines violin_T2D(). Data is already loaded above, so functions.R skips its
# own (guarded) load.

source(here::here("functions/functions.R"))


# Choices derived from the data: ------------------------------------------

genes        <- sort(unique(data_proteins$Gene_name))
layers       <- sort(unique(data_proteins$omic_layer))  # phospho / proteome / transcriptome
default_gene <- if ("IQGAP1" %in% genes) "IQGAP1" else genes[1]

metabolites        <- sort(unique(data_metabolites$CHEMICAL_NAME))
default_metabolite <- metabolites[1]


# Theme + custom CSS: -----------------------------------------------------

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

    /* Card tabs (navset_card_underline): Main effect / Ins.Sensitivity / Sex.
       These render as a classic tabset (ul.nav.nav-underline > li > a) with no
       .nav-link class, so the anchors are targeted directly. Force white on the
       blue card header; active / hovered tab is full-opacity + bold. */
    .card-header .nav > li > a,
    .card-header .nav-link {
      color: #ffffff !important;
      opacity: 0.7;
    }
    .card-header .nav > li.active > a,
    .card-header .nav > li > a:hover,
    .card-header .nav-link.active,
    .card-header .nav-link:hover {
      color: #ffffff !important;
      opacity: 1;
      font-weight: 700;
    }
  ", color_palette[["panel"]], color_palette[["border"]], color_palette[["primary"]]))

# Shared style string for the plot/table panels
panel_style <- "background-color: #FFFFFF; border: 3px solid #1B6CA8; color: black;"


# Graphical user interface: -----------------------------------------------

ui <- page_navbar(
  title = "MAX app",
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
    "Genes, Proteins & Phosphosites",
    sidebarLayout(
      sidebarPanel(
        class = "custom-sidebar",
        div(
          class = "sidebar-section",
          p("Select a gene / protein and an omics layer to view its abundance across the
             exercise time-course."),
          p("Then, toggle between Main effect, Ins.Sensitivity or Sex to view feature abundance
            across different statistical comparisons"),
          selectizeInput("feature", "Gene / protein:", choices = NULL),
          radioButtons(
            "omics_layer", "Omics layer:",
            choices = layers, selected = if ("proteome" %in% layers) "proteome" else layers[1]
          ),
          # Phosphosite picker - only shown for the phosphoproteome layer.
          conditionalPanel(
            condition = "input.omics_layer == 'phosphoproteome'",
            selectizeInput("phosphosite", "Phosphosite:", choices = NULL)
          ),
          p("S, T and Y stand for serine, threonine and tyrosine, respectively"),
          p("M indicates multiplicity of sites: number of concurrent phosphosites in a given peptide")
        ),
        div(
          class = "sidebar-section",
          radioButtons(
            "filetype", "Download the results table as:",
            choices = c("csv", "tsv", "xlsx"), selected = "csv"
          ),
          downloadButton("download_data", "Download table", class = "btn-primary")
        )
      ),
      mainPanel(
        # Report the active card tab to the server. bslib's navset `id` binding
        # can miss Bootstrap 5 tab-change events, so set a Shiny input directly
        # from the native 'shown.bs.tab' event. Generic: any navset with an id
        # gets input$<id>_active (here #comparison and #comparison_met).
        tags$script(HTML(
          "document.addEventListener('shown.bs.tab', function(e) {
             var nav = e.target.closest('ul.nav[id]');
             if (nav) Shiny.setInputValue(nav.id + '_active', e.target.getAttribute('data-value'));
           });"
        )),
        navset_card_underline(
          id = "comparison",   # kept as a fallback source for the active tab
          full_screen = TRUE, title = "Feature Expression",
          nav_panel(
            "Main effect",
            plotly::plotlyOutput("violin_main", height = "45vh")
          ),
          nav_panel(
            "Type 2 Diabetes",
            plotly::plotlyOutput("violin_T2D", height = "45vh")
          ),
          nav_panel(
            "Sex",
            plotly::plotlyOutput("violin_sex", height = "45vh")
          )
        ),
        # Limma results for the selected feature, linked to the active comparison.
        card(
          full_screen = TRUE,
          card_header("Statistics"),
          card_body(
            DT::DTOutput("limma_table"),
            tags$small(
              class = "text-muted",
              "Limma statistics for the selected feature across omics layers in the
               active comparison. Click a transcriptome / proteome row to plot the
               gene, or a phosphosite row to plot that site."
            )
          )
        )
      )
    )
  ),

  # --- Metabolites ---------------------------------------------------------
  nav_panel(
    "Metabolites",
    sidebarLayout(
      sidebarPanel(
        class = "custom-sidebar",
        div(
          class = "sidebar-section",
          p("Select a metabolite to view its abundance across the exercise time-course."),
          p("Then, toggle between Main effect, Type 2 Diabetes or Sex to view it across
             different statistical comparisons."),
          selectizeInput("metabolite", "Metabolite:", choices = NULL)
        ),
        div(
          class = "sidebar-section",
          radioButtons(
            "filetype_met", "Download the results table as:",
            choices = c("csv", "tsv", "xlsx"), selected = "csv"
          ),
          downloadButton("download_met", "Download table", class = "btn-primary")
        )
      ),
      mainPanel(
        navset_card_underline(
          id = "comparison_met",   # input$comparison_met_active via the JS handler
          full_screen = TRUE, title = "Metabolite levels",
          nav_panel("Main effect",     plotly::plotlyOutput("violin_met_main", height = "45vh")),
          nav_panel("Type 2 Diabetes", plotly::plotlyOutput("violin_met_T2D",  height = "45vh")),
          nav_panel("Sex",             plotly::plotlyOutput("violin_met_sex",  height = "45vh"))
        ),
        # Limma results for the selected metabolite, linked to the active comparison.
        card(
          full_screen = TRUE,
          card_header("Statistics"),
          card_body(
            DT::DTOutput("limma_met_table"),
            tags$small(
              class = "text-muted",
              "Limma statistics for the selected metabolite in the active comparison."
            )
          )
        )
      )
    )
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

  # A phosphosite requested by a table row click, consumed once the phosphosite
  # choices have been refreshed for the gene (so it survives a layer switch).
  pending_site <- reactiveVal(NULL)

  # Phosphosites depend on the chosen gene; refresh them when gene / layer change.
  observeEvent(list(input$feature, input$omics_layer), {
    req(input$feature)
    if (identical(input$omics_layer, "phosphoproteome")) {
      sites <- sort(unique(
        data_proteins[.(input$feature, "phosphoproteome"), nomatch = NULL]$PTM_collapse_key
      ))
      # honour a site requested via a row click, else default to the first
      sel <- if (!is.null(pending_site()) && pending_site() %in% sites) pending_site()
             else if (length(sites)) sites[1] else NULL
      pending_site(NULL)
      updateSelectizeInput(session, "phosphosite", choices = sites, selected = sel)
    }
  })

  # Filter the 11M-row data_proteins ONCE per gene + layer change (keyed
  # data.table lookup, ~0 ms) and share it across the three plots + their
  # validate guards. Switching comparison tabs reuses this cached subset instead
  # of re-filtering; the phosphosite is narrowed downstream (in the functions).
  protein_subset <- reactive({
    req(input$feature, input$omics_layer)
    as.data.frame(data_proteins[.(input$feature, input$omics_layer), nomatch = NULL])
  })

  # Main output: the violin from violin_main_effect(), made interactive with ggplotly().
  output$violin_main <- plotly::renderPlotly({
    req(input$feature, input$omics_layer)
    d <- protein_subset()
    validate(need(nrow(d) > 0, "No data for this gene in the selected omics layer."))

    site <- NULL
    if (identical(input$omics_layer, "phosphoproteome")) {
      validate(need(isTRUE(nzchar(input$phosphosite)), "Select a phosphosite."))
      site <- input$phosphosite
    }

    p <- violin_main_effect(d, limma_proteins_T2D,
                    feature = input$feature, omics_layer = input$omics_layer,
                    phosphosite = site)
    plotly::ggplotly(p, tooltip = "text")
  })
  
  
  # Main output: the violin from violin_T2D(), made interactive with ggplotly().
  output$violin_T2D <- plotly::renderPlotly({
    req(input$feature, input$omics_layer)
    d <- protein_subset()
    validate(need(nrow(d) > 0, "No data for this gene in the selected omics layer."))

    site <- NULL
    if (identical(input$omics_layer, "phosphoproteome")) {
      validate(need(isTRUE(nzchar(input$phosphosite)), "Select a phosphosite."))
      site <- input$phosphosite
    }

    p <- violin_T2D(d, limma_proteins_T2D,
                    feature = input$feature, omics_layer = input$omics_layer,
                    phosphosite = site)
    plotly::ggplotly(p, tooltip = "text")
  })
  
  # Main output: the violin from violin_sex(), made interactive with ggplotly().
  output$violin_sex <- plotly::renderPlotly({
    req(input$feature, input$omics_layer)
    d <- protein_subset()
    validate(need(nrow(d) > 0, "No data for this gene in the selected omics layer."))

    site <- NULL
    if (identical(input$omics_layer, "phosphoproteome")) {
      validate(need(isTRUE(nzchar(input$phosphosite)), "Select a phosphosite."))
      site <- input$phosphosite
    }

    p <- violin_sex(d, limma_proteins_sex,
                    feature = input$feature, omics_layer = input$omics_layer,
                    phosphosite = site)
    plotly::ggplotly(p, tooltip = "text")
  })

  # --- Limma results table (linked to the active comparison) ---------------

  # Which comparison tab is showing (defaults before the navset registers).
  active_comparison <- reactive({
    # Prefer the value set by the JS 'shown.bs.tab' handler (reliable on BS5),
    # fall back to the navset id binding, then to the first tab.
    cmp <- input$comparison_active
    if (is.null(cmp) || !nzchar(cmp)) cmp <- input$comparison
    if (is.null(cmp) || !nzchar(cmp)) "Main effect" else cmp
  })

  # Table data: stats for the selected feature in the active comparison. Depends
  # on feature / layer / comparison only (NOT the phosphosite), so clicking a row
  # to change the site does not rebuild the table -> no reactive loop.
  # Table data: the gene's measurements across layers for the active comparison.
  # Depends on feature + comparison only (NOT omics_layer / phosphosite), so
  # clicking a row to change the layer or site does not rebuild it -> no loop.
  limma_view <- reactive({
    req(input$feature)
    limma_results_table(active_comparison(), input$feature,
                        limma_proteins_T2D, limma_proteins_sex)
  })

  # server = FALSE (client-side) so the table re-renders cleanly when the active
  # comparison changes the column set.
  output$limma_table <- DT::renderDT({
    df <- limma_view()
    validate(need(!is.null(df) && nrow(df) > 0, "No limma results for this feature."))
    pcols <- names(df)[grepl("P.Val", names(df), fixed = TRUE)]                  # P.Value + adj.P.Val
    rcols <- names(df)[grepl("logFC", names(df), fixed = TRUE) |
                       grepl("AveExpr", names(df), fixed = TRUE)]
    dt <- DT::datatable(
      df, rownames = FALSE, selection = "single",
      # Bounded, viewport-relative height with internal scroll so the table card
      # fits on screen alongside the plot (no page scrolling). scrollCollapse
      # shrinks it when there are few rows.
      options = list(scrollX = TRUE, scrollY = "28vh", scrollCollapse = TRUE,
                     paging = FALSE, dom = "ft")
    )
    if (length(pcols)) dt <- DT::formatSignif(dt, pcols, 3)
    if (length(rcols)) dt <- DT::formatRound(dt, rcols, 2)
    dt
  }, server = FALSE)

  # Click a row -> switch the plot to that layer (and site). A transcriptome /
  # proteome row sets the omics layer for the gene; a phosphosite row sets the
  # layer + that site.
  observeEvent(input$limma_table_rows_selected, {
    df <- limma_view(); req(df)
    i     <- input$limma_table_rows_selected
    layer <- df$Layer[i]
    feat  <- df$Feature[i]
    if (identical(layer, "phosphoproteome")) {
      if (identical(input$omics_layer, "phosphoproteome")) {
        if (!identical(feat, input$phosphosite)) {
          updateSelectizeInput(session, "phosphosite", selected = feat)
        }
      } else {
        pending_site(feat)                                  # carry site across the layer switch
        updateRadioButtons(session, "omics_layer", selected = "phosphoproteome")
      }
    } else if (!identical(input$omics_layer, layer)) {
      updateRadioButtons(session, "omics_layer", selected = layer)
    }
  })

  # Keep the table's highlighted row in sync with the sidebar selection (the
  # reverse of the row-click linking), so the highlighted row always matches the
  # plotted layer / site. Guarded against re-selecting the current row to avoid
  # ping-ponging with the row-click observer.
  limma_proxy <- DT::dataTableProxy("limma_table")
  observeEvent(list(limma_view(), input$omics_layer, input$phosphosite), {
    df <- limma_view(); req(df)
    target <- if (identical(input$omics_layer, "phosphoproteome")) {
      which(df$Layer == "phosphoproteome" & df$Feature == input$phosphosite)
    } else {
      which(df$Layer == input$omics_layer)
    }
    sel <- isolate(input$limma_table_rows_selected)   # read, don't depend on, the selection
    if (length(target) == 1) {
      if (!identical(as.integer(target), as.integer(sel))) DT::selectRows(limma_proxy, target)
    } else if (!is.null(sel)) {
      DT::selectRows(limma_proxy, NULL)               # no matching row -> clear the highlight
    }
  })

  # Download the limma results table currently on screen (selected feature,
  # active comparison).
  output$download_data <- downloadHandler(
    filename = function() {
      cmp <- gsub("[^A-Za-z0-9]", "", active_comparison())
      paste0("MAX_limma_", input$feature, "_", cmp, ".", input$filetype)
    },
    content = function(file) {
      df <- limma_view()
      if (is.null(df) || !nrow(df)) df <- data.frame(Message = "No results for this feature.")
      write_table_by_type(df, file, input$filetype)
    }
  )

  # --- Metabolites tab ------------------------------------------------------

  # 846 metabolites -> populate the selectize on the server.
  updateSelectizeInput(session, "metabolite", choices = metabolites,
                       selected = default_metabolite, server = TRUE)

  active_comparison_met <- reactive({
    cmp <- input$comparison_met_active
    if (is.null(cmp) || !nzchar(cmp)) "Main effect" else cmp
  })

  # One plotly builder reused for the three comparison tabs.
  met_plot <- function(comparison, limma_tbl) {
    req(input$metabolite)
    validate(need(
      nrow(dplyr::filter(data_metabolites, CHEMICAL_NAME == input$metabolite)) > 0,
      "No data for this metabolite."
    ))
    plotly::ggplotly(
      violin_metabolite(data_metabolites, limma_tbl, input$metabolite, comparison),
      tooltip = "text"
    )
  }
  output$violin_met_main <- plotly::renderPlotly(met_plot("Main effect",     limma_metabolites_TD))
  output$violin_met_T2D  <- plotly::renderPlotly(met_plot("Type 2 Diabetes", limma_metabolites_TD))
  output$violin_met_sex  <- plotly::renderPlotly(met_plot("Sex",             limma_metabolites_sex))

  # Stats table for the selected metabolite, linked to the active comparison.
  limma_met_view <- reactive({
    req(input$metabolite)
    limma_metabolites_table(active_comparison_met(), input$metabolite,
                            limma_metabolites_TD, limma_metabolites_sex)
  })

  output$limma_met_table <- DT::renderDT({
    df <- limma_met_view()
    validate(need(!is.null(df) && nrow(df) > 0, "No limma results for this metabolite."))
    pcols <- names(df)[grepl("P.Val", names(df), fixed = TRUE)]
    rcols <- names(df)[grepl("logFC", names(df), fixed = TRUE) |
                       grepl("AveExpr", names(df), fixed = TRUE)]
    dt <- DT::datatable(df, rownames = FALSE, selection = "none",
                        options = list(scrollX = TRUE, dom = "t", paging = FALSE))
    if (length(pcols)) dt <- DT::formatSignif(dt, pcols, 3)
    if (length(rcols)) dt <- DT::formatRound(dt, rcols, 2)
    dt
  }, server = FALSE)

  output$download_met <- downloadHandler(
    filename = function() {
      cmp <- gsub("[^A-Za-z0-9]", "", active_comparison_met())
      met <- gsub("[^A-Za-z0-9]+", "_", input$metabolite)
      paste0("MAX_metabolite_", met, "_", cmp, ".", input$filetype_met)
    },
    content = function(file) {
      df <- limma_met_view()
      if (is.null(df) || !nrow(df)) df <- data.frame(Message = "No results for this metabolite.")
      write_table_by_type(df, file, input$filetype_met)
    }
  )
}


# Run the application: ----------------------------------------------------

shinyApp(ui = ui, server = server)
