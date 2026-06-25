# SHINY APP FOR <MAX STUDY>: ----------------------------------------------
#
# Skeleton modelled on the shiny_exerome app. It is a *runnable* scaffold:
# it ships with small synthetic demo data so you can launch it immediately,
# then replace the demo data + helper functions with your own.
#
# Structure (mirrors shiny_exerome):
#   - page_navbar + bslib theme + custom CSS
#   - "Home" tab  : study description / citation / contact / logos
#   - data tabs   : sidebarLayout(sidebar = search + indicators + downloads,
#                                  main    = navset_card_underline(plot, table))
#
# Recommended next step: move "Demo data" and "Helper functions" into
#   R/loading_data.R  and  R/functions.R  and source() them (see shiny_exerome),
#   so app.R stays focused on UI + server wiring.
#
# Packages used: shiny, bslib (already in renv.lock) + plotly, ggplot2, dplyr,
#   kableExtra (install these). xlsx download additionally needs openxlsx.
#   Install with, e.g.:
#   renv::install(c("plotly", "ggplot2", "dplyr", "kableExtra", "openxlsx"))


# Libraries: --------------------------------------------------------------

library(shiny)
library(bslib)
# Remaining packages are called with :: (plotly::, ggplot2::, dplyr::, kableExtra::)


# Options / palette: ------------------------------------------------------

# TODO: rebrand to your study's colours.
color_palette <- c(
  bg      = "#F4F7FB",  # app background
  primary = "#1B6CA8",  # navbar / accents
  panel   = "#D6E4F0",  # sidebar background
  border  = "#1B6CA8"   # card / sidebar borders
)


# Demo data: --------------------------------------------------------------
# TODO: DELETE this block and load your real data instead, e.g.:
#   load(here::here("data/abundance.rda"))
#   load(here::here("data/de_results.rda"))
# Keep the column names used below (protein / grouping / abundance, and
# Genes / contrast / logFC / P.Value / adj.P.Val) or update the helpers.

set.seed(1)

demo_groups   <- c("WT", "KO_g2", "KO_g5")
demo_proteins <- c("Ybx3", "Ybx1", "Col1a1", "Tppp3", "Actb",
                   "Gapdh", "Vim", "Fn1", "Hspa8", "Eef2")

# Per-sample protein abundance (3 replicates per group)
demo_abundance <- do.call(rbind, lapply(demo_proteins, function(p) {
  baseline <- runif(1, 18, 24)
  do.call(rbind, lapply(demo_groups, function(g) {
    ko_shift <- if (p == "Ybx3" && grepl("KO", g)) -3 else 0
    data.frame(
      protein   = p,
      grouping  = factor(g, levels = demo_groups),
      sample_id = paste0(g, "_", 1:3),
      abundance = round(rnorm(3, baseline + ko_shift, 0.4), 3),
      stringsAsFactors = FALSE
    )
  }))
}))

# Differential-abundance results for two example contrasts
demo_contrasts <- c("KO_g2 vs WT", "KO_g5 vs WT")
demo_universe  <- c(demo_proteins, paste0("Gene", sprintf("%03d", 1:240)))
demo_de <- do.call(rbind, lapply(demo_contrasts, function(ct) {
  n   <- length(demo_universe)
  lfc <- rnorm(n, 0, 1)
  pv  <- runif(n)^2
  # give a couple of proteins a clear signal so the volcano looks real
  hit <- match(c("Ybx3", "Col1a1"), demo_universe)
  lfc[hit] <- c(-3.2, 1.8)
  pv[hit]  <- c(1e-5, 5e-4)
  data.frame(
    Genes     = demo_universe,
    contrast  = ct,
    logFC     = round(lfc, 3),
    P.Value   = signif(pv, 3),
    adj.P.Val = signif(p.adjust(pv, "BH"), 3),
    stringsAsFactors = FALSE
  )
}))

# Choices for the search boxes
proteins    <- sort(unique(demo_abundance$protein))
de_proteins <- sort(unique(demo_de$Genes))


# Helper functions: -------------------------------------------------------
# TODO: replace plot bodies with your own ggplot code (these mirror the
# boxplot + volcano style in YBX3_KO/R/data_analysis.R).

# Small coloured status circle, reused by the sidebar indicators
create_indicator <- function(label, ok) {
  color <- if (isTRUE(ok)) "green" else "red"
  HTML(sprintf(
    paste0('<div style="display:flex;align-items:center;margin-bottom:5px;">',
           '<div style="width:12px;height:12px;background-color:%s;',
           'border-radius:50%%;margin-right:8px;"></div>',
           '<span>%s</span></div>'),
    color, label
  ))
}

# Styled HTML table (kableExtra) reused by every table output
make_kable <- function(df) {
  kableExtra::kable(df, format = "html") |>
    kableExtra::kable_styling(
      full_width = TRUE,
      bootstrap_options = c("striped", "hover"),
      position = "center"
    ) |>
    kableExtra::row_spec(0, bold = TRUE, font_size = 20) |>
    kableExtra::row_spec(seq_len(nrow(df)), font_size = 16)
}

# Write a data.frame in the format chosen by the radio buttons
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

# Boxplot of one protein's abundance across groups
plot_abundance <- function(.data, .protein) {
  d <- .data |> dplyr::filter(protein == .protein)

  ggplot2::ggplot(d, ggplot2::aes(x = grouping, y = abundance, fill = grouping)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::labs(title = .protein, x = NULL, y = "Abundance (log2)") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
}

# Volcano plot for one contrast, optionally highlighting one protein
plot_volcano <- function(.data, .contrast, .highlight = NULL) {
  d <- .data |>
    dplyr::filter(contrast == .contrast) |>
    dplyr::mutate(
      significant = dplyr::case_when(
        adj.P.Val <= 0.05 & logFC > 0 ~ "Up",
        adj.P.Val <= 0.05 & logFC < 0 ~ "Down",
        TRUE ~ "n.s."
      )
    )
  hl <- d |> dplyr::filter(Genes %in% .highlight)

  ggplot2::ggplot(
    d, ggplot2::aes(x = logFC, y = -log10(P.Value), color = significant, text = Genes)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1.4) +
    ggplot2::geom_point(data = hl, color = "black", size = 2.6) +
    ggplot2::geom_text(
      data = hl, ggplot2::aes(label = Genes),
      color = "black", vjust = -1, size = 4, show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c(Down = "#2c7fb8", n.s. = "grey80", Up = "#d7301f")
    ) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::labs(title = .contrast, x = "log2 fold-change",
                  y = "-log10(P-value)", color = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
}


# Theme + custom CSS: -----------------------------------------------------

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
        p("This app visualises <describe your dataset>."),                       # TODO
        p("The 'Protein abundance' tab shows per-group abundance for a chosen protein."),
        p("The 'Differential abundance' tab shows volcano plots for each contrast."),
        p("Disclaimer: <add any preprint / preliminary-results disclaimer here>.")
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

  # --- Protein abundance ---------------------------------------------------
  nav_panel(
    "Protein abundance",
    sidebarLayout(
      sidebarPanel(
        class = "custom-sidebar",
        div(
          class = "sidebar-section",
          p("Search a protein to view its abundance across groups."),
          selectizeInput(
            "protein", "Browse a protein:",
            choices = proteins, selected = "Ybx3", multiple = FALSE
          ),
          p("Selected protein is:"),
          uiOutput("detect_indicator")
        ),
        div(
          class = "sidebar-section",
          radioButtons(
            "filetype_abundance",
            "Select filetype and download the abundance table:",
            choices = c("csv", "tsv", "xlsx"), selected = "csv"
          ),
          downloadButton("download_abundance", "Download data",
                         class = "btn-primary")
        )
      ),
      mainPanel(
        navset_card_underline(
          full_screen = TRUE, title = "Plot",
          nav_panel("Abundance", plotly::plotlyOutput("abundance_plot"),
                    height = "500px", style = panel_style)
        ),
        navset_card_underline(
          full_screen = TRUE, title = "Table",
          nav_panel("Summary", tableOutput("abundance_table"), style = panel_style)
        )
      )
    )
  ),

  # --- Differential abundance ----------------------------------------------
  nav_panel(
    "Differential abundance",
    sidebarLayout(
      sidebarPanel(
        class = "custom-sidebar",
        div(
          class = "sidebar-section",
          radioButtons(
            "de_contrast", "Contrast:",
            choices = demo_contrasts, selected = demo_contrasts[1]
          )
        ),
        div(
          class = "sidebar-section",
          p("Highlight a protein on the volcano plot:"),
          selectizeInput(
            "de_protein", "Browse a protein:",
            choices = de_proteins, selected = "Ybx3", multiple = FALSE
          ),
          p("In the selected contrast this protein is:"),
          uiOutput("sig_indicator")
        ),
        div(
          class = "sidebar-section",
          radioButtons(
            "filetype_de",
            "Select filetype and download the results table:",
            choices = c("csv", "tsv", "xlsx"), selected = "csv"
          ),
          downloadButton("download_de", "Download data", class = "btn-primary")
        )
      ),
      mainPanel(
        navset_card_underline(
          full_screen = TRUE, title = "Plot",
          nav_panel("Volcano", plotly::plotlyOutput("volcano_plot"),
                    height = "500px", style = panel_style)
        ),
        navset_card_underline(
          full_screen = TRUE, title = "Table",
          nav_panel("Top hits", tableOutput("de_table"), style = panel_style)
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

  # --- Protein abundance ---------------------------------------------------
  output$detect_indicator <- renderUI({
    req(input$protein)
    create_indicator("Detected in dataset", input$protein %in% demo_abundance$protein)
  })

  output$abundance_plot <- plotly::renderPlotly({
    req(input$protein)
    plotly::ggplotly(plot_abundance(demo_abundance, input$protein))
  })

  output$abundance_table <- function() {
    req(input$protein)
    d <- demo_abundance |>
      dplyr::filter(protein == input$protein) |>
      dplyr::group_by(grouping) |>
      dplyr::summarise(
        n    = dplyr::n(),
        mean = round(mean(abundance), 3),
        sd   = round(stats::sd(abundance), 3),
        .groups = "drop"
      )
    make_kable(d)
  }

  output$download_abundance <- downloadHandler(
    filename = function() paste0("MAX_abundance_", input$protein, ".", input$filetype_abundance),
    content  = function(file) {
      d <- demo_abundance |> dplyr::filter(protein == input$protein)
      write_table_by_type(d, file, input$filetype_abundance)
    }
  )

  # --- Differential abundance ----------------------------------------------
  output$sig_indicator <- renderUI({
    req(input$de_protein, input$de_contrast)
    row <- demo_de |>
      dplyr::filter(contrast == input$de_contrast, Genes == input$de_protein)
    sig <- nrow(row) > 0 && row$adj.P.Val[1] <= 0.05
    dir <- if (nrow(row) > 0 && row$logFC[1] > 0) "up" else "down"
    tagList(
      create_indicator(
        if (sig) sprintf("Significant (%s-regulated, adj.P ≤ 0.05)", dir)
        else "Not significant (adj.P > 0.05)",
        sig
      )
    )
  })

  output$volcano_plot <- plotly::renderPlotly({
    req(input$de_contrast)
    plotly::ggplotly(
      plot_volcano(demo_de, input$de_contrast, input$de_protein),
      tooltip = c("text", "x", "y")
    )
  })

  output$de_table <- function() {
    req(input$de_contrast)
    d <- demo_de |>
      dplyr::filter(contrast == input$de_contrast) |>
      dplyr::arrange(P.Value) |>
      head(50)
    make_kable(d)
  }

  output$download_de <- downloadHandler(
    filename = function() {
      paste0("MAX_DE_", gsub(" ", "_", input$de_contrast), ".", input$filetype_de)
    },
    content = function(file) {
      d <- demo_de |>
        dplyr::filter(contrast == input$de_contrast) |>
        dplyr::arrange(P.Value)
      write_table_by_type(d, file, input$filetype_de)
    }
  )
}


# Run the application: ----------------------------------------------------

shinyApp(ui = ui, server = server)
