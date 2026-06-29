

# Data is normally loaded by app.R. These guarded loads let this file also be
# sourced on its own (e.g. interactive testing) without reloading the objects
# when they already exist.
if (!exists("data_proteins")) {
  data_proteins <- vroom::vroom(here::here("data/data_files/data_proteins.txt"))
}
if (!exists("limma_proteins_T2D")) {
  limma_proteins_T2D <- vroom::vroom(here::here("data/limma_outputs/limma_proteins_T2D.txt"))
}
if (!exists("limma_proteins_sex")) {
  limma_proteins_sex <- vroom::vroom(here::here("data/limma_outputs/limma_proteins_sex.txt"))
}
if (!exists("data_metabolites")) {
  data_metabolites <- vroom::vroom(here::here("data/data_files/data_metabolites.txt"))
}
if (!exists("limma_metabolites_TD")) {
  limma_metabolites_TD <- vroom::vroom(here::here("data/limma_outputs/limma_metabolites_TD.txt"))
}
if (!exists("limma_metabolites_sex")) {
  limma_metabolites_sex <- vroom::vroom(here::here("data/limma_outputs/limma_metabolites_sex.txt"))
}

# Helper for downloading data:

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

#' Violins for T2D comparisons
#'
#' @param data Long abundance table (e.g. `data_proteins`) with columns
#'   `Gene_name`, `omic_layer`, `PTM_collapse_key`, `Timepoint`, `Group`, `Value`.
#' @param limma_output limma results table with `Gene_name`, `omics_layer`,
#'   `PTM_collapse_key` and the `P.Value_*` columns.
#' @param feature gene, or protein (matched against `Gene_name`)
#' @param omics_layer one of "transcriptome", "proteome" or "phosphoproteome"
#' @param phosphosite Phosphosite to plot, given as a `PTM_collapse_key` value
#'   (e.g. "AAAS_S174_M1"). Only used when `omics_layer == "phosphoproteome"`
#'   (required there); ignored for the other layers.
#' @param base_size Base font size (pt). The default suits the interactive app;
#'   pass a smaller value (e.g. 6) for small static figures.
#'
#' @returns A ggplot object.
#' @export
#'
#' @examples
#' violin_T2D(data_proteins, limma_proteins_T2D, "IQGAP1", "proteome")
#' violin_T2D(data_proteins, limma_proteins_T2D, "AAAS", "phosphoproteome",
#'            phosphosite = "AAAS_S174_M1")
#'            
violin_T2D <- function(data, limma_output, feature, omics_layer, phosphosite = NULL,
                       base_size = 11) {

# Select the feature and omics layer:

selected_data <- data |> 
  dplyr::filter(
    Gene_name == feature
  ) |> 
  dplyr::filter(
    omic_layer == omics_layer
  )

# Extract the p-values (single limma row for this feature + omics layer).

limma_row <- limma_output |>
  dplyr::filter(Gene_name == feature, omics_layer == .env$omics_layer)

# Phosphoproteome only: narrow the gene down to a single phosphosite. The
# `phosphosite` argument is a PTM_collapse_key value (e.g. "AAAS_S174_M1") and
# is ignored for the proteome / transcriptome layers, where each gene is one
# feature. A gene can carry hundreds of sites, so one must be chosen here.

if (identical(omics_layer, "phosphoproteome")) {
  if (is.null(phosphosite)) {
    stop("`phosphosite` must be supplied when omics_layer == 'phosphoproteome'.",
         call. = FALSE)
  }
  if (!phosphosite %in% selected_data$PTM_collapse_key) {
    stop("Phosphosite '", phosphosite, "' not found for ", feature,
         " in the phosphoproteome data.", call. = FALSE)
  }
  selected_data <- selected_data |> dplyr::filter(PTM_collapse_key == .env$phosphosite)
  limma_row     <- limma_row     |> dplyr::filter(PTM_collapse_key == .env$phosphosite)
}

# Plot title: the gene for proteome / transcriptome, the site for phospho.
plot_title <- if (identical(omics_layer, "phosphoproteome")) {
  paste(phosphosite, "phosphorylation")
} else {
  paste(unique(selected_data$Gene_name), "expression in the", (omics_layer))
}

# Y-axis label: log2 counts for the transcriptome, log2 LFQ intensity for the
# proteome and phosphoproteome (protein / peptide intensities).
y_label <- if (identical(omics_layer, "transcriptome")) {
  "Log2 counts"
} else {
  "Log2 LFQ intensity"
}

# Pull one p-value to 3 significant figures, or NA if the feature/column is
# absent. Using [1] guards against features that occupy several rows (e.g.
# multiple phospho-sites).

pull_p <- function(column) {
  value <- limma_row[[column]]
  if (length(value) == 0) NA_real_ else signif(value[1], 3)
}

p_val_post_NGT <- pull_p("P.Value_post_NGT")
p_val_rec_NGT  <- pull_p("P.Value_rec_NGT")
p_val_post_T2D <- pull_p("P.Value_post_T2D")
p_val_rec_T2D  <- pull_p("P.Value_rec_T2D")

# Create the brackets, keeping ONLY significant comparisons (p <= 0.05). Each
# row carries its facet (Group) and the numeric x positions of the two
# timepoints it spans (base = 1, post = 2, rec = 3 on the discrete x axis).
# Non-significant or missing (NA) comparisons are removed here, so they never
# reach the plot.

y_max <- max(selected_data$Value, na.rm = TRUE)

brackets <- dplyr::bind_rows(
  data.frame(
    xmin = 1, xmax = 2,
    y.position = y_max + 0.5,
    label = c(p_val_post_NGT, p_val_post_T2D),
    Group = c("NGT", "T2D")
  ),
  data.frame(
    xmin = 1, xmax = 3,
    y.position = y_max + 1.2,
    label = c(p_val_rec_NGT, p_val_rec_T2D),
    Group = c("NGT", "T2D")
  )
) |> dplyr::filter(!is.na(label), label <= 0.05)

# Draw each significant bracket as base-ggplot geoms instead of
# ggpubr::geom_bracket: a top bar plus two short downward tips (geom_segment)
# and the p-value (geom_text). These look the same as a ggpubr bracket but,
# unlike GeomBracket, survive plotly::ggplotly() and keep working per facet.
bracket_tip  <- 0.12  # length of the bracket's downward tips, in y (log2) units
label_offset <- 0.15  # gap from the bracket bar up to its p-value label, in y
                      # (log2) units -- increase this to raise the label
bracket_segments <- do.call(rbind, lapply(seq_len(nrow(brackets)), function(i) {
  b <- brackets[i, ]
  data.frame(
    x     = c(b$xmin, b$xmin, b$xmax),
    xend  = c(b$xmax, b$xmin, b$xmax),
    y     = c(b$y.position, b$y.position, b$y.position),
    yend  = c(b$y.position, b$y.position - bracket_tip, b$y.position - bracket_tip),
    Group = b$Group,
    label = b$label
  )
}))

# Box-and-whisker statistics per group x timepoint. We draw the box ourselves
# (geom_rect + geom_segment) instead of geom_boxplot, because ggplotly() renders
# box traces *under* the violin's scatter trace (hiding them); rect/segment geoms
# share the violin's layer and so stay on top.
box_stats <- selected_data |>
  dplyr::filter(!is.na(Value)) |>
  dplyr::group_by(Group, Timepoint) |>
  dplyr::summarise(
    q1   = stats::quantile(Value, 0.25),
    med  = stats::median(Value),
    q3   = stats::quantile(Value, 0.75),
    ymin = min(Value[Value >= stats::quantile(Value, 0.25) - 1.5 * stats::IQR(Value)]),
    ymax = max(Value[Value <= stats::quantile(Value, 0.75) + 1.5 * stats::IQR(Value)]),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    xpos = match(Timepoint, c("base", "post", "rec")),
    xmin = xpos - 0.1,
    xmax = xpos + 0.1,
    # Tooltip text (HTML <br> = line break) shown via tooltip = "text" in ggplotly().
    hover = sprintf(
      "Median: %.2f<br>Q3 (75%%): %.2f<br>Q1 (25%%): %.2f<br>Upper whisker: %.2f<br>Lower whisker: %.2f",
      med, q3, q1, ymax, ymin
    )
  )

# New facet label names for supp variable
  new.labs <- c("non-T2D", "T2D")
  names(new.labs) <- c("NGT", "T2D")

# Attach the same median / quartile summary to the violin data, so hovering the
# violin or the box shows the identical tooltip.
selected_data <- selected_data |>
  dplyr::left_join(dplyr::select(box_stats, Group, Timepoint, hover),
                   by = c("Group", "Timepoint"))

output <- ggplot2::ggplot(selected_data,
                  ggplot2::aes(x = Timepoint, 
                               y = Value)) +
  ggplot2::geom_violin(
    ggplot2::aes(fill = Timepoint, text = hover)) +

  # Boxplot as scatter-layer geoms (whiskers, box, median) so it sits on top of
  # the violin in plotly. Added after the violin -> drawn over it. `text` drives
  # the hover tooltip (median + quartiles) when ggplotly(tooltip = "text").
  ggplot2::geom_segment(
    data = box_stats,
    ggplot2::aes(x = xpos, xend = xpos, y = ymin, yend = ymax, text = hover),
    inherit.aes = FALSE, linewidth = 0.3
  ) +
  ggplot2::geom_rect(
    data = box_stats,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = q1, ymax = q3, text = hover),
    inherit.aes = FALSE, fill = "white", colour = "black", linewidth = 0.3
  ) +
  ggplot2::geom_segment(
    data = box_stats,
    ggplot2::aes(x = xmin, xend = xmax, y = med, yend = med, text = hover),
    inherit.aes = FALSE, linewidth = 0.4
  ) +
    
  # Brackets (bar + tips, then the label) - only when something is significant.
  (if (nrow(brackets) > 0) ggplot2::geom_segment(
    data = bracket_segments,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend, text = paste0("p = ", label)),
    inherit.aes = FALSE, linewidth = 0.3
  ) else NULL) +

  (if (nrow(brackets) > 0) ggplot2::geom_text(
    data = brackets,
    ggplot2::aes(x = (xmin + xmax) / 2, y = y.position + label_offset, label = label,
                 text = paste0("p = ", label)),
    inherit.aes = FALSE, vjust = 0, size = base_size / ggplot2::.pt
  ) else NULL) +
    
  ggplot2::scale_x_discrete(limits = c("base", "post", "rec"),
                   labels = c("Base", "Post", "Rec")) +
  
  ggplot2::scale_fill_manual(
    limits = c("base", "post", "rec"),
    values = c("#6ABFA4", "#F18B64", "#8DA0CB")
  ) +
  
  ggplot2::ylab(y_label) +
  ggplot2::expand_limits(y = c(min(selected_data$Value, na.rm = TRUE) - 1,
                               max(selected_data$Value, na.rm = TRUE) +
                                 if (nrow(brackets) > 0) 1.8 else 0.5)) +
  ggplot2::ggtitle(plot_title) +
  ggplot2::facet_wrap( ~ Group,
              strip.position = "bottom",
              labeller = ggplot2::labeller(Group = new.labs)) +
  ggplot2::theme_light() +
  ggplot2::theme(
    text = ggplot2::element_text(size = base_size),
    legend.position = "none",
    axis.title.x = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(size = base_size + 2, hjust = 0.5),
    strip.text.x = ggplot2::element_text(
      margin = ggplot2::margin(0.01, 0.01, 0.1, 0.01, "cm"),
      colour = "black"
    ),
    strip.background = ggplot2::element_rect(fill = "white"),
    strip.placement = "outside",
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(color = "#B4B4B4"),
    axis.line.y = ggplot2::element_line(color = "#B4B4B4")
  )

return(output)

}

# Example (uncomment to test interactively):
# violin_T2D(data_proteins, limma_proteins_T2D, "IQGAP1", "phosphoproteome", phosphosite = "IQGAP1_S1441_M1")

#' Violin of the overall (main) effect of exercise
#'
#' Like [violin_T2D()], but pools all participants across `Group` (and `Sex`), so
#' each timepoint is a single violin / boxplot with no faceting. Brackets use the
#' overall exercise p-values (`P.Value_post`, `P.Value_rec`).
#'
#' @param data Long abundance table (e.g. `data_proteins`) with columns
#'   `Gene_name`, `omic_layer`, `PTM_collapse_key`, `Timepoint`, `Value`.
#' @param limma_output limma results table with `Gene_name`, `omics_layer`,
#'   `PTM_collapse_key` and the overall `P.Value_post` / `P.Value_rec` columns.
#' @param feature gene, or protein (matched against `Gene_name`)
#' @param omics_layer one of "transcriptome", "proteome" or "phosphoproteome"
#' @param phosphosite Phosphosite to plot, given as a `PTM_collapse_key` value
#'   (e.g. "AAAS_S174_M1"). Only used when `omics_layer == "phosphoproteome"`
#'   (required there); ignored for the other layers.
#' @param base_size Base font size (pt). The default suits the interactive app;
#'   pass a smaller value (e.g. 6) for small static figures.
#'
#' @returns A ggplot object.
#' @export
#'
#' @examples
#' violin_main_effect(data_proteins, limma_proteins_T2D, "IQGAP1", "proteome")
violin_main_effect <- function(data, limma_output, feature, omics_layer, phosphosite = NULL,
                       base_size = 11) {
  
  # Select the feature and omics layer:
  
  selected_data <- data |> 
    dplyr::filter(
      Gene_name == feature
    ) |> 
    dplyr::filter(
      omic_layer == omics_layer
    )
  
  # Extract the p-values (single limma row for this feature + omics layer).
  
  limma_row <- limma_output |>
    dplyr::filter(Gene_name == feature, omics_layer == .env$omics_layer)
  
  # Phosphoproteome only: narrow the gene down to a single phosphosite. The
  # `phosphosite` argument is a PTM_collapse_key value (e.g. "AAAS_S174_M1") and
  # is ignored for the proteome / transcriptome layers, where each gene is one
  # feature. A gene can carry hundreds of sites, so one must be chosen here.
  
  if (identical(omics_layer, "phosphoproteome")) {
    if (is.null(phosphosite)) {
      stop("`phosphosite` must be supplied when omics_layer == 'phosphoproteome'.",
           call. = FALSE)
    }
    if (!phosphosite %in% selected_data$PTM_collapse_key) {
      stop("Phosphosite '", phosphosite, "' not found for ", feature,
           " in the phosphoproteome data.", call. = FALSE)
    }
    selected_data <- selected_data |> dplyr::filter(PTM_collapse_key == .env$phosphosite)
    limma_row     <- limma_row     |> dplyr::filter(PTM_collapse_key == .env$phosphosite)
  }
  
  # Plot title: the gene for proteome / transcriptome, the site for phospho.
  plot_title <- if (identical(omics_layer, "phosphoproteome")) {
    paste(phosphosite, "phosphorylation")
  } else {
    paste(unique(selected_data$Gene_name), "expression in the", (omics_layer))
  }
  
  # Y-axis label: log2 counts for the transcriptome, log2 LFQ intensity for the
  # proteome and phosphoproteome (protein / peptide intensities).
  y_label <- if (identical(omics_layer, "transcriptome")) {
    "Log2 counts"
  } else {
    "Log2 LFQ intensity"
  }
  
  # Pull one p-value to 3 significant figures, or NA if the feature/column is
  # absent. Using [1] guards against features that occupy several rows (e.g.
  # multiple phospho-sites).
  
  pull_p <- function(column) {
    value <- limma_row[[column]]
    if (length(value) == 0) NA_real_ else signif(value[1], 3)
  }
  
  p_val_post <- pull_p("P.Value_post")
  p_val_rec  <- pull_p("P.Value_rec")
  
  # Create the brackets, keeping ONLY significant comparisons (p <= 0.05). Each
  # row carries its facet (Group) and the numeric x positions of the two
  # timepoints it spans (base = 1, post = 2, rec = 3 on the discrete x axis).
  # Non-significant or missing (NA) comparisons are removed here, so they never
  # reach the plot.
  
  y_max <- max(selected_data$Value, na.rm = TRUE)
  
  brackets <- dplyr::bind_rows(
    data.frame(
      xmin = 1, xmax = 2,
      y.position = y_max + 0.5,
      label = c(p_val_post)
    ),
    data.frame(
      xmin = 1, xmax = 3,
      y.position = y_max + 1.2,
      label = c(p_val_rec)
    )
  ) |> dplyr::filter(!is.na(label), label <= 0.05)
  
  # Draw each significant bracket as base-ggplot geoms instead of
  # ggpubr::geom_bracket: a top bar plus two short downward tips (geom_segment)
  # and the p-value (geom_text). These look the same as a ggpubr bracket but,
  # unlike GeomBracket, survive plotly::ggplotly() and keep working per facet.
  bracket_tip  <- 0.12  # length of the bracket's downward tips, in y (log2) units
  label_offset <- 0.15  # gap from the bracket bar up to its p-value label, in y
  # (log2) units -- increase this to raise the label
  bracket_segments <- do.call(rbind, lapply(seq_len(nrow(brackets)), function(i) {
    b <- brackets[i, ]
    data.frame(
      x     = c(b$xmin, b$xmin, b$xmax),
      xend  = c(b$xmax, b$xmin, b$xmax),
      y     = c(b$y.position, b$y.position, b$y.position),
      yend  = c(b$y.position, b$y.position - bracket_tip, b$y.position - bracket_tip),
      label = b$label
    )
  }))
  
  # Box-and-whisker statistics per group x timepoint. We draw the box ourselves
  # (geom_rect + geom_segment) instead of geom_boxplot, because ggplotly() renders
  # box traces *under* the violin's scatter trace (hiding them); rect/segment geoms
  # share the violin's layer and so stay on top.
  box_stats <- selected_data |>
    dplyr::filter(!is.na(Value)) |>
    dplyr::group_by(Timepoint) |>   # main effect: pool across Group (and Sex)
    dplyr::summarise(
      q1   = stats::quantile(Value, 0.25),
      med  = stats::median(Value),
      q3   = stats::quantile(Value, 0.75),
      ymin = min(Value[Value >= stats::quantile(Value, 0.25) - 1.5 * stats::IQR(Value)]),
      ymax = max(Value[Value <= stats::quantile(Value, 0.75) + 1.5 * stats::IQR(Value)]),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      xpos = match(Timepoint, c("base", "post", "rec")),
      xmin = xpos - 0.1,
      xmax = xpos + 0.1,
      # Tooltip text (HTML <br> = line break) shown via tooltip = "text" in ggplotly().
      hover = sprintf(
        "Median: %.2f<br>Q3 (75%%): %.2f<br>Q1 (25%%): %.2f<br>Upper whisker: %.2f<br>Lower whisker: %.2f",
        med, q3, q1, ymax, ymin
      )
    )
  
  
  # Attach the same median / quartile summary to the violin data, so hovering the
  # violin or the box shows the identical tooltip.
  selected_data <- selected_data |>
    dplyr::left_join(dplyr::select(box_stats, Timepoint, hover),
                     by = "Timepoint")
  
  output <- ggplot2::ggplot(selected_data,
                            ggplot2::aes(x = Timepoint, 
                                         y = Value)) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = Timepoint, text = hover)) +
    
    # Boxplot as scatter-layer geoms (whiskers, box, median) so it sits on top of
    # the violin in plotly. Added after the violin -> drawn over it. `text` drives
    # the hover tooltip (median + quartiles) when ggplotly(tooltip = "text").
    ggplot2::geom_segment(
      data = box_stats,
      ggplot2::aes(x = xpos, xend = xpos, y = ymin, yend = ymax, text = hover),
      inherit.aes = FALSE, linewidth = 0.3
    ) +
    ggplot2::geom_rect(
      data = box_stats,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = q1, ymax = q3, text = hover),
      inherit.aes = FALSE, fill = "white", colour = "black", linewidth = 0.3
    ) +
    ggplot2::geom_segment(
      data = box_stats,
      ggplot2::aes(x = xmin, xend = xmax, y = med, yend = med, text = hover),
      inherit.aes = FALSE, linewidth = 0.4
    ) +
    
    # Brackets (bar + tips, then the label) - only when something is significant.
    (if (nrow(brackets) > 0) ggplot2::geom_segment(
      data = bracket_segments,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, text = paste0("p = ", label)),
      inherit.aes = FALSE, linewidth = 0.3
    ) else NULL) +
    
    (if (nrow(brackets) > 0) ggplot2::geom_text(
      data = brackets,
      ggplot2::aes(x = (xmin + xmax) / 2, y = y.position + label_offset, label = label,
                   text = paste0("p = ", label)),
      inherit.aes = FALSE, vjust = 0, size = base_size / ggplot2::.pt
    ) else NULL) +
    
    ggplot2::scale_x_discrete(limits = c("base", "post", "rec"),
                              labels = c("Base", "Post", "Rec")) +
    
    ggplot2::scale_fill_manual(
      limits = c("base", "post", "rec"),
      values = c("#6ABFA4", "#F18B64", "#8DA0CB")
    ) +
    
    ggplot2::ylab(y_label) +
    ggplot2::expand_limits(y = c(min(selected_data$Value, na.rm = TRUE) - 1,
                                 max(selected_data$Value, na.rm = TRUE) +
                                   if (nrow(brackets) > 0) 1.8 else 0.5)) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme_light() +
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size),
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size + 2, hjust = 0.5),
      strip.text.x = ggplot2::element_text(
        margin = ggplot2::margin(0.01, 0.01, 0.1, 0.01, "cm"),
        colour = "black"
      ),
      strip.background = ggplot2::element_rect(fill = "white"),
      strip.placement = "outside",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "#B4B4B4"),
      axis.line.y = ggplot2::element_line(color = "#B4B4B4")
    )
  
  return(output)
  
}

#' Violins for the sex-specific exercise response
#'
#' Like [violin_T2D()], but faceted by `Sex` (Women / Men) instead of disease
#' group, with brackets drawn from the sex-specific timepoint p-values
#' (`P.Value_post_W` / `P.Value_post_M`, `P.Value_rec_W` / `P.Value_rec_M`).
#'
#' @param data Long abundance table (e.g. `data_proteins`) with columns
#'   `Gene_name`, `omic_layer`, `PTM_collapse_key`, `Timepoint`, `Sex`, `Value`.
#' @param limma_output Sex-specific limma table (e.g. `limma_proteins_sex`) with
#'   `Gene_name`, `omics_layer`, `PTM_collapse_key` and the `P.Value_*_W` /
#'   `P.Value_*_M` columns.
#' @param feature gene, or protein (matched against `Gene_name`)
#' @param omics_layer one of "transcriptome", "proteome" or "phosphoproteome"
#' @param phosphosite Phosphosite to plot, given as a `PTM_collapse_key` value
#'   (e.g. "AAAS_S174_M1"). Only used when `omics_layer == "phosphoproteome"`
#'   (required there); ignored for the other layers.
#' @param base_size Base font size (pt). The default suits the interactive app;
#'   pass a smaller value (e.g. 6) for small static figures.
#'
#' @returns A ggplot object.
#' @export
#'
#' @examples
#' violin_sex(data_proteins, limma_proteins_sex, "IQGAP1", "proteome")
violin_sex <- function(data, limma_output, feature, omics_layer, phosphosite = NULL,
                       base_size = 11) {
  
  # Select the feature and omics layer:
  
  selected_data <- data |> 
    dplyr::filter(
      Gene_name == feature
    ) |> 
    dplyr::filter(
      omic_layer == omics_layer
    )
  
  # Extract the p-values (single limma row for this feature + omics layer).
  
  limma_row <- limma_output |>
    dplyr::filter(Gene_name == feature, omics_layer == .env$omics_layer)
  
  # Phosphoproteome only: narrow the gene down to a single phosphosite. The
  # `phosphosite` argument is a PTM_collapse_key value (e.g. "AAAS_S174_M1") and
  # is ignored for the proteome / transcriptome layers, where each gene is one
  # feature. A gene can carry hundreds of sites, so one must be chosen here.
  
  if (identical(omics_layer, "phosphoproteome")) {
    if (is.null(phosphosite)) {
      stop("`phosphosite` must be supplied when omics_layer == 'phosphoproteome'.",
           call. = FALSE)
    }
    if (!phosphosite %in% selected_data$PTM_collapse_key) {
      stop("Phosphosite '", phosphosite, "' not found for ", feature,
           " in the phosphoproteome data.", call. = FALSE)
    }
    selected_data <- selected_data |> dplyr::filter(PTM_collapse_key == .env$phosphosite)
    limma_row     <- limma_row     |> dplyr::filter(PTM_collapse_key == .env$phosphosite)
  }
  
  # Plot title: the gene for proteome / transcriptome, the site for phospho.
  plot_title <- if (identical(omics_layer, "phosphoproteome")) {
    paste(phosphosite, "phosphorylation")
  } else {
    paste(unique(selected_data$Gene_name), "expression in the", (omics_layer))
  }
  
  # Y-axis label: log2 counts for the transcriptome, log2 LFQ intensity for the
  # proteome and phosphoproteome (protein / peptide intensities).
  y_label <- if (identical(omics_layer, "transcriptome")) {
    "Log2 counts"
  } else {
    "Log2 LFQ intensity"
  }
  
  # Pull one p-value to 3 significant figures, or NA if the feature/column is
  # absent. Using [1] guards against features that occupy several rows (e.g.
  # multiple phospho-sites).
  
  pull_p <- function(column) {
    value <- limma_row[[column]]
    if (length(value) == 0) NA_real_ else signif(value[1], 3)
  }
  
  p_val_post_W <- pull_p("P.Value_post_W")
  p_val_rec_W  <- pull_p("P.Value_rec_W")
  p_val_post_M <- pull_p("P.Value_post_M")
  p_val_rec_M  <- pull_p("P.Value_rec_M")
  
  # Create the brackets, keeping ONLY significant comparisons (p <= 0.05). Each
  # row carries its facet (Sex) and the numeric x positions of the two
  # timepoints it spans (base = 1, post = 2, rec = 3 on the discrete x axis).
  # Non-significant or missing (NA) comparisons are removed here, so they never
  # reach the plot.
  
  y_max <- max(selected_data$Value, na.rm = TRUE)
  
  brackets <- dplyr::bind_rows(
    data.frame(
      xmin = 1, xmax = 2,
      y.position = y_max + 0.5,
      label = c(p_val_post_W, p_val_post_M),
      Sex = c("W", "M")
    ),
    data.frame(
      xmin = 1, xmax = 3,
      y.position = y_max + 1.2,
      label = c(p_val_rec_W, p_val_rec_M),
      Sex = c("W", "M")
    )
  ) |> dplyr::filter(!is.na(label), label <= 0.05)
  
  # Draw each significant bracket as base-ggplot geoms instead of
  # ggpubr::geom_bracket: a top bar plus two short downward tips (geom_segment)
  # and the p-value (geom_text). These look the same as a ggpubr bracket but,
  # unlike GeomBracket, survive plotly::ggplotly() and keep working per facet.
  bracket_tip  <- 0.12  # length of the bracket's downward tips, in y (log2) units
  label_offset <- 0.15  # gap from the bracket bar up to its p-value label, in y
  # (log2) units -- increase this to raise the label
  bracket_segments <- do.call(rbind, lapply(seq_len(nrow(brackets)), function(i) {
    b <- brackets[i, ]
    data.frame(
      x     = c(b$xmin, b$xmin, b$xmax),
      xend  = c(b$xmax, b$xmin, b$xmax),
      y     = c(b$y.position, b$y.position, b$y.position),
      yend  = c(b$y.position, b$y.position - bracket_tip, b$y.position - bracket_tip),
      Sex = b$Sex,
      label = b$label
    )
  }))
  
  # Box-and-whisker statistics per Sex x timepoint. We draw the box ourselves
  # (geom_rect + geom_segment) instead of geom_boxplot, because ggplotly() renders
  # box traces *under* the violin's scatter trace (hiding them); rect/segment geoms
  # share the violin's layer and so stay on top.
  box_stats <- selected_data |>
    dplyr::filter(!is.na(Value)) |>
    dplyr::group_by(Sex, Timepoint) |>
    dplyr::summarise(
      q1   = stats::quantile(Value, 0.25),
      med  = stats::median(Value),
      q3   = stats::quantile(Value, 0.75),
      ymin = min(Value[Value >= stats::quantile(Value, 0.25) - 1.5 * stats::IQR(Value)]),
      ymax = max(Value[Value <= stats::quantile(Value, 0.75) + 1.5 * stats::IQR(Value)]),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      xpos = match(Timepoint, c("base", "post", "rec")),
      xmin = xpos - 0.1,
      xmax = xpos + 0.1,
      # Tooltip text (HTML <br> = line break) shown via tooltip = "text" in ggplotly().
      hover = sprintf(
        "Median: %.2f<br>Q3 (75%%): %.2f<br>Q1 (25%%): %.2f<br>Upper whisker: %.2f<br>Lower whisker: %.2f",
        med, q3, q1, ymax, ymin
      )
    )
  
  # New facet label names for supp variable
  new.labs <- c("Women", "Men")
  names(new.labs) <- c("W", "M")
  
  # Attach the same median / quartile summary to the violin data, so hovering the
  # violin or the box shows the identical tooltip.
  selected_data <- selected_data |>
    dplyr::left_join(dplyr::select(box_stats, Sex, Timepoint, hover),
                     by = c("Sex", "Timepoint"))
  
  output <- ggplot2::ggplot(selected_data,
                            ggplot2::aes(x = Timepoint, 
                                         y = Value)) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = Timepoint, text = hover)) +
    
    # Boxplot as scatter-layer geoms (whiskers, box, median) so it sits on top of
    # the violin in plotly. Added after the violin -> drawn over it. `text` drives
    # the hover tooltip (median + quartiles) when ggplotly(tooltip = "text").
    ggplot2::geom_segment(
      data = box_stats,
      ggplot2::aes(x = xpos, xend = xpos, y = ymin, yend = ymax, text = hover),
      inherit.aes = FALSE, linewidth = 0.3
    ) +
    ggplot2::geom_rect(
      data = box_stats,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = q1, ymax = q3, text = hover),
      inherit.aes = FALSE, fill = "white", colour = "black", linewidth = 0.3
    ) +
    ggplot2::geom_segment(
      data = box_stats,
      ggplot2::aes(x = xmin, xend = xmax, y = med, yend = med, text = hover),
      inherit.aes = FALSE, linewidth = 0.4
    ) +
    
    # Brackets (bar + tips, then the label) - only when something is significant.
    (if (nrow(brackets) > 0) ggplot2::geom_segment(
      data = bracket_segments,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, text = paste0("p = ", label)),
      inherit.aes = FALSE, linewidth = 0.3
    ) else NULL) +
    
    (if (nrow(brackets) > 0) ggplot2::geom_text(
      data = brackets,
      ggplot2::aes(x = (xmin + xmax) / 2, y = y.position + label_offset, label = label,
                   text = paste0("p = ", label)),
      inherit.aes = FALSE, vjust = 0, size = base_size / ggplot2::.pt
    ) else NULL) +
    
    ggplot2::scale_x_discrete(limits = c("base", "post", "rec"),
                              labels = c("Base", "Post", "Rec")) +
    
    ggplot2::scale_fill_manual(
      limits = c("base", "post", "rec"),
      values = c("#6ABFA4", "#F18B64", "#8DA0CB")
    ) +
    
    ggplot2::ylab(y_label) +
    ggplot2::expand_limits(y = c(min(selected_data$Value, na.rm = TRUE) - 1,
                                 max(selected_data$Value, na.rm = TRUE) +
                                   if (nrow(brackets) > 0) 1.8 else 0.5)) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::facet_wrap( ~ Sex,
                         strip.position = "bottom",
                         labeller = ggplot2::labeller(Sex = new.labs)) +
    ggplot2::theme_light() +
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size),
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size + 2, hjust = 0.5),
      strip.text.x = ggplot2::element_text(
        margin = ggplot2::margin(0.01, 0.01, 0.1, 0.01, "cm"),
        colour = "black"
      ),
      strip.background = ggplot2::element_rect(fill = "white"),
      strip.placement = "outside",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "#B4B4B4"),
      axis.line.y = ggplot2::element_line(color = "#B4B4B4")
    )
  
  return(output)
  
}

# Example (uncomment to test interactively):
# violin_sex(data_proteins, limma_proteins_sex, "IQGAP1", "phosphoproteome", phosphosite = "IQGAP1_S1441_M1")


#' Limma results table for the selected feature across omics layers
#'
#' Builds a table of limma statistics (logFC, AveExpr, P.Value, adj.P.Val) for
#' the active comparison, with one row per measurement of the gene across the
#' transcriptome, proteome and phosphoproteome (one row per phosphosite). In the
#' app a row can be clicked to switch the plot to that layer / site.
#'
#' @param comparison One of "Main effect", "Type 2 Diabetes" or "Sex" (the
#'   active plot tab). Chooses the limma table and the column scope: the overall
#'   columns, the NGT/T2D pair, or the Women/Men pair respectively.
#' @param feature gene / protein name (matched against `Gene_name`)
#' @param limma_T2D,limma_sex the T2D and sex limma results tables
#'
#' @returns A data.frame with columns `Layer`, `Feature` (gene name, or the
#'   phosphosite id for the phosphoproteome) then the comparison's statistic
#'   columns; or NULL if the feature is absent. Numeric columns are left
#'   unrounded so the caller can format them (e.g. DT::formatSignif).
#' @export
limma_results_table <- function(comparison, feature, limma_T2D, limma_sex) {
  spec <- switch(
    comparison,
    "Main effect"     = list(tbl = limma_T2D, scopes = c(Overall = "")),
    "Type 2 Diabetes" = list(tbl = limma_T2D, scopes = c(NGT = "_NGT", T2D = "_T2D")),
    "Sex"             = list(tbl = limma_sex, scopes = c(Women = "_W", Men = "_M")),
    NULL
  )
  if (is.null(spec)) return(NULL)

  # All measurements of the gene, ordered transcriptome -> proteome -> phospho.
  rows <- spec$tbl |>
    dplyr::filter(Gene_name == feature) |>
    dplyr::mutate(.ord = match(omics_layer, c("transcriptome", "proteome", "phosphoproteome"))) |>
    dplyr::arrange(.ord, PTM_collapse_key)
  if (nrow(rows) == 0) return(NULL)

  out <- data.frame(
    Layer   = rows$omics_layer,
    # gene name for transcriptome/proteome, phosphosite id for the phosphoproteome
    Feature = ifelse(rows$omics_layer == "phosphoproteome", rows$PTM_collapse_key, rows$Gene_name),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  stats <- c("logFC", "AveExpr", "P.Value", "adj.P.Val")
  times <- c(post = "Post", rec = "Rec")
  # Wide layout: one column per (statistic x timepoint x scope), grouped by scope.
  for (sc in seq_along(spec$scopes)) {
    scope_suf <- spec$scopes[[sc]]
    scope_lab <- names(spec$scopes)[sc]
    for (tm in seq_along(times)) {
      t_suf <- names(times)[tm]
      t_lab <- times[[tm]]
      for (st in stats) {
        disp <- if (nzchar(scope_suf)) paste0(st, " (", t_lab, ", ", scope_lab, ")")
                else paste0(st, " (", t_lab, ")")
        out[[disp]] <- rows[[paste0(st, "_", t_suf, scope_suf)]]
      }
    }
  }
  out
}


#' Violin plot for a metabolite (one of three statistical comparisons)
#'
#' The metabolite analogue of [violin_T2D()] / [violin_main_effect()] /
#' [violin_sex()], but for a single feature type (metabolites have no omics
#' layer or phosphosites). `comparison` selects the facet variable and the limma
#' p-value columns:
#'   * "Main effect"     -> all participants pooled, no facet, overall p-values.
#'   * "Type 2 Diabetes" -> faceted by Group (non-T2D / T2D), `*_NGT` / `*_T2D` p-values.
#'   * "Sex"             -> faceted by Sex (Women / Men), `*_W` / `*_M` p-values.
#'
#' Brackets (base->post, base->rec) are shown only for comparisons with p <= 0.05,
#' drawn with base geoms so they survive plotly::ggplotly() (see [violin_T2D()]).
#'
#' @param data Long metabolite table (e.g. `data_metabolites`) with columns
#'   `CHEMICAL_NAME`, `Timepoint`, `Group`, `Sex`, `Value`.
#' @param limma_output limma table for the chosen comparison: `limma_metabolites_TD`
#'   for "Main effect" / "Type 2 Diabetes", `limma_metabolites_sex` for "Sex".
#' @param feature metabolite name (matched against `CHEMICAL_NAME`)
#' @param comparison one of "Main effect", "Type 2 Diabetes" or "Sex"
#' @param base_size Base font size (pt). The default suits the interactive app.
#'
#' @returns A ggplot object.
#' @export
#'
#' @examples
#' violin_metabolite(data_metabolites, limma_metabolites_TD, "glucose", "Main effect")
#' violin_metabolite(data_metabolites, limma_metabolites_sex, "glucose", "Sex")
violin_metabolite <- function(data, limma_output, feature, comparison, base_size = 11) {

  # Comparison-specific config: facet variable, the p-value column suffix per
  # scope, and the facet levels / display labels.
  cfg <- switch(
    comparison,
    "Main effect"     = list(facet = NA_character_, suffixes = c(Overall = ""),
                             levels = NA, labels = NULL),
    "Type 2 Diabetes" = list(facet = "Group", suffixes = c(NGT = "_NGT", T2D = "_T2D"),
                             levels = c("NGT", "T2D"), labels = c(NGT = "non-T2D", T2D = "T2D")),
    "Sex"             = list(facet = "Sex", suffixes = c(W = "_W", M = "_M"),
                             levels = c("W", "M"), labels = c(W = "Women", M = "Men")),
    stop("`comparison` must be 'Main effect', 'Type 2 Diabetes' or 'Sex'.", call. = FALSE)
  )
  faceted <- !is.na(cfg$facet)

  selected_data <- data |> dplyr::filter(CHEMICAL_NAME == feature)
  limma_row     <- limma_output |> dplyr::filter(CHEMICAL_NAME == feature)

  # Pull one p-value to 3 significant figures, or NA if absent.
  pull_p <- function(column) {
    value <- limma_row[[column]]
    if (length(value) == 0) NA_real_ else signif(value[1], 3)
  }

  y_max <- max(selected_data$Value, na.rm = TRUE)

  # Significance brackets per scope (overall, or one per facet level): base->post
  # (x 1->2) and base->rec (x 1->3). Keep only p <= 0.05.
  brackets <- do.call(dplyr::bind_rows, lapply(seq_along(cfg$suffixes), function(s) {
    suf <- cfg$suffixes[[s]]
    b <- dplyr::bind_rows(
      data.frame(xmin = 1, xmax = 2, y.position = y_max + 0.5,
                 label = pull_p(paste0("P.Value_post", suf))),
      data.frame(xmin = 1, xmax = 3, y.position = y_max + 1.2,
                 label = pull_p(paste0("P.Value_rec", suf)))
    )
    if (faceted) b[[cfg$facet]] <- cfg$levels[s]
    b
  }))
  brackets <- brackets |> dplyr::filter(!is.na(label), label <= 0.05)

  # Bracket geometry (bar + two downward tips) + label, base geoms (ggplotly-safe).
  bracket_tip  <- 0.12
  label_offset <- 0.15
  bracket_segments <- do.call(rbind, lapply(seq_len(nrow(brackets)), function(i) {
    b <- brackets[i, ]
    seg <- data.frame(
      x     = c(b$xmin, b$xmin, b$xmax),
      xend  = c(b$xmax, b$xmin, b$xmax),
      y     = c(b$y.position, b$y.position, b$y.position),
      yend  = c(b$y.position, b$y.position - bracket_tip, b$y.position - bracket_tip),
      label = b$label
    )
    if (faceted) seg[[cfg$facet]] <- b[[cfg$facet]]
    seg
  }))

  # Box-and-whisker stats, grouped by the facet variable (if any) + Timepoint, so
  # the box sits on top of the violin in plotly (see violin_T2D() for the why).
  group_vars <- if (faceted) c(cfg$facet, "Timepoint") else "Timepoint"
  box_stats <- selected_data |>
    dplyr::filter(!is.na(Value)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      q1   = stats::quantile(Value, 0.25),
      med  = stats::median(Value),
      q3   = stats::quantile(Value, 0.75),
      ymin = min(Value[Value >= stats::quantile(Value, 0.25) - 1.5 * stats::IQR(Value)]),
      ymax = max(Value[Value <= stats::quantile(Value, 0.75) + 1.5 * stats::IQR(Value)]),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      xpos = match(Timepoint, c("base", "post", "rec")),
      xmin = xpos - 0.1,
      xmax = xpos + 0.1,
      hover = sprintf(
        "Median: %.2f<br>Q3 (75%%): %.2f<br>Q1 (25%%): %.2f<br>Upper whisker: %.2f<br>Lower whisker: %.2f",
        med, q3, q1, ymax, ymin
      )
    )

  # Same median/quartile tooltip on the violin (joined on the grouping vars).
  selected_data <- selected_data |>
    dplyr::left_join(dplyr::select(box_stats, dplyr::all_of(c(group_vars, "hover"))),
                     by = group_vars)

  output <- ggplot2::ggplot(selected_data, ggplot2::aes(x = Timepoint, y = Value)) +
    ggplot2::geom_violin(ggplot2::aes(fill = Timepoint, text = hover)) +
    ggplot2::geom_segment(
      data = box_stats,
      ggplot2::aes(x = xpos, xend = xpos, y = ymin, yend = ymax, text = hover),
      inherit.aes = FALSE, linewidth = 0.3
    ) +
    ggplot2::geom_rect(
      data = box_stats,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = q1, ymax = q3, text = hover),
      inherit.aes = FALSE, fill = "white", colour = "black", linewidth = 0.3
    ) +
    ggplot2::geom_segment(
      data = box_stats,
      ggplot2::aes(x = xmin, xend = xmax, y = med, yend = med, text = hover),
      inherit.aes = FALSE, linewidth = 0.4
    ) +
    (if (nrow(brackets) > 0) ggplot2::geom_segment(
      data = bracket_segments,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, text = paste0("p = ", label)),
      inherit.aes = FALSE, linewidth = 0.3
    ) else NULL) +
    (if (nrow(brackets) > 0) ggplot2::geom_text(
      data = brackets,
      ggplot2::aes(x = (xmin + xmax) / 2, y = y.position + label_offset, label = label,
                   text = paste0("p = ", label)),
      inherit.aes = FALSE, vjust = 0, size = base_size / ggplot2::.pt
    ) else NULL) +
    ggplot2::scale_x_discrete(limits = c("base", "post", "rec"),
                              labels = c("Base", "Post", "Rec")) +
    ggplot2::scale_fill_manual(limits = c("base", "post", "rec"),
                               values = c("#6ABFA4", "#F18B64", "#8DA0CB")) +
    ggplot2::ylab("Log2 relative abundance") +   # TODO: confirm the metabolite unit
    ggplot2::expand_limits(y = c(min(selected_data$Value, na.rm = TRUE) - 1,
                                 max(selected_data$Value, na.rm = TRUE) +
                                   if (nrow(brackets) > 0) 1.8 else 0.5)) +
    ggplot2::ggtitle(paste(feature, "levels")) +
    ggplot2::theme_light() +
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size),
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size + 2, hjust = 0.5),
      strip.text.x = ggplot2::element_text(
        margin = ggplot2::margin(0.01, 0.01, 0.1, 0.01, "cm"),
        colour = "black"
      ),
      strip.background = ggplot2::element_rect(fill = "white"),
      strip.placement = "outside",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "#B4B4B4"),
      axis.line.y = ggplot2::element_line(color = "#B4B4B4")
    )

  # Facet the T2D / Sex comparisons (Main effect is a single pooled panel).
  if (faceted) {
    output <- output + ggplot2::facet_wrap(
      stats::as.formula(paste("~", cfg$facet)),
      strip.position = "bottom",
      labeller = do.call(ggplot2::labeller, stats::setNames(list(cfg$labels), cfg$facet))
    )
  }
  output
}

# Example (uncomment to test interactively):
# violin_metabolite(data_metabolites, limma_metabolites_TD, "glucose", "Type 2 Diabetes")


#' Limma results table for a metabolite
#'
#' The metabolite analogue of [limma_results_table()]: a one-row table of limma
#' statistics (logFC, AveExpr, P.Value, adj.P.Val) for the selected metabolite in
#' the active comparison. (Metabolites are a single feature type, so there is one
#' row and no layer / site columns.)
#'
#' @param comparison One of "Main effect", "Type 2 Diabetes" or "Sex". Chooses
#'   the limma table and the column scope (overall, NGT/T2D, or Women/Men).
#' @param feature metabolite name (matched against `CHEMICAL_NAME`)
#' @param limma_T2D,limma_sex the metabolite T2D and sex limma tables
#'
#' @returns A one-row data.frame: `Metabolite` then the comparison's statistic
#'   columns; or NULL if the metabolite is absent. Numeric columns are unrounded.
#' @export
limma_metabolites_table <- function(comparison, feature, limma_T2D, limma_sex) {
  spec <- switch(
    comparison,
    "Main effect"     = list(tbl = limma_T2D, scopes = c(Overall = "")),
    "Type 2 Diabetes" = list(tbl = limma_T2D, scopes = c(NGT = "_NGT", T2D = "_T2D")),
    "Sex"             = list(tbl = limma_sex, scopes = c(Women = "_W", Men = "_M")),
    NULL
  )
  if (is.null(spec)) return(NULL)

  row <- spec$tbl |> dplyr::filter(CHEMICAL_NAME == feature)
  if (nrow(row) == 0) return(NULL)

  out   <- data.frame(Metabolite = row$CHEMICAL_NAME[1], check.names = FALSE, stringsAsFactors = FALSE)
  stats <- c("logFC", "AveExpr", "P.Value", "adj.P.Val")
  times <- c(post = "Post", rec = "Rec")
  for (sc in seq_along(spec$scopes)) {
    scope_suf <- spec$scopes[[sc]]
    scope_lab <- names(spec$scopes)[sc]
    for (tm in seq_along(times)) {
      t_suf <- names(times)[tm]
      t_lab <- times[[tm]]
      for (st in stats) {
        disp <- if (nzchar(scope_suf)) paste0(st, " (", t_lab, ", ", scope_lab, ")")
                else paste0(st, " (", t_lab, ")")
        out[[disp]] <- row[[paste0(st, "_", t_suf, scope_suf)]][1]
      }
    }
  }
  out
}
