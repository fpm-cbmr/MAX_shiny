

# Data is normally loaded by app.R. These guarded loads let this file also be
# sourced on its own (e.g. interactive testing) without reloading the objects
# when they already exist.
if (!exists("data_proteins")) {
  data_proteins <- vroom::vroom(here::here("data/data_files/data_proteins.txt"))
}
if (!exists("limma_proteins_T2D")) {
  limma_proteins_T2D <- vroom::vroom(here::here("data/limma_outputs/limma_proteins_T2D.txt"))
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
    y.position = y_max + 0.8,
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
                                 if (nrow(brackets) > 0) 1.3 else 0.5)) +
  ggplot2::ggtitle(plot_title) +
  ggplot2::facet_wrap( ~ Group,
              strip.position = "bottom") +
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

