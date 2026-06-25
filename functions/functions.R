

data_proteins <- vroom::vroom(here::here("data/data_files/data_proteins.txt"))
limma_proteins_T2D <- vroom::vroom(here::here("data/limma_outputs/limma_proteins_T2D.txt"))


# Select the feature and omics layer:

selected_data <- data_proteins |> 
  dplyr::filter(
    Gene_name == "IQGAP1"
  ) |> 
  dplyr::filter(
    omic_layer == "proteome"
  )

# Extract the p-values:

p_val_post_NGT <- limma_proteins_T2D |> dplyr::filter(Gene_name == "IQGAP1") |> dplyr::filter(omics_layer == "proteome") |> dplyr::pull(P.Value_post_NGT)
p_val_rec_NGT <- limma_proteins_T2D |> dplyr::filter(Gene_name == "IQGAP1") |> dplyr::filter(omics_layer == "proteome") |> dplyr::pull(P.Value_rec_NGT)
p_val_post_T2D <- limma_proteins_T2D |> dplyr::filter(Gene_name == "IQGAP1") |> dplyr::filter(omics_layer == "proteome") |> dplyr::pull(P.Value_post_T2D)
p_val_rec_T2D <- limma_proteins_T2D |> dplyr::filter(Gene_name == "IQGAP1") |> dplyr::filter(omics_layer == "proteome") |> dplyr::pull(P.Value_rec_T2D)

# Create the brackets:

bracket_post <- data.frame(
  xmin = c(1, 1),
  xmax = c(2, 2),
  y.position = c(max(selected_data$Value + 0.5), max(selected_data$Value + 0.5)),
  label = c(p_val_post_NGT, p_val_rec_NGT),
  Group = c("NGT", "T2D")
)

bracket_rec <- data.frame(
  xmin = c(1, 1),
  xmax = c(3, 3),
  y.position = c(max(selected_data$Value + 0.8), max(selected_data$Value + 0.8)),
  label = c(p_val_rec_NGT, p_val_rec_T2D),
  Group = c("NGT", "T2D")
)

# NR4A3_plot <- 
  
ggplot2::ggplot(selected_data, 
                  ggplot2::aes(x = Timepoint, 
                               y = Value)) +
  ggplot2::geom_violin(
    ggplot2::aes(fill = Timepoint)) +
    
  ggplot2::geom_boxplot(width = 0.2, 
                        outliers = F) +
    
  ggpubr::geom_bracket(
    data = bracket_post,
    xmin = bracket_post$xmin,
    xmax = bracket_post$xmax,
    y.position = bracket_post$y.position,
    label = bracket_post$label,
    label.size = 2
  )  +
    
  ggpubr::geom_bracket(
    data = bracket_rec,
    xmin = bracket_rec$xmin,
    xmax = bracket_rec$xmax,
    y.position = bracket_rec$y.position,
    label = bracket_rec$label,
    label.size = 2
  ) +
    
  ggplot2::scale_x_discrete(limits = c("base", "post", "rec"),
                   labels = c("Base", "Post", "Rec")) +
  
  ggplot2::scale_fill_manual(
    limits = c("base", "post", "rec"),
    values = c("#6ABFA4", "#F18B64", "#8DA0CB")
  ) +
  
  ggplot2::ylab(expression(Log[2] * " count")) +
  ggplot2::expand_limits(y = c(min(selected_data$Value) - 1, max(selected_data$Value) + 1)) +
  ggplot2::ggtitle(paste(unique(selected_data$Gene_name), "expression")) +
  ggplot2::facet_wrap( ~ Group,
              # labeller = ggplot2::labeller(Group = Disease_labels),
              strip.position = "bottom") +
  ggplot2::theme_light() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 6),
    legend.position = "none",
    axis.title.x = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(size = 8, hjust = 0.5),
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
