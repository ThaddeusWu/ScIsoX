#################################################################
#  Visualisation and Analysis Functions for ScIsoX              #
#                                                               #
#  Author: [Siyuan Wu & Ulf Schmitz]                            #
#  Institution: [James Cook University]                         #
#  Date: Jul 29, 2025                                           #
#  Package: ScIsoX V1.1.0                                       #
#################################################################

########################
# Required Libraries   #
########################
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line geom_vline geom_boxplot
#' @importFrom ggplot2 geom_density geom_point theme_minimal labs facet_wrap scale_fill_gradient2
#' @importFrom ggplot2 theme element_text ggtitle xlab ylab
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid textGrob gpar unit
#' @importFrom scales rescale
#' @importFrom stats setNames quantile density kmeans cor as.dist complete.cases reshape
#' @importFrom utils head
#' @importFrom grDevices colorRampPalette dev.list dev.off pdf
#' @importFrom graphics barplot axis segments
#' @importFrom magrittr %>%
#' @importFrom MASS kde2d
#' @importFrom dplyr case_when
NULL

# Define global variables to avoid R CMD check notes
utils::globalVariables(c(
  "Metric", "CellType", "diversity_pattern", "idx", "isoform", "proportion", 
  "quadrant", "short_name", "count", "percentage", "z_norm", 
  "stable_count", "unclassified_count", "Value", "cell_type", "gene",
  "gene_colours"
))

#' Plot multiple threshold visualisations in a grid layout
#' 
#' This function arranges multiple threshold fitting visualisations in a grid layout,
#' which is useful for comparing distributions and threshold choices across
#' different complexity metrics.
#' 
#' @param threshold_plots List of ggplot2 objects showing threshold visualisations
#' @param ncol Number of columns in the grid layout
#' @param title Overall title for the grid
#' 
#' @return A grid arrangement of threshold visualisation plots that can be printed or saved
#' 
#' @examples
#' # Load example data
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' # Create SCHT object with recommended QC parameters
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Calculate complexity metrics
#' tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)
#' 
#' # Plot all threshold visualisations
#' plot_threshold_visualisations(tc_results$threshold_plots, ncol = 2)
#' 
#' # Plot with custom title
#' plot_threshold_visualisations(
#'   tc_results$threshold_plots, 
#'   ncol = 3, 
#'   title = "Complexity Metric Thresholds"
#' )
#' 
#' @export
plot_threshold_visualisations <- function(threshold_plots, ncol = 3, title = "Threshold Fitting Visualisations") {
  # Check that plot list is not empty
  if(length(threshold_plots) == 0) {
    stop("No plots provided to arrange")
  }
  
  # Check if gridExtra is available
  if(!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("This function requires the 'gridExtra' package. Please install it with: install.packages('gridExtra')")
  }
  
  # Optional title
  if(!is.null(title)) {
    title_grob <- grid::textGrob(title, gp = grid::gpar(fontsize = 14, fontface = "bold"))
    
    # Arrange plots with title
    arranged_plot <- gridExtra::grid.arrange(
      gridExtra::arrangeGrob(
        grobs = threshold_plots,
        ncol = ncol
      ),
      top = title_grob
    )
  } else {
    # Arrange plots without title
    arranged_plot <- gridExtra::grid.arrange(
      grobs = threshold_plots,
      ncol = ncol
    )
  }
  
  return(arranged_plot)
}

#' Get colour palette for transcriptome complexity visualisation
#' 
#' Internal function to create a consistent colour palette for visualisations.
#' 
#' @return A list with colour specifications
#' 
#' @keywords internal
.get_tc_palette <- function() {
  # Check if viridis is available
  if (requireNamespace("viridis", quietly = TRUE)) {
    gradient <- viridis::viridis(10)
  } else {
    # Fallback gradient colors similar to viridis
    gradient <- colorRampPalette(c("#440154", "#31688e", "#35b779", "#fde725"))(10)
  }
  
  list(
    main = c("#cc6677", "#40b0a6", "#e1be6a", "#867bb9"),
    highlight = "#f7ac53",  
    gradient = gradient,
    background = "white",
    grid = "#EEEEEE",
    text = "#333333"
  )
}

#' Theme for transcriptomic complexity visualisation
#' 
#' Internal function to create a consistent theme for plots.
#' 
#' @return A ggplot2 theme object
#' 
#' @keywords internal
.tc_theme <- function() {
  palette <- .get_tc_palette()
  
  ggplot2::theme_minimal() +
    ggplot2::theme(
      # Text elements
      text = ggplot2::element_text(family = "sans", color = palette$text),
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0, margin = ggplot2::margin(b = 15)),
      axis.title = ggplot2::element_text(size = 11, face = "bold"),
      axis.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text = ggplot2::element_text(size = 9),
      
      # Panels and backgrounds
      panel.background = ggplot2::element_rect(fill = palette$background, color = NA),
      plot.background = ggplot2::element_rect(fill = palette$background, color = NA),
      panel.grid.major = ggplot2::element_line(color = palette$grid, linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      
      # Legend appearance
      legend.background = ggplot2::element_rect(fill = palette$background, color = NA),
      legend.key = ggplot2::element_blank(),
      
      # Spacing and margins
      plot.margin = ggplot2::margin(15, 15, 15, 15),
      
      # Strip appearance for facets
      strip.background = ggplot2::element_rect(fill = palette$background, color = NA),
      strip.text = ggplot2::element_text(size = 10, face = "bold")
    )
}

#' Get color palette with fallback options
#' 
#' Internal function to get color palettes with fallbacks when RColorBrewer is not available
#' 
#' @param n Number of colors needed
#' @param palette_name Name of the RColorBrewer palette
#' @return Vector of colors
#' 
#' @keywords internal
.get_colour_palette <- function(n, palette_name = "Set3") {
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    if (n <= 12 && palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
      max_colors <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
      return(RColorBrewer::brewer.pal(min(n, max_colors), palette_name)[1:n])
    } else if (n > 12) {
      # For more colors, use colorRampPalette
      base_colours <- RColorBrewer::brewer.pal(min(12, RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]), palette_name)
      return(colorRampPalette(base_colours)(n))
    }
  }
  
  # Fallback palettes when RColorBrewer is not available
  fallback_palettes <- list(
    Set3 = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", 
             "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"),
    Set1 = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", 
             "#a65628", "#f781bf", "#999999"),
    Set2 = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", 
             "#e5c494", "#b3b3b3"),
    Paired = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", 
               "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"),
    Pastel1 = c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6", "#ffffcc", 
                "#e5d8bd", "#fddaec"),
    Pastel2 = c("#b3e2cd", "#fdcdac", "#cbd5e8", "#f4cae4", "#e6f5c9", "#fff2ae", 
                "#f1e2cc", "#cccccc")
  )
  
  base_colours <- if (palette_name %in% names(fallback_palettes)) {
    fallback_palettes[[palette_name]]
  } else {
    fallback_palettes[["Set3"]]  # Default fallback
  }
  
  if (n <= length(base_colours)) {
    return(base_colours[1:n])
  } else {
    return(colorRampPalette(base_colours)(n))
  }
}

#' Get metric-specific classification names
#' 
#' Internal function to return classification names for different metrics.
#' 
#' @param metric Name of the metric
#' @return A vector with two elements: high and low classification names
#' 
#' @keywords internal
.get_metric_names <- function(metric) {
  if(metric == "intra_cellular_isoform_diversity") {
    return(c("Strong Isoform Co-expression", "Weak Isoform Co-expression"))
  } else if(metric == "inter_cellular_isoform_diversity") {
    return(c("High Isoform Diversity", "Low Isoform Diversity"))
  } else if(metric == "intra_cell_type_heterogeneity") {
    return(c("High Cellular Heterogeneity", "Low Cellular Heterogeneity"))
  } else if(metric == "inter_cell_type_specificity") {
    return(c("Cell-Type-Specific Isoform Expression", "Cell-Type-Independent Isoform Expression"))
  } else if(metric == "intra_cell_type_heterogeneity_variability") {
    return(c("Cell-Type-Dependent Heterogeneity", "Cell-Type-Independent Heterogeneity"))
  } else if(metric == "inter_cell_type_difference_variability") {
    return(c("High Cell-Type Distinctions", "Low Cell-Type Distinctions"))
  } else if(metric == "cell_type_coexpression_variability") {
    return(c("Cell-Type-Adaptive Co-Expression", "Cell-Type-Consistent Co-Expression"))
  } else {
    return(c(paste0("High ", metric), paste0("Low ", metric)))
  }
}

#' Prepare data for transcriptomic complexity visualisation
#' 
#' Internal function to prepare data for plotting.
#' 
#' @param tc_results Transcriptomic complexity results object or metrics data frame
#' @param x_metric Name of metric for x-axis
#' @param y_metric Name of metric for y-axis
#' @param highlight_genes Optional vector of gene names to highlight
#' @param label_annotation Column name to use for highlighting/labelling genes
#' @param label_top Number of top genes to label if highlight_genes not provided
#' @param label_direction Direction for selecting genes: "top" (highest values) or "bottom" (lowest values)
#' @param use_thresholds Whether to use thresholds from tc_results
#' @param x_threshold Manual threshold value for x-axis
#' @param y_threshold Manual threshold value for y-axis
#' 
#' @return A list containing prepared data and derived information
#' 
#' @keywords internal
.prepare_tc_data <- function(tc_results, 
                             x_metric = "inter_cellular_isoform_diversity",
                             y_metric = "inter_cell_type_specificity",
                             highlight_genes = NULL,
                             label_annotation = "intra_cell_type_heterogeneity",
                             label_top = 10,
                             label_direction = "top",
                             use_thresholds = TRUE,
                             x_threshold = 0.6,
                             y_threshold = 0.6) {
  
  # Extract metrics data frame
  if(is.data.frame(tc_results)) {
    metrics <- tc_results
  } else if("metrics" %in% names(tc_results)) {
    metrics <- tc_results$metrics
  } else {
    stop("Input must be a transcriptomic complexity results object or a data frame with metrics")
  }
  
  # Get thresholds if available and requested
  if(use_thresholds && "thresholds" %in% names(tc_results)) {
    if(x_metric %in% names(tc_results$thresholds)) {
      x_threshold <- tc_results$thresholds[[x_metric]]
    }
    if(y_metric %in% names(tc_results$thresholds)) {
      y_threshold <- tc_results$thresholds[[y_metric]]
    }
  }
  
  # Get classification names for selected metrics
  x_class_names <- .get_metric_names(x_metric)
  y_class_names <- .get_metric_names(y_metric)
  
  
  all_metrics <- c(
    "intra_cellular_isoform_diversity",
    "inter_cellular_isoform_diversity", 
    "intra_cell_type_heterogeneity",
    "inter_cell_type_specificity",
    "intra_cell_type_heterogeneity_variability",
    "inter_cell_type_difference_variability",
    "cell_type_coexpression_variability"
  )
  
  # Create readable metric names for labels
  readable_metrics <- c(
    "Intra-cellular Isoform Diversity",
    "Inter-cellular Isoform Diversity",
    "Intra-cell-type Heterogeneity",
    "Inter-cell-type Specificity",
    "Intra-cell-type Heterogeneity Variability",
    "Inter-cell-type Difference Variability",
    "Cell-type-specific Co-Expression Variability"
  )
  
  x_label_idx <- match(x_metric, all_metrics)
  x_label <- readable_metrics[x_label_idx]
  
  y_label_idx <- match(y_metric, all_metrics)
  y_label <- readable_metrics[y_label_idx]
  
  # Filter data to remove NA values for the selected metrics
  plot_data <- metrics[!is.na(metrics[[x_metric]]) & !is.na(metrics[[y_metric]]), ]
  
  # Handle case of empty dataset
  if(nrow(plot_data) == 0) {
    stop("No genes have valid values for both selected metrics")
  }
  
  # Define quadrant names
  q1 <- paste0(x_class_names[1], " + ", y_class_names[1])  # High X + High Y (top right)
  q2 <- paste0(x_class_names[2], " + ", y_class_names[1])  # Low X + High Y (top left)
  q3 <- paste0(x_class_names[2], " + ", y_class_names[2])  # Low X + Low Y (bottom left)
  q4 <- paste0(x_class_names[1], " + ", y_class_names[2])  # High X + Low Y (bottom right)
  
  # Assign quadrants
  plot_data$quadrant <- dplyr::case_when(
    plot_data[[x_metric]] > x_threshold & plot_data[[y_metric]] > y_threshold ~ q1,
    plot_data[[x_metric]] <= x_threshold & plot_data[[y_metric]] > y_threshold ~ q2,
    plot_data[[x_metric]] <= x_threshold & plot_data[[y_metric]] <= y_threshold ~ q3,
    TRUE ~ q4
  )
  
  # Create shorthand versions for compact labels
  plot_data$quadrant_short <- dplyr::case_when(
    plot_data[[x_metric]] > x_threshold & plot_data[[y_metric]] > y_threshold ~ "Q1",
    plot_data[[x_metric]] <= x_threshold & plot_data[[y_metric]] > y_threshold ~ "Q2",
    plot_data[[x_metric]] <= x_threshold & plot_data[[y_metric]] <= y_threshold ~ "Q3",
    TRUE ~ "Q4"
  )
  
  # Determine which genes to highlight
  if(is.null(highlight_genes)) {
    # Default: highlight top or bottom genes by label_annotation value
    if(label_annotation %in% colnames(plot_data)) {
      # Determine ordering based on direction
      decreasing_order <- (label_direction == "top")
      ordered_genes <- plot_data$gene[
        order(plot_data[[label_annotation]], decreasing = decreasing_order)
      ]
      top_genes <- ordered_genes[1:min(label_top, nrow(plot_data))]
    } else {
      warning(paste0("Label annotation column '", label_annotation, "' not found. Using gene name for ordering."))
      top_genes <- head(plot_data$gene, label_top)
    }
  } else {
    # Use provided genes that exist in filtered data
    top_genes <- highlight_genes[highlight_genes %in% plot_data$gene]
    if(length(top_genes) == 0) {
      warning("None of the specified highlight genes have valid metric values")
      # Fall back to default behaviour
      top_genes <- head(plot_data$gene, label_top)
    }
  }
  
  # Create dataset of highlighted genes
  highlight_data <- plot_data[plot_data$gene %in% top_genes, ]
  
  # Count genes in each quadrant
  quadrant_counts <- table(plot_data$quadrant)
  
  # Calculate percentages
  quadrant_pct <- round(quadrant_counts / sum(quadrant_counts) * 100, 2)
  
  # Create a data frame for quadrant statistics
  quadrant_df <- data.frame(
    quadrant = names(quadrant_counts),
    count = as.numeric(quadrant_counts),
    percentage = as.numeric(quadrant_pct),
    stringsAsFactors = FALSE
  )
  
  # Add short names to quadrant_df
  quadrant_df$short_name <- dplyr::case_when(
    quadrant_df$quadrant == q1 ~ "Q1",
    quadrant_df$quadrant == q2 ~ "Q2",
    quadrant_df$quadrant == q3 ~ "Q3",
    quadrant_df$quadrant == q4 ~ "Q4"
  )
  
  # Ensure quadrant_df has the correct ordering
  quadrant_df$short_name <- factor(quadrant_df$short_name, levels = c("Q1", "Q2", "Q3", "Q4"))
  quadrant_df <- quadrant_df[order(quadrant_df$short_name), ]
  
  # Create quadrant colours with consistent naming
  palette <- .get_tc_palette()
  quadrant_colours <- palette$main
  names(quadrant_colours) <- c(q1, q2, q3, q4)
  
  # Transform variables with "variability" in name for visualisation only
  need_scaling <- FALSE
  if(grepl("variability", x_metric) || grepl("variability", y_metric)) {
    # Create a copy of the data for visualisation
    plot_data_vis <- plot_data
    
    # Scale the variables that exceed 0-1 range for visualisation purposes only
    if(grepl("variability", x_metric) && max(plot_data[[x_metric]], na.rm = TRUE) > 1) {
      plot_data_vis[[x_metric]] <- (plot_data[[x_metric]] - min(plot_data[[x_metric]], na.rm = TRUE)) / 
        (max(plot_data[[x_metric]], na.rm = TRUE) - min(plot_data[[x_metric]], na.rm = TRUE))
      # Scale threshold for visualisation
      x_threshold_vis <- (x_threshold - min(plot_data[[x_metric]], na.rm = TRUE)) / 
        (max(plot_data[[x_metric]], na.rm = TRUE) - min(plot_data[[x_metric]], na.rm = TRUE))
      need_scaling <- TRUE
    } else {
      x_threshold_vis <- x_threshold
    }
    
    if(grepl("variability", y_metric) && max(plot_data[[y_metric]], na.rm = TRUE) > 1) {
      plot_data_vis[[y_metric]] <- (plot_data[[y_metric]] - min(plot_data[[y_metric]], na.rm = TRUE)) / 
        (max(plot_data[[y_metric]], na.rm = TRUE) - min(plot_data[[y_metric]], na.rm = TRUE))
      # Scale threshold for visualisation
      y_threshold_vis <- (y_threshold - min(plot_data[[y_metric]], na.rm = TRUE)) / 
        (max(plot_data[[y_metric]], na.rm = TRUE) - min(plot_data[[y_metric]], na.rm = TRUE))
      need_scaling <- TRUE
    } else {
      y_threshold_vis <- y_threshold
    }
    
    # Use the transformed data for visualisation only
    vis_data <- plot_data_vis
  } else {
    # If no transformation needed, use original data
    vis_data <- plot_data
    x_threshold_vis <- x_threshold
    y_threshold_vis <- y_threshold
  }
  
  # If we have highlighted genes and transformed data, adjust highlight data too
  if(nrow(highlight_data) > 0 && exists("plot_data_vis")) {
    highlight_data_vis <- highlight_data
    if(grepl("variability", x_metric) && max(plot_data[[x_metric]], na.rm = TRUE) > 1) {
      highlight_data_vis[[x_metric]] <- (highlight_data[[x_metric]] - min(plot_data[[x_metric]], na.rm = TRUE)) / 
        (max(plot_data[[x_metric]], na.rm = TRUE) - min(plot_data[[x_metric]], na.rm = TRUE))
    }
    if(grepl("variability", y_metric) && max(plot_data[[y_metric]], na.rm = TRUE) > 1) {
      highlight_data_vis[[y_metric]] <- (highlight_data[[y_metric]] - min(plot_data[[y_metric]], na.rm = TRUE)) / 
        (max(plot_data[[y_metric]], na.rm = TRUE) - min(plot_data[[y_metric]], na.rm = TRUE))
    }
  } else {
    highlight_data_vis <- highlight_data
  }
  
  # Create label positions for quadrants
  label_positions <- data.frame(
    quadrant = c(q1, q2, q3, q4),
    x = c(0.9, 0.2, 0.2, 0.9),  # Q1, Q2, Q3, Q4 positions
    y = c(0.85, 0.85, 0.15, 0.15)   # Q1, Q2, Q3, Q4 positions
  )
  
  # Merge with statistics
  label_data <- merge(label_positions, quadrant_df, by = "quadrant")
  
  # Return all prepared data structures
  return(list(
    plot_data = plot_data,
    vis_data = vis_data,
    highlight_data = highlight_data,
    highlight_data_vis = highlight_data_vis,
    quadrant_df = quadrant_df,
    label_data = label_data,
    x_threshold = x_threshold,
    y_threshold = y_threshold,
    x_threshold_vis = x_threshold_vis,
    y_threshold_vis = y_threshold_vis,
    x_metric = x_metric,
    y_metric = y_metric,
    x_label = x_label,
    y_label = y_label,
    quadrant_colours = quadrant_colours,
    q1 = q1,
    q2 = q2,
    q3 = q3,
    q4 = q4,
    need_scaling = need_scaling
  ))
}

#' Transcriptomic complexity landscape plot with marginal distributions
#' 
#' Creates a visualisation of transcriptomic complexity landscape with quadrants
#' and marginal distributions. This plot reveals the distribution of genes across
#' different complexity dimensions and helps identify genes with interesting patterns.
#' 
#' @param tc_results Transcriptomic complexity results object or metrics data frame
#' @param x_metric Name of metric for x-axis (default: "inter_cellular_isoform_diversity")
#' @param y_metric Name of metric for y-axis (default: "inter_cell_type_specificity")
#' @param highlight_genes Optional vector of gene names to highlight
#' @param label_annotation Column name to use for highlighting/labelling genes
#' @param n_label Number of genes to label if highlight_genes not provided
#' @param label_direction Direction for selecting genes: "top" (highest values) or "bottom" (lowest values)
#' @param use_thresholds Whether to use thresholds from tc_results
#' @param x_threshold Manual threshold value for x-axis
#' @param y_threshold Manual threshold value for y-axis
#' @param point_transparency Alpha value for points (0-1)
#' 
#' @return A ggplot object with marginal distributions that can be printed or saved
#' 
#' @examples
#' # Load example data and create SCHT
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Calculate complexity metrics
#' tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)
#' 
#' # Create landscape plot with default metrics
#' plot_tc_landscape(tc_results)
#' 
#' # Use different metrics and highlight specific genes
#' # Find some multi-isoform genes to highlight
#' multi_iso_genes <- names(which(tc_results$metrics$n_isoforms > 3))[1:3]
#' plot_tc_landscape(
#'   tc_results, 
#'   x_metric = "intra_cellular_isoform_diversity", 
#'   y_metric = "inter_cellular_isoform_diversity",
#'   highlight_genes = multi_iso_genes
#' )
#'                  
#' # Highlight bottom 10 genes by heterogeneity
#' plot_tc_landscape(
#'   tc_results,
#'   n_label = 10,
#'   label_direction = "bottom"
#' )
#' 
#' @export
plot_tc_landscape <- function(tc_results, 
                              x_metric = "inter_cellular_isoform_diversity",
                              y_metric = "inter_cell_type_specificity",
                              highlight_genes = NULL,
                              label_annotation = "intra_cell_type_heterogeneity",
                              n_label = 10,
                              label_direction = "top",
                              use_thresholds = TRUE,
                              x_threshold = 0.6,
                              y_threshold = 0.6,
                              point_transparency = 0.85) {
  
  # Required packages check with informative error messages
  required_packages <- c("ggplot2", "dplyr", "ggrepel", "ggExtra", "scales")
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if(length(missing_packages) > 0) {
    stop(paste0("The following packages are required but not installed: ", 
                paste(missing_packages, collapse = ", "), 
                ". Please install them with: install.packages(c('", 
                paste(missing_packages, collapse = "', '"), "'))"))
  }
  
  # Prepare data using the shared preparation function
  data_prep <- .prepare_tc_data(
    tc_results = tc_results,
    x_metric = x_metric,
    y_metric = y_metric,
    highlight_genes = highlight_genes,
    label_annotation = label_annotation,
    label_top = n_label,
    label_direction = label_direction,
    use_thresholds = use_thresholds,
    x_threshold = x_threshold,
    y_threshold = y_threshold
  )
  
  # Extract the prepared data components
  vis_data <- data_prep$vis_data
  highlight_data_vis <- data_prep$highlight_data_vis
  label_data <- data_prep$label_data
  quadrant_colours <- data_prep$quadrant_colours
  x_threshold_vis <- data_prep$x_threshold_vis
  y_threshold_vis <- data_prep$y_threshold_vis
  q1 <- data_prep$q1
  q2 <- data_prep$q2
  q3 <- data_prep$q3
  q4 <- data_prep$q4
  
  # Calculate correlation coefficients
  valid_data <- vis_data[!is.na(vis_data[[x_metric]]) & !is.na(vis_data[[y_metric]]), ]
  
  # Calculate correlations
  pearson_cor <- cor(valid_data[[x_metric]], valid_data[[y_metric]], method = "pearson")
  spearman_cor <- cor(valid_data[[x_metric]], valid_data[[y_metric]], method = "spearman")
  kendall_cor <- cor(valid_data[[x_metric]], valid_data[[y_metric]], method = "kendall")
  
  # Format correlation text
  cor_text <- paste0("Pearson r = ", sprintf("%.3f", pearson_cor), 
                     ", Spearman rho = ", sprintf("%.3f", spearman_cor),
                     ", Kendall tau = ", sprintf("%.3f", kendall_cor))
  
  # Create main plot
  p <- ggplot2::ggplot(vis_data, ggplot2::aes(x = .data[[x_metric]], y = .data[[y_metric]]))
  
  # Add background shading for quadrants
  p <- p + ggplot2::annotate("rect", xmin = x_threshold_vis, xmax = 1, ymin = y_threshold_vis, ymax = 1, 
                             fill = quadrant_colours[q1], alpha = 0.1)  # Q1: High X, High Y (top right)
  p <- p + ggplot2::annotate("rect", xmin = 0, xmax = x_threshold_vis, ymin = y_threshold_vis, ymax = 1, 
                             fill = quadrant_colours[q2], alpha = 0.1)  # Q2: Low X, High Y (top left)
  p <- p + ggplot2::annotate("rect", xmin = 0, xmax = x_threshold_vis, ymin = 0, ymax = y_threshold_vis, 
                             fill = quadrant_colours[q3], alpha = 0.1)  # Q3: Low X, Low Y (bottom left)
  p <- p + ggplot2::annotate("rect", xmin = x_threshold_vis, xmax = 1, ymin = 0, ymax = y_threshold_vis, 
                             fill = quadrant_colours[q4], alpha = 0.1)  # Q4: High X, Low Y (bottom right)
  
  # Add threshold lines
  p <- p + ggplot2::geom_vline(xintercept = x_threshold_vis, linetype = "dashed", 
                               colour = "grey60", alpha = 0.7, linewidth = 1)
  p <- p + ggplot2::geom_hline(yintercept = y_threshold_vis, linetype = "dashed", 
                               colour = "grey60", alpha = 0.7, linewidth = 1)
  
  # Add data points
  p <- p + ggplot2::geom_point(
    ggplot2::aes(colour = quadrant),
    size = 2,  # Fixed size for all points
    alpha = point_transparency
  )
  
  # Add quadrant labels with gene counts and percentages
  p <- p + ggplot2::geom_label(
    data = label_data,
    ggplot2::aes(
      x = x, 
      y = y,
      label = paste0(short_name, ": ", count, " genes (", percentage, "%)"),
      fill = quadrant
    ),
    size = 6,
    fontface = "bold",
    color = "white",
    alpha = 0.9,
    label.size = 0,
    show.legend = FALSE
  )
  
  # Add highlighted points if available
  if(nrow(highlight_data_vis) > 0) {
    palette <- .get_tc_palette()
    
    p <- p + ggplot2::geom_point(
      data = highlight_data_vis,
      ggplot2::aes(x = .data[[x_metric]], y = .data[[y_metric]]),
      size = 4, 
      shape = 21,
      fill = palette$highlight,
      colour = "black",
      stroke = 1
    )
    
    # Add gene labels
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = highlight_data_vis,
        ggplot2::aes(label = gene),
        size = 8,
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.4,
        force = 10,
        segment.colour = "grey50",
        segment.size = 0.3,
        min.segment.length = 0,
        max.overlaps = 30
      )
    } else {
      warning("Package 'ggrepel' not installed. Using geom_text() instead. ",
              "Install it with: install.packages('ggrepel')")
      p <- p + ggplot2::geom_text(
        data = highlight_data_vis,
        ggplot2::aes(label = gene),
        size = 3,
        fontface = "bold",
        check_overlap = TRUE
      )
    }
  }
  
  # Add styling
  p <- p + ggplot2::scale_colour_manual(values = quadrant_colours, name = "Classification")
  p <- p + ggplot2::scale_fill_manual(values = quadrant_colours, guide = "none")
  
  # Set axis limits based on data
  if(data_prep$need_scaling) {
    # If using transformed data, use 0-1 range
    x_limits <- c(0, 1)
    y_limits <- c(0, 1)
    
    p <- p + ggplot2::scale_x_continuous(
      limits = x_limits, 
      breaks = seq(0, 1, by = 0.2)
    )
    p <- p + ggplot2::scale_y_continuous(
      limits = y_limits, 
      breaks = seq(0, 1, by = 0.2)
    )
  } else {
    # Calculate actual data ranges with some padding
    x_range <- range(vis_data[[x_metric]], na.rm = TRUE)
    y_range <- range(vis_data[[y_metric]], na.rm = TRUE)
    
    # Add 5% padding to each side
    x_pad <- diff(x_range) * 0.05
    y_pad <- diff(y_range) * 0.05
    
    # Calculate expanded limits
    x_limits <- c(max(0, x_range[1] - x_pad), x_range[2] + x_pad)
    y_limits <- c(max(0, y_range[1] - y_pad), y_range[2] + y_pad)
    
    # Ensure we cover at least 0-1 range regardless of data
    x_limits[1] <- min(x_limits[1], 0)
    x_limits[2] <- max(x_limits[2], 1)
    y_limits[1] <- min(y_limits[1], 0)
    y_limits[2] <- max(y_limits[2], 1)
    
    p <- p + ggplot2::scale_x_continuous(
      limits = x_limits, 
      breaks = seq(floor(10*x_limits[1])/10, ceiling(10*x_limits[2])/10, by = 0.2)
    )
    p <- p + ggplot2::scale_y_continuous(
      limits = y_limits, 
      breaks = seq(floor(10*y_limits[1])/10, ceiling(10*y_limits[2])/10, by = 0.2)
    )
  }
  
  # Use fixed ratio for proper visual representation with square aspect
  p <- p + ggplot2::coord_fixed(ratio = 1)
  
  p <- p + .tc_theme()
  
  # Add appropriate labels and legend setting
  p <- p + ggplot2::labs(
    title = "Transcriptome Complexity Landscape",
    subtitle = paste0("Comparing ", data_prep$x_label, " and ", data_prep$y_label),
    x = data_prep$x_label,
    y = data_prep$y_label,
    caption = cor_text
  ) + 
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 18, face = "bold"), 
      axis.text = ggplot2::element_text(size = 18),
      plot.caption = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5)
    )
    
  
  # Set legend visibility based on parameter
  
  p <- p + ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 12, face = "bold"),
    legend.text = ggplot2::element_text(size = 10)
  ) + 
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = 2,                     
        byrow = TRUE,                  
        override.aes = list(size = 4), 
        title.position = "left"
      )
    )
  
  
  # Create enhanced marginal distributions if ggExtra is available
  if (requireNamespace("ggExtra", quietly = TRUE)) {
    marginal <- ggExtra::ggMarginal(
      p, 
      type = "densigram",  # Combined histogram and density
      xparams = list(
        fill = scales::alpha("#BBD4E9", 0.8),  
        color = scales::alpha("#396CA0", 0.8),
        size = 0.5
      ),
      yparams = list(
        fill = scales::alpha("#BBD4E9", 0.8),  #
        color = scales::alpha("#396CA0", 0.8),
        size = 0.5
      ),
      margins = "both",
      size = 12      
    )
    return(marginal)
  } else {
    warning("Package 'ggExtra' not installed. Marginal distributions not shown. ",
            "Install it with: install.packages('ggExtra')")
    return(p)
  }
}

#' Transcriptomic complexity density plot
#' 
#' Creates a density visualisation of transcriptomic complexity landscape showing
#' gene distribution patterns. This plot provides a smoother representation of
#' gene concentration across the complexity landscape, highlighting hotspots.
#' 
#' @param tc_results Transcriptomic complexity results object or metrics data frame
#' @param x_metric Name of metric for x-axis (default: "inter_cellular_isoform_diversity")
#' @param y_metric Name of metric for y-axis (default: "inter_cell_type_specificity")
#' @param use_thresholds Whether to use thresholds from tc_results
#' @param x_threshold Manual threshold value for x-axis
#' @param y_threshold Manual threshold value for y-axis
#' 
#' @return A ggplot object with density visualisation that can be printed or saved
#' 
#' @examples
#' # Load example data and create SCHT
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Calculate complexity metrics
#' tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)
#' 
#' # Create density plot with default metrics
#' plot_tc_density(tc_results)
#' 
#' # Use different metrics
#' plot_tc_density(
#'   tc_results, 
#'   x_metric = "intra_cellular_isoform_diversity", 
#'   y_metric = "inter_cellular_isoform_diversity"
#' )
#' 
#' @export
plot_tc_density <- function(tc_results, 
                            x_metric = "inter_cellular_isoform_diversity",
                            y_metric = "inter_cell_type_specificity",
                            use_thresholds = TRUE,
                            x_threshold = 0.6,
                            y_threshold = 0.6) {
  
  # Required packages check with informative error messages
  required_packages <- c("ggplot2", "dplyr", "MASS", "viridis", "scales")
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if(length(missing_packages) > 0) {
    stop(paste0("The following packages are required but not installed: ", 
                paste(missing_packages, collapse = ", "), 
                ". Please install them with: install.packages(c('", 
                paste(missing_packages, collapse = "', '"), "'))"))
  }
  
  # Prepare data using the shared preparation function
  data_prep <- .prepare_tc_data(
    tc_results = tc_results,
    x_metric = x_metric,
    y_metric = y_metric,
    use_thresholds = use_thresholds,
    x_threshold = x_threshold,
    y_threshold = y_threshold
  )
  
  # Extract the prepared data components
  vis_data <- data_prep$vis_data
  quadrant_colours <- data_prep$quadrant_colours
  x_threshold_vis <- data_prep$x_threshold_vis
  y_threshold_vis <- data_prep$y_threshold_vis
  q1 <- data_prep$q1
  q2 <- data_prep$q2
  q3 <- data_prep$q3
  q4 <- data_prep$q4
  
  # Calculate expanded data range properly to avoid warnings
  x_range <- range(vis_data[[x_metric]], na.rm = TRUE)
  y_range <- range(vis_data[[y_metric]], na.rm = TRUE)
  
  # Define display limits first (for consistency)
  x_limits <- c(0, 1)
  y_limits <- c(0, 1)
  
  # Generate 2D density estimate 
  dens <- MASS::kde2d(
    vis_data[[x_metric]], 
    vis_data[[y_metric]], 
    n = 100,  
    lims = c(x_limits[1], x_limits[2], y_limits[1], y_limits[2]) 
  )
  
  # Convert to data frame for ggplot 
  dens_df <- expand.grid(  
    x = round(dens$x, 6),  
    y = round(dens$y, 6)
  )
  dens_df$z <- as.vector(dens$z)
  
  # Normalise density for better visualisation
  dens_df$z_norm <- scales::rescale(dens_df$z, to = c(0, 1))
  
  # Filter out NA or zero values 
  dens_df <- dens_df[!is.na(dens_df$z) & is.finite(dens_df$z) & dens_df$z > 0, ]
  
  # Create density plot with contours - using breaks that match the data range
  p_density <- ggplot2::ggplot(dens_df, ggplot2::aes(x = x, y = y, z = z_norm)) +
    ggplot2::geom_tile(ggplot2::aes(fill = z_norm), alpha = 0.7) +
    ggplot2::geom_contour(
      color = "white", 
      alpha = 0.6,
      size = 1,
      breaks = scales::pretty_breaks(n = 6)(range(dens_df$z_norm, na.rm = TRUE))  # Adaptive breaks
    ) +
    ggplot2::geom_vline(
      xintercept = x_threshold_vis, 
      linetype = "dashed", 
      colour = "white", 
      alpha = 0.3,
      size = 1
    ) +
    ggplot2::geom_hline(
      yintercept = y_threshold_vis, 
      linetype = "dashed", 
      colour = "white", 
      alpha = 0.3,
      size = 1
    ) +
    ggplot2::annotate(
      "label",
      x = c(0.9, 0.15, 0.15, 0.9),
      y = c(0.85, 0.85, 0.15, 0.15),
      label = c("Q1", "Q2", "Q3", "Q4"),
      fill = c(quadrant_colours[q1], quadrant_colours[q2], quadrant_colours[q3], quadrant_colours[q4]),
      color = "white",
      size = 4,
      fontface = "bold",
      alpha = 0.9,
      label.size = 0
    ) +
    ggplot2::scale_x_continuous(limits = x_limits, breaks = seq(0, 1, by = 0.2)) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = seq(0, 1, by = 0.2)) +
    ggplot2::scale_fill_gradientn(
      colors = if (requireNamespace("viridis", quietly = TRUE)) {
        viridis::plasma(10, begin = 0, end = 0.9)
      } else {
        colorRampPalette(c("#0d0887", "#6a00a8", "#b12a90", "#e16462", "#fca636", "#f0f921"))(10)
      },
      name = "Density"
    ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(
      title = "Gene Complexity Distribution Density",
      subtitle = paste0("Comparing ", data_prep$x_label, " and ", data_prep$y_label),
      x = data_prep$x_label,
      y = data_prep$y_label
    ) +
    .tc_theme() +
    ggplot2::theme(
      legend.position = "right",
      axis.title = ggplot2::element_text(size = 18, face = "bold"), 
      axis.text = ggplot2::element_text(size = 16)
    )
  
  return(p_density)
}

#' Dual diversity comparison plot
#' 
#' Visualises the relationship between intra-cellular and inter-cellular
#' diversity to identify different isoform usage patterns. This plot helps identify
#' genes with unusual patterns, particularly those with higher intra-cellular
#' than inter-cellular diversity.
#' 
#' @param tc_results Transcriptomic complexity results object or metrics data frame
#' @param label_top Number of top genes (by idi_difference) to label below the diagonal line
#' @param point_transparency Alpha value for points (0-1)
#' @param use_thresholds Whether to use thresholds from tc_results
#' @param x_threshold Threshold value for x-axis (intra-cellular diversity)
#' @param y_threshold Threshold value for y-axis (inter-cellular diversity)
#' 
#' @return A ggplot object with the diversity comparison visualization that can be printed or saved
#' 
#' @examples
#' # Load example data and create SCHT
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Calculate complexity metrics
#' tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)
#' 
#' # Create diversity comparison plot
#' plot_diversity_comparison(tc_results)
#' 
#' # Label more genes
#' plot_diversity_comparison(tc_results, label_top = 20)
#' 
#' # Adjust transparency and thresholds
#' plot_diversity_comparison(
#'   tc_results, 
#'   label_top = 15,
#'   point_transparency = 0.7,
#'   x_threshold = 0.5,
#'   y_threshold = 0.5
#' )
#' 
#' @export
plot_diversity_comparison <- function(tc_results, 
                                      label_top = 10,
                                      point_transparency = 0.85,
                                      use_thresholds = TRUE,
                                      x_threshold = 0.6,
                                      y_threshold = 0.6) {
  
  # Required packages check with informative error messages
  required_packages <- c("ggplot2", "dplyr", "ggrepel", "scales")
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if(length(missing_packages) > 0) {
    stop(paste0("The following packages are required but not installed: ", 
                paste(missing_packages, collapse = ", "), 
                ". Please install them with: install.packages(c('", 
                paste(missing_packages, collapse = "', '"), "'))"))
  }
  
  # Define metrics to use for the dual diversity plot
  x_metric <- "intra_cellular_isoform_diversity"
  y_metric <- "inter_cellular_isoform_diversity"
  
  # Use data preparation function from the landscape plot
  data_prep <- .prepare_tc_data(
    tc_results = tc_results,
    x_metric = x_metric,
    y_metric = y_metric,
    highlight_genes = NULL,
    label_annotation = "intra_cell_type_heterogeneity",
    label_top = label_top * 2, 
    use_thresholds = use_thresholds,
    x_threshold = x_threshold,
    y_threshold = y_threshold
  )
  
  # Extract the prepared data components
  plot_data <- data_prep$plot_data
  vis_data <- data_prep$vis_data
  highlight_data <- data_prep$highlight_data
  highlight_data_vis <- data_prep$highlight_data_vis
  quadrant_colours <- data_prep$quadrant_colours
  q1 <- data_prep$q1
  q2 <- data_prep$q2
  q3 <- data_prep$q3
  q4 <- data_prep$q4
  
  
  # Get the colour palette
  palette <- .get_tc_palette()
  
  # Define diversity pattern categories (based on thresholds)
  vis_data$diversity_pattern <- dplyr::case_when(
    vis_data[[x_metric]] > data_prep$x_threshold & vis_data[[y_metric]] > data_prep$y_threshold ~ 
      "High in both: Multi-isoform cells across population",
    vis_data[[x_metric]] <= data_prep$x_threshold & vis_data[[y_metric]] > data_prep$y_threshold ~ 
      "High inter-cellular, Low intra-cellular: Cell specialisation",
    vis_data[[x_metric]] < data_prep$x_threshold & vis_data[[y_metric]] <= data_prep$y_threshold ~ 
      "Low in both: Single isoform usage",
    TRUE ~ "Low inter-cellular, High intra-cellular: Intra-cellular co-expression"
  ) 
  
  # Create diversity pattern colours
  diversity_colors <- c(
    "High in both: Multi-isoform cells across population" = palette$main[1],
    "High inter-cellular, Low intra-cellular: Cell specialisation" = palette$main[2],
    "Low in both: Single isoform usage" = palette$main[3],
    "Low inter-cellular, High intra-cellular: Intra-cellular co-expression" = palette$main[4]
  )
  
  # Only label genes below diagonal (where inter-cellular < intra-cellular)
  highlight_data_vis <- highlight_data_vis[
    highlight_data_vis[[y_metric]] < highlight_data_vis[[x_metric]], 
  ]
  
  # If we have too many genes, take only the top requested number
  if(nrow(highlight_data_vis) > label_top) {
    highlight_data_vis <- highlight_data_vis[
      order(abs(highlight_data_vis$idi_difference), decreasing = TRUE),
    ][1:label_top, ]
  }
  
  # Create the main plot
  p <- ggplot2::ggplot(vis_data, ggplot2::aes(x = .data[[x_metric]], y = .data[[y_metric]]))
  
  # Add background shading for quadrants
  p <- p + ggplot2::annotate("rect", xmin = data_prep$x_threshold, xmax = 1, ymin = data_prep$y_threshold, ymax = 1, 
                             fill = quadrant_colours[q1], alpha = 0.1)  # Q1: High X, High Y (top right)
  p <- p + ggplot2::annotate("rect", xmin = 0, xmax = data_prep$x_threshold, ymin = data_prep$y_threshold, ymax = 1, 
                             fill = quadrant_colours[q2], alpha = 0.1)  # Q2: Low X, High Y (top left)
  p <- p + ggplot2::annotate("rect", xmin = 0, xmax = data_prep$x_threshold, ymin = 0, ymax = data_prep$y_threshold, 
                             fill = quadrant_colours[q3], alpha = 0.1)  # Q3: Low X, Low Y (bottom left)
  p <- p + ggplot2::annotate("rect", xmin = data_prep$x_threshold, xmax = 1, ymin = 0, ymax = data_prep$y_threshold, 
                             fill = quadrant_colours[q4], alpha = 0.1)  # Q4: High X, Low Y (bottom right)
  
  
  
  # Reference line (where both values are equal)
  p <- p + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                                color = "black", alpha = 0.8, linewidth = 2) +
    
    # Add threshold lines
    ggplot2::geom_vline(xintercept = data_prep$x_threshold, linetype = "dashed", 
                        color = "grey60", alpha = 0.7, linewidth = 1) +
    ggplot2::geom_hline(yintercept = data_prep$y_threshold, linetype = "dashed", 
                        color = "grey60", alpha = 0.7, linewidth = 1) +
    
    # Add data points
    ggplot2::geom_point(
      ggplot2::aes(color = diversity_pattern),
      size = 2,
      alpha = point_transparency
    )
  
  # Find genes below the diagonal to label
  # 1. Extract the original metrics data frame
  if(is.data.frame(tc_results)) {
    metrics <- tc_results
  } else if("metrics" %in% names(tc_results)) {
    metrics <- tc_results$metrics
  } else {
    stop("Input must be a transcriptomic complexity results object or a data frame with metrics")
  }
  
  # 2. Find genes below the diagonal (where y < x)
  below_diagonal <- metrics[metrics[[y_metric]] < metrics[[x_metric]], ]
  
  # 3. Order by idi_difference and take top genes
  below_diagonal <- below_diagonal[order(abs(below_diagonal$idi_difference), decreasing = TRUE), ]
  genes_to_label <- head(below_diagonal, label_top)
  
  # 4. Add highlighted points for these genes
  if(nrow(genes_to_label) > 0) {
    p <- p + ggplot2::geom_point(
      data = genes_to_label,
      ggplot2::aes(x = .data[[x_metric]], y = .data[[y_metric]]),
      size = 4, 
      shape = 21,
      fill = palette$highlight,
      colour = "black",
      stroke = 1.2
    )
    
    # 5. Add gene labels
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = genes_to_label,
        ggplot2::aes(label = gene),
        size = 6,
        fontface = "bold",
        box.padding = 0.6,
        point.padding = 0.5,
        force = 12,
        segment.colour = "black",
        segment.size = 0.4,
        min.segment.length = 0.1,
        max.overlaps = 20
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = genes_to_label,
        ggplot2::aes(label = gene),
        size = 3,
        fontface = "bold",
        check_overlap = TRUE
      )
    }
  }
  
  # Apply styling
  p <- p + ggplot2::scale_color_manual(values = diversity_colors, name = "Diversity Pattern") +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +  # Equal aspect ratio
    .tc_theme() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 18, face = "bold"), 
      axis.text = ggplot2::element_text(size = 16)
    ) +
    # Labels
    ggplot2::labs(
      title = "Dual Diversity Comparison",
      subtitle = "Intra-cellular vs. Inter-cellular Isoform Diversity",
      x = "Intra-cellular Isoform Diversity",
      y = "Inter-cellular Isoform Diversity"
    ) + 
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = 2,                     
        byrow = TRUE,                  
        override.aes = list(size = 4), 
        title.position = "left"
      )
    )
  
  return(p)
}

#' Create ridge plots for visualising complexity metrics distributions
#' 
#' This function creates ridge plots to visualise the distribution of complexity
#' metrics either globally or by cell type. Ridge plots offer a compact way to
#' compare multiple distributions simultaneously.
#' 
#' @param tc_results Results from calculate_isoform_complexity_metrics()
#' @param type Type of plot: "global" (across all metrics) or "cell_type" (comparing cell types)
#' @param cell_types Vector of cell types to display (for cell_type plots, if NULL all will be detected)
#' @param metrics Metrics to display (if NULL, uses all for global, core 3 for cell_type)
#' @param n_celltypes Maximum number of cell types to show if auto-detecting
#' @param label_y_axis Whether to show metric labels (for cell_type plots)
#' 
#' @return A ggplot object that can be printed or saved
#' 
#' @examples
#' # Load example data and create SCHT
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Calculate complexity metrics
#' tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)
#' 
#' # Create global ridge plot for all metrics
#' plot_complexity_ridges(tc_results, type = "global")
#' 
#' # Create cell type-specific ridge plots
#' plot_complexity_ridges(tc_results, type = "cell_type")
#' 
#' # Plot all metrics (default behavior)
#' plot_complexity_ridges(tc_results, type = "global")
#' 
#' @export
plot_complexity_ridges <- function(tc_results, 
                                   type = "global",
                                   cell_types = NULL,
                                   metrics = NULL,
                                   n_celltypes = 10,
                                   label_y_axis = FALSE) {
  
  # Load required packages
  required_packages <- c("ggplot2", "ggridges", "dplyr", "tidyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is needed for this function to work. Please install it."))
    }
  }
  
  # Define default colour palettes
  global_palette <- c(
    "#4C72B0", "#55A868", "#C44E52", "#8172B3", "#CCB974", "#64B5CD", "#BA82B9"
  )
  
  metric_palette <- c(
    "#4878D0", "#EE854A", "#6ACC64", "#D65F5F", "#956CB4", "#8C613C",
    "#DC7EC0", "#82C6E2", "#D5BB67", "#FF9F7F", "#AEC7E8", "#FFD92F",
    "#A6D854", "#E78AC3", "#8DA0CB", "#FC8D62", "#66C2A5", "#E5C494"
  )
  
  # Set up metric names and map to readable labels
  all_metrics <- c(
    "intra_cellular_isoform_diversity",
    "inter_cellular_isoform_diversity", 
    "intra_cell_type_heterogeneity",
    "inter_cell_type_specificity",
    "intra_cell_type_heterogeneity_variability",
    "inter_cell_type_difference_variability",
    "cell_type_coexpression_variability"
  )
  
  readable_metrics <- c(
    "Intra-cellular Isoform Diversity",
    "Inter-cellular Isoform Diversity",
    "Intra-cell-type Heterogeneity",
    "Inter-cell-type Specificity",
    "Intra-cell-type Heterogeneity Variability",
    "Inter-cell-type Difference Variability",
    "Cell-type-specific Co-Expression Variability"
  )
  
  metric_map <- stats::setNames(readable_metrics, all_metrics)
  
  core_cell_metrics <- c(
    "intra_cellular_isoform_diversity",
    "inter_cellular_isoform_diversity",
    "intra_cell_type_heterogeneity"
  )
  
  # Handle global plot type (showing all 7 metrics)
  if (type == "global") {
    
    # Default metrics if not provided
    if (is.null(metrics)) {
      metrics <- all_metrics
    } else {
      # Verify all metrics are valid
      invalid_metrics <- metrics[!metrics %in% all_metrics]
      if (length(invalid_metrics) > 0) {
        stop("Invalid metrics: ", paste(invalid_metrics, collapse = ", "))
      }
    }
    
    # Extract metrics dataframe
    metrics_df <- tc_results$metrics
    
    # Scale the variables that exceed 0-1 range for visualisation purposes only
    hv <- metrics_df$intra_cell_type_heterogeneity_variability
    dv <- metrics_df$inter_cell_type_difference_variability
    sv <- metrics_df$cell_type_coexpression_variability
    
    hv_max <- max(hv, na.rm = TRUE)
    hv_min <- min(hv, na.rm = TRUE)
    dv_max <- max(dv, na.rm = TRUE)
    dv_min <- min(dv, na.rm = TRUE)
    sv_max <- max(sv, na.rm = TRUE)
    sv_min <- min(sv, na.rm = TRUE)
    
    metrics_df$intra_cell_type_heterogeneity_variability <- (hv - hv_min)/(hv_max - hv_min)
    metrics_df$inter_cell_type_difference_variability <- (dv - dv_min)/(dv_max - dv_min)
    metrics_df$cell_type_coexpression_variability <- (sv - sv_min)/(sv_max - sv_min)
    
    # Select only the specified metrics
    metrics_subset <- metrics_df[, metrics, drop = FALSE]
    
    # Gather data for plotting
    if (requireNamespace("tidyr", quietly = TRUE)) {
      plot_data <- metrics_subset %>%
        tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Metric", values_to = "Value") %>%
        dplyr::filter(!is.na(Value))
    } else {
      # Use base R reshape as fallback
      plot_data <- reshape(metrics_subset,
                          direction = "long",
                          varying = list(names(metrics_subset)),
                          v.names = "Value",
                          timevar = "Metric",
                          times = names(metrics_subset))
      plot_data <- plot_data[!is.na(plot_data$Value), ]
      plot_data$Metric <- as.character(plot_data$Metric)
    }
    
    # Create ordered factor with readable metric names
    plot_data$Metric <- factor(plot_data$Metric, 
                               levels = metrics, 
                               labels = metric_map[metrics])
    
    # Check if ggridges is available
    if (!requireNamespace("ggridges", quietly = TRUE)) {
      stop("Package 'ggridges' is required for ridge plots. ",
           "Please install it with: install.packages('ggridges')", call. = FALSE)
    }
    
    # Create the plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Value, y = Metric, fill = Metric)) +
      ggridges::geom_density_ridges2(alpha = 0.7, scale = 2, rel_min_height = 0.01) +
      ggplot2::scale_fill_manual(values = global_palette) +
      ggplot2::labs(
        title = "Distribution of Complexity Metrics",
        x = NULL,
        y = NULL
      ) +
      ggridges::theme_ridges(font_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 12),
        axis.text.x = ggplot2::element_text(size = 10, angle = 0),
        axis.text.y = ggplot2::element_text(face = "bold", size = 12),
        panel.grid.major.x = ggplot2::element_line(linewidth = 0.3, color = "grey85"),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "cm")
      )
    
  } else if (type == "cell_type") {
    # Handle cell type plot
    # Default metrics if not provided
    if (is.null(metrics)) {
      metrics <- core_cell_metrics
    } else {
      # Verify all metrics are valid
      invalid_metrics <- metrics[!metrics %in% all_metrics]
      if (length(invalid_metrics) > 0) {
        stop("Invalid metrics: ", paste(invalid_metrics, collapse = ", "))
      }
    }
    
    # Extract cell type metrics
    ct_metrics <- tc_results$cell_type_metrics
    
    # Create a data frame for plotting
    plot_data <- data.frame(
      Gene = character(),
      CellType = character(),
      Metric = character(),
      Value = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Extract data for each gene and cell type
    for (gene_name in names(ct_metrics)) {
      gene_data <- ct_metrics[[gene_name]]
      
      for (ct_name in names(gene_data)) {
        ct_data <- gene_data[[ct_name]]
        
        for (metric in metrics) {
          if (!is.null(ct_data[[metric]]) && !is.na(ct_data[[metric]])) {
            plot_data <- rbind(plot_data, data.frame(
              Gene = gene_name,
              CellType = ct_name,
              Metric = metric,
              Value = ct_data[[metric]],
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
    
    # If no valid data is found, return an empty plot with a message
    if (nrow(plot_data) == 0) {
      p <- ggplot2::ggplot() + 
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No valid data found") +
        ggplot2::theme_void()
      return(p)
    }
    
    # Auto-detect cell types if not provided
    if (is.null(cell_types)) {
      all_cell_types <- unique(plot_data$CellType)
      
      # If there are too many cell types, limit to top n by data count
      if (length(all_cell_types) > n_celltypes) {
        cell_counts <- table(plot_data$CellType)
        cell_types <- names(sort(cell_counts, decreasing = TRUE)[1:n_celltypes])
      } else {
        cell_types <- all_cell_types
      }
    }
    
    # Filter data for selected cell types
    plot_data <- plot_data %>%
      dplyr::filter(CellType %in% cell_types)
    
    # Map metric codes to readable names
    plot_data$Metric <- factor(plot_data$Metric, 
                               levels = metrics, 
                               labels = metric_map[metrics])
    
    # Create readable cell type labels
    plot_data$CellType <- factor(plot_data$CellType, levels = cell_types)
    
    # Ensure colour palette is long enough
    if (length(metric_palette) < length(metrics)) {
      metric_palette <- colorRampPalette(metric_palette)(length(metrics))
    }
    
    # Create the plot 
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Value, y = Metric, fill = Metric)) +
      ggridges::geom_density_ridges2(alpha = 0.7, scale = 2, rel_min_height = 0.01) +
      ggplot2::scale_fill_manual(values = metric_palette) +
      ggplot2::labs(
        title = "Complexity Metrics by Cell Type",
        x = NULL,
        y = NULL
      ) +
      ggridges::theme_ridges(font_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
        legend.position = "none",
        axis.title.x = ggplot2::element_text(size = 12),
        axis.text.x = ggplot2::element_text(size = 10, angle = 0),
        axis.text.y = ggplot2::element_text(face = "bold", size = 12),
        panel.grid.major.x = ggplot2::element_line(linewidth = 0.3, color = "grey85"),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold", size = 12),
        plot.margin = ggplot2::unit(c(0.5, 0.5, 0.5, 0.5), "cm")
      ) +
      # Facet by cell type instead of metric
      ggplot2::facet_wrap(~ CellType, scales = "free_y", ncol = min(3, length(cell_types)))
    
    # Remove y-axis labels if requested
    if (!label_y_axis) {
      p <- p + ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        legend.position = "right",  
        legend.title = element_blank(), 
        legend.text = element_text(size = 12)
      )
    }
  } else {
    stop("Invalid type. Must be 'global' or 'cell_type'.")
  }
  
  return(p)
}

#' Extract and compare metrics for multiple genes
#'
#' This helper function extracts metrics from stored results for multiple genes,
#' and returns a structured data frame for comparison. It's useful for preparing 
#' data for custom visualisations or statistical analyses.
#'
#' @param tc_metrics Transcriptomic complexity results object
#' @param gene_names Vector of gene names to extract metrics for
#' @param include_mean Whether to include mean values for each gene (default: TRUE)
#' 
#' @return Data frame with all genes' metrics for easy comparison
#' 
#' @examples
#' # Load example data and create SCHT
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Calculate complexity metrics
#' tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)
#' 
#' # Compare metrics for specific genes
#' # Use genes from the middle of the results to ensure they exist
#' available_genes <- rownames(tc_results$metrics)
#' genes_to_compare <- available_genes[c(100, 200, 300, 400)]
#' comparison_data <- compare_gene_metrics(tc_results, genes_to_compare)
#' 
#' # View the comparison
#' print(comparison_data)
#' 
#' @export
compare_gene_metrics <- function(tc_metrics, gene_names, include_mean = TRUE) {
  # Create list to store each gene's metrics
  all_metrics <- list()
  
  # Process each gene
  for(gene in gene_names) {
    # Get pre-calculated cell type metrics
    gene_ct_metrics <- tc_metrics$cell_type_metrics[[gene]]
    
    if(is.null(gene_ct_metrics)) {
      warning("No cell type metrics found for gene: ", gene)
      next
    }
    
    # Extract core metrics from each cell type
    metrics_list <- list()
    for(cell_type in names(gene_ct_metrics)) {
      ct_data <- gene_ct_metrics[[cell_type]]
      
      # Create data frame for this cell type
      df <- data.frame(
        gene = gene,
        cell_type = cell_type,
        intra_cellular_isoform_diversity = ct_data$intra_cellular_isoform_diversity,
        inter_cellular_isoform_diversity = ct_data$inter_cellular_isoform_diversity,
        intra_cell_type_heterogeneity = ct_data$intra_cell_type_heterogeneity,
        simpson_index = ct_data$simpson_index,
        evenness = ct_data$evenness,
        dominant_iso_prop = ct_data$dominant_iso_prop,
        dominant_iso_name = ct_data$dominant_iso_name,
        n_expressed_isoforms = ct_data$n_expressed_isoforms,
        pct_multi_iso_cells = ct_data$pct_multi_iso_cells
      )
      
      metrics_list[[cell_type]] <- df
    }
    
    # Combine all cell types for this gene
    if(length(metrics_list) > 0) {
      gene_metrics <- do.call(rbind, metrics_list)
      
      # Calculate mean values if requested
      if(include_mean) {
        mean_metrics <- data.frame(
          gene = gene,
          cell_type = "MEAN",
          intra_cellular_isoform_diversity = mean(gene_metrics$intra_cellular_isoform_diversity, na.rm = TRUE),
          inter_cellular_isoform_diversity = mean(gene_metrics$inter_cellular_isoform_diversity, na.rm = TRUE),
          intra_cell_type_heterogeneity = mean(gene_metrics$intra_cell_type_heterogeneity, na.rm = TRUE),
          simpson_index = mean(gene_metrics$simpson_index, na.rm = TRUE),
          evenness = mean(gene_metrics$evenness, na.rm = TRUE),
          dominant_iso_prop = mean(gene_metrics$dominant_iso_prop, na.rm = TRUE),
          dominant_iso_name = "MEAN",
          n_expressed_isoforms = mean(gene_metrics$n_expressed_isoforms, na.rm = TRUE),
          pct_multi_iso_cells = mean(gene_metrics$pct_multi_iso_cells, na.rm = TRUE)
        )
        
        # Add mean to the top of the gene's metrics
        gene_metrics <- rbind(mean_metrics, gene_metrics)
      }
      
      # Add to all metrics
      all_metrics[[gene]] <- gene_metrics
    }
  }
  
  # Combine all genes
  if(length(all_metrics) > 0) {
    combined_metrics <- do.call(rbind, all_metrics)
    return(combined_metrics)
  } else {
    warning("No valid metrics found for any of the specified genes")
    return(NULL)
  }
}

#' Plot stacked isoform usage profiles
#' 
#' Creates a stacked bar chart showing proportional isoform usage across
#' different cell types or developmental stages. This visualisation helps identify
#' cell type-specific isoform preferences.
#' 
#' @param scht_obj Single-Cell Hierarchical Tensor object
#' @param gene Gene name to analyse
#' @param cell_type_order Optional vector specifying the order of cell types
#' @param min_prop Minimum proportion to display (minor isoforms grouped as "Other")
#' @param colour_palette Colour palette to use (default: distinct qualitative palette)
#' 
#' @return A ggplot object that can be printed or saved
#' 
#' @examples
#' # Load example data and create SCHT
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Example 1: Plot isoform profile for a gene with multiple isoforms
#' plot_isoform_profile(scht_obj, "Mapk13")
#' 
#' # Example 2: Specify cell type order
#' plot_isoform_profile(
#'   scht_obj, 
#'   "Atl1",
#'   cell_type_order = c("AEC", "HEC", "T1_pre_HSC", "T2_pre_HSC", 
#'                       "E12", "E14", "Adult_HSC")
#' )
#' 
#' # Example 3: Try another gene with error handling
#' tryCatch({
#'   plot_isoform_profile(scht_obj, "Irf8", min_prop = 0.03)
#' }, error = function(e) {
#'   print(paste("Error:", e$message))
#' })
#' 
#' @export
plot_isoform_profile <- function(scht_obj, gene, cell_type_order = NULL, 
                                 min_prop = 0.05, colour_palette = NULL) {
  # Check if required packages are available
  required_packages <- c("ggplot2", "RColorBrewer")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is needed for this function. Please install it."))
    }
  }
  
  # Check if gene exists
  if(!gene %in% names(scht_obj$original_results)) {
    stop(paste("Gene", gene, "not found in SCHT object"))
  }
  
  # Get cell types
  if(is.null(cell_type_order)) {
    cell_types <- names(scht_obj$cell_type_matrices)
  } else {
    # Validate provided cell types
    valid_cell_types <- names(scht_obj$cell_type_matrices)
    if(!all(cell_type_order %in% valid_cell_types)) {
      warning("Some cell types in cell_type_order are not found in the SCHT object")
    }
    cell_types <- cell_type_order[cell_type_order %in% valid_cell_types]
  }
  
  # Collect isoform usage proportions for each cell_type
  usage_data <- list()
  for(cell_type in cell_types) {
    if(gene %in% names(scht_obj$cell_type_matrices[[cell_type]])) {
      iso_mat <- scht_obj$cell_type_matrices[[cell_type]][[gene]]
      
      # Skip cell types with insufficient data
      if(nrow(iso_mat) < 1 || ncol(iso_mat) == 0) next
      
      # Calculate average expression for each isoform in this cell type
      iso_means <- rowMeans(iso_mat, na.rm = TRUE)
      total_expr <- sum(iso_means)
      
      # Skip if no expression
      if(total_expr == 0) next
      
      # Calculate proportions
      props <- iso_means / total_expr
      usage_data[[cell_type]] <- props
    }
  }
  
  # Check if we have any data
  if(length(usage_data) == 0) {
    stop(paste("No cell type data found for gene", gene))
  }
  
  # Get all isoforms
  all_isoforms <- unique(unlist(lapply(usage_data, names)))
  
  # Prepare data for plotting
  plot_data <- data.frame(
    cell_type = character(0),
    isoform = character(0),
    proportion = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for(cell_type in names(usage_data)) {
    props <- usage_data[[cell_type]]
    
    # Group minor isoforms as "Other"
    minor_isos <- props[props < min_prop]
    major_isos <- props[props >= min_prop]
    
    if(length(minor_isos) > 0) {
      other_prop <- sum(minor_isos)
      if(other_prop > 0) {
        major_isos <- c(major_isos, setNames(other_prop, paste0("Other (expression < ", min_prop*100, "%)")))
      }
    }

    for(iso in names(major_isos)) {
      plot_data <- rbind(plot_data, 
                         data.frame(
                           cell_type = cell_type,
                           isoform = iso,
                           proportion = major_isos[iso],
                           stringsAsFactors = FALSE
                         ))
    }
  }
  
  # Order cell types if provided
  plot_data$cell_type <- factor(plot_data$cell_type, levels = names(usage_data))
  
  isoforms <- unique(plot_data$isoform)
  plot_data$isoform <- factor(plot_data$isoform, levels = isoforms)
  
  # Prepare colours for better differentiation
  isoforms <- unique(plot_data$isoform)
  n_colours <- length(isoforms)
  
  if(is.null(colour_palette)) {
    # Generate a colour palette
    if(n_colours <= 12) {
      # For 12 or fewer colours, use qualitative palette
      colours <- .get_colour_palette(n_colours, "Set3")
    } else {
      # For more colours, create a more differentiated palette
      # Combine palettes for more distinct colours
      pal1 <- .get_colour_palette(min(12, n_colours), "Paired")
      pal2 <- .get_colour_palette(min(12, max(3, n_colours-12)), "Set3")
      pal3 <- .get_colour_palette(min(8, max(3, n_colours-24)), "Set2")
      
      # Add viridis colours if available
      pal4 <- if (requireNamespace("viridis", quietly = TRUE)) {
        viridis::viridis(max(0, n_colours - 32), begin = 0.1, end = 0.9)
      } else {
        colorRampPalette(c("#440154", "#31688e", "#35b779", "#fde725"))(max(0, n_colours - 32))
      }
      
      colours <- c(pal1, pal2, pal3, pal4)
      colours <- colours[1:n_colours]
    }
  } else {
    # Use provided palette
    colours <- colour_palette
    if(length(colours) < n_colours) {
      warning("Provided colour palette has fewer colours than isoforms. Recycling colours.")
      colours <- rep(colours, length.out = n_colours)
    }
  }
  
  # Ensure "Other" is grey if present
  other_name <- paste0("Other (expression < ", min_prop*100, "%)")
  
  if (other_name %in% isoforms) {
    other_idx <- which(isoforms == other_name)
    colours[other_idx] <- "grey90"
  }
  
  # Create the stacked bar chart
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = cell_type, y = proportion, fill = isoform)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black", size = 0.2) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.5),
      axis.title = ggplot2::element_text(face = "bold", size = 14),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 14, hjust = 0.5),
      legend.title = ggplot2::element_text(face = "bold", size = 12),
      legend.text = ggplot2::element_text(size = 10),
      legend.position = "right",
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key.size = ggplot2::unit(0.8, "cm"),
      axis.text = ggplot2::element_text(colour = "black", size = 12),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.spacing.y = ggplot2::unit(0.2, "cm"),
      legend.key = ggplot2::element_rect(color = "black", size = 0.2)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(byrow = TRUE, ncol = 1, keywidth = ggplot2::unit(0.8, "cm"), 
                                                 keyheight = ggplot2::unit(0.5, "cm"))) +
    ggplot2::labs(title = paste("Isoform Usage Profile for", gene),
                  x = "Cell Type",
                  y = "Proportion",
                  fill = "Isoform")
  
  return(p)
}

#' Plot isoform usage transitions across cell types
#' 
#' Creates a line plot showing how isoform proportions change across
#' different cell types. This visualisation is particularly useful for
#' developmental trajectories or other ordered cell type progressions.
#' 
#' @param scht_obj Single-Cell Hierarchical Tensor object
#' @param gene Gene name to analyse
#' @param cell_type_order Vector specifying the order of cell types (required for this plot)
#' @param selected_isoforms Optional vector of isoform names to plot. If provided, only these isoforms will be shown (ignoring min_prop)
#' @param min_prop Minimum proportion to display (default: 0.05). Ignored if selected_isoforms is provided
#' @param smooth Apply smoothing to lines (default: TRUE)
#' @param colour_palette Colour palette to use (default: viridis plasma)
#' 
#' @return A ggplot object that can be printed or saved
#' 
#' @examples
#' # Load example data and create SCHT
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Example 1: Create isoform transition plot
#' plot_isoform_transitions(
#'   scht_obj, 
#'   gene = "Atl1",
#'   cell_type_order = c("AEC", "HEC", "T1_pre_HSC", "T2_pre_HSC", 
#'                       "E12", "E14", "Adult_HSC")
#' )
#' 
#' # Example 2: Another gene with custom colors
#' plot_isoform_transitions(
#'   scht_obj, 
#'   gene = "Mapk13",
#'   cell_type_order = c("AEC", "HEC", "E12", "E14", "Adult_HSC"),
#'   colour_palette = RColorBrewer::brewer.pal(5, "Set2")
#' )
#' 
#' @export
plot_isoform_transitions <- function(scht_obj, gene, cell_type_order, 
                                     selected_isoforms = NULL,
                                     min_prop = 0.05, smooth = TRUE, colour_palette = NULL) {
  # Check if required packages are available
  required_packages <- c("ggplot2", "RColorBrewer")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is needed for this function. Please install it."))
    }
  }
  
  # Check if gene exists
  if(!gene %in% names(scht_obj$original_results)) {
    stop(paste("Gene", gene, "not found in SCHT object"))
  }
  
  # Cell type order is required for this plot
  if(is.null(cell_type_order)) {
    stop("cell_type_order must be provided for transition plot")
  }
  
  # Validate provided cell types
  valid_cell_types <- names(scht_obj$cell_type_matrices)
  if(!all(cell_type_order %in% valid_cell_types)) {
    warning("Some cell types in cell_type_order are not found in the SCHT object")
  }
  cell_types <- cell_type_order[cell_type_order %in% valid_cell_types]
  
  if(length(cell_types) < 2) {
    stop("Need at least 2 valid cell types for transition plot")
  }
  
  # Collect isoform usage proportions for each cell type
  usage_data <- list()
  for(cell_type in cell_types) {
    if(gene %in% names(scht_obj$cell_type_matrices[[cell_type]])) {
      iso_mat <- scht_obj$cell_type_matrices[[cell_type]][[gene]]
      
      # Skip cell types with insufficient data
      if(nrow(iso_mat) < 1 || ncol(iso_mat) == 0) next
      
      # Calculate average expression for each isoform in this cell type
      iso_means <- rowMeans(iso_mat, na.rm = TRUE)
      total_expr <- sum(iso_means)
      
      # Skip if no expression
      if(total_expr == 0) next
      
      # Calculate proportions
      props <- iso_means / total_expr
      usage_data[[cell_type]] <- props
    }
  }
  
  # Check if we have enough data
  if(length(usage_data) < 2) {
    stop(paste("Not enough cell type data found for gene", gene))
  }
  
  # Get all isoforms and determine which to plot
  all_isoforms <- unique(unlist(lapply(usage_data, names)))
  
  if (!is.null(selected_isoforms)) {
    # Use user-specified isoforms
    major_isoforms <- selected_isoforms[selected_isoforms %in% all_isoforms]
    if (length(major_isoforms) == 0) {
      stop("None of the selected isoforms are found in the data")
    }
  } else {
    # Filter by min_prop
    major_isoforms <- c()
    
    for(iso in all_isoforms) {
      # Check if isoform exceeds min_prop in any cell type
      for(cell_type in names(usage_data)) {
        props <- usage_data[[cell_type]]
        if(iso %in% names(props) && props[iso] >= min_prop) {
          major_isoforms <- c(major_isoforms, iso)
          break
        }
      }
    }
    
    major_isoforms <- unique(major_isoforms)
  }
  
  # Prepare data for plotting
  plot_data <- data.frame(
    cell_type = character(0),
    isoform = character(0),
    proportion = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for(cell_type in names(usage_data)) {
    props <- usage_data[[cell_type]]
    
    # Add major isoforms
    for(iso in major_isoforms) {
      prop <- ifelse(iso %in% names(props), props[iso], 0)
      plot_data <- rbind(plot_data, 
                         data.frame(
                           cell_type = cell_type,
                           isoform = iso,
                           proportion = prop,
                           stringsAsFactors = FALSE
                         ))
    }
    
    # Add "Other" if needed (only when not using selected_isoforms)
    if (is.null(selected_isoforms)) {
      other_prop <- 1 - sum(props[names(props) %in% major_isoforms])
      if(other_prop > 0.01) { 
        plot_data <- rbind(plot_data, 
                           data.frame(
                             cell_type = cell_type,
                             isoform = "Other",
                             proportion = other_prop,
                             stringsAsFactors = FALSE
                           ))
      }
    }
  }
  
  # Order cell types
  plot_data$cell_type <- factor(plot_data$cell_type, levels = cell_types)
  
  # Create cell type index for x-axis
  cell_type_idx <- data.frame(
    cell_type = cell_types,
    idx = 1:length(cell_types),
    stringsAsFactors = FALSE
  )
  plot_data <- merge(plot_data, cell_type_idx, by = "cell_type")
  
  isoforms <- unique(plot_data$isoform)
  n_colours <- length(isoforms)
  
 if(is.null(colour_palette)) {
    if(n_colours <= 8) {
      # Palette for up to 8 colours
      colours <- c(
        "#4878D0",  # Soft blue
        "#EE854A",  # Soft orange
        "#6ACC64",  # Soft green
        "#D65F5F",  # Soft red
        "#956CB4",  # Soft purple
        "#8C613C",  # Soft brown
        "#DC7EC0",  # Soft pink
        "#797979"   # Soft grey
      )[1:min(8, n_colours)]
    } else if(n_colours <= 12) {
      # Extended palette for up to 12 colours
      colours <- c(
        "#4878D0",  # Soft blue
        "#EE854A",  # Soft orange
        "#6ACC64",  # Soft green
        "#D65F5F",  # Soft red
        "#956CB4",  # Soft purple
        "#8C613C",  # Soft brown
        "#DC7EC0",  # Soft pink
        "#797979",  # Soft grey
        "#82C6E2",  # Light blue
        "#FFB481",  # Light orange
        "#ACCF91",  # Light green
        "#B5A5D5"   # Light purple
      )[1:min(12, n_colours)]
    } else {
      # For more colours, use a combination of softer palettes
      base_colours <- c(
        "#4878D0", "#EE854A", "#6ACC64", "#D65F5F", 
        "#956CB4", "#8C613C", "#DC7EC0", "#797979",
        "#82C6E2", "#FFB481", "#ACCF91", "#B5A5D5"
      )
      
      # Add pastel colours
      pastel_colours <- .get_colour_palette(8, "Pastel1")
      
      # Add some colours from Pastel2
      pastel2_colours <- .get_colour_palette(8, "Pastel2")
      
      # For very many colours, add some muted viridis colours
      viridis_muted <- if (requireNamespace("viridis", quietly = TRUE)) {
        viridis::viridis(max(0, n_colours - 28), begin = 0.1, end = 0.9, alpha = 0.7)
      } else {
        colorRampPalette(c("#440154", "#31688e", "#35b779", "#fde725"))(max(0, n_colours - 28))
      }
      
      # Combine and make sure we have enough colours
      colours <- c(base_colours, pastel_colours, pastel2_colours, viridis_muted)
      # Ensure we have enough colours and take only what we need
      colours <- colours[1:n_colours]
    }
  } else {
    # Use provided palette
    colours <- colour_palette
    if(length(colours) < n_colours) {
      warning("Provided colour palette has fewer colours than isoforms. Recycling colours.")
      colours <- rep(colours, length.out = n_colours)
    }
  }
  
  # Ensure "Other" is grey if present
  if("Other" %in% isoforms) {
    other_idx <- which(isoforms == "Other")
    colours[other_idx] <- "#7F7F7F" 
  }
  

  # Check if we have enough data points for smoothing
  n_timepoints <- length(unique(plot_data$idx))
  use_smooth <- smooth && n_timepoints >= 4  # Need at least 4 points for reliable smoothing
  
  if (smooth && !use_smooth) {
    message("Note: Smoothing disabled due to insufficient data points (< 4 cell types)")
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = idx, y = proportion, color = isoform, group = isoform)) +
    ggplot2::geom_point(size = 3, alpha = 0.9, stroke = 0.5) +
    {if(use_smooth) ggplot2::geom_smooth(se = FALSE, span = 0.7, method = "loess", size = 1.2) else ggplot2::geom_line(size = 1.2)}
  
  # Add labels and styling
  p <- p +
    ggplot2::scale_x_continuous(breaks = cell_type_idx$idx, labels = cell_type_idx$cell_type) +
    ggplot2::scale_color_manual(values = colours) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.5),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90", linetype = "dashed", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold", size = 18),
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 14, hjust = 0.5),
      legend.title = ggplot2::element_text(face = "bold", size = 12),
      legend.text = ggplot2::element_text(size = 10),
      legend.position = "right",
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key.size = ggplot2::unit(0.8, "cm"),
      axis.text = ggplot2::element_text(colour = "black", size = 14),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    ggplot2::labs(title = paste("Isoform Usage Transitions for", gene),
                  x = "Cell Type",
                  y = "Proportion",
                  colour = "Isoform")
  
  return(p)
}

#' Create gene complexity radar chart
#' 
#' Creates a radar chart showing multiple complexity metrics for selected genes.
#' This visualisation provides a compact way to compare different genes across
#' multiple complexity dimensions simultaneously.
#' 
#' @param tc_metrics Data frame from calculate_isoform_complexity_metrics() or the results object
#' @param genes Vector of gene names to include
#' @param scale_type How to scale metrics: "global" (across all genes), "per_metric" (each metric independently), or "none"
#' 
#' @return A ggplot object that can be printed or saved
#' 
#' @examples
#' # Load example data and create SCHT
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Calculate complexity metrics
#' tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)
#' 
#' # Create radar plot for specific genes
#' # Use genes from the middle of the results
#' available_genes <- rownames(tc_results$metrics)
#' genes_to_plot <- available_genes[c(100, 200, 300, 400, 500)]
#' 
#' # Note: This example requires the 'ggradar' package
#' if(requireNamespace("ggradar", quietly = TRUE)) {
#'   plot_complexity_radar(tc_results, genes_to_plot)
#' }
#' 
#' @export
plot_complexity_radar <- function(tc_metrics, genes, scale_type = "global") {
  # Check for required packages
  if (!requireNamespace("ggradar", quietly = TRUE)) {
    stop("Package 'ggradar' is needed for this function. Please install it with: install.packages('ggradar')")
  }
  
  # Extract complexity metrics from results object if needed
  if(is.data.frame(tc_metrics)) {
    metrics <- tc_metrics
  } else if("metrics" %in% names(tc_metrics)) {
    metrics <- tc_metrics$metrics
  } else {
    stop("Input must be a transcriptomic complexity results object or a data frame with metrics")
  }
  
  # Ensure provided genes exist
  valid_genes <- genes[genes %in% metrics$gene]
  if(length(valid_genes) == 0) {
    stop("None of the specified genes found in metrics")
  }
  
  # If some genes weren't found, print a warning
  if(length(valid_genes) < length(genes)) {
    missing_genes <- setdiff(genes, valid_genes)
    warning(paste("The following genes were not found:", paste(missing_genes, collapse = ", ")))
  }
  
  # Select valid genes
  gene_data <- metrics[metrics$gene %in% valid_genes, ]
  
  # Define metrics to include in radar chart
  radar_metrics <- c(
    "intra_cellular_isoform_diversity",
    "inter_cellular_isoform_diversity", 
    "intra_cell_type_heterogeneity",
    "inter_cell_type_specificity",
    "intra_cell_type_heterogeneity_variability",
    "inter_cell_type_difference_variability",
    "cell_type_coexpression_variability"
  )
  
  # Create user-friendly labels for metrics
  metric_labels <- c(
    "Intra-cellular\nIsoform Diversity",
    "Inter-cellular\nIsoform Diversity",
    "Intra-cell-type\nHeterogeneity",
    "Inter-cell-type\nSpecificity",
    "Heterogeneity\nVariability",
    "Difference\nVariability",
    "Co-expression\nVariability"
  )
  
  # Check if all required metrics are in the data
  missing_metrics <- setdiff(radar_metrics, colnames(gene_data))
  if(length(missing_metrics) > 0) {
    stop(paste("The following metrics are missing from the data:", 
               paste(missing_metrics, collapse = ", ")))
  }
  
  # Create radar data with only gene column and radar metrics
  radar_data <- gene_data[, c("gene", radar_metrics)]
  
  # Process each metric
  for(i in seq_along(radar_metrics)) {
    col <- radar_metrics[i]
    
    # Replace NAs with 0
    radar_data[[col]][is.na(radar_data[[col]])] <- 0
    
    # Check if all values are 0 or NA
    if(max(radar_data[[col]], na.rm = TRUE) <= 0) {
      warning(paste("Metric", col, "has all zero or NA values. Setting to small default value."))
      radar_data[[col]] <- 0.001  # Small default value
    }
  }
  
  # Scale values based on scaling type
  if(scale_type == "global") {
    # Find global min and max across all metrics
    all_values <- unlist(radar_data[, radar_metrics])
    min_val <- min(all_values, na.rm = TRUE)
    max_val <- max(all_values, na.rm = TRUE)
    
    # Apply global scaling to all metrics
    for(col in radar_metrics) {
      if(max_val > min_val) {
        radar_data[[col]] <- (radar_data[[col]] - min_val) / (max_val - min_val)
      } else {
        radar_data[[col]] <- 0.5  # Default if no variation
      }
    }
  } else if(scale_type == "per_metric") {
    # Scale each metric independently to 0-1
    for(col in radar_metrics) {
      min_val <- min(radar_data[[col]], na.rm = TRUE)
      max_val <- max(radar_data[[col]], na.rm = TRUE)
      
      if(max_val > min_val) {
        radar_data[[col]] <- (radar_data[[col]] - min_val) / (max_val - min_val)
      } else {
        radar_data[[col]] <- 0.5  # Default if no variation
      }
    }
  }
  # If scale_type is "none", we leave the values as they are
  
  # Create a new dataframe with renamed columns for the plot
  plot_data <- radar_data
  colnames(plot_data) <- c("gene", metric_labels)
  
  # Remove rows with any NA values to prevent ggradar errors
  complete_rows <- complete.cases(plot_data[, -1])  # Exclude gene column
  if(!all(complete_rows)) {
    removed_genes <- plot_data$gene[!complete_rows]
    warning(paste("Removing genes with NA values:", paste(removed_genes, collapse = ", ")))
    plot_data <- plot_data[complete_rows, ]
  }
  
  # Check if any data remains
  if(nrow(plot_data) == 0) {
    stop("No complete data available for plotting after removing NA values")
  }
  
  # Create custom colours for genes
  n_genes <- nrow(plot_data)
  valid_genes <- plot_data$gene
  if(n_genes <= 8) {
    gene_colours <- .get_colour_palette(n_genes, "Set1")
  } else {
    # For more genes, use colorRampPalette
    base_colours <- .get_colour_palette(8, "Set1")
    gene_colours <- colorRampPalette(base_colours)(n_genes)
  }
  
  # Create the radar chart
  radar_plot <- ggradar::ggradar(
    plot_data,
    grid.min = 0, 
    grid.mid = 0.5, 
    grid.max = 1,
    axis.label.size = 5,
    grid.label.size = 4,
    group.point.size = 3,
    group.line.width = 1.2,       
    background.circle.colour = "grey95",
    legend.position = "right",
    legend.text.size = 9,
    base.size = 15,  
    axis.label.offset = 1.15,  
    gridline.mid.colour = "grey65",
    group.colours = gene_colours  
  ) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 18, hjust = 0.5, face = "bold"),
      plot.margin = ggplot2::margin(5, 10, 5, 10)
    ) +
    ggplot2::labs(colour = "Gene",
                  title = "Gene Complexity Comparison")  
  
  return(radar_plot)
}

#' Create radar chart for a single gene across cell types with enhanced visualization
#' 
#' An enhanced version of the original plot_cell_type_complexity_radar function with
#' better visual elements and higher resolution output. This function helps
#' visualise how complexity metrics vary across different cell types for a single gene.
#' 
#' @param tc_results Transcriptomic complexity results object
#' @param gene_name Gene name to analyse
#' @param metrics Vector of metric names to include in radar chart
#' @param scale_values Whether to scale values from 0-1 for better comparability
#' 
#' @return A ggplot object that can be printed or saved
#' 
#' @examples
#' # Load example data and create SCHT
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Calculate complexity metrics
#' tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)
#' 
#' # Find a gene with cell type-specific patterns
#' high_specificity_genes <- names(tc_results$metrics$inter_cell_type_specificity)[
#'   tc_results$metrics$inter_cell_type_specificity > 0.5
#' ]
#' 
#' # Note: This example requires the 'ggradar' package
#' if(requireNamespace("ggradar", quietly = TRUE) && length(high_specificity_genes) > 0) {
#'   plot_single_gene_radar_cell_type(tc_results, high_specificity_genes[1])
#' }
#' 
#' @export
plot_single_gene_radar_cell_type <- function(tc_results, gene_name, 
                                             metrics = NULL, scale_values = TRUE) {
  # Check if ggradar is installed
  if (!requireNamespace("ggradar", quietly = TRUE)) {
    stop("Package 'ggradar' is needed for this function. Please install it with: install.packages('ggradar')")
  }
  
  # Set default metrics if not provided
  if(is.null(metrics)) {
    metrics <- c("intra_cellular_isoform_diversity", "inter_cellular_isoform_diversity", 
                 "intra_cell_type_heterogeneity", "dominant_iso_prop", "n_expressed_isoforms")
  }
  
  # Validate gene name
  if(!(gene_name %in% names(tc_results$cell_type_metrics))) {
    stop(paste("Gene", gene_name, "not found in the data"))
  }
  
  # Extract gene data
  gene_data <- tc_results$cell_type_metrics[[gene_name]]
  
  metric_labels <- c(
    "intra_cellular_isoform_diversity" = "Intra-cellular\nisoform diversity",
    "inter_cellular_isoform_diversity" = "Inter-cellular\nisoform diversity",
    "intra_cell_type_heterogeneity" = "Intra-cell-type\nheterogeneity",
    "dominant_iso_prop" = "Dominant\nisoform prop",
    "n_expressed_isoforms" = "Isoform\ncount"
  )
  
  # Use provided labels if available, otherwise create default ones
  readable_metrics <- sapply(metrics, function(m) {
    if(m %in% names(metric_labels)) {
      return(metric_labels[[m]])
    } else {
      # Create a default multiline label
      label <- gsub("_", " ", m)
      label <- paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label)))
      
      if(nchar(label) > 10) {
        words <- strsplit(label, " ")[[1]]
        if(length(words) > 1) {
          mid_point <- ceiling(length(words)/2)
          label <- paste(paste(words[1:mid_point], collapse = " "), 
                         paste(words[(mid_point+1):length(words)], collapse = " "),
                         sep = "\n")
        } else {
          mid_char <- ceiling(nchar(label)/2)
          label <- paste(substr(label, 1, mid_char),
                         substr(label, mid_char+1, nchar(label)),
                         sep = "\n")
        }
      }
      return(label)
    }
  })
  
  # Initialize empty data frame for radar data
  radar_data <- data.frame(cell_type = character(), stringsAsFactors = FALSE)
  
  # Collect data for each cell type
  for(cell_type in names(gene_data)) {
    ct_metrics <- gene_data[[cell_type]]
    
    # Create a row for this cell type
    cell_type_row <- data.frame(cell_type = cell_type, stringsAsFactors = FALSE)
    
    # Extract values for specified metrics
    for(i in seq_along(metrics)) {
      metric <- metrics[i]
      readable_metric <- readable_metrics[i]
      
      if(metric %in% names(ct_metrics)) {
        cell_type_row[[readable_metric]] <- ct_metrics[[metric]]
      } else {
        cell_type_row[[readable_metric]] <- NA
      }
    }
    
    # Add to radar data
    radar_data <- rbind(radar_data, cell_type_row)
  }
  
  # Check if we have sufficient data
  if(nrow(radar_data) == 0) {
    stop("No valid metric data found for gene ", gene_name)
  }
  
  # Scale values if requested
  if(scale_values) {
    for(metric in readable_metrics) {
      metric_data <- radar_data[[metric]]
      min_val <- min(metric_data, na.rm = TRUE)
      max_val <- max(metric_data, na.rm = TRUE)
      
      # Only scale if we have variation
      if(max_val > min_val) {
        radar_data[[metric]] <- (metric_data - min_val) / (max_val - min_val)
      } else {
        # If no variation
        radar_data[[metric]] <- 0.5
      }
    }
  }
  
  # Remove rows with any NA values to prevent ggradar errors
  complete_rows <- complete.cases(radar_data[, readable_metrics])
  if(!all(complete_rows)) {
    removed_cell_types <- radar_data$cell_type[!complete_rows]
    warning(paste("Removing cell types with NA values:", paste(removed_cell_types, collapse = ", ")))
    radar_data <- radar_data[complete_rows, ]
  }
  
  # Check if any data remains
  if(nrow(radar_data) == 0) {
    stop("No complete data available for plotting after removing NA values")
  }
  
  # Determine colour palette based on number of cell types
  n_cell_types <- nrow(radar_data)
  if(n_cell_types <= 8) {
    colour_palette <- .get_colour_palette(n_cell_types, "Set2")
  } else {
    # For more cell types, use colorRampPalette
    base_colours <- .get_colour_palette(8, "Set2")
    colour_palette <- colorRampPalette(base_colours)(n_cell_types)
  }
  
  # Create the radar chart with improved visual elements
  p <- ggradar::ggradar(
    radar_data,
    base.size = 15,
    font.radar = "sans",
    grid.min = 0,
    grid.mid = 0.5,
    grid.max = 1,
    axis.label.size = 6,
    axis.label.offset = 1.15,
    label.centre.y = FALSE,
    axis.line.colour = "grey60",
    grid.line.width = 0.8,
    group.line.width = 1.2,
    group.point.size = 3,
    group.colours = colour_palette,
    background.circle.colour = "white",
    legend.position = "right",
    legend.text.size = 11
  ) +
    ggplot2::labs(
      title = paste0("Complexity Profile for ", gene_name, " Across Cell Types"),
      subtitle = "Radar chart showing complexity metrics across cell types"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 14, hjust = 0.5),
      axis.text = ggplot2::element_text(face = "bold", size = 12),
      legend.text = ggplot2::element_text(size = 11)
    ) 
  
  return(p)
}

#' Create radar charts by cell type with better label handling
#' 
#' Creates a separate radar chart for each cell type comparing multiple genes,
#' with improved label handling and a single combined legend for all plots.
#' This function is useful for comparing multiple genes across different 
#' cell types in a structured grid format.
#' 
#' @param tc_results Transcriptomic complexity results object
#' @param gene_names Vector of gene names to compare
#' @param cell_types Vector of cell types to include (if NULL, all detected cell types are used)
#' @param metrics Vector of metric names to include in radar chart
#' @param scale_type How to scale metrics: "global" (across all cell types) or "per_cell_type" (each cell type independently). Default is "per_cell_type"
#' @param ncol Number of columns for arranging the plots in a grid
#' 
#' @return A combined plot grid with informative legend and clean labels
#' 
#' @examples
#' # Load example data and create SCHT
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Calculate complexity metrics
#' \dontrun{
#' tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)
#' 
#' # Select genes that exist in the complexity results
#' available_genes <- names(tc_results$metrics$inter_cellular_isoform_diversity)
#' # Pick genes at different complexity levels
#' n_genes <- length(available_genes)
#' genes_of_interest <- available_genes[c(1, 
#'                                       round(n_genes/4), 
#'                                       round(n_genes/2))]
#' 
#' # Note: This example requires the 'ggradar' package
#' if(requireNamespace("ggradar", quietly = TRUE) && 
#'    length(genes_of_interest) > 0) {
#'   # Check that genes exist in cell type data
#'   genes_to_plot <- genes_of_interest[genes_of_interest %in% 
#'                                     rownames(tc_results$cell_type_metrics[[1]])]
#'   
#'   if(length(genes_to_plot) > 0) {
#'     # Create multi-gene radar plots with per-cell-type scaling
#'     plot_compare_multiple_genes_radar_cell_type(tc_results, genes_to_plot)
#'     
#'     # Use global scaling for direct comparison
#'     plot_compare_multiple_genes_radar_cell_type(tc_results, genes_to_plot, 
#'                                                scale_type = "global")
#'   } else {
#'     cat("No suitable genes found for plotting\n")
#'   }
#' }
#' }
#' 
#' @export
plot_compare_multiple_genes_radar_cell_type <- function(tc_results, gene_names, cell_types = NULL,
                                                   metrics = NULL, scale_type = "per_cell_type",
                                                   ncol = 3) {
  # Check for required packages
  if (!requireNamespace("ggradar", quietly = TRUE)) {
    stop("Package 'ggradar' is needed for this function. Please install it with: install.packages('ggradar')")
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is needed for this function. Please install it with: install.packages('cowplot')")
  }
  
  # Validate scale_type parameter
  valid_scale_types <- c("global", "per_cell_type")
  if (!scale_type %in% valid_scale_types) {
    warning(paste("Invalid scale_type '", scale_type, "'. Using 'per_cell_type' instead.", sep = ""))
    scale_type <- "per_cell_type"
  }
  
  # Set default metrics if not provided
  if(is.null(metrics)) {
    metrics <- c("intra_cellular_isoform_diversity", "inter_cellular_isoform_diversity", 
                 "intra_cell_type_heterogeneity", "dominant_iso_prop", "n_expressed_isoforms")
  }
  
  # Validate gene names
  valid_genes <- gene_names[gene_names %in% names(tc_results$cell_type_metrics)]
  if(length(valid_genes) == 0) {
    stop("None of the provided gene names were found in the data")
  }
  if(length(valid_genes) < length(gene_names)) {
    missing_genes <- setdiff(gene_names, valid_genes)
    warning(paste("The following genes were not found in the data:", 
                  paste(missing_genes, collapse = ", ")))
  }
  
  # Get all unique cell types if not specified
  if(is.null(cell_types)) {
    all_cell_types <- unique(unlist(lapply(valid_genes, function(gene) {
      names(tc_results$cell_type_metrics[[gene]])
    })))
    cell_types <- all_cell_types
  }
  
  # Make short but descriptive metric names for radar plots
  metric_labels <- c(
    "intra_cellular_isoform_diversity" = "Intra-cellular\nisoform diversity",
    "inter_cellular_isoform_diversity" = "Inter-cellular\nisoform\ndiversity",
    "intra_cell_type_heterogeneity" = "Intra-cell-type\nheterogeneity",
    "dominant_iso_prop" = "Dominant\nisoform prop",
    "n_expressed_isoforms" = "Isoform\ncount"
  )
  
  # Use provided labels if available, otherwise create default ones
  readable_metrics <- sapply(metrics, function(m) {
    if(m %in% names(metric_labels)) {
      return(metric_labels[[m]])
    } else {
      # Create a default multiline label
      label <- gsub("_", " ", m)
      label <- paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label)))
      
      if(nchar(label) > 10) {
        words <- strsplit(label, " ")[[1]]
        if(length(words) > 1) {
          mid_point <- ceiling(length(words)/2)
          label <- paste(paste(words[1:mid_point], collapse = " "), 
                         paste(words[(mid_point+1):length(words)], collapse = " "),
                         sep = "\n")
        } else {
          mid_char <- ceiling(nchar(label)/2)
          label <- paste(substr(label, 1, mid_char),
                         substr(label, mid_char+1, nchar(label)),
                         sep = "\n")
        }
      }
      return(label)
    }
  })
  
  # Collect all data first to prepare for scaling
  all_data <- list()
  
  # Process each cell type and collect data
  for(cell_type in cell_types) {
    # Create data frame for this cell type
    cell_data <- data.frame(gene = character(), stringsAsFactors = FALSE)
    
    # Add metrics as columns
    for(i in seq_along(metrics)) {
      cell_data[[readable_metrics[i]]] <- numeric()
    }
    
    # Collect data for each gene in this cell type
    for(gene in valid_genes) {
      if(!gene %in% names(tc_results$cell_type_metrics)) {
        next
      }
      
      gene_data <- tc_results$cell_type_metrics[[gene]]
      
      # Check if this cell type exists for this gene
      if(cell_type %in% names(gene_data)) {
        ct_metrics <- gene_data[[cell_type]]
        
        # Create a row for this gene
        gene_row <- data.frame(gene = gene, stringsAsFactors = FALSE)
        
        # Extract values for each metric
        has_data <- FALSE
        for(i in seq_along(metrics)) {
          metric <- metrics[i]
          if(metric %in% names(ct_metrics) && !is.na(ct_metrics[[metric]])) {
            gene_row[[readable_metrics[i]]] <- ct_metrics[[metric]]
            has_data <- TRUE
          } else {
            gene_row[[readable_metrics[i]]] <- NA
          }
        }
        
        # Only add if we have at least one valid metric
        if(has_data) {
          cell_data <- rbind(cell_data, gene_row)
        }
      }
    }
    
    # Skip cell types with no data
    if(nrow(cell_data) > 0) {
      all_data[[cell_type]] <- cell_data
    } else {
      warning(paste("No data available for cell type:", cell_type))
    }
  }
  
  # Create colour palette for genes (consistent across all plots)
  n_genes <- length(valid_genes)
  if(n_genes <= 8) {
    gene_colours <- .get_colour_palette(n_genes, "Set1")
  } else {
    # For more genes, use colorRampPalette
    base_colours <- .get_colour_palette(8, "Set1")
    gene_colours <- colorRampPalette(base_colours)(n_genes)
  }
  names(gene_colours) <- valid_genes
  
  # Apply scaling based on scale_type
  if(scale_type == "global") {
    # For global scaling, find global min/max across all metrics and cell types
    all_values <- unlist(lapply(all_data, function(cell_data) {
      unlist(cell_data[, readable_metrics])
    }))
    
    global_min <- min(all_values, na.rm = TRUE)
    global_max <- max(all_values, na.rm = TRUE)
    
    # Apply global scaling to all data
    for(cell_type in names(all_data)) {
      for(metric in readable_metrics) {
        if(global_max > global_min) {
          all_data[[cell_type]][[metric]] <- (all_data[[cell_type]][[metric]] - global_min) / (global_max - global_min)
        } else {
          all_data[[cell_type]][[metric]] <- 0.5
        }
      }
    }
  } else {
    # For per_cell_type scaling (default), each cell type is scaled independently
    for(cell_type in names(all_data)) {
      for(metric in readable_metrics) {
        metric_values <- all_data[[cell_type]][[metric]]
        metric_min <- min(metric_values, na.rm = TRUE)
        metric_max <- max(metric_values, na.rm = TRUE)
        
        if(metric_max > metric_min) {
          all_data[[cell_type]][[metric]] <- (metric_values - metric_min) / (metric_max - metric_min)
        } else {
          all_data[[cell_type]][[metric]] <- 0.5
        }
      }
    }
  }
  
  # Create a list to store plots for each cell type
  cell_type_plots <- list()
  
  # Create a separate label guide plot
  label_data <- data.frame(
    gene = "Guide",
    stringsAsFactors = FALSE
  )
  
  # Add placeholder values for each metric
  for(i in seq_along(metrics)) {
    label_data[[readable_metrics[i]]] <- 0.5
  }
  
  # Create the label guide plot with proper grid values
  guide_plot <- ggradar::ggradar(
    label_data,
    base.size = 16,
    font.radar = "sans",
    values.radar = c("", "", ""),
    grid.min = 0,
    grid.mid = 0.5,
    grid.max = 1,
    axis.label.size = 4,
    axis.label.offset = 0.8,
    label.centre.y = FALSE, 
    axis.line.colour = "grey60",
    grid.line.width = 0.5,
    group.line.width = 0,
    group.point.size = 0,
    group.colours = "transparent",
    background.circle.colour = "white",
    legend.position = "none",
    plot.legend = FALSE
  ) +
    ggplot2::labs(title = "Metric Guide") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  # Create radar plots for each cell type
  for(cell_type in names(all_data)) {
    cell_data <- all_data[[cell_type]]
    
    # Remove rows with any NA values to prevent ggradar errors
    complete_rows <- complete.cases(cell_data[, readable_metrics])
    if(!all(complete_rows)) {
      removed_genes <- cell_data$gene[!complete_rows]
      warning(paste("Removing genes with NA values in cell type", cell_type, ":", 
                    paste(removed_genes, collapse = ", ")))
      cell_data <- cell_data[complete_rows, ]
    }
    
    # Skip if no data remains after NA removal
    if(nrow(cell_data) == 0) {
      warning(paste("No complete data for cell type:", cell_type, "- skipping"))
      next
    }
    
    # Ensure gene colours match the valid genes in this cell type
    genes_in_cell_type <- cell_data$gene
    cell_gene_colours <- gene_colours[names(gene_colours) %in% genes_in_cell_type]
    
    # Create radar plot without axis labels
    p <- ggradar::ggradar(
      cell_data,
      base.size = 15,
      font.radar = "sans",
      grid.min = 0,
      grid.mid = 0.5,
      grid.max = 1,
      axis.label.size = 0,  
      axis.label.offset = 1.15,
      label.centre.y = FALSE,
      axis.line.colour = "grey60",
      grid.line.width = 0.9,
      group.line.width = 1.2,
      group.point.size = 3,
      group.colours = cell_gene_colours,
      background.circle.colour = "white",
      legend.position = "none",
      legend.text.size = 10
    ) +
      ggplot2::labs(
        title = paste0(cell_type)
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5)
      )
    
    # Store the plot
    cell_type_plots[[cell_type]] <- p
  }
  
  # Check if we have any plots
  if(length(cell_type_plots) == 0) {
    stop("No plots were generated. Please check your input data.")
  }
  
  # Create a common legend for all plots
  legend_data <- data.frame(
    gene = valid_genes,
    stringsAsFactors = FALSE
  )
  
  # Add dummy values for each metric
  for(i in seq_along(metrics)) {
    legend_data[[readable_metrics[i]]] <- 0.5
  }
  
  # Create a plot just to extract the legend
  legend_plot <- ggradar::ggradar(
    legend_data,
    base.size = 15,
    group.colours = gene_colours[valid_genes],
    axis.label.offset = 1,
    legend.position = "right",
    legend.text.size = 12
  ) +
    ggplot2::labs(colour = "Gene") +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 10)
    )
  
  # Extract the legend
  tryCatch({
    legend <- cowplot::get_legend(legend_plot)
  }, error = function(e) {
    warning("Failed to extract legend, using a default guide area instead.")
    return(NULL)
  })
  
  # Add the guide plot to the list of plots
  all_plots <- c(list(guide = guide_plot), cell_type_plots)
  
  # Create the main plot grid
  main_grid <- cowplot::plot_grid(
    plotlist = all_plots,
    ncol = ncol,
    align = "hv"
  )
  
  # Create title and subtitle
  title <- cowplot::ggdraw() + 
    cowplot::draw_label(
      "Gene Comparison Across Cell Types",
      fontface = "bold",
      size = 18,
      x = 0.5,
      hjust = 0.5
    )
  
  subtitle <- cowplot::ggdraw() + 
    cowplot::draw_label(
      "Each panel shows a different cell type",
      size = 14,
      x = 0.5,
      hjust = 0.5
    )
  
  # Combine title, subtitle, and main plot grid
  titled_grid <- cowplot::plot_grid(
    title, subtitle, main_grid,
    ncol = 1,
    rel_heights = c(0.1, 0.05, 1)
  )
  
  # Add the legend to the right if available
  if(!is.null(legend)) {
    final_plot <- cowplot::plot_grid(
      titled_grid, legend,
      rel_widths = c(0.85, 0.15),
      ncol = 2
    )
  } else {
    final_plot <- titled_grid
  }
  
  return(final_plot)
}

#' Compare transcriptomic complexity density differences between groups
#' 
#' This function calculates and visualises the density differences of genes
#' in a 2D metric space between different groups or conditions. It helps identify
#' regions where gene distributions shift between groups.
#' 
#' @param tc_results_list List of transcriptomic complexity results objects
#' @param group_names Vector of names for each group (if NULL, generic names are used)
#' @param pair_indices List of integer pairs for which groups to compare (if NULL, consecutive pairs are used)
#' @param x_metric Name of metric for x-axis (default: "inter_cellular_isoform_diversity")
#' @param y_metric Name of metric for y-axis (default: "inter_cell_type_specificity")
#' @param grid_size Size of the density estimation grid (higher = more detailed but slower)
#' @param main_title Main title for the combined plot
#' 
#' @return A ggplot object (single comparison) or patchwork combined plot (multiple comparisons)
#' 
#' @examples
#' \dontrun{
#' # This example requires multiple conditions/datasets
#' # Here we demonstrate with subsampled data to simulate conditions
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' # Create SCHT for full dataset
#' scht_full <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Calculate metrics for demonstration
#' tc_results_full <- calculate_isoform_complexity_metrics(scht_full, verbose = FALSE)
#' 
#' # For demonstration, use the same results as different "conditions"
#' tc_results_list <- list(
#'   baseline = tc_results_full,
#'   treatment = tc_results_full  # In practice, these would be different
#' )
#'                         
#' # Compare the conditions
#' plot_compare_tc_density_difference(
#'   tc_results_list, 
#'   group_names = c("Baseline", "Treatment")
#' )
#' }
#' 
#' @export
plot_compare_tc_density_difference <- function(tc_results_list,
                                               group_names = NULL,
                                               pair_indices = NULL,
                                               x_metric = "inter_cellular_isoform_diversity",
                                               y_metric = "inter_cell_type_specificity",
                                               grid_size = 100,
                                               main_title = "Density Differences Between Groups") {
  
  # Required packages check with informative error messages
  required_packages <- c("ggplot2", "patchwork", "MASS")
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if(length(missing_packages) > 0) {
    stop(paste0("The following packages are required but not installed: ", 
                paste(missing_packages, collapse = ", "), 
                ". Please install them with: install.packages(c('", 
                paste(missing_packages, collapse = "', '"), "'))"))
  }
  
  # Check inputs
  if(!is.list(tc_results_list)) {
    stop("tc_results_list must be a list of complexity metrics results")
  }
  
  # If group_names not provided, create generic names
  if(is.null(group_names)) {
    group_names <- paste0("Group_", seq_along(tc_results_list))
  }
  
  # Ensure group_names has the correct length
  if(length(group_names) != length(tc_results_list)) {
    warning("Length of group_names doesn't match length of tc_results_list. Using generic names.")
    group_names <- paste0("Group_", seq_along(tc_results_list))
  }
  
  all_metrics <- c(
    "intra_cellular_isoform_diversity",
    "inter_cellular_isoform_diversity", 
    "intra_cell_type_heterogeneity",
    "inter_cell_type_specificity",
    "intra_cell_type_heterogeneity_variability",
    "inter_cell_type_difference_variability",
    "cell_type_coexpression_variability"
  )
  
  # Create readable metric names for labels
  readable_metrics <- c(
    "Intra-cellular Isoform Diversity",
    "Inter-cellular Isoform Diversity",
    "Intra-cell-type Heterogeneity",
    "Inter-cell-type Specificity",
    "Intra-cell-type Heterogeneity Variability",
    "Inter-cell-type Difference Variability",
    "Cell-type-specific Co-Expression Variability"
  )
  
  x_label_idx <- match(x_metric, all_metrics)
  x_label <- readable_metrics[x_label_idx]
  
  y_label_idx <- match(y_metric, all_metrics)
  y_label <- readable_metrics[y_label_idx]
  
  # Define custom colour palette
  colour_palette <- colorRampPalette(c("#780522", "#A81428", "#C6403D", "#F5AC8B", "#FAD4BF", "#FBE3D6","#F8EAE1",
                                      "#EDF2F6", "#C1DDE9", "#84BDDA", "#74B0D2", "#3685BB", "#256CAE","#134B87"))
  
  # Process each group's data and calculate densities
  densities <- list()
  for (i in seq_along(tc_results_list)) {
    # Check if .prepare_tc_data exists
    if(exists(".prepare_tc_data", mode = "function")) {
      data_prep <- .prepare_tc_data(tc_results_list[[i]], x_metric, y_metric)
      data_i <- data_prep$vis_data
    } else {
      # Fallback if .prepare_tc_data doesn't exist
      if(x_metric %in% names(tc_results_list[[i]]) && y_metric %in% names(tc_results_list[[i]])) {
        data_i <- data.frame(
          x_metric = tc_results_list[[i]][[x_metric]],
          y_metric = tc_results_list[[i]][[y_metric]]
        )
        names(data_i) <- c(x_metric, y_metric)
      } else {
        # Check if data is in a nested structure
        if("tc_data" %in% names(tc_results_list[[i]])) {
          if(x_metric %in% names(tc_results_list[[i]]$tc_data) && y_metric %in% names(tc_results_list[[i]]$tc_data)) {
            data_i <- data.frame(
              x_metric = tc_results_list[[i]]$tc_data[[x_metric]],
              y_metric = tc_results_list[[i]]$tc_data[[y_metric]]
            )
            names(data_i) <- c(x_metric, y_metric)
          } else {
            stop(paste("Metrics", x_metric, "and", y_metric, "not found in tc_results_list[[", i, "]]$tc_data"))
          }
        } else {
          stop(paste("Metrics", x_metric, "and", y_metric, "not found in tc_results_list[[", i, "]]"))
        }
      }
    }
    
    # Calculate density
    density_i <- MASS::kde2d(
      data_i[[x_metric]], 
      data_i[[y_metric]], 
      n = grid_size,
      lims = c(0, 1, 0, 1)
    )
    
    # Store density
    densities[[i]] <- density_i
  }
  
  # Define which pairs to compare
  if (is.null(pair_indices)) {
    # If not specified, compare all consecutive pairs
    pairs <- list()
    for (i in 1:(length(group_names)-1)) {
      pairs[[i]] <- c(i, i+1)
    }
  } else {
    # Use specified pairs
    pairs <- pair_indices
  }
  
  # Calculate all density differences to find global maximum for unified scale
  diff_data_list <- list()
  
  for (i in seq_along(pairs)) {
    pair <- pairs[[i]]
    idx1 <- pair[1]
    idx2 <- pair[2]
    
    # Ensure indices are valid
    if(idx1 > length(densities) || idx2 > length(densities)) {
      warning(paste("Invalid pair indices:", idx1, idx2, "- skipping"))
      next
    }
    
    # Calculate density difference
    diff_z <- densities[[idx2]]$z - densities[[idx1]]$z
    diff_data_list[[i]] <- diff_z
  }
  
  # Find global maximum absolute value for unified scale
  global_max_abs <- max(sapply(diff_data_list, function(diff) {
    if(!is.null(diff)) max(abs(diff), na.rm = TRUE) else 0
  }))
  
  # Create plots for each pair
  diff_plots <- list()
  
  for (i in seq_along(pairs)) {
    pair <- pairs[[i]]
    idx1 <- pair[1]
    idx2 <- pair[2]
    
    if(is.null(diff_data_list[[i]])) next
    
    # Convert to data frame for plotting
    diff_df <- expand.grid(
      x = densities[[idx1]]$x,
      y = densities[[idx1]]$y
    )
    diff_df$diff <- as.vector(diff_data_list[[i]])
    
    # Create plot
    plot_title <- paste(group_names[idx1], "\u2192", group_names[idx2])
    
    p <- ggplot2::ggplot(diff_df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = diff)) +
      ggplot2::scale_fill_gradientn(
        colours = colour_palette(100),
        limits = c(-global_max_abs, global_max_abs),  
        name = "Density\nDifference"
      ) +
      ggplot2::stat_contour(ggplot2::aes(z = diff), colour = "grey20", alpha = 0.7, breaks = 0, size = 0.5) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::labs(
        title = plot_title,
        x = x_label,
        y = y_label
      ) +
      ggplot2::theme_minimal(base_size = 16) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.title = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(size = 16, face = "bold"), 
        axis.text = ggplot2::element_text(size = 14)
      )

    diff_plots[[i]] <- p
  }
  
  # For a single plot, return it directly
  if (length(diff_plots) == 1) {
    return(diff_plots[[1]])
  }
  
  # Arrange multiple plots in a grid
  n_plots <- length(diff_plots)
  if(n_plots <= 2) {
    n_cols <- n_plots
  } else {
    n_cols <- ceiling(sqrt(n_plots))
  }
  
  # For multiple plots, arrange in a grid
  combined_plot <- patchwork::wrap_plots(diff_plots, ncol = n_cols) +
    patchwork::plot_annotation(
      title = main_title,
      subtitle = paste("Comparing", x_label, "and", y_label),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 14)
      )
    )
  
  return(combined_plot)
}

#' Create heatmaps for transcriptome complexity metrics across groups using ComplexHeatmap
#' 
#' This function creates heatmaps for visualising multiple complexity metrics
#' across different groups or conditions. It can also show changes in metrics
#' between consecutive groups. The function supports different methods for
#' selecting genes to display.
#' 
#' @param tc_results_list List of transcriptomic complexity results objects
#' @param groups Vector of names for each group (if NULL, generic names are used)
#' @param metrics Vector of metric names to include in heatmaps
#' @param n_top_genes Number of top genes to include in heatmaps
#' @param selection_method Method for selecting genes ("variance", "magnitude", or "custom")
#' @param custom_genes Vector of custom gene names to use when selection_method is "custom"
#' @param cluster_genes Whether to cluster genes in heatmaps
#' @param show_changes Whether to create separate heatmaps showing changes between groups
#' 
#' @return A list containing:
#'   \item{heatmaps}{List of ComplexHeatmap objects for each metric}
#'   \item{change_heatmaps}{List of ComplexHeatmap objects showing changes between groups}
#'   \item{metric_matrices}{List of matrices with metric values}
#'   \item{change_matrices}{List of matrices with change values}
#'   \item{top_genes}{Vector of selected gene names}
#'   \item{selection_method}{Method used for gene selection}
#' 
#' @examples
#' \dontrun{
#' # This example requires multiple conditions/datasets  
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' # Create SCHT
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#' 
#' # Calculate metrics
#' tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)
#' 
#' # For demonstration, create a list with the same results
#' tc_results_list <- list(
#'   baseline = tc_results,
#'   treatment = tc_results  # In practice, these would be different
#' )
#'                         
#' # Create heatmaps with default settings
#' heatmap_results <- plot_compare_tc_complexity_heatmap(
#'   tc_results_list, 
#'   groups = c("Baseline", "Treatment")
#' )
#'                   
#' # Display the first heatmap if ComplexHeatmap is available
#' if(requireNamespace("ComplexHeatmap", quietly = TRUE)) {
#'   print(heatmap_results$heatmaps$intra_cellular_isoform_diversity)
#' }
#' }
#' 
#' @export
plot_compare_tc_complexity_heatmap <- function(tc_results_list, 
                                               groups = NULL,
                                               metrics = c(
                                                 "intra_cellular_isoform_diversity",
                                                 "inter_cellular_isoform_diversity", 
                                                 "intra_cell_type_heterogeneity",
                                                 "inter_cell_type_specificity",
                                                 "intra_cell_type_heterogeneity_variability",
                                                 "inter_cell_type_difference_variability",
                                                 "cell_type_coexpression_variability"
                                               ),
                                               n_top_genes = 50,
                                               selection_method = "variance",
                                               custom_genes = NULL, 
                                               cluster_genes = FALSE,
                                               show_changes = TRUE) {
  
  # Required packages check with informative error messages
  required_packages <- c("ComplexHeatmap", "circlize", "tidyr", "dplyr", "grid")
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if(length(missing_packages) > 0) {
    stop(paste0("The following packages are required but not installed: ", 
                paste(missing_packages, collapse = ", "), 
                ". Please install the Bioconductor packages with BiocManager::install() and CRAN packages with install.packages()"))
  }
  
  # If groups not provided, use generic names
  if (is.null(groups)) {
    groups <- paste0("Group_", seq_along(tc_results_list))
  }
  
  # Ensure tc_results_list and groups have the same length
  if (length(tc_results_list) != length(groups)) {
    stop("Length of tc_results_list must match length of groups")
  }
  
  # Find genes present in all groups (this part stays the same)
  genes_by_group <- lapply(tc_results_list, function(tc) {
    if ("metrics" %in% names(tc) && "gene" %in% names(tc$metrics)) {
      return(tc$metrics$gene)
    } else if ("gene" %in% names(tc)) {
      return(tc$gene)
    } else {
      # Try to find genes in other common structures
      for (field in names(tc)) {
        if (is.data.frame(tc[[field]]) && "gene" %in% names(tc[[field]])) {
          return(tc[[field]]$gene)
        }
      }
      stop("Could not find gene information in tc_results structure")
    }
  })
  
  common_genes <- Reduce(intersect, genes_by_group)
  
  if (length(common_genes) == 0) {
    stop("No genes common to all groups found")
  }
  
  # Create list to store matrices for each metric (this stays the same)
  metric_matrices <- list()
  
  # This extraction function stays the same
  .extract_metric_value <- function(tc_results, gene, metric) {
    # Try different possible structures
    if ("metrics" %in% names(tc_results) && "gene" %in% names(tc_results$metrics)) {
      idx <- which(tc_results$metrics$gene == gene)
      if (length(idx) > 0 && metric %in% names(tc_results$metrics)) {
        return(tc_results$metrics[idx, metric])
      }
    } else if (all(c("gene", metric) %in% names(tc_results))) {
      idx <- which(tc_results$gene == gene)
      if (length(idx) > 0) {
        return(tc_results[idx, metric])
      }
    } else {
      # Try to find in nested structures
      for (field in names(tc_results)) {
        if (is.data.frame(tc_results[[field]]) && 
            all(c("gene", metric) %in% names(tc_results[[field]]))) {
          idx <- which(tc_results[[field]]$gene == gene)
          if (length(idx) > 0) {
            return(tc_results[[field]][idx, metric])
          }
        }
      }
    }
    return(NA)
  }
  
  # Fill metric matrices (unchanged)
  for (metric in metrics) {
    # Create matrix for this metric - rows are genes, columns are groups
    metric_matrix <- matrix(NA, nrow = length(common_genes), ncol = length(groups))
    rownames(metric_matrix) <- common_genes
    colnames(metric_matrix) <- groups
    
    # Fill in values
    for (i in seq_along(tc_results_list)) {
      tc_results <- tc_results_list[[i]]
      group <- groups[i]
      
      for (gene in common_genes) {
        value <- .extract_metric_value(tc_results, gene, metric)
        if (!is.null(value) && !is.na(value)) {
          metric_matrix[gene, group] <- value
        }
      }
    }
    
    # Store matrix
    metric_matrices[[metric]] <- metric_matrix
  }
  
  # For change matrices, calculate changes between consecutive groups (unchanged)
  change_matrices <- list()
  if (show_changes && length(groups) > 1) {
    for (metric in metrics) {
      metric_matrix <- metric_matrices[[metric]]
      change_matrix <- matrix(NA, nrow = length(common_genes), ncol = length(groups) - 1)
      rownames(change_matrix) <- common_genes
      
      # Create column names for changes
      change_colnames <- character(length(groups) - 1)
      for (i in 1:(length(groups) - 1)) {
        change_colnames[i] <- paste0(groups[i], "\u2192", groups[i+1])
      }
      colnames(change_matrix) <- change_colnames
      
      # Calculate changes - keeping NAs as they are
      for (i in 1:(length(groups) - 1)) {
        # Only calculate difference if both values are non-NA
        for (j in 1:nrow(metric_matrix)) {
          if (!is.na(metric_matrix[j, i]) && !is.na(metric_matrix[j, i+1])) {
            change_matrix[j, i] <- metric_matrix[j, i+1] - metric_matrix[j, i]
          }
          # else leave as NA
        }
      }
      
      # Store change matrix
      change_matrices[[metric]] <- change_matrix
    }
  }
  
  # Select top genes for visualisation (unchanged)
  if (selection_method == "custom" && !is.null(custom_genes)) {
    # Use user-provided gene list
    valid_custom_genes <- intersect(custom_genes, common_genes)
    
    if (length(valid_custom_genes) == 0) {
      stop("None of the custom genes provided are present in the dataset")
    }
    
    if (length(valid_custom_genes) < length(custom_genes)) {
      warning(paste("Only", length(valid_custom_genes), "out of", 
                    length(custom_genes), "custom genes were found in the dataset"))
    }
    
    top_genes <- valid_custom_genes
  } else if (selection_method == "variance") {
    # Method 1: Calculate variance across groups for each metric
    gene_variances <- matrix(0, nrow = length(common_genes), ncol = length(metrics))
    rownames(gene_variances) <- common_genes
    colnames(gene_variances) <- metrics
    
    for (i in seq_along(metrics)) {
      metric <- metrics[i]
      metric_matrix <- metric_matrices[[metric]]
      
      # Calculate row variances (safely omit NA values)
      gene_variances[, metric] <- apply(metric_matrix, 1, function(x) {
        if (all(is.na(x))) return(0)
        return(var(x, na.rm = TRUE))
      })
    }
    
    # Calculate overall variance score
    overall_variance <- rowSums(gene_variances, na.rm = TRUE)
    
    # Find top genes by variance
    top_genes <- names(sort(overall_variance, decreasing = TRUE)[1:min(n_top_genes, length(common_genes))])
  } else if (selection_method == "magnitude") {
    # Method 2: Calculate absolute magnitude of changes
    if (!show_changes) {
      warning("selection_method='magnitude' works best with show_changes=TRUE")
    }
    
    gene_magnitudes <- matrix(0, nrow = length(common_genes), ncol = length(metrics))
    rownames(gene_magnitudes) <- common_genes
    colnames(gene_magnitudes) <- metrics
    
    for (i in seq_along(metrics)) {
      metric <- metrics[i]
      if (show_changes && !is.null(change_matrices[[metric]])) {
        change_matrix <- change_matrices[[metric]]
        # Sum of absolute changes, omitting NAs
        gene_magnitudes[, metric] <- apply(change_matrix, 1, function(x) {
          if (all(is.na(x))) return(0)
          return(sum(abs(x), na.rm = TRUE))
        })
      } else {
        # Alternative: range of values, omitting NAs
        metric_matrix <- metric_matrices[[metric]]
        gene_magnitudes[, metric] <- apply(metric_matrix, 1, function(x) {
          if (all(is.na(x))) return(0)
          return(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
        })
      }
    }
    
    # Calculate overall magnitude score
    overall_magnitude <- rowSums(gene_magnitudes, na.rm = TRUE)
    
    # Find top genes by magnitude
    top_genes <- names(sort(overall_magnitude, decreasing = TRUE)[1:min(n_top_genes, length(common_genes))])
  } else {
    stop("Invalid selection_method. Use 'variance', 'magnitude', or 'custom'")
  }
  
  # Set up colorblind-friendly palette function for ComplexHeatmap
  .get_colour_palette <- function(is_diverging = FALSE, matrix) {
    if (is_diverging) {
      # Colorblind-friendly diverging palette (blue-white-red)
      colors <- c("#B2182BFF", "#D6604DFF", "#F4A582FF", "#FDDBC7FF", 
                  "#D1E5F0FF", "#92C5DEFF", "#4393C3FF", "#2166ACFF")
      
      # Get max absolute value for symmetric color scale
      max_abs <- max(abs(matrix), na.rm = TRUE)
      
      # Create a color mapping function with circlize
      return(circlize::colorRamp2(
        c(-max_abs, -max_abs/2, 0, max_abs/2, max_abs),
        c(colors[1], colors[3], "white", colors[6], colors[8])
      ))
    } else {
      # Colorblind-friendly sequential palette (purples)
      colors <- c("#F9DDDAFF", "#F3BEC7FF", "#E8A0BCFF", "#D785B5FF", "#BF6DB0FF", 
                  "#A159A9FF", "#7C489CFF", "#573B88FF")
      
      # Get min and max values
      min_val <- min(matrix, na.rm = TRUE)
      max_val <- max(matrix, na.rm = TRUE)
      
      # Create color mapping function
      return(circlize::colorRamp2(
        seq(min_val, max_val, length.out = length(colors)),
        colors
      ))
    }
  }
  
  # Create readable metric names
  readable_metrics <- c(
    "Intra-cellular Isoform Diversity",
    "Inter-cellular Isoform Diversity",
    "Intra-cell-type Heterogeneity",
    "Inter-cell-type Specificity",
    "Intra-cell-type Heterogeneity Variability",
    "Inter-cell-type Difference Variability",
    "Cell-type-specific Co-Expression Variability"
  )
  names(readable_metrics) <- metrics
  
  # New function to create a ComplexHeatmap
  .create_complex_heatmap <- function(matrix, title, is_change = FALSE, cluster_rows = cluster_genes) {
    # Select appropriate colour palette
    is_diverging <- is_change || any(matrix < 0, na.rm = TRUE)
    col_fun <- .get_colour_palette(is_diverging, matrix)
    
    # Handle NA values
    na_col <- "grey90"
    
    # Set up clustering parameters
    if (cluster_rows) {
      # Define a custom distance function that considers NA patterns
      cluster_rows_custom <- function(mat) {
        # Create a distance matrix considering NA patterns
        n <- nrow(mat)
        dist_mat <- matrix(0, n, n)
        na_pattern <- is.na(mat)
        na_distance_weight <- 0.5
        
        for (i in 1:(n-1)) {
          for (j in (i+1):n) {
            # Pattern distance: proportion of columns where NA pattern differs
            pattern_diff <- sum(na_pattern[i,] != na_pattern[j,]) / ncol(mat)
            
            # Value distance: only for columns where both have values
            valid_idx <- which(!na_pattern[i,] & !na_pattern[j,])
            value_dist <- 0
            if (length(valid_idx) > 0) {
              # Normalise by number of valid comparisons
              value_dist <- sqrt(sum((mat[i,valid_idx] - mat[j,valid_idx])^2)) / sqrt(length(valid_idx))
            }
            
            # Combine distances with weighting
            if (length(valid_idx) > 0) {
              combined_dist <- (1 - na_distance_weight) * value_dist + 
                na_distance_weight * pattern_diff
            } else {
              combined_dist <- pattern_diff
            }
            
            dist_mat[i,j] <- dist_mat[j,i] <- combined_dist
          }
        }
        
        # Convert to dist object and perform hierarchical clustering
        hclust_result <- stats::hclust(stats::as.dist(dist_mat), method = "complete")
        return(hclust_result)
      }
      
      # Use tryCatch for robust error handling
      ht <- tryCatch({
        # Try custom distance-based clustering
        ComplexHeatmap::Heatmap(
          matrix,
          name = if (is_diverging) "Change" else "Value",
          col = col_fun,
          
          # Clustering parameters
          cluster_rows = cluster_rows_custom(matrix),
          cluster_columns = FALSE,
          
          # Cell appearance
          rect_gp = grid::gpar(col = NA),  # No cell borders
          na_col = na_col,
          
          # Row and column names
          row_names_gp = grid::gpar(fontsize = 10),
          column_names_gp = grid::gpar(fontsize = 12),
          column_names_rot = 45,
          
          # Title
          column_title = title,
          column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
          
          # Legend parameters
          heatmap_legend_param = list(
            title = if (is_diverging) "Change" else "Value",
            title_position = "topcenter",
            title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = 10),
            legend_width = grid::unit(2, "cm")
          )
        )
      }, error = function(e) {
        # If custom clustering fails, try standard ComplexHeatmap clustering
        message("Custom clustering with NA patterns failed: ", e$message, 
                ". Trying standard method...")
        
        tryCatch({
          ComplexHeatmap::Heatmap(
            matrix,
            name = if (is_diverging) "Change" else "Value",
            col = col_fun,
            
            # Standard clustering
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            clustering_distance_rows = "euclidean",
            clustering_method_rows = "complete",
            
            # Cell appearance
            rect_gp = grid::gpar(col = NA),
            na_col = na_col,
            
            # Row and column names
            row_names_gp = grid::gpar(fontsize = 10),
            column_names_gp = grid::gpar(fontsize = 12),
            column_names_rot = 45,
            
            # Title
            column_title = title,
            column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
            
            # Legend parameters
            heatmap_legend_param = list(
              title = if (is_diverging) "Change" else "Value",
              title_position = "topcenter",
              title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
              labels_gp = grid::gpar(fontsize = 10),
              legend_width = grid::unit(2, "cm")
            )
          )
        }, error = function(e2) {
          # If all clustering methods fail, fall back to sorting rows
          message("Standard clustering also failed: ", e2$message, ". Using row sorting instead.")
          
          # Sort rows by mean value (ignoring NAs)
          row_means <- rowMeans(matrix, na.rm = TRUE)
          # Handle completely NA rows
          row_means[is.na(row_means)] <- min(row_means, na.rm = TRUE) - 1
          
          # Sort in descending order
          order_idx <- order(row_means, decreasing = TRUE)
          
          ComplexHeatmap::Heatmap(
            matrix[order_idx, , drop = FALSE],
            name = if (is_diverging) "Change" else "Value",
            col = col_fun,
            
            # No clustering
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            
            # Cell appearance
            rect_gp = grid::gpar(col = NA),
            na_col = na_col,
            
            # Row and column names
            row_names_gp = grid::gpar(fontsize = 10),
            column_names_gp = grid::gpar(fontsize = 12),
            column_names_rot = 45,
            
            # Title
            column_title = paste0(title, " (rows sorted)"),
            column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
            
            # Legend parameters
            heatmap_legend_param = list(
              title = if (is_diverging) "Change" else "Value",
              title_position = "topcenter",
              title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
              labels_gp = grid::gpar(fontsize = 10),
              legend_width = grid::unit(2, "cm")
            )
          )
        })
      })
    } else {
      # Heatmap without clustering
      ht <- ComplexHeatmap::Heatmap(
        matrix,
        name = if (is_diverging) "Change" else "Value",
        col = col_fun,
        
        # No clustering
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        
        # Cell appearance
        rect_gp = grid::gpar(col = NA),
        na_col = na_col,
        
        # Row and column names
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 12),
        column_names_rot = 45,
        
        # Title
        column_title = title,
        column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
        
        # Legend parameters
        heatmap_legend_param = list(
          title = if (is_diverging) "Change" else "Value",
          title_position = "topcenter",
          title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
          labels_gp = grid::gpar(fontsize = 10),
          legend_width = grid::unit(2, "cm")
        )
      )
    }
    
    return(ht)
  }
  
  # Create standard heatmaps using ComplexHeatmap
  heatmaps <- list()
  change_heatmaps <- list()
  
  for (metric in metrics) {
    # Get matrix for top genes only
    matrix_subset <- metric_matrices[[metric]][top_genes, , drop = FALSE]
    
    # Create heatmap
    heatmap_title <- paste0(readable_metrics[metric], " Across Groups")
    heatmaps[[metric]] <- .create_complex_heatmap(matrix_subset, 
                                                  title = heatmap_title, 
                                                  is_change = FALSE)
  }
  
  # Create change heatmaps if requested
  if (show_changes && length(groups) > 1) {
    for (metric in metrics) {
      # Get change matrix for top genes only
      change_matrix_subset <- change_matrices[[metric]][top_genes, , drop = FALSE]
      
      # Create heatmap
      change_heatmap_title <- paste0("Changes in ", readable_metrics[metric], " Between Groups")
      change_heatmaps[[metric]] <- .create_complex_heatmap(change_matrix_subset, 
                                                           title = change_heatmap_title,
                                                           is_change = TRUE)
    }
  }
  
  # Return results
  return(list(
    heatmaps = heatmaps,
    change_heatmaps = change_heatmaps,
    metric_matrices = metric_matrices,
    change_matrices = if(show_changes) change_matrices else NULL,
    top_genes = top_genes,
    selection_method = selection_method
  ))
}
