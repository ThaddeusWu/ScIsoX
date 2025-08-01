##################################################################
#  Complexity Metrics Calculation Functions for Single-cell      #
#  Hierarchical Tensor (SCHT) structure                          #
#                                                                #
#  Author: [Siyuan Wu & Ulf Schmitz]                             #
#  Institution: [James Cook University]                          #
#  Date: Jul 29, 2025                                            #
#  Package: ScIsoX V1.1.0                                        #
##################################################################

########################
# Required Libraries   #
########################
#' @importFrom stats na.omit var sd quantile density kmeans dnorm IQR approx cor
#' @importFrom utils txtProgressBar setTxtProgressBar head
#' @importFrom dplyr case_when
#' @importFrom mclust Mclust densityMclust mclustBIC
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line geom_vline
#' @importFrom ggplot2 geom_density geom_point theme_minimal labs annotate
#' @importFrom progress progress_bar
#' @importFrom knitr kable
#' @importFrom car powerTransform
#' @importFrom diptest dip.test
#' @importFrom moments skewness
NULL
# Define global variables
utils::globalVariables(c("Value", "x", "y", "Component"))

#' Calculate intra-cellular isoform diversity index (IDI)
#'
#' This function calculates the diversity of isoform usage within individual cells,
#' weighting by the expression level of each cell. This approach captures
#' cell-level diversity (intra-cellular) and measures the average tendency
#' for cells to co-express multiple isoforms of a gene.
#'
#' @param iso_mat Matrix of isoform expression, rows = isoforms, columns = cells
#' @return Named list containing intra_cellular_isoform_diversity and individual cell IDI values
#' @keywords internal
.calculate_intra_cellular_isoform_diversity <- function(iso_mat) {
  # Get cell sums
  cell_sums <- colSums(iso_mat)
  
  # Pre-allocate result vector
  cell_entropies <- numeric(ncol(iso_mat))
  
  # Process only non-zero cells
  nonzero_cells <- which(cell_sums > 0)
  
  # Calculate entropy for each non-zero cell
  for(j in nonzero_cells) {
    # Get non-zero expression values for this cell
    cell_expr <- iso_mat[, j]
    nonzero_indices <- which(cell_expr > 0)
    
    # If only one isoform is expressed, entropy is 0
    if(length(nonzero_indices) <= 1) {
      cell_entropies[j] <- 0
      next
    }
    
    # Calculate proportions for non-zero indices
    props <- cell_expr[nonzero_indices] / cell_sums[j]
    
    # Calculate Shannon entropy on non-zero values only
    shannon <- -sum(props * log2(props), na.rm = TRUE)
    
    # Normalise by maximum possible entropy
    expressed_isoforms <- length(nonzero_indices)
    cell_entropies[j] <- shannon / log2(expressed_isoforms)
  }
  
  # Calculate weighted average entropy
  if(sum(cell_sums) == 0) {
    intra_cellular_isoform_diversity <- 0
  } else {
    intra_cellular_isoform_diversity <- sum(cell_entropies * cell_sums) / sum(cell_sums)
  }
  
  # Return both the overall metric and individual cell values
  return(list(
    intra_cellular_isoform_diversity = intra_cellular_isoform_diversity,
    cell_entropies = cell_entropies
  ))
}

#' Calculate inter-cellular isoform diversity index (IDI)
#'
#' This function calculates the diversity of isoform usage across the entire
#' cell population, using average expression levels. This captures
#' population-level diversity (inter-cellular) and measures the overall
#' diversity of isoforms used across the entire cell population.
#'
#' @param iso_mat Matrix of isoform expression, rows = isoforms, columns = cells
#' @return Inter-cellular isoform diversity index
#' @keywords internal
.calculate_inter_cellular_isoform_diversity <- function(iso_mat) {
  # Calculate average expression for each isoform
  iso_means <- rowMeans(iso_mat, na.rm = TRUE)
  
  # Skip if no expression
  if(sum(iso_means) == 0) return(0)
  
  # Calculate total expression
  total_expr <- sum(iso_means)
  
  # Calculate relative proportions for ALL isoforms
  iso_props <- iso_means / total_expr
  
  # Calculate Shannon index, skipping log(0) issues with na.rm=TRUE
  shannon <- -sum(iso_props * log2(iso_props), na.rm = TRUE)
  
  # Normalise by maximum possible entropy using TOTAL isoform count
  n_isoforms <- length(iso_props)
  inter_cellular_isoform_diversity <- ifelse(n_isoforms > 1, shannon / log2(n_isoforms), 0)
  
  return(inter_cellular_isoform_diversity)
}

#' Calculate Jensen-Shannon distance between two probability distributions
#'
#' Calculates Jensen-Shannon distance between two probability vectors,
#' taking advantage of efficiency by operating only on non-zero elements.
#'
#' @param p First probability vector
#' @param q Second probability vector
#' @return Jensen-Shannon distance
#' @keywords internal
.calculate_js_distance <- function(p, q) {
  # Get indices where either p or q is non-zero
  nonzero_indices <- unique(c(which(p > 0), which(q > 0)))
  
  # If either distribution is empty, return 0
  if(length(nonzero_indices) == 0) return(0)
  
  # Subset to non-zero indices only
  p_sub <- p[nonzero_indices]
  q_sub <- q[nonzero_indices]
  
  # Calculate mixture distribution
  m <- (p_sub + q_sub) / 2
  
  # Calculate KL divergences, carefully avoiding log(0)
  p_pos <- p_sub > 0
  q_pos <- q_sub > 0
  
  # Calculate KL(p||m) for non-zero p
  kl_pm <- 0
  if(any(p_pos)) {
    p_indices <- which(p_pos)
    kl_pm <- sum(p_sub[p_indices] * log2(p_sub[p_indices] / m[p_indices]))
  }
  
  # Calculate KL(q||m) for non-zero q
  kl_qm <- 0
  if(any(q_pos)) {
    q_indices <- which(q_pos)
    kl_qm <- sum(q_sub[q_indices] * log2(q_sub[q_indices] / m[q_indices]))
  }
  
  # Calculate JS distance (square root of JS divergence)
  js_dist <- sqrt(0.5 * kl_pm + 0.5 * kl_qm)
  
  return(js_dist)
}

#' Calculate intra-cell type cellular heterogeneity
#'
#' This function measures the cell-to-cell variation in isoform usage patterns
#' within each cell type. High values indicate substantial cellular
#' heterogeneity within a cell type, suggesting subpopulations with distinct
#' isoform preferences.
#'
#' @param cell_type_mat Matrix of isoform expression in a specific cell type
#' @return Average Jensen-Shannon distance between cells within the cell type
#' @keywords internal
.calculate_intra_cell_type_heterogeneity <- function(cell_type_mat) {
  # Need at least 2 cells
  if(ncol(cell_type_mat) < 2) return(0)
  
  # Get cell sums for normalisation
  cell_sums <- colSums(cell_type_mat)
  
  # Skip cells with no expression
  valid_cells <- which(cell_sums > 0)
  if(length(valid_cells) < 2) return(0)  # Need at least 2 valid cells
  
  # Calculate proportions for each cell
  n_valid <- length(valid_cells)
  prop_list <- vector("list", n_valid)
  
  # Pre-calculate all cell proportion vectors
  for(i in 1:n_valid) {
    j <- valid_cells[i]
    cell_expr <- cell_type_mat[, j]
    nonzero_indices <- which(cell_expr > 0)
    props <- numeric(nrow(cell_type_mat))
    props[nonzero_indices] <- cell_expr[nonzero_indices] / cell_sums[j]
    prop_list[[i]] <- props
  }
  
  # Calculate all pairwise distances
  distances <- numeric(n_valid * (n_valid - 1) / 2)
  counter <- 0
  
  for(i in 1:(n_valid-1)) {
    for(j in (i+1):n_valid) {
      counter <- counter + 1
      distances[counter] <- .calculate_js_distance(prop_list[[i]], prop_list[[j]])
    }
  }
  
  # Return average JS distance
  return(mean(distances))
}

#' Calculate inter-cell type specificity (cell type specificity)
#'
#' This function calculates how differently a gene uses its isoforms
#' across different cell types using Jensen-Shannon distance.
#' It measures the cell type-specificity of isoform usage patterns.
#'
#' @param scht_obj Single-Cell Hierarchical Tensor object
#' @param gene Gene name to analyse
#' @return List containing cell type specificity and related information
#' @keywords internal
.calculate_inter_cell_type_specificity <- function(scht_obj, gene) {
  # Get all cell types
  cell_types <- names(scht_obj$cell_type_matrices)
  if(length(cell_types) < 2) {
    return(list(
      inter_cell_type_specificity = NA,
      single_cell_type = NA,
      cell_type_name = NA,
      expressed_cell_types = NA,
      cell_type_props = list(),
      pairwise_distances = c()
    ))
  }
  
  # Collect isoform usage proportions for each cell type
  cell_type_props <- list()
  expressed_cell_types <- c()
  
  for(cell_type in cell_types) {
    if(gene %in% names(scht_obj$cell_type_matrices[[cell_type]])) {
      iso_mat <- scht_obj$cell_type_matrices[[cell_type]][[gene]]
      
      # Skip cell types with insufficient data
      if(nrow(iso_mat) < 2 || ncol(iso_mat) == 0) next
      
      # Calculate average expression for each isoform in this cell type
      iso_means <- rowMeans(iso_mat)
      total_expr <- sum(iso_means)
      
      # Skip if no expression
      if(total_expr == 0) next
      
      # Get only non-zero means for efficiency
      nonzero_indices <- which(iso_means > 0)
      if(length(nonzero_indices) == 0) next
      
      # Calculate proportions for non-zero indices
      cell_type_props_vec <- numeric(length(iso_means))
      cell_type_props_vec[nonzero_indices] <- iso_means[nonzero_indices] / total_expr
      cell_type_props[[cell_type]] <- cell_type_props_vec
      expressed_cell_types <- c(expressed_cell_types, cell_type)
    }
  }
  
  # Check if only expressed in a single cell type
  if(length(cell_type_props) == 1) {
    return(list(
      inter_cell_type_specificity = 1.0,  # Maximum specificity for single-cell_type genes
      single_cell_type = TRUE,
      cell_type_name = names(cell_type_props)[1],
      expressed_cell_types = expressed_cell_types,
      cell_type_props = cell_type_props,
      pairwise_distances = c()
    ))
  }
  
  # Need at least two cell types to calculate specificity
  if(length(cell_type_props) < 2) {
    return(list(
      inter_cell_type_specificity = NA,
      single_cell_type = NA,
      cell_type_name = NA,
      expressed_cell_types = expressed_cell_types,
      cell_type_props = cell_type_props,
      pairwise_distances = c()
    ))
  }
  
  # Calculate pairwise JS distances
  cell_types_list <- names(cell_type_props)
  n_cell_types <- length(cell_types_list)
  distances <- numeric(n_cell_types * (n_cell_types - 1) / 2)
  counter <- 0
  
  for(i in 1:(n_cell_types-1)) {
    for(j in (i+1):n_cell_types) {
      counter <- counter + 1
      distances[counter] <- .calculate_js_distance(
        cell_type_props[[cell_types_list[i]]],
        cell_type_props[[cell_types_list[j]]]
      )
    }
  }
  
  # Return all information
  return(list(
    inter_cell_type_specificity = mean(distances, na.rm = TRUE),
    single_cell_type = FALSE,
    cell_type_name = NA,
    expressed_cell_types = expressed_cell_types,
    cell_type_props = cell_type_props,
    pairwise_distances = distances
  ))
}

#' Calculate CV of intra-cell type heterogeneity
#'
#' This function calculates the coefficient of variation of intra-cell_type
#' heterogeneity values across different cell types. High values indicate that
#' certain cell types exhibit dramatically higher cellular heterogeneity
#' than others, suggesting functionally distinct subpopulations within specific cell types.
#'
#' @param scht_obj Single-Cell Hierarchical Tensor object
#' @param gene Gene name to analyse
#' @return List containing CV and cell type-specific heterogeneity values
#' @keywords internal
.calculate_intra_cell_type_heterogeneity_variability <- function(scht_obj, gene) {
  # Get all cell types
  cell_types <- names(scht_obj$cell_type_matrices)
  cell_type_heterogeneity <- numeric(length(cell_types))
  names(cell_type_heterogeneity) <- cell_types
  valid_cell_types <- 0
  
  # Calculate heterogeneity for each cell type
  for(cell_type in cell_types) {
    if(gene %in% names(scht_obj$cell_type_matrices[[cell_type]])) {
      cell_type_mat <- scht_obj$cell_type_matrices[[cell_type]][[gene]]
      
      # Skip cell types with insufficient data
      if(nrow(cell_type_mat) < 2 || ncol(cell_type_mat) < 3) {
        cell_type_heterogeneity[cell_type] <- NA
        next
      }
      
      # Calculate intra-cell_type heterogeneity
      cell_type_heterogeneity[cell_type] <- .calculate_intra_cell_type_heterogeneity(cell_type_mat)
      valid_cell_types <- valid_cell_types + 1
    } else {
      cell_type_heterogeneity[cell_type] <- NA
    }
  }
  
  # Calculate coefficient of variation
  if(valid_cell_types >= 2) {
    mean_heterogeneity <- mean(cell_type_heterogeneity, na.rm = TRUE)
    sd_heterogeneity <- sd(cell_type_heterogeneity, na.rm = TRUE)
    cv <- ifelse(mean_heterogeneity > 0, sd_heterogeneity / mean_heterogeneity, 0)
  } else {
    cv <- NA
  }
  
  return(list(
    heterogeneity_variability = cv,
    cell_type_heterogeneity = cell_type_heterogeneity,
    valid_cell_types = valid_cell_types
  ))
}

#' Calculate CV of inter-cell type difference
#'
#' This function calculates the coefficient of variation of pairwise
#' Jensen-Shannon distances between cell types. High values indicate that
#' certain cell type pairs exhibit much more dramatic differences in isoform usage
#' than others, suggesting specialised functional distinctions between specific cell populations.
#'
#' @param specificity_result Result from .calculate_inter_cell_type_specificity()
#' @return Coefficient of variation of cell type difference distances
#' @keywords internal
.calculate_inter_cell_type_difference_variability <- function(specificity_result) {
  # Extract pairwise distances
  distances <- specificity_result$pairwise_distances
  
  # Need at least 2 pairwise distances
  if(length(distances) < 2) return(NA)
  
  # Calculate coefficient of variation
  mean_dist <- mean(distances, na.rm = TRUE)
  sd_dist <- sd(distances, na.rm = TRUE)
  cv <- ifelse(mean_dist > 0, sd_dist / mean_dist, 0)
  
  return(cv)
}

#' Calculate CV of cell-type-specific co-expression varaibility
#'
#' This function calculates the coefficient of variation of mean intra-cellular
#' diversity across cell types. High values indicate that a gene employs dramatically
#' different isoform co-expression patterns in different cell types.
#'
#' @param scht_obj Single-Cell Hierarchical Tensor object
#' @param gene Gene name to analyse
#' @return List containing CV and cell type-specific mean IDI values
#' @keywords internal
.calculate_cell_type_coexpression_variability <- function(scht_obj, gene) {
  # Get all cell types
  cell_types <- names(scht_obj$cell_type_matrices)
  cell_type_mean_idi <- numeric(length(cell_types))
  names(cell_type_mean_idi) <- cell_types
  valid_cell_types <- 0
  
  # Calculate mean IDI for each cell type
  for(cell_type in cell_types) {
    if(gene %in% names(scht_obj$cell_type_matrices[[cell_type]])) {
      cell_type_mat <- scht_obj$cell_type_matrices[[cell_type]][[gene]]
      
      # Skip cell types with insufficient data
      if(nrow(cell_type_mat) < 2 || ncol(cell_type_mat) < 2) {
        cell_type_mean_idi[cell_type] <- NA
        next
      }
      
      # Calculate cell IDI values
      cell_results <- .calculate_intra_cellular_isoform_diversity(cell_type_mat)
      
      # Calculate mean for this cell type
      cell_type_mean_idi[cell_type] <- mean(cell_results$cell_entropies, na.rm = TRUE)
      valid_cell_types <- valid_cell_types + 1
    } else {
      cell_type_mean_idi[cell_type] <- NA
    }
  }
  
  # Calculate coefficient of variation
  if(valid_cell_types >= 2) {
    mean_idi <- mean(cell_type_mean_idi, na.rm = TRUE)
    sd_idi <- sd(cell_type_mean_idi, na.rm = TRUE)
    cv <- ifelse(mean_idi > 0, sd_idi / mean_idi, 0)
  } else {
    cv <- NA
  }
  
  return(list(
    coexpression_variability = cv,
    cell_type_mean_idi = cell_type_mean_idi,
    valid_cell_types = valid_cell_types
  ))
}

#' Enhanced data preprocessing for threshold detection
#'
#' This function provides comprehensive preprocessing including outlier detection,
#' transformation, and distribution classification, with standardised handling for
#' zero-inflation, tail inflation, extreme skewness, and multimodality.
#'
#' @param values Numeric vector of values
#' @param handle_outliers Whether to temporarily exclude extreme outliers
#' @param feature_name Name of the feature (for messages)
#' @param verbose Whether to print detailed messages
#' @return List containing processed data and preprocessing information
#' @keywords internal
.preprocess_data <- function(values, handle_outliers = TRUE, feature_name = "feature", verbose = FALSE) {
  # Remove missing values
  values <- values[is.finite(values)]
  values <- values[!is.na(values)]
  
  if (length(values) < 10) {
    return(list(
      processed_values = values,
      original_values = values,
      method = "none",
      explanation = "Insufficient data for preprocessing",
      distribution_type = "insufficient_data",
      skewness = NA,
      is_zero_inflated = FALSE,
      is_right_tail_inflated = FALSE,
      is_multimodal = FALSE,
      is_extremely_skewed = FALSE
    ))
  }
  
  # Calculate basic statistics
  stats <- list(
    min = min(values),
    max = max(values),
    mean = mean(values),
    median = median(values),
    q25 = quantile(values, 0.25),
    q75 = quantile(values, 0.75),
    sd = sd(values)
  )
  stats$iqr <- stats$q75 - stats$q25
  stats$range <- stats$max - stats$min
  
  # Calculate skewness 
  skewness <- moments::skewness(values)
  
  # Flag for extreme skewness
  is_extremely_skewed <- abs(skewness) >= 1.2
  is_skewed <- abs(skewness) >= 0.5
  skew_severity <- abs(skewness)
  
  # Detect zero-inflation/tail-inflation using histogram-based approach
  # Determine optimal bin width using Freedman-Diaconis rule
  bw <- 2 * stats$iqr / length(values)^(1/3)
  if (!is.finite(bw) || bw <= 0) {
    # Fallback if FD rule gives bad result
    bw <- (stats$max - stats$min) / min(30, sqrt(length(values)))
  }
  
  # Calculate optimal number of bins
  data_range <- stats$max - stats$min
  n_bins <- max(10, min(50, ceiling(data_range / bw)))
  
  # Create breaks that ensure the first bin starts exactly at zero
  breaks <- c(stats$min, seq(bw, stats$max * 1.001, length.out = n_bins))
  
  # Generate histogram
  hist_result <- hist(values, breaks = breaks, plot = FALSE)
  
  # Determine tails boundaries
  near_zero_threshold <- hist_result$breaks[2]
  near_right_tail_threshold <- hist_result$breaks[length(hist_result$breaks)-1]
  
  # Analyse first and second bins
  first_bin_count <- hist_result$counts[1]
  first_bin_proportion <- first_bin_count/sum(hist_result$counts)
  if (length(hist_result$counts) > 1) {
    second_bin_count <- hist_result$counts[2]
    bin_ratio <- first_bin_count / max(second_bin_count, 1)  # Avoid division by zero
  } else {
    # Handle case with only one bin
    bin_ratio <- Inf
  }
  
  is_zero_inflated <- (bin_ratio > 2) & (first_bin_proportion > 0.05)
  
  # Analyse last and second last bins
  last_bin_count <- hist_result$counts[length(hist_result$counts)]
  last_bin_proportion <- last_bin_count/sum(hist_result$counts)
  if (length(hist_result$counts) > 1) {
    second_last_bin_count <- hist_result$counts[length(hist_result$counts)-1]
    right_tail_bin_ratio <- last_bin_count / max(second_last_bin_count, 1)  # Avoid division by zero
  } else {
    # Handle case with only one bin
    right_tail_bin_ratio <- Inf
  }
  
  is_right_tail_inflated <- right_tail_bin_ratio > 2 & (last_bin_proportion > 0.05)
  
  # Detect tail inflation using robust methods
  iqr_factor <- max(1.5, min(3.0, 3.0 - (abs(skewness) - 1) * 0.5))
  upper_fence <- stats$q75 + iqr_factor * stats$iqr
  lower_fence <- stats$q25 - iqr_factor * stats$iqr
  
  is_upper_outlier <- values > upper_fence
  is_lower_outlier <- values < lower_fence
  is_outlier <- is_upper_outlier | is_lower_outlier
  
  outlier_count <- sum(is_outlier)
  outlier_proportion <- outlier_count / length(values)
  
  
  # Improved multimodality detection that's more sensitive
  is_multimodal <- .detect_multimodality(values)
  
  # Determine distribution type for method selection
  # 1. Extreme skewness, 2. Zero-inflation, 3. Multimodality, 4. Skewness, 5. Unimodal
  if (is_extremely_skewed) {
    if (skewness > 0) {
      distribution_type <- "extreme_right_skewed"
    } else {
      distribution_type <- "extreme_left_skewed"
    }
  } else if (is_zero_inflated && is_right_tail_inflated) {
    distribution_type <- "both_tails_inflated"
  } else if (is_zero_inflated || is_right_tail_inflated) {
    if (is_zero_inflated) {
      distribution_type <- "zero_inflated"
    } else if (is_right_tail_inflated) {
      distribution_type <- "right_tail_inflated"
    }
  } else if (is_multimodal) {
    distribution_type <- "multimodal"
  } else if (is_skewed) {
    if (skewness > 0) {
      distribution_type <- "right_skewed"
    } else {
      distribution_type <- "left_skewed"
    }
  } else {
    distribution_type <- "unimodal"
  }
  
  # Decision logic for preprocessing approach
  if (is_extremely_skewed) {
    # For extremely skewed data, use more aggressive transformation
    lambda <- ifelse(skewness > 0, 0, 2) 
    transformed_values <- .yeo_johnson_transform(values, lambda)
    
    # Remove any remaining outliers after transformation
    t_stats <- list(
      mean = mean(transformed_values),
      median = median(transformed_values),
      q25 = quantile(transformed_values, 0.25),
      q75 = quantile(transformed_values, 0.75)
    )
    t_stats$iqr <- t_stats$q75 - t_stats$q25
    
    t_upper_fence <- t_stats$q75 + 2 * t_stats$iqr
    t_lower_fence <- t_stats$q25 - 2 * t_stats$iqr
    
    t_is_outlier <- transformed_values > t_upper_fence | transformed_values < t_lower_fence
    
    return(list(
      processed_values = transformed_values[!t_is_outlier],
      original_values = values,
      method = "extreme_skew_transform",
      explanation = sprintf("Applied aggressive transformation for extreme skewness (%.2f)", skewness),
      distribution_type = distribution_type,
      lambda = lambda,
      transform_func = function(x) .yeo_johnson_transform(x, lambda),
      inverse_func = function(x) .yeo_johnson_inverse(x, lambda),
      skewness = skewness,
      is_extremely_skewed = is_extremely_skewed,
      is_zero_inflated = is_zero_inflated,
      is_multimodal = is_multimodal,
      is_right_tail_inflated = is_right_tail_inflated,
      safe_range = list(lower = 0.4, upper = 0.7)
    ))
  } else if (is_zero_inflated && is_right_tail_inflated) {
    # Extract non-tails component
    non_tails_values <- values[(values > near_zero_threshold) & (values < near_right_tail_threshold)]
    
    if (length(non_tails_values) < 10) {
      # Too few non-tails values
      return(list(
        processed_values = values,
        original_values = values,
        method = "preserve_both_tails_inflated",
        explanation = "Preserved both tails-inflated data (insufficient non-zero values)",
        distribution_type = distribution_type,
        is_zero_inflated = TRUE,
        is_extremely_skewed = is_extremely_skewed,
        right_tail_threshold = near_right_tail_threshold,
        zero_threshold = near_zero_threshold,
        skewness = skewness,
        is_right_tail_inflated = TRUE
      ))
    }
    
    # Check if non-tails component is multimodal
    non_tails_multimodal <- .detect_multimodality(non_tails_values)
    
    if (non_tails_multimodal) {
      #  non-tails component is multimodal, treat accordingly
      return(list(
        processed_values = non_tails_values,
        original_values = values,
        method = "both_tails_inflated_multimodal",
        explanation = "Extracted non-tails values which show multimodal distribution",
        distribution_type = distribution_type, 
        is_zero_inflated = TRUE,
        is_right_tail_inflated = TRUE,
        is_multimodal = TRUE,
        is_extremely_skewed = is_extremely_skewed,
        right_tail_threshold = near_right_tail_threshold,
        zero_threshold = near_zero_threshold,
        skewness = skewness
      ))
    }
    
    # Process the non-tails component
    non_tails_skewness <- moments::skewness(non_tails_values)
    
    if (abs(non_tails_skewness) > 0.8) {
      # Skewed non-tails component, apply transformation
      lambda <- .find_optimal_lambda(non_tails_values)
      transformed_non_tails <- .yeo_johnson_transform(non_tails_values, lambda)
      
      return(list(
        processed_values = transformed_non_tails,
        original_values = values,
        non_tails_values = non_tails_values,
        method = "both_tails_inflated_with_transform",
        explanation = sprintf("Extracted non-tails values (%.1f%%) and applied transformation",
                              100 * length(non_tails_values) / length(values)),
        distribution_type = distribution_type,
        is_zero_inflated = TRUE,
        is_right_tail_inflated = TRUE,
        is_extremely_skewed = is_extremely_skewed,
        right_tail_threshold = near_right_tail_threshold,
        zero_threshold = near_zero_threshold,
        lambda = lambda,
        transform_func = function(x) .yeo_johnson_transform(x, lambda),
        inverse_func = function(x) .yeo_johnson_inverse(x, lambda),
        skewness = skewness,
        non_tails_skewness = non_tails_skewness
      ))
    } else {
      # non-tails component not highly skewed
      return(list(
        processed_values = non_tails_values,
        original_values = values,
        method = "both_tails_inflated_extract",
        explanation = sprintf("Extracted non-tails values (%.1f%%) for threshold determination",
                              100 * length(non_tails_values) / length(values)),
        distribution_type = distribution_type,
        is_zero_inflated = TRUE,
        is_right_tail_inflated = TRUE,
        is_extremely_skewed = is_extremely_skewed,
        right_tail_threshold = near_right_tail_threshold,
        zero_threshold = near_zero_threshold,
        skewness = skewness,
        non_tails_skewness = non_tails_skewness
      ))
    }
  } else if (is_zero_inflated || is_right_tail_inflated) {
    if (is_zero_inflated) {
      # Extract non-zero component
      non_zero_values <- values[values > near_zero_threshold]
      
      if (length(non_zero_values) < 10) {
        # Too few non-zero values
        return(list(
          processed_values = values,
          original_values = values,
          method = "preserve_zero_inflated",
          explanation = "Preserved zero-inflated data (insufficient non-zero values)",
          distribution_type = distribution_type,
          is_zero_inflated = TRUE,
          is_extremely_skewed = is_extremely_skewed,
          zero_threshold = near_zero_threshold,
          skewness = skewness,
          is_right_tail_inflated = FALSE
        ))
      }
      
      # Check if non-zero component is multimodal
      non_zero_multimodal <- .detect_multimodality(non_zero_values)
      
      if (non_zero_multimodal) {
        # Non-zero component is multimodal, treat accordingly
        return(list(
          processed_values = non_zero_values,
          original_values = values,
          method = "zero_inflated_multimodal",
          explanation = "Extracted non-zero values which show multimodal distribution",
          distribution_type = distribution_type, 
          is_zero_inflated = TRUE,
          is_multimodal = TRUE,
          is_right_tail_inflated = FALSE,
          is_extremely_skewed = is_extremely_skewed,
          zero_threshold = near_zero_threshold,
          skewness = skewness
        ))
      }
      
      # Process the non-zero component
      non_zero_skewness <- moments::skewness(non_zero_values)
      
      if (abs(non_zero_skewness) > 0.8) {
        # Skewed non-zero component, apply transformation
        lambda <- .find_optimal_lambda(non_zero_values)
        transformed_non_zero <- .yeo_johnson_transform(non_zero_values, lambda)
        
        return(list(
          processed_values = transformed_non_zero,
          original_values = values,
          non_zero_values = non_zero_values,
          method = "zero_inflated_with_transform",
          explanation = sprintf("Extracted non-zero values (%.1f%%) and applied transformation",
                                100 * length(non_zero_values) / length(values)),
          distribution_type = distribution_type,
          is_zero_inflated = TRUE,
          is_right_tail_inflated = FALSE,
          is_extremely_skewed = is_extremely_skewed,
          zero_threshold = near_zero_threshold,
          lambda = lambda,
          transform_func = function(x) .yeo_johnson_transform(x, lambda),
          inverse_func = function(x) .yeo_johnson_inverse(x, lambda),
          skewness = skewness,
          non_zero_skewness = non_zero_skewness
        ))
      } else {
        # Non-zero component not highly skewed
        return(list(
          processed_values = non_zero_values,
          original_values = values,
          method = "zero_inflated_extract",
          explanation = sprintf("Extracted non-zero values (%.1f%%) for threshold determination",
                                100 * length(non_zero_values) / length(values)),
          distribution_type = distribution_type,
          is_zero_inflated = TRUE,
          is_right_tail_inflated = FALSE,
          is_extremely_skewed = is_extremely_skewed,
          zero_threshold = near_zero_threshold,
          skewness = skewness,
          non_zero_skewness = non_zero_skewness
        ))
      }
    } else if (is_right_tail_inflated) {
      # Extract non-right tail component
      non_right_tail_values <- values[values < near_right_tail_threshold]
      
      if (length(non_right_tail_values) < 10) {
        # Too few non-right tail values
        return(list(
          processed_values = values,
          original_values = values,
          method = "preserve_right_tail_inflated",
          explanation = "Preserved right tail-inflated data (insufficient non-zero values)",
          distribution_type = distribution_type,
          is_zero_inflated = FALSE,
          is_extremely_skewed = is_extremely_skewed,
          right_tail_threshold = near_right_tail_threshold,
          skewness = skewness,
          is_right_tail_inflated = TRUE
        ))
      }
      
      # Check if non-right tail component is multimodal
      non_right_tail_multimodal <- .detect_multimodality(non_right_tail_values)
      
      if (non_right_tail_multimodal) {
        #  non-right tail component is multimodal, treat accordingly
        return(list(
          processed_values = non_right_tail_values,
          original_values = values,
          method = "right_tail_inflated_multimodal",
          explanation = "Extracted non-right tail values which show multimodal distribution",
          distribution_type = distribution_type, 
          is_zero_inflated = FALSE,
          is_right_tail_inflated = TRUE,
          is_multimodal = TRUE,
          is_extremely_skewed = is_extremely_skewed,
          right_tail_threshold = near_right_tail_threshold,
          skewness = skewness
        ))
      }
      
      # Process the non-right component
      non_right_tail_skewness <- moments::skewness(non_right_tail_values)
      
      if (abs(non_right_tail_skewness) > 0.8) {
        # Skewed non-right tail component, apply transformation
        lambda <- .find_optimal_lambda(non_right_tail_values)
        transformed_non_right_tail <- .yeo_johnson_transform(non_right_tail_values, lambda)
        
        return(list(
          processed_values = transformed_non_right_tail,
          original_values = values,
          non_right_tail_values = non_right_tail_values,
          method = "right_tail_inflated_with_transform",
          explanation = sprintf("Extracted non-right tail values (%.1f%%) and applied transformation",
                                100 * length(non_right_tail_values) / length(values)),
          distribution_type = distribution_type,
          is_zero_inflated = FALSE,
          is_right_tail_inflated = TRUE,
          is_extremely_skewed = is_extremely_skewed,
          right_tail_threshold = near_right_tail_threshold,
          lambda = lambda,
          transform_func = function(x) .yeo_johnson_transform(x, lambda),
          inverse_func = function(x) .yeo_johnson_inverse(x, lambda),
          skewness = skewness,
          non_right_tail_skewness = non_right_tail_skewness
        ))
      } else {
        # non-right tail component not highly skewed
        return(list(
          processed_values = non_right_tail_values,
          original_values = values,
          method = "right_tail_inflated_extract",
          explanation = sprintf("Extracted non-right tail values (%.1f%%) for threshold determination",
                                100 * length(non_right_tail_values) / length(values)),
          distribution_type = distribution_type,
          is_zero_inflated = FALSE,
          is_right_tail_inflated = TRUE,
          is_extremely_skewed = is_extremely_skewed,
          right_tail_threshold = near_right_tail_threshold,
          skewness = skewness,
          non_right_tail_skewness = non_right_tail_skewness
        ))
      }
    }
  } else if (is_multimodal) {
    # For multimodal data, preserve original distribution
    return(list(
      processed_values = values,
      original_values = values,
      method = "preserve_multimodal",
      explanation = "Preserved multimodal structure",
      distribution_type = distribution_type,
      is_multimodal = TRUE,
      is_zero_inflated = is_zero_inflated,
      is_extremely_skewed = is_extremely_skewed,
      is_right_tail_inflated = is_right_tail_inflated,
      skewness = skewness
    ))
  } else if (is_skewed && abs(skewness) > 0.8) {
    # For significantly skewed (but not extreme) data
    lambda <- .find_optimal_lambda(values)
    transformed_values <- .yeo_johnson_transform(values, lambda)
    
    return(list(
      processed_values = transformed_values,
      original_values = values,
      method = "yeo_johnson",
      explanation = sprintf("Applied Yeo-Johnson transformation (lambda = %.2f) for skewness", lambda),
      distribution_type = distribution_type,
      lambda = lambda,
      transform_func = function(x) .yeo_johnson_transform(x, lambda),
      inverse_func = function(x) .yeo_johnson_inverse(x, lambda),
      is_extremely_skewed = is_extremely_skewed,
      skewness = skewness
    ))
  } else if (handle_outliers && outlier_proportion < 0.1 && outlier_count > 0) {
    # Moderate outliers, temporarily exclude
    temp_values <- values[!is_outlier]
    return(list(
      processed_values = temp_values,
      original_values = values,
      method = "temp_exclude_outliers",
      explanation = sprintf("Temporarily excluded %d outliers (%.1f%%) for threshold determination",
                            outlier_count, 100 * outlier_proportion),
      distribution_type = distribution_type,
      is_extremely_skewed = is_extremely_skewed,
      excluded_indices = which(is_outlier),
      outlier_count = outlier_count,
      skewness = skewness
    ))
  } else {
    # No major issues - use original data
    return(list(
      processed_values = values,
      original_values = values,
      method = "no_transformation",
      explanation = "No preprocessing needed (distribution suitable for analysis)",
      distribution_type = distribution_type,
      is_extremely_skewed = is_extremely_skewed,
      skewness = skewness,
      is_zero_inflated = is_zero_inflated,
      is_right_tail_inflated = is_right_tail_inflated
    ))
  }
}

#' Improved multimodality detection
#'
#' Uses multiple methods to detect multimodality with balanced sensitivity and specificity
#'
#' @param values Numeric vector of values
#' @param verbose Logical indicating whether to print detailed information
#' @return Logical indicating if distribution is multimodal
#' @keywords internal
.detect_multimodality <- function(values, verbose = FALSE) {
  # Initialise all logical flags to ensure proper return type
  is_multimodal_dip <- FALSE
  is_multimodal_kde <- FALSE
  is_multimodal_gmm <- FALSE
  
  # Log information if verbose is TRUE
  if(verbose) {
    message("Checking for multimodality in data with ", length(values), " values")
  }
  
  # Method 1: Use Hartigan's dip test if available
  if (requireNamespace("diptest", quietly = TRUE)) {
    # Use try() to safely handle any errors during test execution
    dip_result <- try(diptest::dip.test(values), silent = TRUE)
    if (!inherits(dip_result, "try-error")) {
      is_multimodal_dip <- dip_result$p.value < 0.05
      if(verbose && is_multimodal_dip) {
        message("Dip test indicates multimodality (p-value: ", 
                round(dip_result$p.value, 4), ")")
      }
    }
  }
  
  # Method 2: Kernel Density Estimation (KDE) peak detection
  kde_result <- try({
    # Calculate basic statistics for reference
    median_val <- median(values)
    q25 <- quantile(values, 0.25)
    q75 <- quantile(values, 0.75)
    iqr <- q75 - q25
    data_range <- max(values) - min(values)
    
    # Calculate skewness for later use
    skewness <- moments::skewness(values)
    
    # Try slightly different bandwidths for better peak detection
    bandwidths <- c(0.8, 1.0, 1.2)
    
    for (bw in bandwidths) {
      # Calculate density with more aggressive smoothing for stability
      kde <- density(values, adjust = bw, n = min(1024, max(256, length(values))))
      
      # Find peaks in density using sign changes in the first derivative
      peaks <- which(diff(sign(diff(kde$y))) < 0) + 1
      
      # Find valleys for later analysis
      valleys <- which(diff(sign(diff(kde$y))) > 0) + 1
      
      # Filter out isolated peaks with limited data support
      filtered_peaks <- c()
      for (i in seq_along(peaks)) {
        peak_x <- kde$x[peaks[i]]
        
        # Calculate data density around peak
        window_width <- data_range * 0.1  
        in_window <- values >= (peak_x - window_width) & values <= (peak_x + window_width)
        data_proportion <- sum(in_window) / length(values)
        
        # requirements for tail peaks
        is_in_tail <- peak_x < quantile(values, 0.10) || peak_x > quantile(values, 0.90)
        min_support <- ifelse(is_in_tail, 0.15, 0.05)  
        
        # Calculate percentile position of peak
        percentile_position <- mean(values <= kde$x[peaks[i]])
        
        # For right skew, be very strict with right tail peaks
        if (skewness > 0.5 && percentile_position > 0.90) {
          min_support <- 0.20  
        }
        # For left skew, be very strict with left tail peaks
        else if (skewness < -0.5 && percentile_position < 0.15) {
          min_support <- 0.20  
        }
        
        # Only keep peaks with sufficient data support
        if (data_proportion >= min_support) {
          # For tail peaks, also check relative prominence
          if (is_in_tail) {
            # Find local background density
            tail_region <- which(if (peak_x < quantile(values, 0.10)) {
              kde$x <= quantile(values, 0.3)
            } else {
              kde$x >= quantile(values, 0.7)
            })
            
            if (length(tail_region) > 0) {
              # Calculate average density in the tail region
              avg_tail_density <- mean(kde$y[tail_region])
              
              # Require peak to be significantly higher than background
              peak_prominence_ratio <- kde$y[peaks[i]] / avg_tail_density
              
              # Skip this peak if not prominent enough
              if (peak_prominence_ratio < 1.0) {
                next
              }
            }
          }
          filtered_peaks <- c(filtered_peaks, peaks[i])
        }
      }
      
      # Merge peaks that are too close to each other
      merged_peaks <- c()
      if (length(filtered_peaks) > 0) {
        # Sort peaks by position
        sorted_indices <- order(kde$x[filtered_peaks])
        sorted_peaks <- filtered_peaks[sorted_indices]
        
        # Define minimum distance as a percentage of data range
        min_peak_distance <- data_range * 0.05
        
        # Start with first peak
        current_group <- c(sorted_peaks[1])
        
        # Group nearby peaks
        if (length(sorted_peaks) > 1) {
          for (i in 2:length(sorted_peaks)) {
            if ((kde$x[sorted_peaks[i]] - kde$x[sorted_peaks[i-1]]) < min_peak_distance) {
              # Add to current group
              current_group <- c(current_group, sorted_peaks[i])
            } else {
              # Select peak with highest density from current group
              highest_idx <- current_group[which.max(kde$y[current_group])]
              merged_peaks <- c(merged_peaks, highest_idx)
              # Start new group
              current_group <- c(sorted_peaks[i])
            }
          }
          
          # Add the last group
          highest_idx <- current_group[which.max(kde$y[current_group])]
          merged_peaks <- c(merged_peaks, highest_idx)
        } else {
          merged_peaks <- sorted_peaks
        }
        
        # Replace filtered_peaks with merged_peaks
        filtered_peaks <- merged_peaks
      }
      
      # Validate peaks using prominence
      validated_peaks <- c()
      for (p in filtered_peaks) {
        # Calculate peak prominence
        left_points <- which(kde$x < kde$x[p])
        right_points <- which(kde$x > kde$x[p])
        
        if (length(left_points) > 0 && length(right_points) > 0) {
          # Find lowest points between this peak and adjacent peaks
          left_valley <- max(left_points[kde$y[left_points] == min(kde$y[left_points[kde$x[left_points] > quantile(kde$x[left_points], 0.5)]])])
          right_valley <- min(right_points[kde$y[right_points] == min(kde$y[right_points[kde$x[right_points] < quantile(kde$x[right_points], 0.5)]])])
          
          # Calculate prominence (height above higher saddle)
          left_min_y <- kde$y[left_valley]
          right_min_y <- kde$y[right_valley]
          prominence <- kde$y[p] - max(left_min_y, right_min_y)
          rel_prominence <- prominence / kde$y[p]
          
          # Only keep peaks with sufficient prominence
          if (rel_prominence > 0.05) {
            validated_peaks <- c(validated_peaks, p)
          }
        }
      }
      
      # Use validated peaks
      if (length(validated_peaks) > 0) {
        filtered_peaks <- validated_peaks
      }
      
      # Replace original peaks with filtered ones
      peaks <- filtered_peaks
      
      if (length(peaks) >= 2) {
        # Extract peak heights and positions
        peak_heights <- kde$y[peaks]
        peak_positions <- kde$x[peaks]
        
        # Filter to retain only significant peaks 
        peak_height_threshold <- quantile(peak_heights, 0.25) 
        significant_peaks <- peaks[peak_heights > peak_height_threshold]
        
        if(verbose) {
          message("KDE analysis (bandwidth ", bw, ") found ", length(peaks), 
                  " peaks with ", length(significant_peaks), " significant peaks")
        }
        
        if (length(significant_peaks) >= 2) {
          # Sort peak positions for examining separation
          ordered_peaks <- order(kde$x[significant_peaks])
          ordered_peak_x <- kde$x[significant_peaks[ordered_peaks]]
          ordered_peak_heights <- peak_heights[significant_peaks[ordered_peaks] - min(peaks) + 1]
          
          # Calculate peak gaps
          peak_gaps <- diff(ordered_peak_x)
          
          # Consider distribution potentially multimodal if largest gap is substantial
          max_gap_relative <- max(peak_gaps) / data_range
          
          if (max_gap_relative > 0.05) {
            is_multimodal_kde <- TRUE
            if (verbose) {
              message("Distribution determined to be multimodal based on ", 
                      length(significant_peaks), " significant peaks with max gap ratio: ", 
                      round(max_gap_relative, 3))
            }
          }
        }
      }
    }
    
    # Return the calculated flag value
    is_multimodal_kde
  }, silent = TRUE)
  
  # Handle any errors that may have occurred during KDE analysis
  if (inherits(kde_result, "try-error") || !is.logical(kde_result)) {
    if(verbose && inherits(kde_result, "try-error")) {
      message("KDE analysis failed: ", as.character(attributes(kde_result)$condition))
    }
    is_multimodal_kde <- FALSE
  }
  
  # Method 3: Test for multimodality using Gaussian mixture models (GMM)
  if (requireNamespace("mclust", quietly = TRUE)) {
    gmm_result <- try({
      # Fit mixture models with 1-3 components
      model <- suppressWarnings(mclust::Mclust(values, G = 1:3, verbose = FALSE))
      
      # Initialise local GMM flag
      local_is_multimodal_gmm <- FALSE
      
      # Check if multiple components fit better
      if (!is.null(model) && model$G > 1) {
        # Extract component means
        means <- model$parameters$mean
        
        if(verbose) {
          message("GMM fitted with ", model$G, " components")
        }
        
        # Handle different variance structures in mclust models
        if (is.list(model$parameters$variance)) {
          variances <- model$parameters$variance$sigmasq
        } else {
          variances <- rep(model$parameters$variance, length(means))
        }
        
        # Calculate standard deviations
        sds <- sqrt(variances)
        
        # Sort components by mean value
        means_order <- order(means)
        sorted_means <- means[means_order]
        sorted_sds <- sds[means_order]
        
        # Calculate standardised distances between adjacent component means
        distances <- numeric(length(means) - 1)
        for (i in 1:(length(means) - 1)) {
          distances[i] <- abs(sorted_means[i+1] - sorted_means[i]) /
            sqrt((sorted_sds[i]^2 + sorted_sds[i+1]^2) / 2)
        }
        
        max_distance <- max(distances)
        bic_diff <- model$BIC[model$G] - model$BIC[1]
        adaptive_threshold <- max(0.5, min(1.5, exp(-bic_diff/200)))
        
        if (verbose) {
          message("Max standardised distance: ", round(max_distance, 4), 
                  ", Adaptive threshold: ", round(adaptive_threshold, 4), 
                  " (BIC diff: ", round(bic_diff, 2), ")")
        }
        
        if(!is.na(max_distance) && is.finite(max_distance) && max_distance > adaptive_threshold) {
          local_is_multimodal_gmm <- TRUE
          if (verbose) {
            message("GMM indicates multimodality")
          }
        }
      }
      
      # Return the calculated GMM flag
      local_is_multimodal_gmm
    }, silent = TRUE)
    
    # Handle errors and ensure logical return type
    if (inherits(gmm_result, "try-error")) {
      if(verbose) {
        message("GMM analysis failed: ", as.character(attributes(gmm_result)$condition))
      }
      is_multimodal_gmm <- FALSE
    } else if (!is.logical(gmm_result)) {
      is_multimodal_gmm <- FALSE
    } else {
      is_multimodal_gmm <- gmm_result
    }
  }
  
  # Combine results - determine if any method detected multimodality
  # Apply as.logical to ensure proper type conversion
  # Handle potential NA values by providing a default
  method_count <- sum(c(is_multimodal_dip, is_multimodal_kde, is_multimodal_gmm))
  result <- as.logical(method_count >= 2)
  
  if (is.na(result)) {
    result <- FALSE
  }
  
  if(verbose) {
    message("Final multimodality determination: ", 
            ifelse(result, "multimodal", "not multimodal"),
            " (dip: ", is_multimodal_dip, 
            ", kde: ", is_multimodal_kde, 
            ", gmm: ", is_multimodal_gmm, ")")
  }
  
  return(result)
}

#' Apply Yeo-Johnson transformation
#'
#' @param x Numeric vector to transform
#' @param lambda Transformation parameter
#' @return Transformed values
#' @keywords internal
.yeo_johnson_transform <- function(x, lambda) {
  if (is.null(lambda)) lambda <- 1  # No transformation
  
  pos_idx <- x >= 0
  neg_idx <- !pos_idx
  
  result <- x
  
  # Transform positive values
  if (any(pos_idx)) {
    if (abs(lambda) < 1e-6) {
      # Near zero lambda case (log transformation)
      result[pos_idx] <- log(x[pos_idx] + 1)
    } else {
      # Standard case
      result[pos_idx] <- ((x[pos_idx] + 1)^lambda - 1) / lambda
    }
  }
  
  # Transform negative values
  if (any(neg_idx)) {
    if (abs(lambda - 2) < 1e-6) {
      # Special case near lambda = 2
      result[neg_idx] <- -log(-x[neg_idx] + 1)
    } else {
      result[neg_idx] <- -((-x[neg_idx] + 1)^(2 - lambda) - 1) / (2 - lambda)
    }
  }
  
  return(result)
}

#' Inverse Yeo-Johnson transformation
#'
#' @param x Transformed values
#' @param lambda Transformation parameter
#' @return Original scale values
#' @keywords internal
.yeo_johnson_inverse <- function(x, lambda) {
  if (is.null(lambda)) return(x)  # No transformation
  
  pos_idx <- x >= 0
  neg_idx <- !pos_idx
  
  result <- x
  
  # Inverse transform positive values
  if (any(pos_idx)) {
    if (abs(lambda) < 1e-6) {
      # Near zero lambda case (exp transformation)
      result[pos_idx] <- exp(x[pos_idx]) - 1
    } else {
      # Standard case
      result[pos_idx] <- (lambda * x[pos_idx] + 1)^(1/lambda) - 1
    }
  }
  
  # Inverse transform negative values
  if (any(neg_idx)) {
    if (abs(lambda - 2) < 1e-6) {
      # Special case near lambda = 2
      result[neg_idx] <- -(exp(-x[neg_idx]) - 1)
    } else {
      result[neg_idx] <- -(((2 - lambda) * (-x[neg_idx]) + 1)^(1/(2 - lambda)) - 1)
    }
  }
  
  return(result)
}

#' Determine optimal Yeo-Johnson transformation parameter
#'
#' @param values Numeric vector to transform
#' @return Optimal lambda parameter
#' @keywords internal
.find_optimal_lambda <- function(values) {
  if (requireNamespace("car", quietly = TRUE)) {
    # Use car package's powerTransform function if available
    lambda <- try(car::powerTransform(values, family = "yjPower")$lambda, silent = TRUE)
    
    if (!inherits(lambda, "try-error") && is.finite(lambda)) {
      return(lambda)
    }
  }
  
  # Simplified lambda selection using grid search
  lambda_candidates <- seq(-2, 2, by = 0.1)
  best_lambda <- 1  # Default to no transformation
  best_normality <- -Inf
  
  for (l in lambda_candidates) {
    transformed <- .yeo_johnson_transform(values, l)
    
    # Assess normality
    q25 <- quantile(transformed, 0.25)
    q50 <- quantile(transformed, 0.50)
    q75 <- quantile(transformed, 0.75)
    # Symmetry measure
    normality_score <- -abs((q75 - q50) - (q50 - q25))
    
    if (normality_score > best_normality) {
      best_normality <- normality_score
      best_lambda <- l
    }
  }
  
  return(best_lambda)
}

#' Adaptive curve inflection detection with skewness-based search strategy
#'
#' @param analysis_values Numeric vector of values
#' @param verbose Logical indicating whether to print detailed information
#' @return List with inflection detection results
#' @keywords internal
.detect_curve_inflection <- function(analysis_values, verbose = FALSE) {
  # Generate high-resolution density estimation
  kde <- density(analysis_values, adjust = 1.0, n = 1024)
  
  # Get basic statistics
  min_val <- min(analysis_values)
  max_val <- max(analysis_values)
  range_val <- max_val - min_val
  median_val <- median(analysis_values)
  peak_idx <- which.max(kde$y)
  peak_x <- kde$x[peak_idx]
  max_density <- max(kde$y)
  
  # Calculate peak position (normalised)
  peak_position <- (peak_x - min_val) / range_val
  
  # Calculate skewness to determine distribution shape
  skewness <- moments::skewness(analysis_values)
  is_right_skewed <- skewness > 0.5
  is_left_skewed <- skewness < -0.5
  
  # Calculate first derivative (slope)
  first_deriv <- diff(kde$y) / diff(kde$x)
  
  # Smooth derivatives for improved stability
  smooth_deriv <- stats::filter(first_deriv, rep(1/5, 5), sides = 2)
  smooth_deriv <- smooth_deriv[!is.na(smooth_deriv)]
  
  # Calculate second derivative (curvature)
  second_deriv <- diff(smooth_deriv) / diff(kde$x[-1])
  second_deriv <- stats::filter(second_deriv, rep(1/5, 5), sides = 2)
  second_deriv <- second_deriv[!is.na(second_deriv)]
  
  if (verbose) {
    message("Searching for inflection points in distribution with skewness ", 
            round(skewness, 2), " and peak position at ", round(peak_position, 2))
  }
  
  # Determine search ranges based on skewness
  left_search_range <- NULL
  right_search_range <- NULL
  
  # Left search range (before peak)
  left_search_start <- max(1, peak_idx - round(peak_idx * 0.8))
  left_search_end <- max(1, peak_idx - 3)  # Allow a small gap from peak
  
  # Right search range (after peak)
  right_search_start <- min(length(kde$x), peak_idx + 3)  # Allow a small gap from peak
  right_search_end <- min(length(second_deriv), peak_idx + round((length(kde$x) - peak_idx) * 0.8))
  
  # Adjust search strategy based on distribution shape
  if (is_right_skewed) {
    # For right-skewed distributions, focus more on right side
    if (verbose) message("Right-skewed distribution detected, focusing on right side")
    
    # Prioritise right search with extended range
    right_search_end <- min(length(second_deriv), peak_idx + round((length(kde$x) - peak_idx) * 0.9))
    right_search_range <- right_search_start:right_search_end
    
    # Still include left side but with reduced priority
    left_search_range <- left_search_start:left_search_end
    
  } else if (is_left_skewed) {
    # For left-skewed distributions, focus more on left side
    if (verbose) message("Left-skewed distribution detected, focusing on left side")
    
    # Prioritise left search with extended range
    left_search_start <- max(1, peak_idx - round(peak_idx * 0.9))
    left_search_range <- left_search_start:left_search_end
    
    # Still include right side but with reduced priority
    right_search_range <- right_search_start:right_search_end
    
  } else {
    # For symmetric distributions, search both sides equally
    if (verbose) message("Symmetric distribution detected, searching both sides equally")
    left_search_range <- left_search_start:left_search_end
    right_search_range <- right_search_start:right_search_end
  }
  
  if (verbose) {
    message("Left search range: ", length(left_search_range), " points, Right search range: ", 
            length(right_search_range), " points")
  }
  
  # Ensure search ranges are valid
  valid_left_search <- !is.null(left_search_range) && length(left_search_range) >= 5
  valid_right_search <- !is.null(right_search_range) && length(right_search_range) >= 5
  
  if (!valid_left_search && !valid_right_search) {
    return(list(has_inflection = FALSE, reason = "insufficient_data_for_search"))
  }
  
  # Find significant negative curvature points in both ranges
  curvature_threshold <- -0.1 * max(abs(second_deriv))
  
  # Store inflection points from both sides
  left_neg_curve_points <- c()
  right_neg_curve_points <- c()
  
  # Search left side if valid
  if (valid_left_search) {
    left_neg_curve_points <- left_search_range[second_deriv[left_search_range] < curvature_threshold]
  }
  
  # Search right side if valid
  if (valid_right_search) {
    right_neg_curve_points <- right_search_range[second_deriv[right_search_range] < curvature_threshold]
  }
  
  # Combine and filter all negative curvature points
  neg_curve_points <- c(left_neg_curve_points, right_neg_curve_points)
  filtered_neg_curve_points <- c()
  
  if (length(neg_curve_points) > 0) {
    for (point_idx in neg_curve_points) {
      curve_x <- kde$x[point_idx + 1]  # +1 due to diff operation
      
      # Check if this point is in the tail region
      percentile_position <- mean(analysis_values <= curve_x)
      is_in_tail <- percentile_position < 0.15 || percentile_position > 0.85
      
      if (is_in_tail) {
        # Apply stricter filtering for tail points
        window_width <- range_val * 0.1
        in_window <- analysis_values >= (curve_x - window_width) & 
          analysis_values <= (curve_x + window_width)
        data_proportion <- sum(in_window) / length(analysis_values)
        
        # Skip if insufficient data support in the tail
        if (data_proportion < 0.1) {
          next
        }
        
        # Check prominence - relative height of inflection point
        tail_region <- if (percentile_position < 0.15) {
          which(kde$x <= quantile(kde$x, 0.3))
        } else {
          which(kde$x >= quantile(kde$x, 0.7))
        }
        
        # Skip if tail region is too small
        if (length(tail_region) < 5) {
          next
        }
        
        # Calculate average density in tail region
        avg_tail_density <- mean(kde$y[tail_region])
        density_at_point <- kde$y[point_idx]
        
        # Skip tail points without sufficient prominence
        if (density_at_point / avg_tail_density < 1.5) {
          next
        }
      }
      
      # Add to filtered points
      filtered_neg_curve_points <- c(filtered_neg_curve_points, point_idx)
    }
    
    # If we have filtered points, replace the original ones
    if (length(filtered_neg_curve_points) > 0) {
      neg_curve_points <- filtered_neg_curve_points
    }
  }
  
  if (length(neg_curve_points) > 0) {
    # Apply skewness-specific selection strategy
    if (is_right_skewed) {
      # Prioritise points on the right side
      right_side_points <- neg_curve_points[neg_curve_points > peak_idx]
      if (length(right_side_points) > 0) {
        neg_curve_points <- right_side_points
      }
    } else if (is_left_skewed) {
      # Prioritise points on the left side
      left_side_points <- neg_curve_points[neg_curve_points < peak_idx]
      if (length(left_side_points) > 0) {
        neg_curve_points <- left_side_points
      }
    }
    
    # Find the strongest negative curvature point
    strongest_curve_idx <- neg_curve_points[which.min(second_deriv[neg_curve_points])]
    curve_x <- kde$x[strongest_curve_idx + 1]  # +1 due to diff operation
    curve_pos <- (curve_x - min_val) / range_val
    
    neg_curvatures <- sort(abs(second_deriv[second_deriv < 0]))
    percentile_pos <- which(neg_curvatures >= abs(second_deriv[strongest_curve_idx]))[1] / length(neg_curvatures)
    curve_strength <- min(0.95, percentile_pos)
    
    
    # Determine which side the selected point came from
    side <- ifelse(strongest_curve_idx < peak_idx, "left", "right")
    
    # Calculate data balance
    prop_below <- mean(analysis_values <= curve_x)
    prop_above <- 1 - prop_below
    
    # Basic check: ensure inflection point is in a reasonable position (avoid extreme edges)
    if (curve_pos > 0.2 && curve_pos < 0.8) {
      # Post check: ensure there's sufficient data on both sides of the inflection
      min_proportion <- 0.05  # Default minimum proportion
      
      # Adjust minimum proportion based on skewness
      if (is_right_skewed && side == "right") {
        min_proportion <- 0.05  # Allow more imbalanced splits for right-side points in right-skewed data
      } else if (is_left_skewed && side == "left") {
        min_proportion <- 0.05  # Allow more imbalanced splits for left-side points in left-skewed data
      }
      
      # Check data balance with adjusted minimum
      if (min(prop_below, prop_above) > min_proportion) {
        # Calculate reliability - based on curvature strength and position
        # Position preference depends on distribution shape
        position_score <- 0
        
        if (is_right_skewed) {
          if (side == "right") {
            position_score <- 1 - min(abs(curve_pos - 0.3) * 1.5, 1)
          } else {
            position_score <- 0.5  # Lower score for left-side points
          }
        } else if (is_left_skewed) {
          if (side == "left") {
            position_score <- 1 - min(abs(curve_pos - 0.7) * 1.5, 1)
          } else {
            position_score <- 0.5  # Lower score for right-side points
          }
        } else {
          # For symmetric, prefer central points
          position_score <- 1 - min(abs(curve_pos - 0.5) * 1.5, 1)
        }
        
        
        # Final reliability calculation
        reliability <- min(0.95, curve_strength * 0.6 + position_score * 0.4)
        
        # Adjust reliability based on distribution-specific factors
        if (is_right_skewed && side == "right" && curve_pos > 0.5) {
          reliability <- min(0.95, reliability + 0.05)  # Slight boost for ideal points
        } else if (is_left_skewed && side == "left" && curve_pos < 0.5) {
          reliability <- min(0.95, reliability + 0.05)  # Slight boost for ideal points
        }
        
        if (verbose) {
          message(sprintf("Detected inflection point: %.3f (normalised: %.2f), strength: %.2f, reliability: %.2f (%s of peak)", 
                          curve_x, curve_pos, curve_strength, reliability, side))
        }
        
        return(list(
          has_inflection = TRUE,
          inflection_position = curve_x,
          inflection_strength = curve_strength,
          normalised_position = curve_pos,
          reliability = reliability,
          side_of_peak = side,
          skewness = skewness
        ))
      }
    }
  }
  
  # No suitable inflection point found
  return(list(has_inflection = FALSE, reason = "no_significant_inflection"))
}

#' Improved multimodal threshold detection with enhanced mixture model logic
#'
#' @param analysis_values Processed values to analyse
#' @param preprocessing Preprocessing result
#' @param verbose Whether to print detailed messages
#' @return List containing threshold and method details
#' @keywords internal
.multimodal_threshold_method <- function(analysis_values, preprocessing, verbose = FALSE) {
  
  # Get distribution statistics
  min_val <- min(analysis_values)
  max_val <- max(analysis_values)
  range_val <- max_val - min_val
  median_val <- median(analysis_values)
  
  # Improved mixture modelling approach with better intersection selection
  mixture_result <- try({
    # Define verbose if not provided by parent function
    if (!exists("verbose")) verbose <- FALSE
    
    if (requireNamespace("mclust", quietly = TRUE)) {
      # Fit mixture models with 3 components 
      model <- suppressWarnings(mclust::Mclust(
        analysis_values,
        G = 3,  
        verbose = FALSE
      ))
      
      if (!is.null(model) && !is.null(model$parameters)) {
        # Extract model parameters
        means <- model$parameters$mean
        
        # Define number of components based on the selected model
        n_components <- length(means)
        
        if (verbose) message(sprintf("Mclust selected %d components", n_components))
        
        # Get component weights
        weights <- model$parameters$pro
        
        # Extract standard deviations
        if ("sigmasq" %in% names(model$parameters$variance)) {
          sds <- sqrt(model$parameters$variance$sigmasq)
        } else if (is.numeric(model$parameters$variance)) {
          sds <- rep(sqrt(model$parameters$variance), n_components)
        } else {
          # Extract from variance matrices
          sds <- numeric(n_components)
          for (i in 1:n_components) {
            sds[i] <- sqrt(model$parameters$variance$sigma[1, 1, i])
          }
        }
        
        # Sort components by mean for easier interpretation
        sorted_idx <- order(means)
        means <- means[sorted_idx]
        sds <- sds[sorted_idx]
        weights <- weights[sorted_idx]
        
        if (verbose) {
          message("Component means: ", paste(round(means, 3), collapse=", "))
          message("Component SDs: ", paste(round(sds, 3), collapse=", "))
          message("Component weights: ", paste(round(weights, 3), collapse=", "))
        }
        
        # Compute KDE for reference and density evaluation
        kde <- density(analysis_values, adjust = 0.8, n = 512)
        kde_max <- max(kde$y)
        
        # Consider all possible component pairs
        threshold_candidates <- list()
        
        for (i in 1:(n_components-1)) {
          for (j in (i+1):n_components) {
            m1 <- means[i]
            m2 <- means[j]
            s1 <- sds[i]
            s2 <- sds[j]
            w1 <- weights[i]
            w2 <- weights[j]
            
            # Skip component pairs that are too close (likely part of same mode)
            if (abs(m2 - m1) < (s1 + s2) * 0.5) {
              if (verbose) message(sprintf("Skipping components %d and %d as they are too close", i, j))
              next
            }
            
            # Find intersection of the two density curves
            threshold <- NULL
            
            # Compute weighted density functions
            weighted_density1 <- function(x) w1 * dnorm(x, m1, s1)
            weighted_density2 <- function(x) w2 * dnorm(x, m2, s2)
            
            # Calculate the distance between means in terms of SDs
            distance_in_sds <- abs(m2 - m1) / sqrt((s1^2 + s2^2)/2)
            
            if (abs(s1 - s2) < 1e-6) {
              # Equal variances - analytical solution
              if (w1 == w2) {
                # With equal weights, intersection is just the midpoint
                threshold <- (m1 + m2) / 2
              } else {
                # With unequal weights, weighted midpoint
                # This is the analytical solution for equal variance Gaussians
                threshold <- ((m1 * w2) + (m2 * w1)) / (w1 + w2)
              }
              if (verbose) message("Equal variances, analytical solution")
            } else {
              # Different variances - find intersections
              
              # Create a search grid focused between the means
              grid_points <- seq(
                min(m1 - 2*s1, m2 - 2*s2), 
                max(m1 + 2*s1, m2 + 2*s2), 
                length.out = 1000
              )
              
              # Calculate density difference at each point
              density_diffs <- weighted_density1(grid_points) - weighted_density2(grid_points)
              
              # Find sign changes (intersections)
              sign_changes <- which(diff(sign(density_diffs)) != 0)
              
              if (length(sign_changes) > 0) {
                intersections <- numeric(length(sign_changes))
                
                # For each sign change, find the precise intersection
                for (k in 1:length(sign_changes)) {
                  idx <- sign_changes[k]
                  x0 <- grid_points[idx]
                  x1 <- grid_points[idx + 1]
                  
                  # Refine intersection using linear interpolation
                  y0 <- density_diffs[idx]
                  y1 <- density_diffs[idx + 1]
                  
                  # Linear interpolation for more accurate intersection
                  intersections[k] <- x0 - y0 * (x1 - x0) / (y1 - y0)
                }
                
                # Calculate densities at intersections for both components
                intersection_densities <- sapply(intersections, function(x) {
                  d1 <- weighted_density1(x)
                  d2 <- weighted_density2(x)
                  # Return average of the two densities (they should be equal at intersection)
                  return((d1 + d2) / 2)
                })
                
                # Filter out intersections that occur at very low densities
                # This helps avoid choosing trivial intersections where the curves barely meet
                min_density_threshold <- 0.05 * max(
                  sapply(means, function(mu) weighted_density1(mu)),
                  sapply(means, function(mu) weighted_density2(mu))
                )
                
                valid_intersections <- intersections[intersection_densities >= min_density_threshold]
                valid_densities <- intersection_densities[intersection_densities >= min_density_threshold]
                
                if (length(valid_intersections) > 0) {
                  # Among valid intersections, prioritise those between the means
                  between_means <- valid_intersections >= min(m1, m2) & valid_intersections <= max(m1, m2)
                  
                  if (any(between_means)) {
                    # If there are intersections between means, prefer the one with highest density
                    between_idx <- which(between_means)
                    max_density_idx <- between_idx[which.max(valid_densities[between_means])]
                    threshold <- valid_intersections[max_density_idx]
                  } else {
                    # If no intersections between means, select the closest to weighted mean
                    weighted_mean <- (m1 * w1 + m2 * w2) / (w1 + w2)
                    threshold <- valid_intersections[which.min(abs(valid_intersections - weighted_mean))]
                  }
                } else {
                  # Fallback to weighted midpoint if no valid intersections found
                  threshold <- (m1 * w1 + m2 * w2) / (w1 + w2)
                  if (verbose) message("No valid intersections with sufficient density, using weighted midpoint")
                }
              } else {
                # No intersections found - use weighted midpoint
                threshold <- (m1 * w1 + m2 * w2) / (w1 + w2)
                if (verbose) message("No intersections found, using weighted midpoint")
              }
            }
            
            # Ensure threshold is within data range
            data_range <- range(analysis_values, na.rm = TRUE)
            threshold <- min(max(threshold, data_range[1]), data_range[2])
            
            # Calculate separation quality
            separation <- distance_in_sds
            
            
            # Evaluate threshold quality
            # Get density at threshold using KDE (from actual data)
            threshold_kde_idx <- which.min(abs(kde$x - threshold))
            threshold_density <- kde$y[threshold_kde_idx] / kde_max  # Normalise to [0,1]
            
            # Calculate mode prominence for both components
            mode1_idx <- which.min(abs(kde$x - m1))
            mode2_idx <- which.min(abs(kde$x - m2))
            mode1_height <- kde$y[mode1_idx] / kde_max
            mode2_height <- kde$y[mode2_idx] / kde_max
            
            # Skip if one mode is insignificant (< 10% of max density)
            if (mode1_height < 0.1 || mode2_height < 0.1) {
              if (verbose) message(sprintf("Skipping components %d and %d as one mode is insignificant", i, j))
              next
            }
            
            
            # Store candidate
            threshold_candidates[[length(threshold_candidates) + 1]] <- list(
              threshold = threshold,
              components = c(i, j),
              separation = separation,
              component_means = c(m1, m2),
              component_sds = c(s1, s2),
              component_weights = c(w1, w2)
            )
            
            if (verbose) {
              message(sprintf("Candidate threshold %.3f between components %d and %d", 
                              threshold, i, j))
              message(sprintf("  Separation: %.2f", separation))
            }
          }
        }
        
        # Filter for candidates with reasonable separation
        valid_candidates <- threshold_candidates[sapply(threshold_candidates, 
                                                        function(x) x$separation > 0.8)]
        
        # If no candidates with sufficient separation, try less strict criteria
        if (length(valid_candidates) == 0) {
          valid_candidates <- threshold_candidates[sapply(threshold_candidates, 
                                                          function(x) x$separation > 0.5)]
          if (verbose && length(valid_candidates) > 0) {
            message("Using less strict separation criteria")
          }
        }
        
        if (length(valid_candidates) > 0) {
          # Sort by separation
          separations <- sapply(valid_candidates, function(x) x$separation)
          best_idx <- which.max(separations)
          best_candidate <- valid_candidates[[best_idx]]
          
          method_type <- "multimodal.mixture.model"
          
          # Adjust reliability based on separation
          reliability <- min(0.95, 0.5 + best_candidate$separation * 0.15 )
          
          if (verbose) {
            message(sprintf("Selected threshold at %.3f with separation %.3f and reliability %.2f", 
                            best_candidate$threshold, best_candidate$separation, reliability))
            message("SUCCESS: Mixture model method successfully found a threshold")
          }
          
          
          return(list(
            threshold = best_candidate$threshold,
            method = method_type,
            reliability = reliability,
            separation = best_candidate$separation,
            selected_components = best_candidate$components
          ))
          
        } else if (verbose) {
          message("No mixture component pairs with sufficient separation found")
          message("WARNING: Mixture model method failed to find a suitable threshold, trying other methods")
        }
      }
    }
  }, silent = TRUE)
  
  if (!inherits(mixture_result, "try-error") && is.list(mixture_result)) {
    return(mixture_result)
  } else if (inherits(mixture_result, "try-error")) {
    if (verbose) {
      message("ERROR in mixture model method: ", conditionMessage(attr(mixture_result, "condition")))
      message("Falling back to alternative methods")
    }
  } else if (!is.list(mixture_result)) {
    if (verbose) {
      message("WARNING: Mixture model method did not return usable results, trying other methods")
    }
  }
  
  # Try KDE approach
  # Try multiple bandwidths to better identify peaks and valleys
  best_threshold <- NULL
  best_score <- -Inf
  best_method <- ""
  best_peaks <- NULL  
  
  # Try different bandwidths for more robust peak detection
  bandwidths <- c(0.8, 1.0, 1.2)
  
  for (bw in bandwidths) {
    # Generate density estimation
    kde <- density(analysis_values, adjust = bw, n = 512)
    
    # Find peaks (maxima)
    peaks <- which(diff(sign(diff(kde$y))) < 0) + 1
    peak_heights <- kde$y[peaks]
    peak_x <- kde$x[peaks]
    
    # Find valleys (minima)
    valleys <- which(diff(sign(diff(kde$y))) > 0) + 1
    valley_heights <- kde$y[valleys]
    valley_x <- kde$x[valleys]
    
    # Need at least 2 peaks and 1 valley
    if (length(peaks) >= 2 && length(valleys) >= 1) {
      # Filter out isolated peaks with limited data support
      filtered_peaks <- c()
      for (i in seq_along(peaks)) {
        peak_x_val <- kde$x[peaks[i]]
        
        # Calculate data density around peak
        window_width <- range_val * 0.08  
        in_window <- analysis_values >= (peak_x_val - window_width) & 
          analysis_values <= (peak_x_val + window_width)
        data_proportion <- sum(in_window) / length(analysis_values)
        
        # More stringent requirements for tail peaks
        is_in_tail <- peak_x_val < quantile(analysis_values, 0.15) || 
          peak_x_val > quantile(analysis_values, 0.85)
        min_support <- ifelse(is_in_tail, 0.2, 0.05)  
        
        # Only keep peaks with sufficient data support
        if (data_proportion >= min_support) {
          # For tail peaks, also check relative prominence
          if (is_in_tail) {
            # Calculate percentile position
            percentile_pos <- mean(analysis_values <= peak_x_val)
            
            # Find local background density
            tail_region <- if (percentile_pos < 0.15) {
              kde$x <= quantile(analysis_values, 0.3)
            } else {
              kde$x >= quantile(analysis_values, 0.7)
            }
            
            # Calculate average density in region
            avg_density <- mean(kde$y[tail_region])
            
            # Calculate prominence ratio
            prominence_ratio <- kde$y[peaks[i]] / avg_density
            
            # Skip if not prominent enough
            if (prominence_ratio < 1.8) {
              next
            }
          }
          
          filtered_peaks <- c(filtered_peaks, peaks[i])
        }
      }
      
      # Merge peaks that are too close to each other
      merged_peaks <- c()
      if (length(filtered_peaks) > 0) {
        # Sort peaks by position
        sorted_indices <- order(kde$x[filtered_peaks])
        sorted_peaks <- filtered_peaks[sorted_indices]
        
        # Define minimum distance as a percentage of data range
        min_peak_distance <- range_val * 0.1
        
        # Start with first peak
        current_group <- c(sorted_peaks[1])
        
        # Group nearby peaks
        if (length(sorted_peaks) > 1) {
          for (i in 2:length(sorted_peaks)) {
            if ((kde$x[sorted_peaks[i]] - kde$x[sorted_peaks[i-1]]) < min_peak_distance) {
              # Add to current group
              current_group <- c(current_group, sorted_peaks[i])
            } else {
              # Select peak with highest density from current group
              highest_idx <- current_group[which.max(kde$y[current_group])]
              merged_peaks <- c(merged_peaks, highest_idx)
              # Start new group
              current_group <- c(sorted_peaks[i])
            }
          }
          
          # Add the last group
          highest_idx <- current_group[which.max(kde$y[current_group])]
          merged_peaks <- c(merged_peaks, highest_idx)
        } else {
          merged_peaks <- sorted_peaks
        }
        
        # Replace filtered_peaks with merged_peaks
        filtered_peaks <- merged_peaks
      }
      
      # Keep the original code below, but use our filtered and merged peaks
      
      # More stringent peak significance threshold
      significant_peaks <- filtered_peaks[kde$y[filtered_peaks] > 0.1 * max(kde$y)]
      significant_peak_x <- kde$x[significant_peaks]
      significant_peak_heights <- kde$y[significant_peaks]
      
      if (length(significant_peaks) >= 2) {
        # Sort peaks by position
        ordered_peaks <- order(significant_peak_x)
        ordered_peak_x <- significant_peak_x[ordered_peaks]
        ordered_peak_heights <- significant_peak_heights[ordered_peaks]
        
        # Calculate data density around each peak
        peak_data_proportions <- numeric(length(ordered_peak_x))
        for (i in seq_along(ordered_peak_x)) {
          # Calculate window around peak
          window_width <- range_val * 0.1
          in_window <- analysis_values >= (ordered_peak_x[i] - window_width) &
            analysis_values <= (ordered_peak_x[i] + window_width)
          peak_data_proportions[i] <- sum(in_window) / length(analysis_values)
        }
        
        # Evaluate each valley between peaks
        valley_scores <- numeric()
        valley_thresholds <- numeric()
        
        for (i in 1:(length(ordered_peak_x)-1)) {
          # Find valleys between adjacent peaks
          current_peak_x <- ordered_peak_x[i]
          next_peak_x <- ordered_peak_x[i+1]
          
          between_valleys <- valleys[valley_x > current_peak_x &
                                       valley_x < next_peak_x]
          
          if (length(between_valleys) > 0) {
            for (v_idx in seq_along(between_valleys)) {
              v <- between_valleys[v_idx]
              potential_threshold <- kde$x[v]
              
              # Better valley depth calculation
              left_peak_height <- ordered_peak_heights[i]
              right_peak_height <- ordered_peak_heights[i+1]
              valley_height <- kde$y[v]
              
              # Depth relative to lower peak 
              min_peak_height <- min(left_peak_height, right_peak_height)
              depth_ratio <- (min_peak_height - valley_height) / min_peak_height
              depth_score <- min(1, depth_ratio * 3)  # Cap at 1
              
              # Increase requirement for valley depth
              if (depth_ratio < 0.4) {  
                next  # Skip valleys that aren't deep enough
              }
              
              # Position scoring to avoid left-bias
              # Calculate normalised position
              position <- (potential_threshold - min_val) / range_val
              
              # Strong penalty for extreme left positions, moderate for extreme right
              if (position < 0.25) {
                # Quadratic penalty for far left positions
                position_score <- position^2 * 16 
              } else if (position > 0.75) {
                # Linear penalty for far right positions 
                position_score <- (1 - position) * 4  
              } else {
                # Full score for central positions
                position_score <- 1.0
              }
              
              # Better balance scoring
              # Proportion of data on each side of threshold
              prop_below <- mean(analysis_values <= potential_threshold)
              prop_above <- 1 - prop_below
              
              # Balance ratio that punishes extreme imbalance
              balance_ratio <- min(prop_below, prop_above) / max(prop_below, prop_above)
              
              # Apply non-linear transformation to heavily penalise extreme imbalance
              balance_score <- balance_ratio^0.7  # Power < 1 makes moderate imbalance ok
              
              # Peak importance based on data distribution
              # Peaks with more data should be more important
              left_peak_importance <- peak_data_proportions[i]
              right_peak_importance <- peak_data_proportions[i+1]
              peaks_importance_score <- (left_peak_importance + right_peak_importance)
              
              # Adjust scoring weights for better valley selection
              if (position < 0.3) {
                # Far left - very strong bias correction
                valley_score <- 0.4 * depth_score + 0.35 * balance_score +
                  0.15 * position_score + 0.1 * peaks_importance_score
              } else if (position < 0.4) {
                # Left side - strong bias correction
                valley_score <- 0.45 * depth_score + 0.25 * balance_score +
                  0.2 * position_score + 0.1 * peaks_importance_score
              } else {
                # Centre/right - balanced weights
                valley_score <- 0.55 * depth_score + 0.15 * balance_score +
                  0.2 * position_score + 0.1 * peaks_importance_score
              }
              
              # Save scores and thresholds
              valley_scores <- c(valley_scores, valley_score)
              valley_thresholds <- c(valley_thresholds, potential_threshold)
              
              if (verbose) {
                message(sprintf("Valley at %.3f (pos=%.2f): depth=%.2f, balance=%.2f (%.1f%%/%.1f%%), score=%.3f",
                                potential_threshold, position, depth_score,
                                balance_score, 100*prop_below, 100*prop_above, valley_score))
              }
            }
          }
        }
        
        # Choose the best valley
        if (length(valley_scores) > 0) {
          best_valley_idx <- which.max(valley_scores)
          potential_threshold <- valley_thresholds[best_valley_idx]
          potential_score <- valley_scores[best_valley_idx]
          
          # More conservative update criterion
          # Only update if significantly better than previous
          if (potential_score > best_score + 0.05) {
            best_score <- potential_score
            best_threshold <- potential_threshold
            best_peaks <- significant_peaks  # Store the significant peaks for this bandwidth
            
            # More descriptive method naming
            valley_quality <- ifelse(potential_score > 0.8, "deep",
                                     ifelse(potential_score > 0.5, "moderate", "shallow"))
            best_method <- paste0("multimodal.kde.valley.", valley_quality)
          }
        }
      }
    }
  }
  
  # If KDE method found a threshold, return it
  if (!is.null(best_threshold)) {
    # Determine reliability based on score
    if (best_score > 0.8) {
      reliability <- 0.95
    } else if (best_score > 0.6) {
      reliability <- 0.75
    } else if (best_score > 0.4) {
      reliability <- 0.65
    } else {
      reliability <- 0.5
    }
    
    return(list(
      threshold = best_threshold,
      method = best_method,
      reliability = reliability,
      peaks = best_peaks  # Return peak positions for visualisation
    ))
  }
  
  # Try curve inflection method with consistent error handling
  inflection_result <- try({
    # Check for shoulder pattern in the distribution
    inflection_info <- .detect_curve_inflection(analysis_values, verbose)
    
    if (!inflection_info$has_inflection) {
      stop("No inflection point detected")
    }
    
    if (verbose) message("Detected curve inflection")
    
    # Calculate reliability based on distribution type
    min_reliability <- 0.7  # Default threshold
    
    # Check reliability threshold
    if (inflection_info$reliability <= min_reliability) {
      stop("Inflection detection reliability is low (", 
           round(inflection_info$reliability, 2), ")")
    }
    
    # Return the result
    list(
      threshold = inflection_info$inflection_position,
      method = "curve.inflection",
      reliability = inflection_info$reliability,
      peaks = NULL  # No peaks for inflection method
    )
  }, silent = TRUE)
  
  # Return inflection result if successful
  if (!inherits(inflection_result, "try-error") && is.list(inflection_result)) {
    return(inflection_result)
  } else if (verbose && inherits(inflection_result, "try-error")) {
    message("Curve inflection method failed: ", as.character(attributes(inflection_result)$condition))
  }
  
  
  # Existing k-means fallback logic
  kmeans_result <- try({
    km <- kmeans(analysis_values, centers = 2, nstart = 25)
    centres <- km$centers
    
    # Ensure centres are ordered
    if (centres[1] > centres[2]) {
      clusters <- km$cluster
      clusters[clusters == 1] <- 3
      clusters[clusters == 2] <- 1
      clusters[clusters == 3] <- 2
      centres <- centres[c(2, 1)]
    } else {
      clusters <- km$cluster
    }
    
    # Get values for each cluster
    c1_vals <- analysis_values[clusters == 1]
    c2_vals <- analysis_values[clusters == 2]
    
    # Robust threshold calculation using quantiles
    c1_upper <- quantile(c1_vals, 0.95)
    c2_lower <- quantile(c2_vals, 0.05)
    
    threshold <- (c1_upper + c2_lower) / 2
    
    # Calculate separation quality
    c1_mean <- mean(c1_vals)
    c2_mean <- mean(c2_vals)
    pooled_sd <- sqrt((var(c1_vals) + var(c2_vals)) / 2)
    separation <- abs(c2_mean - c1_mean) / pooled_sd
    
    # Create simulated peaks for visualisation (at cluster centres)
    kde <- density(analysis_values, adjust = 1.0)
    kmeans_peaks <- c()
    
    # Find the closest points on density curve to the cluster centres
    for (centre in c(c1_mean, c2_mean)) {
      closest_idx <- which.min(abs(kde$x - centre))
      kmeans_peaks <- c(kmeans_peaks, closest_idx)
    }
    
    list(
      threshold = threshold,
      method = paste0("multimodal.kmeans.",
                      ifelse(separation > 5, "excellent",
                             ifelse(separation > 3, "good", "moderate"))),
      reliability = min(0.8, 0.5 + (separation / 10)),
      peaks = kmeans_peaks  # Store simulated peaks for visualisation
    )
  }, silent = TRUE)
  
  if (!inherits(kmeans_result, "try-error")) {
    return(kmeans_result)
  }
  
  # Ultimate fallback
  return(list(
    threshold = quantile(analysis_values, 0.66),
    method = "multimodal.fallback.q66",
    reliability = 0.5,
    peaks = NULL  # No peaks for fallback method
  ))
}

#' Extreme skewness threshold detection method
#'
#' Specialised implementation for extremely skewed distributions
#' that focuses on finding a reasonable threshold in the body of
#' the distribution rather than the extreme tail.
#'
#' @param analysis_values Processed values to analyse
#' @param preprocessing Preprocessing result
#' @param verbose Whether to print detailed messages
#' @return List containing threshold and method details
#' @keywords internal
.extreme_skew_threshold_method <- function(analysis_values, preprocessing, verbose = FALSE) {
  # Default safe range
  safe_range <- c(0.4, 0.7)
  
  
  # Get key percentiles within safe range
  q_lower <- quantile(analysis_values, safe_range[1])
  q_upper <- quantile(analysis_values, safe_range[2])
  
  # Calculate kernel density in safe range with extra smoothing for stability
  kde <- try({
    density(analysis_values, adjust = 1.5)  # Extra smoothing
  }, silent = TRUE)
  
  if (!inherits(kde, "try-error")) {
    # Extract points in safe range
    in_range_idx <- which(kde$x >= q_lower & kde$x <= q_upper)
    
    if (length(in_range_idx) >= 5) {
      x_in_range <- kde$x[in_range_idx]
      y_in_range <- kde$y[in_range_idx]
      
      # Find inflection point (where second derivative changes sign)
      # Calculate approximation to second derivative
      if (length(y_in_range) >= 5) {
        # First derivative
        dy <- diff(y_in_range) / diff(x_in_range)
        
        # Second derivative
        d2y <- diff(dy) / diff(x_in_range[-length(x_in_range)])
        
        # Find where second derivative crosses zero (inflection point)
        # This often represents a good threshold in skewed distributions
        inflection_points <- which(diff(sign(d2y)) != 0)
        
        if (length(inflection_points) > 0) {
          # Take the most significant inflection point (largest change)
          idx <- inflection_points[which.max(abs(d2y[inflection_points]))]
          
          if (idx < length(x_in_range) - 1) {
            # Use the inflection point as threshold
            threshold <- x_in_range[idx + 1]  # +1 because of diff operations
            
            # Check if this point makes sense (should be in safe range)
            if (threshold >= q_lower && threshold <= q_upper) {
              return(list(
                threshold = threshold,
                method = "extreme_skew.inflection_point",
                reliability = 0.85
              ))
            }
          }
        }
      }
      
      # If inflection point method fails, find the "elbow" in the density curve
      # This is often a good threshold in skewed distributions
      if (preprocessing$distribution_type == "extreme_right_skewed") {
        # For right skew, find point of maximum curvature from left side
        # This is approximated as maximum second derivative
        if (length(y_in_range) >= 5) {
          # Find where density curve starts to flatten
          idx <- which.max(abs(d2y))
          if (idx < length(x_in_range) - 1) {
            threshold <- x_in_range[idx + 1]
            
            # Check if this point makes sense
            if (threshold >= q_lower && threshold <= q_upper) {
              return(list(
                threshold = threshold,
                method = "extreme_skew.max_curvature",
                reliability = 0.8
              ))
            }
          }
        }
      } else {
        # For left skew, similar approach but from right side
        if (length(y_in_range) >= 5) {
          idx <- which.max(abs(d2y))
          if (idx < length(x_in_range) - 1) {
            threshold <- x_in_range[idx + 1]
            
            # Check if this point makes sense
            if (threshold >= q_lower && threshold <= q_upper) {
              return(list(
                threshold = threshold,
                method = "extreme_skew.max_curvature",
                reliability = 0.8
              ))
            }
          }
        }
      }
    }
  }
  
  # If all else fails, use a reasonable percentile
  if (preprocessing$distribution_type == "extreme_right_skewed") {
    threshold <- quantile(analysis_values, 0.45)  # Lower percentile for right skew
    method <- "extreme_skew.percentile_45"
  } else {
    threshold <- quantile(analysis_values, 0.55)  # Higher percentile for left skew
    method <- "extreme_skew.percentile_55"
  }
  
  return(list(
    threshold = threshold,
    method = method,
    reliability = 0.75
  ))
}

#' Both tails-inflated threshold detection method
#'
#' Specialised implementation for both tails-inflated distributions that
#' properly handles the non-tails component to find an appropriate threshold.
#'
#' @param analysis_values Processed values to analyse
#' @param preprocessing Preprocessing result
#' @param verbose Whether to print detailed messages
#' @return List containing threshold and method details
#' @keywords internal
.both_tails_inflated_threshold_method <- function(analysis_values, preprocessing, verbose = FALSE) {
  # For both tails-inflated data, we've already extracted the non-tails component
  # in preprocessing, so analysis_values should contain only the non-tails part
  
  # Check if the non-tails component shows multimodality
  if ("is_multimodal" %in% names(preprocessing) &&
      preprocessing$is_multimodal) {
    
    # Use multimodal method on the non-tails component
    multimodal_result <- .multimodal_threshold_method(analysis_values, preprocessing, verbose)
    
    # Return result but mark as both tails-inflated multimodal
    multimodal_result$method <- paste0("both_tails_inflated.", multimodal_result$method)
    
    return(multimodal_result)
    
  } else {
    # Use unimodal method on the non-tails component
    unimodal_result <- .unimodal_threshold_method(analysis_values, preprocessing, verbose)
    
    # Return result but mark as both tails-inflated unimodal
    unimodal_result$method <- paste0("both_tails_inflated.", unimodal_result$method)
    
    return(unimodal_result)
  }
  
  # Fallback method using percentiles when other methods fail
  # Higher percentile for zero-inflated data to focus on non-zero distribution
  threshold <- quantile(analysis_values, 0.75)
  
  return(list(
    threshold = threshold,
    method = "both_tails_inflated.percentile_75",
    reliability = 0.7
  ))
}

#' Zero-inflated threshold detection method
#'
#' Specialised implementation for zero-inflated distributions that
#' properly handles the non-zero component to find an appropriate threshold.
#'
#' @param analysis_values Processed values to analyse
#' @param preprocessing Preprocessing result
#' @param verbose Whether to print detailed messages
#' @return List containing threshold and method details
#' @keywords internal
.zero_inflated_threshold_method <- function(analysis_values, preprocessing, verbose = FALSE) {
  # For zero-inflated data, we've already extracted the non-zero component
  # in preprocessing, so analysis_values should contain only the non-zero part
  
  # Check if the non-zero component shows multimodality
  if ("is_multimodal" %in% names(preprocessing) &&
      preprocessing$is_multimodal) {
    
    # Use multimodal method on the non-zero component
    multimodal_result <- .multimodal_threshold_method(analysis_values, preprocessing, verbose)
    
    # Return result but mark as zero-inflated multimodal
    multimodal_result$method <- paste0("zero_inflated.", multimodal_result$method)
    
    return(multimodal_result)
    
  } else {
    # Use unimodal method on the non-zero component
    unimodal_result <- .unimodal_threshold_method(analysis_values, preprocessing, verbose)
    
    # Return result but mark as zero-inflated unimodal
    unimodal_result$method <- paste0("zero_inflated.", unimodal_result$method)
    
    return(unimodal_result)
  }
  
  # Fallback method using percentiles when other methods fail
  # Higher percentile for zero-inflated data to focus on non-zero distribution
  threshold <- quantile(analysis_values, 0.75)
  
  return(list(
    threshold = threshold,
    method = "zero_inflated.percentile_75",
    reliability = 0.7
  ))
}

#' Right tail-inflated threshold detection method
#'
#' Specialised implementation for right tail-inflated distributions that
#' properly handles the non-right tail component to find an appropriate threshold.
#'
#' @param analysis_values Processed values to analyse
#' @param preprocessing Preprocessing result
#' @param verbose Whether to print detailed messages
#' @return List containing threshold and method details
#' @keywords internal
.right_tail_inflated_threshold_method <- function(analysis_values, preprocessing, verbose = FALSE) {
  # For right tail-inflated data, we've already extracted the non-right tail component
  # in preprocessing, so analysis_values should contain only the non-right tail part
  
  # Check if the non-zero component shows multimodality
  if ("is_multimodal" %in% names(preprocessing) &&
      preprocessing$is_multimodal) {
    
    # Use multimodal method on the non-zero component
    multimodal_result <- .multimodal_threshold_method(analysis_values, preprocessing, verbose)
    
    # Return result but mark as right tail-inflated multimodal
    multimodal_result$method <- paste0("right_tail_inflated.", multimodal_result$method)
    
    return(multimodal_result)
    
  } else {
    # Use multimodal method on the non-zero component
    unimodal_result <- .unimodal_threshold_method(analysis_values, preprocessing, verbose)
    
    # Return result but mark as right tail-inflated multimodal
    unimodal_result$method <- paste0("right_tail_inflated.", unimodal_result$method)
    
    return(unimodal_result)
  }
  

  # Fallback method using percentiles when other methods fail
  # Higher percentile for zero-inflated data to focus on non-zero distribution
  threshold <- quantile(analysis_values, 0.75)
  
  return(list(
    threshold = threshold,
    method = "right_tail_inflated.percentile_75",
    reliability = 0.7
  ))
}

#' Unimodal distribution threshold detection method
#'
#' Implementation for unimodal and mildly skewed distributions.
#'
#' @param analysis_values Processed values to analyse
#' @param preprocessing Preprocessing result
#' @param verbose Whether to print detailed messages
#' @return List containing threshold and method details
#' @keywords internal
.unimodal_threshold_method <- function(analysis_values, preprocessing, verbose = FALSE) {
  # Calculate basic statistics for determining the optimal cutoff
  mean_val <- mean(analysis_values)
  median_val <- median(analysis_values)
  sd_val <- sd(analysis_values)
  q25 <- quantile(analysis_values, 0.25)
  q75 <- quantile(analysis_values, 0.75)
  iqr <- q75 - q25
  q10 <- quantile(analysis_values, 0.10)
  q90 <- quantile(analysis_values, 0.90)
  
  # Get skewness from preprocessing
  skewness <- moments::skewness(analysis_values)
  if (is.null(skewness) || !is.finite(skewness)) {
    skewness <- preprocessing$skewness
  }
  
  # 1. Optimal cutoff detection based on density curve analysis
  optimal_density_point <- try({
    # Calculate density with moderate smoothing
    kde <- density(analysis_values, adjust = 1.0)
    x_values <- kde$x
    y_values <- kde$y
    
    # Determine distribution shape pattern
    peak_idx <- which.max(y_values)
    peak_x <- x_values[peak_idx]
    peak_y <- y_values[peak_idx]
    
    # Calculate position-normalised distance from median
    median_position <- mean(analysis_values <= median_val)
    peak_position <- mean(analysis_values <= peak_x)
    position_shift <- abs(peak_position - 0.5)
    
    # Calculate density derivatives
    first_derivative <- diff(y_values) / diff(x_values)
    second_derivative <- diff(first_derivative) / diff(x_values[-length(x_values)])
    
    if (abs(skewness) < 0.5) {
      # Nearly symmetric distribution - use median + 0.5 * SD
      threshold <- median_val + 0.5 * sd_val
      
      # Check if this point makes sense within distribution
      threshold_quantile <- mean(analysis_values <= threshold)
      if (threshold_quantile < 0.6 || threshold_quantile > 0.8) {
        # Fallback to 75th percentile if result seems off
        threshold <- q75
      }
      
      return(list(
        threshold = threshold,
        method = "unimodal.balanced_threshold",
        reliability = 0.85
      ))
    } else if (abs(skewness) > 0.5) {
      
      # Try curve inflection method with consistent error handling
      inflection_result <- try({
        # Check for shoulder pattern in the distribution
        inflection_info <- .detect_curve_inflection(analysis_values, verbose)
        
        if (!inflection_info$has_inflection) {
          stop("No inflection point detected")
        }
        
        if (verbose) message("Detected curve inflection")
        
        # Calculate reliability based on distribution type
        min_reliability <- 0.7  # Default threshold
        
        # Check reliability threshold
        if (inflection_info$reliability <= min_reliability) {
          stop("Inflection detection reliability is low (", 
               round(inflection_info$reliability, 2), ")")
        }
        
        # Return the result
        list(
          threshold = inflection_info$inflection_position,
          method = "curve.inflection",
          reliability = inflection_info$reliability,
          peaks = NULL  # No peaks for inflection method
        )
      }, silent = TRUE)
      
      # Return inflection result if successful
      if (!inherits(inflection_result, "try-error") && is.list(inflection_result)) {
        return(inflection_result)
      } else if (verbose && inherits(inflection_result, "try-error")) {
        message("Curve inflection method failed: ", as.character(attributes(inflection_result)$condition))
      } 
    }
    
    if (skewness > 0) {  
      # Right-skewed distribution - find curve "elbow" on right side
      # Start from peak and move right
      right_segment <- peak_idx:length(x_values)
      right_x <- x_values[right_segment]
      right_y <- y_values[right_segment]
      
      # Find region where curve starts flattening
      # Calculate curvature in right segment
      if (length(right_segment) > 20) {
        right_first_deriv <- diff(right_y) / diff(right_x)
        right_second_deriv <- diff(right_first_deriv) / diff(right_x[-length(right_x)])
        
        # Find points where curvature changes significantly
        curve_change_points <- which(abs(right_second_deriv) > 0.1 * max(abs(right_second_deriv)))
        
        if (length(curve_change_points) > 0) {
          # Filter to points in reasonable range (between median and 90th percentile)
          median_idx_in_segment <- which.min(abs(right_x - median_val))
          q90_idx_in_segment <- which.min(abs(right_x - q90))
          
          valid_points <- curve_change_points[curve_change_points > median_idx_in_segment & 
                                                curve_change_points < q90_idx_in_segment]
          
          if (length(valid_points) > 0) {
            # Find point of maximum curve change in valid range
            best_idx <- valid_points[which.max(abs(right_second_deriv[valid_points]))]
            threshold <- right_x[best_idx]
            
            return(list(
              threshold = threshold,
              method = "unimodal.right_skew_elbow",
              reliability = 0.75
            ))
          }
        }
      }
      
      # Fallback for right-skewed: use adaptive percentile
      # More skewed -> lower percentile to avoid extreme tail
      percentile <- max(0.6, min(0.8, 0.8 - skewness * 0.05))
      threshold <- quantile(analysis_values, percentile)
      
      return(list(
        threshold = threshold,
        method = paste0("unimodal.right_skew_p", round(percentile * 100)),
        reliability = 0.65
      ))
    } else {
      # Left-skewed distribution - find curve "elbow" on left side
      # Similar approach but reversed
      left_segment <- 1:peak_idx
      left_x <- x_values[left_segment]
      left_y <- y_values[left_segment]
      
      if (length(left_segment) > 20) {
        left_first_deriv <- diff(left_y) / diff(left_x)
        left_second_deriv <- diff(left_first_deriv) / diff(left_x[-length(left_x)])
        
        # Find points where curvature changes significantly
        curve_change_points <- which(abs(left_second_deriv) > 0.1 * max(abs(left_second_deriv)))
        
        if (length(curve_change_points) > 0) {
          # Filter to points in reasonable range (between 10th percentile and median)
          q10_idx_in_segment <- which.min(abs(left_x - q10))
          median_idx_in_segment <- which.min(abs(left_x - median_val))
          
          valid_points <- curve_change_points[curve_change_points > q10_idx_in_segment & 
                                                curve_change_points < median_idx_in_segment]
          
          if (length(valid_points) > 0) {
            # Find point of maximum curve change in valid range
            best_idx <- valid_points[which.max(abs(left_second_deriv[valid_points]))]
            threshold <- left_x[best_idx]
            
            return(list(
              threshold = threshold,
              method = "unimodal.left_skew_elbow",
              reliability = 0.75
            ))
          }
        }
      }
      
      # Fallback for left-skewed: use adaptive percentile
      # More skewed -> higher percentile to avoid extreme tail
      percentile <- min(0.8, max(0.6, 0.6 - skewness * 0.05))
      threshold <- quantile(analysis_values, percentile)
      
      return(list(
        threshold = threshold,
        method = paste0("unimodal.left_skew_p", round(percentile * 100)),
        reliability = 0.65
      ))
    }
  }, silent = TRUE)
  
  if (!inherits(optimal_density_point, "try-error") && is.list(optimal_density_point)) {
    return(optimal_density_point)
  }
  
  # Fallback: Try mixture modelling with optimal component selection
  try_mixture <- try({
    if (requireNamespace("mclust", quietly = TRUE)) {
      # Fit mixture models with 1-3 components
      model <- suppressWarnings(mclust::Mclust(
        analysis_values,
        G = 1:3,  # Try 1-3 components
        verbose = FALSE
      ))
      
      if (!is.null(model) && !is.null(model$parameters)) {
        if (model$G == 1) {
          # Single component - use mean + 0.5 * standard deviation
          mean_val <- model$parameters$mean
          
          if (is.numeric(model$parameters$variance)) {
            sd_val <- sqrt(model$parameters$variance)
          } else if ("sigmasq" %in% names(model$parameters$variance)) {
            sd_val <- sqrt(model$parameters$variance$sigmasq)
          } else {
            sd_val <- sqrt(model$parameters$variance$sigma[1, 1, 1])
          }
          
          threshold <- mean_val + 0.5 * sd_val
          
          # Ensure threshold is in reasonable range (60-80th percentile)
          threshold_quantile <- mean(analysis_values <= threshold)
          if (threshold_quantile < 0.6 || threshold_quantile > 0.8) {
            threshold <- quantile(analysis_values, 0.75)
          }
          
          return(list(
            threshold = threshold,
            method = "unimodal.single_component",
            reliability = 0.7
          ))
        } else {
          # Extract model parameters
          n_components <- model$G
          means <- model$parameters$mean
          weights <- model$parameters$pro
          
          # Extract standard deviations
          if ("sigmasq" %in% names(model$parameters$variance)) {
            sds <- sqrt(model$parameters$variance$sigmasq)
          } else if (is.numeric(model$parameters$variance)) {
            sds <- rep(sqrt(model$parameters$variance), length(means))
          } else {
            # Extract from variance matrices
            sds <- numeric(length(means))
            for (i in 1:length(means)) {
              sds[i] <- sqrt(model$parameters$variance$sigma[1, 1, i])
            }
          }
          
          # Consider all possible component pairs
          threshold_candidates <- list()
          
          for (i in 1:(n_components-1)) {
            for (j in (i+1):n_components) {
              m1 <- means[i]
              m2 <- means[j]
              s1 <- sds[i]
              s2 <- sds[j]
              w1 <- weights[i]
              w2 <- weights[j]
              
              # Calculate weighted mean for reference
              weighted_mean <- (m1 * w1 + m2 * w2) / (w1 + w2)
              
              # Find intersection of density curves
              if (abs(s1 - s2) < 1e-10) {
                # Equal variances - weighted midpoint
                threshold <- ((m1 * w2) + (m2 * w1)) / (w1 + w2)
              } else {
                # Different variances - solve quadratic equation with weight correction
                # Equation accounts for mixing weights
                c_weight_term <- log((s2 * w1) / (s1 * w2))  # Weight correction term
                a <- 1/(2*s1^2) - 1/(2*s2^2)
                b <- m2/s2^2 - m1/s1^2
                c <- m1^2/(2*s1^2) - m2^2/(2*s2^2) - c_weight_term
                
                discriminant <- b^2 - 4*a*c
                if (discriminant >= 0) {
                  # Real solutions exist
                  x1 <- (-b + sqrt(discriminant)) / (2*a)
                  x2 <- (-b - sqrt(discriminant)) / (2*a)
                  
                  # Choose solution between means or closest to weighted mean
                  thresholds <- c()
                  if ((x1 >= m1 && x1 <= m2) || (x1 <= m1 && x1 >= m2)) {
                    thresholds <- c(thresholds, x1)
                  }
                  if ((x2 >= m1 && x2 <= m2) || (x2 <= m1 && x2 >= m2)) {
                    thresholds <- c(thresholds, x2)
                  }
                  
                  if (length(thresholds) > 0) {
                    # Choose threshold closest to weighted mean
                    threshold <- thresholds[which.min(abs(thresholds - weighted_mean))]
                  } else {
                    # No valid intersection found - use weighted mean
                    threshold <- weighted_mean
                  }
                } else {
                  # No real solution - use weighted mean
                  threshold <- weighted_mean
                }
              }
              
              # Calculate separation quality for reliability assessment
              separation <- abs(m2 - m1) / sqrt((s1^2 + s2^2) / 2)
              
              # Calculate combined weight of these components
              combined_weight <- w1 + w2
              
              # Quality score combining separation and component weights
              quality_score <- separation * combined_weight
              
              # Store this candidate
              threshold_candidates[[length(threshold_candidates) + 1]] <- list(
                threshold = threshold,
                components = c(i, j),
                separation = separation,
                component_means = c(m1, m2),
                component_sds = c(s1, s2),
                component_weights = c(w1, w2),
                quality_score = quality_score,
                weighted_mean = weighted_mean
              )
            }
          }
          
          # Filter for candidates with reasonable separation
          valid_candidates <- threshold_candidates[sapply(threshold_candidates, 
                                                          function(x) x$separation > 1.2)]
          
          if (length(valid_candidates) > 0) {
            # Sort by quality score
            quality_scores <- sapply(valid_candidates, function(x) x$quality_score)
            best_idx <- which.max(quality_scores)
            best_candidate <- valid_candidates[[best_idx]]
            
            threshold <- best_candidate$threshold
            
            # Ensure threshold is reasonable by checking against data distribution
            data_median <- median(analysis_values)
            data_iqr <- stats::IQR(analysis_values, na.rm = TRUE)
            data_range <- max(analysis_values) - min(analysis_values)
            
            # If threshold is excessively far from median, make adjustment
            if (abs(threshold - data_median) > min(data_iqr, data_range/4)) {
              # Calculate adjustment toward median (compromise between model and data centre)
              adjustment_factor <- 0.4  # 40% adjustment toward median
              original_threshold <- threshold
              threshold <- threshold * (1 - adjustment_factor) + data_median * adjustment_factor
              
              if (verbose) {
                message("Adjusting threshold from ", round(original_threshold, 3), 
                        " toward median (", round(data_median, 3), 
                        ") to ", round(threshold, 3))
              }
              method_suffix <- ".adjusted"
            } else {
              method_suffix <- ""
            }
            
            return(list(
              threshold = threshold,
              method = paste0("unimodal.mixture_model", method_suffix),
              reliability = min(0.9, 0.7 + best_candidate$separation * 0.1),
              separation = best_candidate$separation,
              selected_components = best_candidate$components
            ))
          }
        }
      }
    }
  }, silent = TRUE)
  
  # 3. Try Otsu's method (robust for more complex distributions)
  otsu_result <- try({
    # Determine appropriate number of bins
    n_bins <- min(256, max(20, ceiling(length(analysis_values) / 10)))
    
    # Create histogram
    hist_data <- hist(analysis_values, breaks = n_bins, plot = FALSE)
    counts <- hist_data$counts
    bin_centers <- hist_data$mids
    
    # Calculate cumulative distributions
    p <- counts / sum(counts)
    omega <- cumsum(p)
    mu <- cumsum(p * bin_centers)
    mu_t <- mu[length(mu)]
    
    # Calculate between-class variance
    sigma_b_squared <- numeric(length(omega))
    
    # Avoid division by zero
    valid_indices <- omega > 0 & omega < 1
    sigma_b_squared[valid_indices] <- ((mu_t * omega[valid_indices] - mu[valid_indices])^2) /
      (omega[valid_indices] * (1 - omega[valid_indices]))
    
    # Find optimal threshold
    max_idx <- which.max(sigma_b_squared)
    
    if (length(max_idx) > 0 && sigma_b_squared[max_idx] > 0) {
      threshold <- bin_centers[max_idx]
      
      # Check threshold quality
      below_mean <- mean(analysis_values[analysis_values <= threshold])
      above_mean <- mean(analysis_values[analysis_values > threshold])
      separation <- (above_mean - below_mean) / sd(analysis_values)
      
      # Check if threshold is in a reasonable range
      threshold_quantile <- mean(analysis_values <= threshold)
      if (threshold_quantile < 0.25 || threshold_quantile > 0.85) {
        # If outside reasonable range, adjust toward median
        threshold <- 0.6 * threshold + 0.4 * median_val
      }
      
      return(list(
        threshold = threshold,
        method = "unimodal.otsu",
        reliability = min(0.85, 0.6 + separation * 0.1)
      ))
    }
  }, silent = TRUE)
  
  if (!inherits(otsu_result, "try-error") && is.list(otsu_result)) {
    return(otsu_result)
  }
  
  # 4. Fallback to adaptive percentile method based on skewness
  if (skewness > 0.5) {
    # For moderately right-skewed, use lower percentile (65th)
    percentile <- 0.65
  } else if (skewness < -0.5) {
    # For moderately left-skewed, use higher percentile (80th)
    percentile <- 0.8
  } else if (abs(skewness) < 0.5) {
    # For nearly symmetric, use 75th percentile
    percentile <- 0.75
  } else {
    percentile <- 0.75
  }
  
  threshold <- quantile(analysis_values, percentile)
  
  return(list(
    threshold = threshold,
    method = paste0("unimodal.percentile_", round(percentile * 100)),
    reliability = 0.7
  ))
}

#' Advanced Statistical Threshold Detection System
#'
#' This system implements a comprehensive approach to threshold detection with
#' consistent preprocessing and optimised handling for challenging distributions
#' including multimodality, zero-inflation, and extreme skewness.
#'
#' @param values Vector of numeric values to analyse
#' @param feature_name Name of the feature being analysed
#' @param min_threshold Minimum acceptable threshold value
#' @param max_threshold Maximum acceptable threshold value (NULL = use max(values))
#' @param default_threshold Default threshold to use if all methods fail
#' @param n_bootstrap Number of bootstrap iterations for stability assessment
#' @param cross_validate Whether to use cross-validation for reliability assessment
#' @param handle_outliers Whether to temporarily exclude extreme outliers
#' @param verbose Whether to print detailed progress messages
#' @return A list containing threshold, reliability metrics, and method information
#' @keywords internal
.detect_optimal_threshold <- function(values,
                                      feature_name = "feature",
                                      min_threshold = 0.1,
                                      max_threshold = NULL,
                                      default_threshold = 0.6,
                                      n_bootstrap = 100,
                                      cross_validate = TRUE,
                                      handle_outliers = TRUE,
                                      verbose = FALSE) {
  
  # Remove missing values and verify sufficient data
  values <- values[is.finite(values)]
  values <- values[!is.na(values)]
  
  if (length(values) < 20) {
    if (verbose) message("Insufficient data for reliable threshold determination.")
    return(list(
      threshold = default_threshold,
      reliability = 0,
      method = "default.insufficient_data",
      explanation = "Insufficient data for reliable analysis"
    ))
  }
  
  # Set max threshold if not provided
  if (is.null(max_threshold)) {
    max_threshold <- max(values, na.rm = TRUE)
  }
  
  # Calculate basic statistics for reference
  stats <- list(
    n = length(values),
    min = min(values),
    max = max(values),
    mean = mean(values),
    median = median(values),
    q25 = quantile(values, 0.25),
    q75 = quantile(values, 0.75),
    sd = sd(values)
  )
  
  stats$iqr <- stats$q75 - stats$q25
  stats$range <- stats$max - stats$min
  
  #-----------------------------------------------------
  # 1. PREPROCESSING AND DISTRIBUTION CLASSIFICATION
  #-----------------------------------------------------
  
  # Apply enhanced preprocessing to handle all distribution types
  preprocessing <- .preprocess_data(values, handle_outliers, feature_name, verbose)
  processed_values <- preprocessing$processed_values
  
  # For threshold calculations, use the processed values
  analysis_values <- processed_values
  
  # Update stats for processed values if they differ from original
  if (!identical(processed_values, values)) {
    proc_stats <- list(
      n = length(processed_values),
      min = min(processed_values),
      max = max(processed_values),
      mean = mean(processed_values),
      median = median(processed_values),
      q25 = quantile(processed_values, 0.25),
      q75 = quantile(processed_values, 0.75),
      sd = sd(processed_values)
    )
    proc_stats$iqr <- proc_stats$q75 - proc_stats$q25
    proc_stats$range <- proc_stats$max - proc_stats$min
  } else {
    proc_stats <- stats
  }
  
  # Determine if this is a difference/variability feature
  is_diff_var <- grepl("difference|variability", feature_name, ignore.case = TRUE)
  
  # Log distribution characteristics if verbose
  if (verbose) {
    message("\nDistribution characteristics for ", feature_name, ":")
    message("- Sample size: ", proc_stats$n)
    message("- Range: [", round(proc_stats$min, 3), ", ", round(proc_stats$max, 3), "]")
    message("- Mean: ", round(proc_stats$mean, 3), ", Median: ", round(proc_stats$median, 3))
    message("- Preprocessing: ", preprocessing$explanation)
    message("- Distribution type: ", preprocessing$distribution_type)
    
    if ("is_zero_inflated" %in% names(preprocessing) && preprocessing$is_zero_inflated) {
      message("- Zero-inflated: Yes")
    }
    
    if ("is_right_tail_inflated" %in% names(preprocessing) && preprocessing$is_right_tail_inflated) {
      message("- Right tail-inflated: Yes")
    }
    
    if ("is_extremely_skewed" %in% names(preprocessing) && preprocessing$is_extremely_skewed) {
      message("- Extreme skewness detected (", round(preprocessing$skewness, 2), ")")
    } else if ("skewness" %in% names(preprocessing) && abs(preprocessing$skewness) > 0.5) {
      message("- Skewed: ", ifelse(preprocessing$skewness > 0, "right", "left"),
              " (", round(preprocessing$skewness, 2), ")")
    }
    
    if ("is_multimodal" %in% names(preprocessing) && preprocessing$is_multimodal) {
      message("- Multimodal: Yes")
    }
    
    if (is_diff_var) {
      message("- Feature type: Difference/variability metric")
    }
  }
  
  #-----------------------------------------------------
  # 2. SPECIALISED METHOD SELECTION
  #-----------------------------------------------------
  
  # Initialise storage for results
  single_method_result <- NULL
  
  # Dispatch to specialised method based on distribution type
  if (grepl("extreme", preprocessing$distribution_type)) {
    # For extremely skewed distributions
    single_method_result <- .extreme_skew_threshold_method(analysis_values, preprocessing, verbose)
  } else if (preprocessing$distribution_type == "multimodal") {
    # For multimodal distributions
    single_method_result <- .multimodal_threshold_method(analysis_values, preprocessing, verbose)
  } else if (preprocessing$distribution_type == "both_tails_inflated") {
    # For both tails-inflated distributions
    single_method_result <- .both_tails_inflated_threshold_method(analysis_values, preprocessing, verbose)
  } else if (preprocessing$distribution_type == "zero_inflated") {
    # For zero-inflated distributions
    single_method_result <- .zero_inflated_threshold_method(analysis_values, preprocessing, verbose)
  } else if (preprocessing$distribution_type == "right_tail_inflated") {
    # For right tail-inflated distributions
    single_method_result <- .right_tail_inflated_threshold_method(analysis_values, preprocessing, verbose)
  } else {
    # For unimodal or mildly skewed distributions
    single_method_result <- .unimodal_threshold_method(analysis_values, preprocessing, verbose)
  }
  
  # Store results from specialised method
  thresholds <- list()
  methods <- list()
  reliability_scores <- list()
  
  if (!is.null(single_method_result)) {
    # Convert threshold back to original scale if necessary
    if ("inverse_func" %in% names(preprocessing)) {
      single_method_result$threshold <- preprocessing$inverse_func(single_method_result$threshold)
    }
    
    # Store results
    thresholds$specialised <- single_method_result$threshold
    methods$specialised <- single_method_result$method
    reliability_scores$specialised <- single_method_result$reliability
  }
  
  #-----------------------------------------------------
  # 3. VALIDATION USING BOOTSTRAP (if requested)
  #-----------------------------------------------------
  
  bootstrap_ci <- NULL
  
  if (n_bootstrap > 0 && length(thresholds) > 0) {
    # Initialise bootstrap storage
    bootstrap_thresholds <- numeric(n_bootstrap)
    
    # Run bootstrap iterations
    if (verbose) {
      message(paste("Running", n_bootstrap, "bootstrap iterations to assess threshold stability..."))
    }
    
    # Main bootstrap loop
    for (i in 1:n_bootstrap) {
      # Sample with replacement
      bootstrap_sample <- sample(values, replace = TRUE)
      
      # Apply comparable preprocessing to bootstrap sample
      bootstrap_prep <- .preprocess_data(bootstrap_sample, handle_outliers,
                                         paste0(feature_name, "_bootstrap"), FALSE)
      bootstrap_values <- bootstrap_prep$processed_values
      
      # Apply same type of specialised method
      bootstrap_result <- try({
        if (grepl("extreme", bootstrap_prep$distribution_type)) {
          result <- .extreme_skew_threshold_method(bootstrap_values, bootstrap_prep, FALSE)
        } else if (bootstrap_prep$distribution_type == "multimodal") {
          result <- .multimodal_threshold_method(bootstrap_values, bootstrap_prep, FALSE)
        } else if (bootstrap_prep$distribution_type == "both_tails_inflated") {
          result <- .both_tails_inflated_threshold_method(bootstrap_values, bootstrap_prep, FALSE)
        } else if (bootstrap_prep$distribution_type == "zero_inflated") {
          result <- .zero_inflated_threshold_method(bootstrap_values, bootstrap_prep, FALSE)
        } else if (bootstrap_prep$distribution_type == "right_tail_inflated") {
          result <- .right_tail_inflated_threshold_method(bootstrap_values, bootstrap_prep, FALSE)
        } else {
          result <- .unimodal_threshold_method(bootstrap_values, bootstrap_prep, FALSE)
        }
        
        # Convert back to original scale if needed
        if ("inverse_func" %in% names(bootstrap_prep)) {
          result$threshold <- bootstrap_prep$inverse_func(result$threshold)
        }
        
        bootstrap_thresholds[i] <- result$threshold
      }, silent = TRUE)
      
      if (inherits(bootstrap_result, "try-error")) {
        bootstrap_thresholds[i] <- NA
      }
    }
    
    # Remove NA values
    bootstrap_thresholds <- bootstrap_thresholds[!is.na(bootstrap_thresholds)]
    
    if (length(bootstrap_thresholds) >= n_bootstrap * 0.5) {
      # Calculate bootstrap statistics
      bootstrap_mean <- mean(bootstrap_thresholds)
      bootstrap_sd <- sd(bootstrap_thresholds)
      bootstrap_cv <- bootstrap_sd / stats$range
      bootstrap_ci <- quantile(bootstrap_thresholds, c(0.025, 0.975))
      
      # Update reliability based on bootstrap stability
      # More stable -> higher reliability
      cv_reliability_factor <- 1.0  # Default no change
      
      if (bootstrap_cv < 0.05) {
        # Very stable
        cv_reliability_factor <- 1.1
      } else if (bootstrap_cv < 0.1) {
        # Stable
        cv_reliability_factor <- 1.0
      } else if (bootstrap_cv < 0.15) {
        # Moderately stable
        cv_reliability_factor <- 0.9
      } else if (bootstrap_cv < 0.2) {
        # Somewhat unstable
        cv_reliability_factor <- 0.8
      } else {
        # Highly unstable
        cv_reliability_factor <- 0.6
      }
      
      # Apply reliability adjustment
      reliability_scores$specialised <- reliability_scores$specialised * cv_reliability_factor
      
      # Cap reliability at 1.0
      reliability_scores$specialised <- min(1.0, reliability_scores$specialised)
      
      # Store bootstrap statistics
      thresholds$bootstrap_mean <- bootstrap_mean
      thresholds$bootstrap_sd <- bootstrap_sd
      thresholds$bootstrap_cv <- bootstrap_cv
      thresholds$bootstrap_ci_lower <- bootstrap_ci[1]
      thresholds$bootstrap_ci_upper <- bootstrap_ci[2]
    } else {
      # Not enough valid bootstrap results - reduce reliability
      reliability_scores$specialised <- reliability_scores$specialised * 0.8
    }
  }
  
  #-----------------------------------------------------
  # 4. CROSS-VALIDATION (if enabled)
  #-----------------------------------------------------
  
  if (cross_validate && length(values) >= 100) {
    # Determine appropriate number of folds
    if (length(values) >= 500) {
      k_folds <- 10
    } else if (length(values) >= 200) {
      k_folds <- 5
    } else {
      k_folds <- 3
    }
    
    fold_size <- length(values) %/% k_folds
    
    if (verbose) {
      message(paste("Performing", k_folds, "fold cross-validation for threshold consistency..."))
    }
    
    # Create stratified folds to better preserve distribution
    sorted_indices <- order(values)
    folds <- vector("list", k_folds)
    
    for (i in 1:length(sorted_indices)) {
      fold_idx <- (i %% k_folds) + 1
      folds[[fold_idx]] <- c(folds[[fold_idx]], sorted_indices[i])
    }
    
    # Initialise storage for cross-validation results
    cv_thresholds <- numeric(k_folds)
    
    # Run cross-validation
    for (i in 1:k_folds) {
      # Create training set (all data except current fold)
      train_indices <- unlist(folds[-i])
      train_values <- values[train_indices]
      
      # Skip if insufficient training data
      if (length(train_values) < 20) {
        if (verbose) message("Fold ", i, " has insufficient training data. Skipping.")
        cv_thresholds[i] <- NA
        next
      }
      
      fold_stats <- list(
        n = length(train_values),
        min = min(train_values, na.rm = TRUE),
        max = max(train_values, na.rm = TRUE),
        mean = mean(train_values, na.rm = TRUE),
        median = median(train_values, na.rm = TRUE),
        q10 = quantile(train_values, 0.1, na.rm = TRUE),
        q25 = quantile(train_values, 0.25, na.rm = TRUE),
        q75 = quantile(train_values, 0.75, na.rm = TRUE),
        q90 = quantile(train_values, 0.9, na.rm = TRUE)
      )
      
      fold_stats$iqr <- fold_stats$q75 - fold_stats$q25
      fold_stats$range <- fold_stats$max - fold_stats$min
      
      # Preprocess the training set
      fold_prep <- .preprocess_data(train_values, handle_outliers,
                                    paste0(feature_name, "_cv_fold", i), FALSE)
      fold_values <- fold_prep$processed_values
      
      # Apply specialised method to this fold
      fold_result <- try({
        if (grepl("extreme", fold_prep$distribution_type)) {
          result <- .extreme_skew_threshold_method(fold_values, fold_prep, FALSE)
        } else if (fold_prep$distribution_type == "multimodal") {
          result <- .multimodal_threshold_method(fold_values, fold_prep, FALSE)
        } else if (fold_prep$distribution_type == "both_tails_inflated") {
          result <- .both_tails_inflated_threshold_method(fold_values, fold_prep, FALSE)
        } else if (fold_prep$distribution_type == "zero_inflated") {
          result <- .zero_inflated_threshold_method(fold_values, fold_prep, FALSE)
        } else if (fold_prep$distribution_type == "right_tail_inflated") {
          result <- .right_tail_inflated_threshold_method(fold_values, fold_prep, FALSE)
        } else {
          result <- .unimodal_threshold_method(fold_values, fold_prep, FALSE)
        }
        
        # Convert back to original scale if needed
        if ("inverse_func" %in% names(fold_prep)) {
          result$threshold <- fold_prep$inverse_func(result$threshold)
        }
        
        if (!is.na(result$threshold)) {
          if (preprocessing$distribution_type == "multimodal") {
            if (!is.na(fold_stats$q10) && !is.na(fold_stats$q90) &&
                (result$threshold < fold_stats$q10 || result$threshold > fold_stats$q90)) {
              result$threshold <- (fold_stats$q25 + fold_stats$q75) / 2
            }
          } else if (grepl("extreme_right", preprocessing$distribution_type)) {
            if (!is.na(fold_stats$q75) && result$threshold > fold_stats$q75) {
              result$threshold <- fold_stats$q50
            }
          } else if (grepl("extreme_left", preprocessing$distribution_type)) {
            if (!is.na(fold_stats$q25) && result$threshold < fold_stats$q25) {
              result$threshold <- fold_stats$q50
            }
          }
        } else {
          result$threshold <- fold_stats$median
        }
        
        cv_thresholds[i] <- result$threshold
      }, silent = TRUE)
      
      if (inherits(fold_result, "try-error")) {
        cv_thresholds[i] <- NA
        if (verbose) {
          message("Fold ", i, " failed: ", fold_result[1])
        }
      }
    }
    
    # Remove NA values
    cv_thresholds <- cv_thresholds[!is.na(cv_thresholds)]
    
    if (length(cv_thresholds) >= max(2, k_folds * 0.6)) {
      # Calculate cross-validation statistics
      cv_mean <- mean(cv_thresholds)
      cv_sd <- sd(cv_thresholds)
      
      # Scale by data range for consistent comparison
      if (stats$range > 0) {
        cv_var <- cv_sd / stats$range
      } else {
        cv_var <- cv_sd / (mean(values) + 1e-10)
      }
      
      # Calculate min-max range
      cv_range <- (max(cv_thresholds) - min(cv_thresholds)) / stats$range
      
      # Adjust reliability based on cross-validation consistency
      if (cv_var < 0.05) {
        # Very consistent
        reliability_scores$specialised <- reliability_scores$specialised * 1.1
      } else if (cv_var < 0.1) {
        # Consistent
        reliability_scores$specialised <- reliability_scores$specialised * 1.0
      } else if (cv_var < 0.15) {
        # Moderately consistent
        reliability_scores$specialised <- reliability_scores$specialised * 0.9
      } else if (cv_var < 0.2) {
        # Somewhat inconsistent
        reliability_scores$specialised <- reliability_scores$specialised * 0.8
      } else {
        # Highly inconsistent
        reliability_scores$specialised <- reliability_scores$specialised * 0.6
      }
      
      # Cap reliability at 1.0
      reliability_scores$specialised <- min(1.0, reliability_scores$specialised)
      
      # Store cross-validation statistics
      thresholds$cv_mean <- cv_mean
      thresholds$cv_sd <- cv_sd
      thresholds$cv_var <- cv_var
      thresholds$cv_range <- cv_range
      thresholds$cv_min <- min(cv_thresholds)
      thresholds$cv_max <- max(cv_thresholds)
    } else {
      # Not enough valid cross-validation results - reduce reliability
      reliability_scores$specialised <- reliability_scores$specialised * 0.8
    }
  }
  
  #-----------------------------------------------------
  # 5. FINAL THRESHOLD SELECTION AND VALIDATION
  #-----------------------------------------------------
  
  # Final threshold is the specialised method's threshold
  final_threshold <- thresholds$specialised
  final_reliability <- reliability_scores$specialised
  final_method <- methods$specialised
  
  # Generate explanation
  explanation <- paste0("Using ", methods$specialised, " method for ",
                        preprocessing$distribution_type, " distribution")
  
  #-----------------------------------------------------
  # 6. SANITY CHECK - ensure threshold is reasonable
  #-----------------------------------------------------
  if (!exists("q10", where = stats)) {
    stats$q10 <- tryCatch(quantile(values, 0.1, na.rm = TRUE), error = function(e) NA)
    stats$q90 <- tryCatch(quantile(values, 0.9, na.rm = TRUE), error = function(e) NA)
  }
  
  # Define reasonable limits based on distribution type
  if (grepl("extreme_right", preprocessing$distribution_type)) {
    # For extreme right skew, threshold shouldn't be too high
    if (!is.na(final_threshold) && !is.na(stats$q75) && final_threshold > stats$q75) {
      safe_threshold <- stats$median + 0.5 * stats$iqr
      if (verbose) {
        message("Threshold too high for extreme right skew. Adjusting from ",
                round(final_threshold, 3), " to ", round(safe_threshold, 3))
      }
      final_threshold <- safe_threshold
      explanation <- paste0(explanation, " (adjusted for extreme skew)")
    }
  } else if (grepl("extreme_left", preprocessing$distribution_type)) {
    # For extreme left skew, threshold shouldn't be too low
    if (!is.na(final_threshold) && !is.na(stats$q25) && final_threshold < stats$q25) {
      safe_threshold <- stats$median - 0.5 * stats$iqr
      if (verbose) {
        message("Threshold too low for extreme left skew. Adjusting from ",
                round(final_threshold, 3), " to ", round(safe_threshold, 3))
      }
      final_threshold <- safe_threshold
      explanation <- paste0(explanation, " (adjusted for extreme skew)")
    }
  } else if (preprocessing$distribution_type == "multimodal") {
    low_condition <- !is.na(final_threshold) && !is.na(stats$q10) && final_threshold < stats$q10
    high_condition <- !is.na(final_threshold) && !is.na(stats$q90) && final_threshold > stats$q90
    
    if (low_condition || high_condition) {
      safe_threshold <- (stats$q25 + stats$q75) / 2
      if (verbose) {
        message("Threshold outside reasonable range for multimodal distribution. Adjusting from ",
                round(final_threshold, 3), " to ", round(safe_threshold, 3))
      }
      final_threshold <- safe_threshold
      explanation <- paste0(explanation, " (adjusted to reasonable range)")
    }
  }
  
  if (is.na(final_threshold)) {
    if (verbose) message("Final threshold is NA, using median instead")
    final_threshold <- stats$median
    explanation <- paste0(explanation, " (NA value corrected)")
  }
  
  # Ensure threshold is within allowed bounds
  final_threshold <- max(min_threshold, min(max_threshold, final_threshold))
  
  # Create final result
  final_result <- list(
    threshold = final_threshold,
    reliability = final_reliability,
    method = final_method,
    explanation = explanation,
    preprocessing = preprocessing,
    all_thresholds = thresholds,
    all_methods = methods,
    all_reliabilities = reliability_scores
  )
  
  # Add confidence interval if available
  if (!is.null(bootstrap_ci)) {
    final_result$confidence_interval <- bootstrap_ci
  } else {
    # Create approximate confidence interval
    ci_width <- 0.1 * stats$iqr
    final_result$confidence_interval <- c(final_threshold - ci_width, final_threshold + ci_width)
  }
  
  # Add diagnostic information
  final_result$diagnostics <- list(
    n = stats$n,
    data_min = stats$min,
    data_max = stats$max,
    data_mean = stats$mean,
    data_median = stats$median,
    data_q25 = stats$q25,
    data_q75 = stats$q75,
    data_sd = stats$sd,
    data_iqr = stats$iqr,
    distribution_type = preprocessing$distribution_type,
    is_zero_inflated = "is_zero_inflated" %in% names(preprocessing) && preprocessing$is_zero_inflated,
    is_extremely_skewed = "is_extremely_skewed" %in% names(preprocessing) && preprocessing$is_extremely_skewed,
    is_multimodal = "is_multimodal" %in% names(preprocessing) && preprocessing$is_multimodal,
    preprocessing_method = preprocessing$method,
    preprocessing_explanation = preprocessing$explanation
  )
  
  # Log final decision
  if (verbose) {
    message(paste("Final threshold for", feature_name, ":", round(final_result$threshold, 3)))
    message(final_result$explanation)
    message(paste("Reliability score:", round(final_result$reliability, 2)))
  }
  
  return(final_result)
}

#' Determine optimal thresholds for all metrics
#'
#' @param metrics_df Data frame containing calculated metrics
#' @param default_thresholds List of default thresholds for each metric
#' @param visualise Whether to create visualisations
#' @param min_samples Minimum number of samples required for reliable modelling
#' @param verbose Whether to print detailed progress messages
#' @return List with thresholds, visualisations, and statistics
#' @keywords internal
.determine_optimal_thresholds <- function(metrics_df,
                                          default_thresholds = list(
                                            intra_cellular_isoform_diversity = 0.6,
                                            inter_cellular_isoform_diversity = 0.6,
                                            intra_cell_type_heterogeneity = 0.4,
                                            inter_cell_type_specificity = 0.6,
                                            intra_cell_type_heterogeneity_variability = 0.5,
                                            inter_cell_type_difference_variability = 0.3,
                                            cell_type_coexpression_variability = 0.4
                                          ),
                                          visualise = TRUE,
                                          min_samples = 20,
                                          verbose = TRUE) {
  
  message("Determining optimal classification thresholds for all seven core metrics...")
  
  # Calculate and report NA statistics
  core_metrics <- c(
    "intra_cellular_isoform_diversity",
    "inter_cellular_isoform_diversity",
    "intra_cell_type_heterogeneity",
    "inter_cell_type_specificity",
    "intra_cell_type_heterogeneity_variability",
    "inter_cell_type_difference_variability",
    "cell_type_coexpression_variability"
  )
  
  core_metrics_names <-  c("Intra-cellular Isoform Diversity",
                           "Inter-cellular Isoform Diversity",
                           "Intra-cell-type Heterogeneity",
                           "Inter-cell-type Specificity",
                           "Intra-cell-type Heterogeneity Variability",
                           "Inter-cell-type Difference Variability",
                           "Cell-type-specific Co-expression Variability")
  
  na_counts <- sapply(metrics_df[core_metrics], function(x) sum(is.na(x)))
  na_percentages <- round(100 * na_counts / nrow(metrics_df), 1)
  
  message("\nNA value proportions for core metrics:")
  for(i in seq_along(core_metrics)) {
    metric_name <- core_metrics_names[i]
    message(sprintf("  - %s: %d NA values (%.1f%%)",
                    metric_name, na_counts[i], na_percentages[i]))
  }
  message("\n")
  
  # List to store all thresholds
  thresholds <- list()
  
  # List to store visualisations
  plot_list <- list()
  
  # List to store NA statistics
  na_stats <- list(
    counts = na_counts,
    percentages = na_percentages
  )
  
  # Process each metric
  for (i in seq_along(core_metrics)) {
    metric <- core_metrics[i]
    values <- metrics_df[[metric]]
    metric_name <- core_metrics_names[i]
    
    # Use the advanced threshold detection with enhanced algorithms
    result <- .detect_optimal_threshold(
      values = values,
      feature_name = metric,
      min_threshold = 0.1,
      default_threshold = default_thresholds[[metric]],
      n_bootstrap = 100,
      cross_validate = length(values) >= 100,
      handle_outliers = TRUE,
      verbose = verbose
    )
    
    # Store threshold
    thresholds[[metric]] <- result$threshold
    
    # Create visualisation if requested
    if (visualise) {
      plot_list[[metric]] <- .visualise_threshold_fitting(
        values,
        result$threshold,
        metric_name,
        model = list(
          method = result$method,
          reliability = result$reliability,
          explanation = result$explanation,
          diagnostics = result$diagnostics
        )
      )
    }
  }
  
  # Summary message
  message(paste0(
    "\nOptimal thresholds determined:",
    "\n- Intra-cellular Isoform Diversity: ", round(thresholds$intra_cellular_isoform_diversity, 3),
    "\n- Inter-cellular Isoform Diversity: ", round(thresholds$inter_cellular_isoform_diversity, 3),
    "\n- Intra-cell-type Heterogeneity: ", round(thresholds$intra_cell_type_heterogeneity, 3),
    "\n- Inter-cell-type Specificity: ", round(thresholds$inter_cell_type_specificity, 3),
    "\n- Intra-cell-type Heterogeneity Variability: ", round(thresholds$intra_cell_type_heterogeneity_variability, 3),
    "\n- Inter-cell-type Difference Variability: ", round(thresholds$inter_cell_type_difference_variability, 3),
    "\n- Cell-type-specific Co-expression Variability: ", round(thresholds$cell_type_coexpression_variability, 3)
  ))
  
  # Return thresholds, visualisations, and NA statistics
  return(list(
    thresholds = thresholds,
    plots = plot_list,
    na_stats = na_stats
  ))
}

#' Enhanced visualisation of threshold determination
#'
#' @param values Vector of numeric values
#' @param threshold Determined threshold value
#' @param metric_name Name of the metric (for title)
#' @param model Optional model or diagnostics information
#' @return ggplot2 object with visualisation
#' @import ggplot2
#' @keywords internal
.visualise_threshold_fitting <- function(values, threshold, metric_name, model = NULL) {
  # Remove NA and infinite values
  values <- values[is.finite(values)]
  values <- na.omit(values)
  
  # Check if we have enough data
  if(length(values) < 5) {
    # Return empty plot with message
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = paste("Insufficient data for", metric_name, "visualisation")) +
             ggplot2::theme_void())
  }
  
  # Format metric name for title
 # readable_name <- gsub("_", " ", metric_name)
 # readable_name <- paste0(toupper(substr(readable_name, 1, 1)),
 #                         substr(readable_name, 2, nchar(readable_name)))
  
  # Calculate basic statistics for reference
  stats <- list(
    n = length(values),
    min = min(values),
    max = max(values),
    mean = mean(values),
    median = median(values),
    q25 = quantile(values, 0.25),
    q75 = quantile(values, 0.75),
    sd = sd(values)
  )
  
  stats$iqr <- stats$q75 - stats$q25
  
  # Create base dataframe for histogram
  df <- data.frame(Value = values)
  
  # Determine optimal bin width using Freedman-Diaconis rule
  bw <- 2 * stats$iqr / length(values)^(1/3)
  if (!is.finite(bw) || bw <= 0) {
    # Fallback if FD rule gives bad result
    bw <- (stats$max - stats$min) / min(30, sqrt(length(values)))
  }
  
  # Count optimal number of bins
  data_range <- stats$max - stats$min
  n_bins <- max(10, min(50, ceiling(data_range / bw)))
  
  # Create base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Value)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = n_bins,
      fill = "#BBD4E9",
      color = "#396CA0",  
      alpha = 0.8
    ) +
    ggplot2::geom_density(
      fill = "#396CA0",
      alpha = 0.3,
      color = "#233E5B",
      linewidth = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = threshold,
      color = "#D64550",
      linewidth = 1.2,
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = paste(metric_name),
      x = metric_name,
      y = "Density"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10, lineheight = 1.2),  # Increased line height for subtitle
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "gray92"),
      plot.margin = ggplot2::margin(t = 15, r = 10, b = 15, l = 10)  # Increased top and bottom margins
    )
  
  # Add threshold label with improved visibility
  p <- p + ggplot2::annotate(
    "label",
    x = threshold,
    y = max(density(values)$y) * 1.1,
    label = paste0("Threshold = ", round(threshold, 3)),
    fill = "#FDF6E3",
    color = "#D64550",
    size = 3.5,
    fontface = "bold"
  )
  
  
  # Prepare subtitle elements as separate lines
  dist_type_text <- NULL
  method_text <- NULL
  
  if (is.list(model)) {
    # Extract distribution type information
    if (is.list(model$diagnostics) && !is.null(model$diagnostics$distribution_type)) {
      dist_type <- model$diagnostics$distribution_type
      # Make distribution type more readable
      dist_type <- gsub("_", " ", dist_type)
      dist_type_text <- paste0("Distribution type: ", dist_type)
    }
    
    # Extract method information
    if (!is.null(model$method)) {
      method <- gsub("_", ".", model$method)
      method_text <- paste0("Using ", method, " method")
    }
    
    # Combine subtitle lines with a newline character
    subtitle_lines <- c()
    if (!is.null(dist_type_text)) subtitle_lines <- c(subtitle_lines, dist_type_text)
    if (!is.null(method_text)) subtitle_lines <- c(subtitle_lines, method_text)
    
    if (length(subtitle_lines) > 0) {
      p <- p + ggplot2::labs(subtitle = paste(subtitle_lines, collapse = "\n"))
    }
  }
  
  return(p)
}

#' Create classification labels for each core metric
#'
#' This function classifies genes based on their complexity metrics.
#' It creates meaningful category labels for each dimension of complexity.
#' NA values are preserved and represent biologically meaningful cases
#' such as cell type-specific expression or single isoform expression.
#'
#' @param metrics_df Data frame containing calculated metrics
#' @param thresholds List of threshold values for classification
#' @return Data frame with added classification columns
#' @keywords internal
.classify_genes <- function(metrics_df, thresholds) {
  # 1. Intra-cellular Isoform Diversity Classification
  metrics_df$intra_cellular_isoform_diversity_class <- ifelse(
    is.na(metrics_df$intra_cellular_isoform_diversity), "Unclassified",
    ifelse(metrics_df$intra_cellular_isoform_diversity > thresholds$intra_cellular_isoform_diversity,
           "Strong Isoform Co-expression", "Weak Isoform Co-expression")
  )
  
  # 2. Inter-cellular Isoform Diversity Classification
  metrics_df$inter_cellular_isoform_diversity_class <- ifelse(
    is.na(metrics_df$inter_cellular_isoform_diversity), "Unclassified",
    ifelse(metrics_df$inter_cellular_isoform_diversity > thresholds$inter_cellular_isoform_diversity,
           "High Isoform Diversity", "Low Isoform Diversity")
  )
  
  # 3. Intra-cell-type Heterogeneity Classification
  metrics_df$intra_cell_type_heterogeneity_class <- ifelse(
    is.na(metrics_df$intra_cell_type_heterogeneity), "Unclassified",
    ifelse(metrics_df$intra_cell_type_heterogeneity > thresholds$intra_cell_type_heterogeneity,
           "High Cellular Heterogeneity", "Low Cellular Heterogeneity")
  )
  
  # 4. Inter-cell-type Specificity Classification
  metrics_df$inter_cell_type_specificity_class <- ifelse(
    is.na(metrics_df$inter_cell_type_specificity), "Single-Cell Type Expression",
    ifelse(metrics_df$inter_cell_type_specificity > thresholds$inter_cell_type_specificity,
           "Cell-Type-Specific Isoform Expression", "Cell-Type-Independent Isoform Expression")
  )
  
  # 5. Intra-cell-type Heterogeneity Variability Classification
  metrics_df$intra_cell_type_heterogeneity_variability_class <- ifelse(
    is.na(metrics_df$intra_cell_type_heterogeneity_variability), "Insufficient Cell Type Data",
    ifelse(metrics_df$intra_cell_type_heterogeneity_variability > thresholds$intra_cell_type_heterogeneity_variability,
           "Variable Heterogeneity Across Cell Types", "Consistent Heterogeneity Across Cell Types")
  )
  
  # 6. Inter-cell-type Difference Variability Classification
  metrics_df$inter_cell_type_difference_variability_class <- ifelse(
    is.na(metrics_df$inter_cell_type_difference_variability), "Insufficient Difference Data",
    ifelse(metrics_df$inter_cell_type_difference_variability > thresholds$inter_cell_type_difference_variability,
           "High Cell-Type Distinctions", "Low Cell-Type Distinctions")
  )
  
  # 7. Cell-type-specific Co-expression Variability Classification
  metrics_df$cell_type_coexpression_variability_class <- ifelse(
    is.na(metrics_df$cell_type_coexpression_variability), "Insufficient Data",
    ifelse(metrics_df$cell_type_coexpression_variability > thresholds$cell_type_coexpression_variability,
           "Cell-Type-Adaptive Co-expression", "Cell-Type-Consistent Co-expression")
  )
  
  # Traditional complexity category (based on diversity and specificity)
  metrics_df$complexity_category <- case_when(
    is.na(metrics_df$inter_cellular_isoform_diversity) | is.na(metrics_df$inter_cell_type_specificity) ~
      "Unclassified",
    
    metrics_df$inter_cellular_isoform_diversity > thresholds$inter_cellular_isoform_diversity &
      metrics_df$inter_cell_type_specificity > thresholds$inter_cell_type_specificity ~
      "High Diversity + High Specificity",
    
    metrics_df$inter_cellular_isoform_diversity <= thresholds$inter_cellular_isoform_diversity &
      metrics_df$inter_cell_type_specificity > thresholds$inter_cell_type_specificity ~
      "Low Diversity + High Specificity",
    
    metrics_df$inter_cellular_isoform_diversity > thresholds$inter_cellular_isoform_diversity &
      metrics_df$inter_cell_type_specificity <= thresholds$inter_cell_type_specificity ~
      "High Diversity + Low Specificity",
    
    TRUE ~ "Low Diversity + Low Specificity"
  )
  
  return(metrics_df)
}

#' Process a single gene for complexity metrics
#'
#' Helper function that processes a single gene to calculate all complexity metrics.
#'
#' @param scht_obj Single-Cell Hierarchical Tensor object
#' @param gene Gene name to analyse
#' @return Data frame row with complexity metrics for the gene, or NULL if insufficient data
#' @keywords internal
.process_gene_complexity <- function(scht_obj, gene) {
  #------------ 1. Get Expression Data ------------#
  
  # Get global expression matrix
  iso_mat <- scht_obj$original_results[[gene]]
  
  # Skip genes with insufficient data
  if(nrow(iso_mat) < 2 || ncol(iso_mat) == 0) return(NULL)
  
  #------------ 2. Calculate Core Metrics ------------#
  
  # Calculate intra-cellular isoform diversity (Core Metric 1)
  cell_idi_results <- .calculate_intra_cellular_isoform_diversity(iso_mat)
  intra_cellular_isoform_diversity <- cell_idi_results$intra_cellular_isoform_diversity
  
  # Calculate inter-cellular isoform diversity (Core Metric 2)
  # Using the specified diversity function
  inter_cellular_isoform_diversity <- .calculate_inter_cellular_isoform_diversity(iso_mat)
  
  # Calculate the difference between inter and intra-cellular diversity
  idi_difference <- inter_cellular_isoform_diversity - intra_cellular_isoform_diversity
  
  # Calculate inter-cell-type specificity (Core Metric 3)
  specificity_result <- .calculate_inter_cell_type_specificity(scht_obj, gene)
  inter_cell_type_specificity <- specificity_result$inter_cell_type_specificity
  is_single_cell_type <- specificity_result$single_cell_type
  specific_cell_type <- specificity_result$cell_type_name
  expressed_cell_types <- specificity_result$expressed_cell_types
  
  # Calculate intra-cell-type heterogeneity (Core Metric 4)
  # and its variability (Core Metric 5)
  heterogeneity_result <- .calculate_intra_cell_type_heterogeneity_variability(scht_obj, gene)
  intra_cell_type_heterogeneity <- mean(heterogeneity_result$cell_type_heterogeneity, na.rm = TRUE)
  intra_cell_type_heterogeneity_variability <- heterogeneity_result$heterogeneity_variability
  
  # Calculate inter-cell-type difference variability (Core Metric 6)
  inter_cell_type_difference_variability <- .calculate_inter_cell_type_difference_variability(specificity_result)
  
  # Calculate cell-type-specific co-expression variability (Core Metric 7)
  coexpression_result <- .calculate_cell_type_coexpression_variability(scht_obj, gene)
  cell_type_coexpression_variability <- coexpression_result$coexpression_variability
  
  #------------ 3. Calculate Additional Metrics ------------#
  
  # Additional metrics for comprehensive characterisation
  iso_means <- rowMeans(iso_mat)
  nonzero_indices <- which(iso_means > 0)
  
  # Calculate proportion statistics
  total_expr <- sum(iso_means)
  iso_props <- rep(0, length(iso_means))
  iso_props[nonzero_indices] <- iso_means[nonzero_indices] / total_expr
  
  # Find dominant isoform
  if(length(nonzero_indices) > 0) {
    max_idx <- which.max(iso_props)
    dominant_iso_prop <- iso_props[max_idx]
    dominant_iso_name <- rownames(iso_mat)[max_idx]
  } else {
    dominant_iso_prop <- 0
    dominant_iso_name <- ""
  }
  
  n_isoforms <- nrow(iso_mat)
  n_expressed_isoforms <- length(nonzero_indices)
  
  # Calculate Simpson diversity index (complementary diversity measure)
  simpson_index <- 1 - sum(iso_props^2)
  
  # Calculate evenness (how equally expressed the isoforms are)
  evenness <- ifelse(n_expressed_isoforms > 1,
                     inter_cellular_isoform_diversity / log2(n_expressed_isoforms),
                     0)
  
  # Calculate cell-level metrics
  cell_sums <- colSums(iso_mat)
  cells_expressing <- sum(cell_sums > 0)
  pct_cells_expressing <- cells_expressing / ncol(iso_mat) * 100
  
  # Calculate percentage of cells expressing multiple isoforms
  if(cells_expressing > 0) {
    multi_iso_count <- 0
    for(j in which(cell_sums > 0)) {
      if(sum(iso_mat[, j] > 0) > 1) {
        multi_iso_count <- multi_iso_count + 1
      }
    }
    pct_multi_iso_cells <- multi_iso_count / cells_expressing * 100
  } else {
    pct_multi_iso_cells <- 0
  }
  
  #------------ 4. Compile Results ------------#
  
  gene_results <- data.frame(
    gene = gene,
    
    # Core metrics
    intra_cellular_isoform_diversity = intra_cellular_isoform_diversity,
    inter_cellular_isoform_diversity = inter_cellular_isoform_diversity,
    intra_cell_type_heterogeneity = intra_cell_type_heterogeneity,
    inter_cell_type_specificity = inter_cell_type_specificity,
    intra_cell_type_heterogeneity_variability = intra_cell_type_heterogeneity_variability,
    inter_cell_type_difference_variability = inter_cell_type_difference_variability,
    cell_type_coexpression_variability = cell_type_coexpression_variability,
    
    # Additional characterisation metrics
    idi_difference = idi_difference,
    dominant_iso_prop = dominant_iso_prop,
    dominant_iso_name = as.character(dominant_iso_name),
    n_isoforms = n_isoforms,
    n_expressed_isoforms = n_expressed_isoforms,
    simpson_index = simpson_index,
    evenness = evenness,
    
    # Single-cell_type information
    is_single_cell_type = is_single_cell_type,
    specific_cell_type = ifelse(is.na(specific_cell_type), "", as.character(specific_cell_type)),
    expressed_cell_types_count = length(expressed_cell_types),
    
    # Cell-level metrics
    cells_expressing = cells_expressing,
    pct_cells_expressing = pct_cells_expressing,
    pct_multi_iso_cells = pct_multi_iso_cells,
    
    stringsAsFactors = FALSE
  )
  
  return(gene_results)
}

#' Calculate cell type-specific metrics for a gene
#'
#' This function calculates a comprehensive set of diversity and heterogeneity metrics
#' for each cell type. It enables detailed comparison of isoform usage patterns
#' across different cell populations, revealing tissue-specific regulatory mechanisms.
#'
#' @param scht_obj Single-Cell Hierarchical Tensor object
#' @param gene Gene name to analyse
#' @return List of cell type-specific metrics, or NULL if insufficient data
#' @keywords internal
.calculate_cell_type_metrics <- function(scht_obj, gene) {
  # Initialise result list
  cell_type_metrics <- list()
  
  # Get all cell types
  cell_types <- names(scht_obj$cell_type_matrices)
  
  # Calculate metrics for each cell type
  for(cell_type in cell_types) {
    if(gene %in% names(scht_obj$cell_type_matrices[[cell_type]])) {
      # Obtain the cell type-specific expression matrix
      ct_iso_mat <- scht_obj$cell_type_matrices[[cell_type]][[gene]]
      
      # Skip cell types with insufficient data
      if(nrow(ct_iso_mat) < 2 || ncol(ct_iso_mat) == 0) next
      
      #------------ 1. Core Metric Calculations ------------#
      
      # Calculate cell type-specific intra-cellular isoform diversity
      # Measures how cells within this cell type co-express multiple isoforms
      ct_intra_results <- .calculate_intra_cellular_isoform_diversity(ct_iso_mat)
      
      # Calculate cell type-specific inter-cellular isoform diversity
      # Measures the overall diversity of isoforms used by this cell type
      ct_inter_diversity <- .calculate_inter_cellular_isoform_diversity(ct_iso_mat)
      
      # Calculate cell type heterogeneity
      # Measures cell-to-cell variation in isoform usage within this cell type
      ct_heterogeneity <- .calculate_intra_cell_type_heterogeneity(ct_iso_mat)
      
      #------------ 2. Additional Metric Calculations ------------#
      
      # Calculate average expression for each isoform in this cell type
      iso_means <- rowMeans(ct_iso_mat)
      nonzero_indices <- which(iso_means > 0)
      
      # Calculate proportion statistics
      total_expr <- sum(iso_means)
      iso_props <- rep(0, length(iso_means))
      if(total_expr > 0) {
        iso_props[nonzero_indices] <- iso_means[nonzero_indices] / total_expr
      }
      
      # Find dominant isoform in this cell type
      # Identifies the isoform with highest expression and its relative proportion
      if(length(nonzero_indices) > 0) {
        max_idx <- which.max(iso_props)
        dominant_iso_prop <- iso_props[max_idx]
        dominant_iso_name <- rownames(ct_iso_mat)[max_idx]
      } else {
        dominant_iso_prop <- 0
        dominant_iso_name <- ""
      }
      
      # Count expressed isoforms in this cell type
      n_expressed_isoforms <- length(nonzero_indices)
      
      # Calculate Simpson diversity index
      # Alternative measure of diversity that's more sensitive to dominant species
      simpson_index <- 1 - sum(iso_props^2)
      
      # Calculate evenness (how equally the isoforms are expressed)
      # Normalised diversity that controls for the number of expressed isoforms
      evenness <- ifelse(n_expressed_isoforms > 1,
                         ct_inter_diversity / log2(n_expressed_isoforms),
                         0)
      
      # Calculate cell-level metrics specific to this cell type
      cell_sums <- colSums(ct_iso_mat)
      cells_expressing <- sum(cell_sums > 0)
      pct_cells_expressing <- cells_expressing / ncol(ct_iso_mat) * 100
      
      # Calculate percentage of cells expressing multiple isoforms
      # Measures the prevalence of co-expression at the single-cell level
      if(cells_expressing > 0) {
        multi_iso_count <- 0
        for(j in which(cell_sums > 0)) {
          if(sum(ct_iso_mat[, j] > 0) > 1) {
            multi_iso_count <- multi_iso_count + 1
          }
        }
        pct_multi_iso_cells <- multi_iso_count / cells_expressing * 100
      } else {
        pct_multi_iso_cells <- 0
      }
      
      #------------ 3. Compile Results for This Cell Type ------------#
      
      # Store comprehensive set of metrics for this cell type
      cell_type_metrics[[cell_type]] <- list(
        # Core metrics
        intra_cellular_isoform_diversity = ct_intra_results$intra_cellular_isoform_diversity,
        inter_cellular_isoform_diversity = ct_inter_diversity,
        intra_cell_type_heterogeneity = ct_heterogeneity,
        
        # Isoform dominance metrics
        dominant_iso_prop = dominant_iso_prop,
        dominant_iso_name = as.character(dominant_iso_name),
        
        # Diversity characterisation metrics
        n_expressed_isoforms = n_expressed_isoforms,
        simpson_index = simpson_index,
        evenness = evenness,
        
        # Cell-level metrics
        cells_expressing = cells_expressing,
        pct_cells_expressing = pct_cells_expressing,
        pct_multi_iso_cells = pct_multi_iso_cells
      )
    }
  }
  
  # Return NULL if no valid cell types found
  if(length(cell_type_metrics) == 0) return(NULL)
  
  return(cell_type_metrics)
}

#' Calculate all complexity metrics for all genes in SCHT object
#'
#' @description
#' This function calculates a comprehensive set of metrics to quantify
#' transcriptomic complexity at multiple levels: within cells, across cells,
#' within cell types, and across cell types. It analyses isoform diversity,
#' cell type specificity, and expression patterns to provide a multidimensional
#' characterisation of transcriptomic complexity.
#'
#' @details
#' The function calculates seven core complexity metrics:
#'
#' 1. Intra-cellular Isoform Diversity: Measures the tendency for cells to co-express multiple
#'    isoforms of a gene
#'
#' 2. Inter-cellular Isoform Diversity: Measures the overall diversity of isoforms used
#'    across the entire cell population
#'
#' 3. Intra-cell type Heterogeneity: Measures the cell-to-cell variation in isoform
#'    usage within each cell type
#'
#' 4. Inter-cell type Specificity: Measures how differently a gene uses its isoforms
#'    across different cell types
#'
#' 5. Intra-cell type Heterogeneity Variability: Measures whether certain cell types show
#'    particularly high cellular heterogeneity
#'
#' 6. Inter-cell type Difference Variability: Measures whether certain cell type pairs
#'    show particularly significant differences in isoform usage patterns
#'
#' 7. Cell-type-specific Co-expression Variability: Measures whether a gene employs different
#'    isoform co-expression patterns in different cell types
#'
#' For each metric, a classification column is created to categorise genes based on
#' statistically determined thresholds. This provides a comprehensive complexity
#' profile for each gene across multiple dimensions.
#'
#' Some metrics may contain NA values, which have specific biological meanings:
#' - Inter-cell type metrics are NA for genes expressed in only one cell type
#' - Heterogeneity metrics are NA for genes with insufficient cells per cell type
#' - Diversity metrics are NA for genes with only a single expressed isoform
#' These NA values are properly handled during threshold determination and classification.
#'
#' @param scht_obj Single-Cell Hierarchical Tensor object
#' @param default_thresholds List of default thresholds for each metric (optional)
#' @param data_driven_thresholds Whether to use data-driven threshold detection (default: TRUE)
#' @param visualise Whether to create visualisations of threshold determination (default: TRUE)
#' @param batch_size Number of genes to process at once (default: 500)
#' @param min_samples Minimum number of samples required for reliable threshold modelling
#' @param verbose Whether to print detailed progress messages
#' @return A transcriptomic_complexity object (list) containing:
#'   \itemize{
#'     \item metrics: Data frame with complexity metrics and classifications for each gene
#'     \item single_cell_type_genes: List of genes expressed only in specific cell types
#'     \item thresholds: The thresholds used for classification
#'     \item threshold_plots: Visualisations of threshold fitting for each metric
#'     \item na_statistics: Statistics about NA value proportions for each metric
#'   }
#'   With performance attribute containing:
#'   \itemize{
#'     \item \code{total_time_sec}: Total processing time in seconds
#'     \item \code{memory_used_mb}: Memory utilised in megabytes
#'   }
#'
#' @examples
#' # Load example data
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' # Create SCHT object first
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   n_hvg = 3000,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   verbose = FALSE
#' )
#' 
#' # Example 1: Calculate complexity metrics with default settings
#' tc_results <- calculate_isoform_complexity_metrics(
#'   scht_obj,
#'   verbose = TRUE
#' )
#' 
#' # Examine the results
#' print(names(tc_results))
#' print(head(tc_results$metrics))
#' 
#' # Check classification distribution
#' print(table(tc_results$metrics$intra_cellular_isoform_diversity_class))
#' 
#' # Example 2: Use custom thresholds
#' custom_thresholds <- list(
#'   intra_cellular_isoform_diversity = 0.7,
#'   inter_cellular_isoform_diversity = 0.7,
#'   intra_cell_type_heterogeneity = 0.5,
#'   inter_cell_type_specificity = 0.5,
#'   intra_cell_type_heterogeneity_variability = 0.8,
#'   inter_cell_type_difference_variability = 0.6,
#'   cell_type_coexpression_variability = 1.0
#' )
#' 
#' \donttest{
#' tc_custom <- calculate_isoform_complexity_metrics(
#'   scht_obj,
#'   default_thresholds = custom_thresholds,
#'   data_driven_thresholds = FALSE,  # Use only custom thresholds
#'   verbose = FALSE
#' )
#' }
#' 
#' # Example 3: Process in smaller batches for memory efficiency
#' \donttest{
#' tc_batched <- calculate_isoform_complexity_metrics(
#'   scht_obj,
#'   batch_size = 100,  # Process 100 genes at a time
#'   verbose = TRUE
#' )
#' }
#' 
#' # Example 4: Access specific metrics for genes of interest
#' genes_of_interest <- head(rownames(tc_results$metrics), 10)
#' gene_metrics <- tc_results$metrics[
#'   tc_results$metrics$gene %in% genes_of_interest,
#' ]
#' print(gene_metrics[, c("gene", "intra_cellular_isoform_diversity", 
#'                        "inter_cellular_isoform_diversity")])
#' 
#' # Example 5: Identify highly complex genes
#' complex_genes <- tc_results$metrics[
#'   tc_results$metrics$intra_cellular_isoform_diversity_class == "high" &
#'   tc_results$metrics$inter_cellular_isoform_diversity_class == "high" &
#'   !is.na(tc_results$metrics$intra_cellular_isoform_diversity),
#' ]
#' print(paste("Number of highly complex genes:", nrow(complex_genes)))
#' 
#' @export
calculate_isoform_complexity_metrics <- function(scht_obj,
                                                 default_thresholds = list(
                                                   intra_cellular_isoform_diversity = 0.6,
                                                   inter_cellular_isoform_diversity = 0.6,
                                                   intra_cell_type_heterogeneity = 0.4,
                                                   inter_cell_type_specificity = 0.6,
                                                   intra_cell_type_heterogeneity_variability = 0.5,
                                                   inter_cell_type_difference_variability = 0.3,
                                                   cell_type_coexpression_variability = 0.4
                                                 ),
                                                 data_driven_thresholds = TRUE,
                                                 visualise = TRUE,
                                                 batch_size = 500,
                                                 min_samples = 20,
                                                 verbose = TRUE) {
  
  # Start timing and memory tracking
  start_time <- Sys.time()
  gc_start <- gc(full = TRUE)
  mem_start <- sum(gc_start[, "used"])
  
  # Handle both SCHT and IntegratedSCHT objects
  if (inherits(scht_obj, "IntegratedSCHT")) {
    # For IntegratedSCHT objects, use original_results
    gene_list <- scht_obj$original_results
  } else if (inherits(scht_obj, "SCHT")) {
    # For regular SCHT objects, use the object directly as the gene list
    gene_list <- scht_obj
  } else {
    stop("Input must be a SCHT or IntegratedSCHT object")
  }
  
  # Count total genes to process
  total_genes <- length(names(gene_list))
  message(paste("Calculating complexity metrics for", total_genes, "genes..."))
  
  # Check if progress package is available
  has_progress <- requireNamespace("progress", quietly = TRUE)
  if(!has_progress) {
    warning("progress package not available. Installing basic progress display. Consider installing the 'progress' package for better progress reporting.")
  }
  
  # Initialise results data frame
  results <- data.frame()
  
  # Initialise list to store single-cell_type genes
  single_cell_type_genes <- list()
  
  # Initialise list for cell type-specific metrics
  cell_type_metrics <- list()
  
  # Process genes in batches for better memory management
  gene_batches <- split(names(gene_list),
                        ceiling(seq_along(names(gene_list)) / batch_size))
  
  # Create a progress bar
  if(has_progress) {
    pb <- progress::progress_bar$new(
      format = "  Processing genes [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
      total = total_genes,
      clear = F,
      width = 100
    )
  } else {
    # Fallback to simpler progress indication without txtProgressBar
    message(paste("Processing", total_genes, "genes in", length(gene_batches), "batches"))
    progress_step <- max(1, floor(total_genes / 20))  # Report progress every ~5%
  }
  
  gene_counter <- 0
  
  # Process each batch
  for(batch_idx in seq_along(gene_batches)) {
    batch_genes <- gene_batches[[batch_idx]]
    batch_results <- list()
    
    # Process each gene in the batch
    for(gene in batch_genes) {
      gene_counter <- gene_counter + 1
      
      # Update progress
      if(has_progress) {
        pb$tick()
      } else if(gene_counter %% progress_step == 0 || gene_counter == total_genes) {
        message(sprintf("  Completed %d of %d genes (%.1f%%)",
                        gene_counter, total_genes, 100 * gene_counter / total_genes))
      }
      
      # Process this gene with appropriate diversity function
      gene_result <- .process_gene_complexity(scht_obj, gene)
      
      # Skip if null result
      if(is.null(gene_result)) next
      
      # Add to batch results
      batch_results[[gene]] <- gene_result
      
      # Analyse individual cell type metrics independently
      gene_ct_metrics <- .calculate_cell_type_metrics(scht_obj, gene)
      
      # Store if valid cell type metrics are available
      if(!is.null(gene_ct_metrics) && length(gene_ct_metrics) > 0) {
        cell_type_metrics[[gene]] <- gene_ct_metrics
      }
      
      # Store single-cell_type genes by cell type
      if(gene_result$is_single_cell_type && gene_result$specific_cell_type != "") {
        specific_cell_type <- gene_result$specific_cell_type
        if(!specific_cell_type %in% names(single_cell_type_genes)) {
          single_cell_type_genes[[specific_cell_type]] <- c()
        }
        single_cell_type_genes[[specific_cell_type]] <- c(single_cell_type_genes[[specific_cell_type]], gene)
      }
    }
    
    # Combine batch results
    if(length(batch_results) > 0) {
      batch_df <- do.call(rbind, batch_results)
      results <- rbind(results, batch_df)
    }
    
    # Clean memory after each batch
    gc(verbose = FALSE)
  }
  
  message(paste("Computed metrics for", nrow(results), "genes with sufficient data"))
  
  # Initialise thresholds, plots and NA stats lists
  threshold_values <- default_thresholds
  threshold_plots <- list()
  na_stats <- NULL
  
  # Determine optimal thresholds if enabled and sufficient data available
  if(data_driven_thresholds && nrow(results) >= min_samples) {
    optimal_results <- .determine_optimal_thresholds(
      results,
      default_thresholds = default_thresholds,
      visualise = visualise,
      min_samples = min_samples,
      verbose = verbose
    )
    threshold_values <- optimal_results$thresholds
    threshold_plots <- optimal_results$plots
    na_stats <- optimal_results$na_stats
  } else {
    message(paste("Using default thresholds for classification"))
    
    # Calculate basic NA statistics
    core_metrics <- c(
      "intra_cellular_isoform_diversity",
      "inter_cellular_isoform_diversity",
      "intra_cell_type_heterogeneity",
      "inter_cell_type_specificity",
      "intra_cell_type_heterogeneity_variability",
      "inter_cell_type_difference_variability",
      "cell_type_coexpression_variability"
    )
    
    core_metrics_names <-  c("Intra-cellular Isoform Diversity",
                             "Inter-cellular Isoform Diversity",
                             "Intra-cell-type Heterogeneity",
                             "Inter-cell-type Specificity",
                             "Intra-cell-type Heterogeneity Variability",
                             "Inter-cell-type Difference Variability",
                             "Cell-type-specific Co-expression Variability")
    
    na_counts <- sapply(results[core_metrics], function(x) sum(is.na(x)))
    na_percentages <- round(100 * na_counts / nrow(results), 1)
    
    na_stats <- list(
      counts = na_counts,
      percentages = na_percentages
    )
    
    # Report NA statistics
    message("\nNA value proportions for core metrics:")
    for(i in seq_along(core_metrics)) {
      metric_name <- core_metrics_names[i]
      message(sprintf("  - %s: %d NA values (%.1f%%)",
                      metric_name, na_counts[i], na_percentages[i]))
    }
    
    # If visualisation requested but not using data-driven thresholds,
    # create basic visualisations
    if(visualise) {
      for(metric in names(default_thresholds)) {
        if(metric %in% colnames(results)) {
          threshold_plots[[metric]] <- .visualise_threshold_fitting(
            results[[metric]],
            default_thresholds[[metric]],
            metric_name
          )
        }
      }
    }
  }
  
  # Add classification columns for each dimension
  results <- .classify_genes(results, threshold_values)
  
  # Sort by inter_cell_type_specificity for easier analysis
  results <- results[order(results$inter_cell_type_specificity, decreasing = TRUE), ]
  
  # Calculate performance metrics
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  gc_end <- gc(full = TRUE)
  mem_end <- sum(gc_end[, "used"])
  memory_used <- (mem_end - mem_start) / 1024  # Convert to MB
  
  if (verbose) {
    message(sprintf("\nComplexity metric calculation completed in %.2f seconds (%.2f minutes)", 
                   total_time, total_time/60))
    message(sprintf("Memory utilised: %.2f MB", memory_used))
  }
  
  # Return the metrics, single_cell_type_genes list, thresholds and visualisations
  TC_results <- list(
    metrics = results,
    cell_type_metrics = cell_type_metrics,
    single_cell_type_genes = single_cell_type_genes,
    thresholds = threshold_values,
    threshold_plots = threshold_plots,
    na_statistics = na_stats
  )
  
  # Add performance as an attribute
  attr(TC_results, "performance") <- list(
    total_time_sec = total_time,
    memory_used_mb = memory_used
  )
  
  class(TC_results) <- "transcriptomic_complexity"
  
  return(TC_results)
}

########################################
# Transcriptomic Complexity S3 Methods #
########################################

#' Print method for transcriptomic_complexity objects
#'
#' @param x A transcriptomic_complexity object
#' @param ... Additional arguments (not used)
#' @return Invisibly returns x
#' @examples
#' # Create SCHT and calculate complexity metrics
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
#' tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)
#' 
#' # Print the results
#' print(tc_results)
#' @export
print.transcriptomic_complexity <- function(x, ...) {
  cat("Isoform Complexity Analysis Result\n")
  cat("=================================\n\n")
  
  cat(sprintf("Total genes analysed: %d\n", nrow(x$metrics)))
  cat(sprintf("Classification thresholds determined by: %s\n",
              ifelse(is.null(x$na_statistics), "default values", "data-driven modelling")))
  
  cat("\nTop complexity classes:\n")
  
  # Get top complexity categories
  complexity_counts <- table(x$metrics$complexity_category)
  complexity_pcts <- round(100 * complexity_counts / sum(complexity_counts), 1)
  top_categories <- names(sort(complexity_counts, decreasing = TRUE))[1:length(complexity_counts)]
  
  for(i in seq_along(top_categories)) {
    cat(sprintf("  %d. %s: %d genes (%.1f%%)\n",
                i, top_categories[i], complexity_counts[top_categories[i]],
                complexity_pcts[top_categories[i]]))
  }
  
  cat("\nCall summary() for detailed statistics\n")
  
  invisible(x)
}

#' Summary method for transcriptomic_complexity objects
#'
#' @param object A transcriptomic_complexity object
#' @param ... Additional arguments (not used)
#' @return A list with summary statistics
#' @examples
#' # Using the same complexity results
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
#' tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)
#' 
#' # Get summary statistics
#' summary_stats <- summary(tc_results)
#' print(summary_stats$metric_summaries)
#' @export
summary.transcriptomic_complexity <- function(object, ...) {
  
  metrics_df <- object$metrics
  thresholds <- object$thresholds
  na_statistics <- object$na_statistics
  
  # Basic statistics
  gene_count <- nrow(metrics_df)
  
  # Create summary text with better formatting
  cat("Isoform Complexity Analysis Summary:\n")
  cat(sprintf("Total genes analysed: %d\n\n", gene_count))
  
  # Add NA value statistics if provided
  if(!is.null(na_statistics)) {
    cat("NA Value Statistics (representing biologically meaningful cases):\n")
    
    core_metrics <- c(
      "intra_cellular_isoform_diversity",
      "inter_cellular_isoform_diversity",
      "intra_cell_type_heterogeneity",
      "inter_cell_type_specificity",
      "intra_cell_type_heterogeneity_variability",
      "inter_cell_type_difference_variability",
      "cell_type_coexpression_variability"
    )
    
    core_metrics_names <-  c("Intra-cellular Isoform Diversity",
                             "Inter-cellular Isoform Diversity",
                             "Intra-cell-type Heterogeneity",
                             "Inter-cell-type Specificity",
                             "Intra-cell-type Heterogeneity Variability",
                             "Inter-cell-type Difference Variability",
                             "Cell-type-specific Co-expression Variability")
    
    na_explanations <- c(
      "Occurs when genes have only one expressed isoform, indicating lack of alternative splicing",
      "Occurs when only one isoform is expressed across the entire cell population",
      "Occurs when a cell type has insufficient cells expressing the gene (< 3 cells)",
      "Occurs for genes expressed in only a single cell type (high cell type-specificity)",
      "Occurs when a gene is expressed in fewer than 2 cell types with sufficient data",
      "Occurs when there are insufficient pairwise differences between cell types to calculate variability",
      "Occurs when a gene lacks sufficient cell type-specific expression data to compare co-expression patterns"
    )
    
    for(i in seq_along(core_metrics)) {
      metric_name <- core_metrics_names[i]
      count <- na_statistics$counts[i]
      percentage <- na_statistics$percentages[i]
      explanation <- na_explanations[i]
      
      cat(sprintf("  - %-35s: %4d genes (%.1f%%) - %s\n",
                  metric_name, count, percentage, explanation))
    }
    
    cat("\n")
  }
  
  cat("Classification Distribution Across Complexity Dimensions:\n")
  
  # List of classification columns
  class_columns <- c(
    "intra_cellular_isoform_diversity_class",
    "inter_cellular_isoform_diversity_class",
    "intra_cell_type_heterogeneity_class",
    "inter_cell_type_specificity_class",
    "intra_cell_type_heterogeneity_variability_class",
    "inter_cell_type_difference_variability_class",
    "cell_type_coexpression_variability_class"
  )
  
  core_class_columns_names <-  c("Intra-cellular Isoform Diversity",
                                 "Inter-cellular Isoform Diversity",
                                 "Intra-cell-type Heterogeneity",
                                 "Inter-cell-type Specificity",
                                 "Intra-cell-type Heterogeneity Variability",
                                 "Inter-cell-type Difference Variability",
                                 "Cell-type-specific Co-expression Variability")
  
  
  # Add distribution for each dimension with better formatting
  for(col in class_columns) {
    # Skip if column doesn't exist
    if(!col %in% colnames(metrics_df)) next
    
    # Get counts
    counts <- table(metrics_df[[col]])
    percentages <- round(100 * counts / gene_count, 1)
    
    # Add to summary with proper formatting
    idx <- match(col, class_columns)
    metrics_name <- core_class_columns_names[idx]
    
    cat(sprintf("\n%s:\n", metrics_name))
    
    for(category in names(counts)) {
      cat(sprintf("  - %-40s: %4d genes (%.1f%%)\n",
                  category, counts[category], percentages[category]))
    }
  }
  
  # Add traditional complexity categories 
  cat("\nTraditional Complexity Categories:\n")
  
  complexity_counts <- table(metrics_df$complexity_category)
  complexity_percentages <- round(100 * complexity_counts / gene_count, 1)
  
  for(category in names(complexity_counts)) {
    cat(sprintf("  - %-40s: %4d genes (%.1f%%)\n",
                category, complexity_counts[category], complexity_percentages[category]))
  }
  
  # Single-cell_type genes
  single_cell_type_count <- sum(metrics_df$is_single_cell_type, na.rm = TRUE)
  cat("\nCell-Type-Specific Genes:\n")
  cat(sprintf("  - %-40s: %4d genes (%.1f%%)\n",
              "Single cell-type genes", single_cell_type_count, round(100*single_cell_type_count/gene_count, 1)))
  
  # Print metric statistics as a formatted table
  cat("\nCore Metrics Statistics:\n")
  
  metric_stats <- data.frame(
    Metric = c("Intra-cellular Isoform Diversity",
               "Inter-cellular Isoform Diversity",
               "Intra-cell-type Heterogeneity",
               "Inter-cell-type Specificity",
               "Intra-cell-type Heterogeneity Variability",
               "Inter-cell-type Difference Variability",
               "Cell-type-specific Co-expression Variability"),
    Mean = c(mean(metrics_df$intra_cellular_isoform_diversity, na.rm = TRUE),
             mean(metrics_df$inter_cellular_isoform_diversity, na.rm = TRUE),
             mean(metrics_df$intra_cell_type_heterogeneity, na.rm = TRUE),
             mean(metrics_df$inter_cell_type_specificity, na.rm = TRUE),
             mean(metrics_df$intra_cell_type_heterogeneity_variability, na.rm = TRUE),
             mean(metrics_df$inter_cell_type_difference_variability, na.rm = TRUE),
             mean(metrics_df$cell_type_coexpression_variability, na.rm = TRUE)),
    Median = c(median(metrics_df$intra_cellular_isoform_diversity, na.rm = TRUE),
               median(metrics_df$inter_cellular_isoform_diversity, na.rm = TRUE),
               median(metrics_df$intra_cell_type_heterogeneity, na.rm = TRUE),
               median(metrics_df$inter_cell_type_specificity, na.rm = TRUE),
               median(metrics_df$intra_cell_type_heterogeneity_variability, na.rm = TRUE),
               median(metrics_df$inter_cell_type_difference_variability, na.rm = TRUE),
               median(metrics_df$cell_type_coexpression_variability, na.rm = TRUE)),
    SD = c(sd(metrics_df$intra_cellular_isoform_diversity, na.rm = TRUE),
           sd(metrics_df$inter_cellular_isoform_diversity, na.rm = TRUE),
           sd(metrics_df$intra_cell_type_heterogeneity, na.rm = TRUE),
           sd(metrics_df$inter_cell_type_specificity, na.rm = TRUE),
           sd(metrics_df$intra_cell_type_heterogeneity_variability, na.rm = TRUE),
           sd(metrics_df$inter_cell_type_difference_variability, na.rm = TRUE),
           sd(metrics_df$cell_type_coexpression_variability, na.rm = TRUE)),
    Min = c(min(metrics_df$intra_cellular_isoform_diversity, na.rm = TRUE),
            min(metrics_df$inter_cellular_isoform_diversity, na.rm = TRUE),
            min(metrics_df$intra_cell_type_heterogeneity, na.rm = TRUE),
            min(metrics_df$inter_cell_type_specificity, na.rm = TRUE),
            min(metrics_df$intra_cell_type_heterogeneity_variability, na.rm = TRUE),
            min(metrics_df$inter_cell_type_difference_variability, na.rm = TRUE),
            min(metrics_df$cell_type_coexpression_variability, na.rm = TRUE)),
    Max = c(max(metrics_df$intra_cellular_isoform_diversity, na.rm = TRUE),
            max(metrics_df$inter_cellular_isoform_diversity, na.rm = TRUE),
            max(metrics_df$intra_cell_type_heterogeneity, na.rm = TRUE),
            max(metrics_df$inter_cell_type_specificity, na.rm = TRUE),
            max(metrics_df$intra_cell_type_heterogeneity_variability, na.rm = TRUE),
            max(metrics_df$inter_cell_type_difference_variability, na.rm = TRUE),
            max(metrics_df$cell_type_coexpression_variability, na.rm = TRUE)),
    Threshold = c(thresholds$intra_cellular_isoform_diversity,
                  thresholds$inter_cellular_isoform_diversity,
                  thresholds$intra_cell_type_heterogeneity,
                  thresholds$inter_cell_type_specificity,
                  thresholds$intra_cell_type_heterogeneity_variability,
                  thresholds$inter_cell_type_difference_variability,
                  thresholds$cell_type_coexpression_variability),
    NA_Percent = c(
      mean(is.na(metrics_df$intra_cellular_isoform_diversity)) * 100,
      mean(is.na(metrics_df$inter_cellular_isoform_diversity)) * 100,
      mean(is.na(metrics_df$intra_cell_type_heterogeneity)) * 100,
      mean(is.na(metrics_df$inter_cell_type_specificity)) * 100,
      mean(is.na(metrics_df$intra_cell_type_heterogeneity_variability)) * 100,
      mean(is.na(metrics_df$inter_cell_type_difference_variability)) * 100,
      mean(is.na(metrics_df$cell_type_coexpression_variability)) * 100
    )
  )
  
  # Format metrics statistics
  metric_stats[, 2:8] <- round(metric_stats[, 2:8], 3)
  
  # Print table
  metric_df_formatted <- metric_stats
  if (requireNamespace("knitr", quietly = TRUE)) {
    print(knitr::kable(metric_df_formatted, row.names = FALSE))
  } else {
    # Otherwise use built-in print
    print(metric_df_formatted)
  }
  
  # Display performance metrics if available
  perf <- attr(object, "performance")
  if (!is.null(perf)) {
    cat("\nPerformance metrics:\n")
    cat(sprintf("  Processing time: %.2f seconds (%.2f minutes)\n", 
                perf$total_time_sec, perf$total_time_sec/60))
    cat(sprintf("  Memory utilised: %.2f MB\n", perf$memory_used_mb))
  }
  
  # Instead of printing everything, return the results as a list
  summary_results <- list(
    gene_count = gene_count,
    metric_stats = metric_stats,
    complexity_counts = complexity_counts,
    single_cell_type_count = single_cell_type_count,
    thresholds = thresholds,
    na_statistics = na_statistics,
    performance = perf
  )
  
  invisible(summary_results)
}

##############################
# Gene Selection Functions   #
##############################

#' Select genes for further investigation based on complexity classification
#'
#' @param metrics_df Data frame containing complexity metrics
#' @param category Classification of interest (e.g., "Strong Isoform Co-expression", "Cell Type-Specific")
#' @param column Classification column to filter on
#' @param top_n Number of top genes to select (by relevance)
#' @param sort_by Column to sort by (default depends on classification)
#' @return Vector of selected gene names
#' @export
select_genes_of_interest <- function(metrics_df,
                                     category,
                                     column = NULL,
                                     top_n = 20,
                                     sort_by = NULL) {
  # Handle column selection
  if(is.null(column)) {
    # Try to find matching column based on category
    if(grepl("Isoform Co-expression", category)) {
      column <- "intra_cellular_isoform_diversity_class"
    } else if(grepl("Isoform Diversity", category)) {
      column <- "inter_cellular_isoform_diversity_class"
    } else if(grepl("Cellular Heterogeneity", category) && !grepl("Variability", category)) {
      column <- "intra_cell_type_heterogeneity_class"
    } else if(grepl("Cell-Type-Specific", category) || grepl("Cell-Type-Independent", category)) {
      column <- "inter_cell_type_specificity_class"
    } else if(grepl("Variable Heterogeneity Across Cell Types", category) || grepl("Consistent Heterogeneity Across Cell Types", category)) {
      column <- "intra_cell_type_heterogeneity_variability_class"
    } else if(grepl("Cell Type Distinctions", category)) {
      column <- "inter_cell_type_difference_variability_class"
    } else if(grepl("Cell-Type-Adaptive Co-expression", category) || grepl("Cell-Type-Consistent Co-expression", category)) {
      column <- "cell_type_coexpression_variability_class"
    } else if(grepl("Diversity", category) && grepl("Specificity", category)) {
      column <- "complexity_category"
    } else {
      stop("Could not identify appropriate column for category: ", category)
    }
  }
  
  # Check if column exists
  if(!column %in% colnames(metrics_df)) {
    stop("Column not found in metrics data frame: ", column)
  }
  
  # Filter genes by category
  category_genes <- metrics_df[metrics_df[[column]] == category, ]
  
  # If no genes found, return empty vector
  if(nrow(category_genes) == 0) {
    message("No genes found in category: ", category)
    return(character(0))
  }
  
  # Determine column to sort by if not specified
  if(is.null(sort_by)) {
    if(column == "intra_cellular_isoform_diversity_class") {
      sort_by <- "intra_cellular_isoform_diversity"
    } else if(column == "inter_cellular_isoform_diversity_class") {
      sort_by <- "inter_cellular_isoform_diversity"
    } else if(column == "intra_cell_type_heterogeneity_class") {
      sort_by <- "intra_cell_type_heterogeneity"
    } else if(column == "inter_cell_type_specificity_class") {
      sort_by <- "inter_cell_type_specificity"
    } else if(column == "intra_cell_type_heterogeneity_variability_class") {
      sort_by <- "intra_cell_type_heterogeneity_variability"
    } else if(column == "inter_cell_type_difference_variability_class") {
      sort_by <- "inter_cell_type_difference_variability"
    } else if(column == "cell_type_coexpression_variability_class") {
      sort_by <- "cell_type_coexpression_variability"
    } else {
      sort_by <- "inter_cell_type_specificity"  # Default
    }
  }
  
  # Check if sort column exists
  if(!sort_by %in% colnames(category_genes)) {
    warning("Sort column not found, using gene name instead: ", sort_by)
    sort_by <- "gene"
  }
  
  # Sort genes
  if(sort_by == "gene") {
    category_genes <- category_genes[order(category_genes$gene), ]
  } else {
    # Most metrics are better when higher
    category_genes <- category_genes[order(category_genes[[sort_by]], decreasing = TRUE), ]
  }
  
  # Select top N genes
  top_genes <- head(category_genes$gene, top_n)
  
  return(top_genes)
}

#' Find multi-dimensional complexity patterns
#'
#' This function identifies genes that match specific combinations of complexity
#' classifications across multiple dimensions. It enables discovery of genes with
#' interesting or unusual complexity profiles.
#'
#' @param metrics_df Data frame containing complexity metrics
#' @param pattern Named list of category patterns to match
#' @param top_n Number of top genes to select
#' @param sort_by Column to sort by (default: "inter_cell_type_specificity")
#' @return Vector of selected gene names
#' @export
find_complexity_pattern <- function(metrics_df,
                                    pattern,
                                    top_n = 20,
                                    sort_by = "inter_cell_type_specificity") {
  # Check that pattern is a named list
  if(!is.list(pattern) || is.null(names(pattern))) {
    stop("Pattern must be a named list with classification column names as keys")
  }
  
  # Initial filter
  genes_idx <- rep(TRUE, nrow(metrics_df))
  
  # Apply each filter
  for(col in names(pattern)) {
    # Check if column exists
    if(!col %in% colnames(metrics_df)) {
      warning("Column not found in metrics data frame: ", col)
      next
    }
    
    # Apply filter
    genes_idx <- genes_idx & (metrics_df[[col]] == pattern[[col]])
  }
  
  # Extract matching genes
  matching_genes <- metrics_df[genes_idx, ]
  
  # If no genes found, return empty vector
  if(nrow(matching_genes) == 0) {
    message("No genes found matching the specified pattern")
    return(character(0))
  }
  
  # Check if sort column exists
  if(!sort_by %in% colnames(matching_genes)) {
    warning("Sort column not found, using gene name instead: ", sort_by)
    sort_by <- "gene"
  }
  
  # Sort genes
  if(sort_by == "gene") {
    matching_genes <- matching_genes[order(matching_genes$gene), ]
  } else {
    # Most metrics are better when higher
    matching_genes <- matching_genes[order(matching_genes[[sort_by]], decreasing = TRUE), ]
  }
  
  # Select top N genes
  top_genes <- head(matching_genes$gene, top_n)
  
  message(paste0("Found ", nrow(matching_genes), " genes matching the pattern. Returning top ",
                 length(top_genes), " sorted by ", sort_by))
  
  return(top_genes)
}
