##################################################################
#  Co-expression Analysis Functions for ScIsoX                   #
#                                                                #
#  Comprehensive co-expression analysis including correlation    #
#  calculation, isoform switching detection, and conservation    #
#  analysis across cell types                                    #
#                                                                #
#  Author: Siyuan Wu & Ulf Schmitz                               #
#  Date: Jul 29, 2025                                            #
#  Package: ScIsoX V1.1.0                                        #
##################################################################

# Define global variables to avoid R CMD check notes
utils::globalVariables(c("correlation", "pair"))

#' Calculate isoform co-expression correlation matrix
#'
#' @description
#' Calculates correlation matrix between isoforms of a gene using
#' the specified correlation method.
#'
#' @param scht_obj SCHT or IntegratedSCHT object
#' @param gene Gene name to analyse
#' @param method Correlation method: "pearson", "spearman", or "kendall"
#' @param min_cells Minimum number of cells required for analysis
#' @param min_expression Minimum mean expression threshold
#'
#' @return A list containing:
#'   \itemize{
#'     \item cor_matrix: Correlation matrix between isoforms
#'     \item n_isoforms: Number of isoforms
#'     \item n_cells: Number of cells
#'     \item method: Correlation method used
#'   }
#'
#' @export
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
#' # Example 1: Find genes with multiple isoforms in the filtered data
#' # Get the actual gene list from the SCHT object
#' if (inherits(scht_obj, "IntegratedSCHT")) {
#'   gene_list <- names(scht_obj$original_results)
#' } else {
#'   gene_list <- names(scht_obj)
#' }
#' 
#' # Find genes with multiple isoforms
#' multi_iso_genes <- character()
#' for (g in gene_list[1:min(50, length(gene_list))]) {
#'   if (inherits(scht_obj, "IntegratedSCHT")) {
#'     n_iso <- nrow(scht_obj$original_results[[g]])
#'   } else {
#'     n_iso <- nrow(scht_obj[[g]])
#'   }
#'   if (n_iso > 2) multi_iso_genes <- c(multi_iso_genes, g)
#' }
#' 
#' # Analyze the first multi-isoform gene found
#' if (length(multi_iso_genes) > 0) {
#'   cor_result <- calculate_isoform_coexpression(scht_obj, multi_iso_genes[1])
#'   cat("Gene:", multi_iso_genes[1], "\n")
#'   cat("Number of isoforms:", cor_result$n_isoforms, "\n")
#'   print(cor_result$cor_matrix)
#'   
#'   # Check switching patterns
#'   switching <- detect_isoform_switching(cor_result)
#'   if (switching$n_switching_pairs > 0) {
#'     cat("Switching pairs found:", switching$n_switching_pairs, "\n")
#'   }
#' }
calculate_isoform_coexpression <- function(scht_obj, 
                                         gene, 
                                         method = "pearson",
                                         min_cells = 10,
                                         min_expression = 0) {
  
  # Handle both SCHT and IntegratedSCHT objects
  if (inherits(scht_obj, "IntegratedSCHT")) {
    gene_list <- scht_obj$original_results
  } else if (inherits(scht_obj, "SCHT")) {
    gene_list <- scht_obj
  } else {
    stop("Input must be a SCHT or IntegratedSCHT object")
  }
  
  # Check if gene exists
  if(!gene %in% names(gene_list)) {
    stop(paste("Gene", gene, "not found in SCHT object"))
  }
  
  # Get isoform expression matrix
  iso_mat <- gene_list[[gene]]
  
  # Skip if insufficient data
  if(nrow(iso_mat) < 2 || ncol(iso_mat) < min_cells) {
    stop(paste("Gene", gene, "has insufficient data (needs multiple isoforms and at least", 
               min_cells, "cells)"))
  }
  
  # Filter by expression if specified
  if (min_expression > 0) {
    mean_expr <- rowMeans(iso_mat)
    keep_isoforms <- mean_expr >= min_expression
    
    if (sum(keep_isoforms) < 2) {
      stop(paste("Gene", gene, "has fewer than 2 isoforms after expression filtering"))
    }
    
    iso_mat <- iso_mat[keep_isoforms, , drop = FALSE]
  }
  
  # Calculate correlation matrix between isoforms
  cor_mat <- cor(t(iso_mat), method = method, use = "pairwise.complete.obs")
  
  return(list(
    cor_matrix = cor_mat,
    n_isoforms = nrow(iso_mat),
    n_cells = ncol(iso_mat),
    method = method,
    gene = gene
  ))
}

#' Calculate gene co-expression across all cell types
#'
#' @description
#' Calculates co-expression patterns for a gene across all cell types,
#' returning both overall and cell type-specific correlations.
#'
#' @param scht_obj IntegratedSCHT object with cell type information
#' @param gene Gene name to analyse
#' @param method Correlation method: "pearson", "spearman", or "kendall"
#' @param min_cells Minimum cells per cell type
#' @param min_expression Minimum mean expression threshold
#'
#' @return A list containing overall and cell type-specific correlations
#'
#' @examples
#' # Load example data
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' # Create IntegratedSCHT object with cell type information
#' integrated_scht <- create_scht(
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
#'   require_cell_type = TRUE,  # Creates IntegratedSCHT
#'   verbose = FALSE
#' )
#' 
#' # Use specific genes known to have multiple isoforms
#' test_gene <- "Mapk13"  # Known multi-isoform gene in the example data
#' 
#' tryCatch({
#'   # Calculate coexpression across all cell types
#'   ct_cor_results <- calculate_gene_coexpression_all_celltypes(
#'     integrated_scht,
#'     gene = test_gene,
#'     method = "pearson",
#'     min_cells = 10
#'   )
#'   
#'   # Examine overall correlation
#'   cat("\nOverall correlation matrix for", test_gene, ":\n")
#'   print(round(ct_cor_results$overall$cor_matrix, 2))
#'   
#'   # Examine cell type-specific patterns
#'   cat("\nFound correlations for", 
#'       length(ct_cor_results$cell_types), "cell types\n")
#'   
#'   # Compare correlations across cell types
#'   for (ct in names(ct_cor_results$cell_types)) {
#'     ct_data <- ct_cor_results$cell_types[[ct]]
#'     cat("\nCell type:", ct, "\n")
#'     cat("Number of cells:", ct_data$n_cells, "\n")
#'     cat("Number of isoforms:", ct_data$n_isoforms, "\n")
#'     
#'     # Show correlation matrix
#'     if (!is.null(ct_data$cor_matrix)) {
#'       cat("Correlation matrix:\n")
#'       print(round(ct_data$cor_matrix, 2))
#'     }
#'   }
#'   
#'   # Identify cell type-specific switching
#'   switching_summary <- lapply(names(ct_cor_results$cell_types), function(ct) {
#'     ct_data <- ct_cor_results$cell_types[[ct]]
#'     if (!is.null(ct_data$cor_matrix)) {
#'       switching <- detect_isoform_switching(ct_data)
#'       return(data.frame(
#'         cell_type = ct,
#'         n_switching_pairs = switching$n_switching_pairs,
#'         n_strong_switching = switching$n_strong_switching
#'       ))
#'     }
#'     return(NULL)
#'   })
#'   
#'   switching_df <- do.call(rbind, switching_summary[!sapply(switching_summary, is.null)])
#'   if (!is.null(switching_df) && nrow(switching_df) > 0) {
#'     cat("\nSwitching summary by cell type:\n")
#'     print(switching_df)
#'   }
#'   
#'   # Use different correlation methods
#'   ct_cor_spearman <- calculate_gene_coexpression_all_celltypes(
#'     integrated_scht,
#'     gene = test_gene,
#'     method = "spearman"
#'   )
#'   
#'   # Compare methods
#'   cat("\nComparing Pearson vs Spearman correlations (overall):\n")
#'   pearson_vals <- ct_cor_results$overall$cor_matrix[upper.tri(ct_cor_results$overall$cor_matrix)]
#'   spearman_vals <- ct_cor_spearman$overall$cor_matrix[upper.tri(ct_cor_spearman$overall$cor_matrix)]
#'   cat("Mean absolute difference:", 
#'       round(mean(abs(pearson_vals - spearman_vals)), 3), "\n")
#' }, error = function(e) {
#'   cat("Error analyzing gene", test_gene, ":", e$message, "\n")
#' })
#' 
#' @export
calculate_gene_coexpression_all_celltypes <- function(scht_obj,
                                                    gene,
                                                    method = "pearson",
                                                    min_cells = 10,
                                                    min_expression = 0) {
  
  if (!inherits(scht_obj, "IntegratedSCHT")) {
    stop("This function requires an IntegratedSCHT object with cell type information")
  }
  
  # Calculate overall correlation
  overall_cor <- calculate_isoform_coexpression(
    scht_obj, gene, method, min_cells, min_expression
  )
  
  # Calculate cell type-specific correlations
  ct_matrices <- scht_obj$cell_type_matrices
  celltype_correlations <- list()
  
  for (ct in names(ct_matrices)) {
    if (!gene %in% names(ct_matrices[[ct]])) next
    
    tryCatch({
      # Get the expression matrix for this cell type
      expr_matrix <- ct_matrices[[ct]][[gene]]
      
      # Check if we have enough cells
      if (ncol(expr_matrix) < min_cells) next
      
      # Filter by expression
      mean_expr <- rowMeans(expr_matrix)
      keep_isoforms <- mean_expr >= min_expression
      if (sum(keep_isoforms) < 2) next
      
      # Calculate correlation for this cell type
      expr_matrix_filtered <- expr_matrix[keep_isoforms, , drop = FALSE]
      
      # Remove isoforms with zero variance
      iso_vars <- apply(expr_matrix_filtered, 1, var)
      keep_variable <- iso_vars > 0
      if (sum(keep_variable) < 2) next
      
      expr_matrix_filtered <- expr_matrix_filtered[keep_variable, , drop = FALSE]
      cor_matrix <- cor(t(expr_matrix_filtered), method = method, use = "pairwise.complete.obs")
      
      # Create result in the same format as calculate_isoform_coexpression
      ct_cor <- list(
        cor_matrix = cor_matrix,
        n_isoforms = nrow(expr_matrix_filtered),
        n_cells = ncol(expr_matrix),
        method = method,
        gene = gene,
        cell_type = ct
      )
      
      celltype_correlations[[ct]] <- ct_cor
    }, error = function(e) {
      # Skip cell types with insufficient data
      NULL
    })
  }
  
  return(list(
    overall = overall_cor,
    cell_types = celltype_correlations,
    gene = gene,
    method = method
  ))
}

#' Detect isoform switching patterns
#'
#' @description
#' Identifies potential isoform switching events based on negative
#' correlations between isoforms.
#'
#' @param cor_result Result from calculate_isoform_coexpression
#' @param threshold Correlation threshold for switching detection (default: -0.3)
#' @param strong_threshold Threshold for strong switching (default: -0.5)
#'
#' @return A list containing switching pairs and statistics
#'
#' @examples
#' # Load example data
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' # Create SCHT object
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
#' # Use specific genes known to have multiple isoforms
#' test_genes <- c("Mapk13", "Atl1", "Irf8")
#' 
#' for (gene in test_genes) {
#'   tryCatch({
#'     # Calculate coexpression for the gene
#'     cor_result <- calculate_isoform_coexpression(
#'       scht_obj, 
#'       gene = gene,
#'       method = "pearson"
#'     )
#'     
#'     # Detect switching with default threshold (-0.3)
#'     switching_default <- detect_isoform_switching(cor_result)
#'     cat("\nGene:", gene, "\n")
#'     cat("Found", switching_default$n_switching_pairs, "switching pairs\n")
#'     
#'     # Detect switching with more stringent threshold
#'     switching_strict <- detect_isoform_switching(
#'       cor_result, 
#'       threshold = -0.4,
#'       strong_threshold = -0.6
#'     )
#'     
#'     # Examine switching pairs
#'     if (switching_default$n_switching_pairs > 0) {
#'       cat("Switching pairs detected:\n")
#'       print(switching_default$switching_pairs)
#'       
#'       # Check specific pair types
#'       strong_pairs <- switching_default$switching_pairs[
#'         switching_default$switching_pairs$strength == "strong", 
#'       ]
#'       if (nrow(strong_pairs) > 0) {
#'         cat("Strong switching pairs (r <", 
#'             switching_default$strong_threshold, "):\n")
#'         print(strong_pairs)
#'       }
#'     }
#'     
#'     # For cell type-specific analysis (only if IntegratedSCHT)
#'     if (inherits(scht_obj, "IntegratedSCHT")) {
#'       ct_cor <- calculate_gene_coexpression_all_celltypes(
#'         scht_obj, 
#'         gene = gene
#'       )
#'       
#'       # Check switching in specific cell types
#'       for (ct in names(ct_cor$cell_types)) {
#'         ct_switching <- detect_isoform_switching(ct_cor$cell_types[[ct]])
#'         if (ct_switching$n_switching_pairs > 0) {
#'           cat("Cell type", ct, "has", 
#'               ct_switching$n_switching_pairs, "switching pairs\n")
#'         }
#'       }
#'     }
#'   }, error = function(e) {
#'     cat("Could not analyze gene", gene, ":", e$message, "\n")
#'   })
#' }
#' 
#' @export
detect_isoform_switching <- function(cor_result, 
                                   threshold = -0.3,
                                   strong_threshold = -0.5) {
  
  if (is.list(cor_result) && !is.null(cor_result$cor_matrix)) {
    cor_matrix <- cor_result$cor_matrix
  } else if (is.matrix(cor_result)) {
    cor_matrix <- cor_result
  } else {
    stop("Input must be a correlation matrix or result from calculate_isoform_coexpression")
  }
  
  # Find negative correlations
  negative_pairs <- which(cor_matrix < threshold & upper.tri(cor_matrix), arr.ind = TRUE)
  
  if (nrow(negative_pairs) == 0) {
    return(list(
      switching_pairs = data.frame(),
      n_switching_pairs = 0,
      n_strong_switching = 0,
      summary = "No significant negative correlations found"
    ))
  }
  
  # Create pairs data frame
  switching_pairs <- data.frame(
    isoform1 = rownames(cor_matrix)[negative_pairs[, 1]],
    isoform2 = rownames(cor_matrix)[negative_pairs[, 2]],
    correlation = apply(negative_pairs, 1, function(idx) cor_matrix[idx[1], idx[2]]),
    stringsAsFactors = FALSE
  )
  
  # Add switching strength
  switching_pairs$strength <- ifelse(
    switching_pairs$correlation < strong_threshold, 
    "strong", 
    "moderate"
  )
  
  # Sort by correlation strength
  switching_pairs <- switching_pairs[order(switching_pairs$correlation), ]
  
  return(list(
    switching_pairs = switching_pairs,
    n_switching_pairs = nrow(switching_pairs),
    n_strong_switching = sum(switching_pairs$strength == "strong"),
    threshold = threshold,
    strong_threshold = strong_threshold
  ))
}

#' Calculate co-expression statistics
#'
#' @description
#' Generates comprehensive statistics for co-expression analysis results.
#'
#' @param cor_result Result from calculate_isoform_coexpression
#' @param include_pairwise Include pairwise correlation details
#'
#' @return A list of statistical summaries
#'
#' @export
calculate_coexpression_stats <- function(cor_result, include_pairwise = FALSE) {
  
  if (is.list(cor_result) && !is.null(cor_result$cor_matrix)) {
    cor_matrix <- cor_result$cor_matrix
  } else if (is.matrix(cor_result)) {
    cor_matrix <- cor_result
  } else {
    stop("Input must be a correlation matrix or result from calculate_isoform_coexpression")
  }
  
  # Extract upper triangle
  upper_tri <- cor_matrix[upper.tri(cor_matrix)]
  
  # Basic statistics
  stats <- list(
    n_isoforms = nrow(cor_matrix),
    n_pairs = length(upper_tri),
    mean_correlation = mean(upper_tri, na.rm = TRUE),
    median_correlation = median(upper_tri, na.rm = TRUE),
    sd_correlation = sd(upper_tri, na.rm = TRUE),
    min_correlation = min(upper_tri, na.rm = TRUE),
    max_correlation = max(upper_tri, na.rm = TRUE),
    
    # Categorized counts
    n_strong_positive = sum(upper_tri > 0.5, na.rm = TRUE),
    n_moderate_positive = sum(upper_tri > 0.3 & upper_tri <= 0.5, na.rm = TRUE),
    n_weak_positive = sum(upper_tri > 0 & upper_tri <= 0.3, na.rm = TRUE),
    n_weak_negative = sum(upper_tri < 0 & upper_tri >= -0.3, na.rm = TRUE),
    n_moderate_negative = sum(upper_tri < -0.3 & upper_tri >= -0.5, na.rm = TRUE),
    n_strong_negative = sum(upper_tri < -0.5, na.rm = TRUE)
  )
  
  # Add percentages
  stats$pct_positive = sum(upper_tri > 0, na.rm = TRUE) / length(upper_tri) * 100
  stats$pct_negative = sum(upper_tri < 0, na.rm = TRUE) / length(upper_tri) * 100
  
  # Pairwise details if requested
  if (include_pairwise && length(upper_tri) > 0) {
    pairs_idx <- which(upper.tri(cor_matrix), arr.ind = TRUE)
    stats$pairwise_details <- data.frame(
      isoform1 = rownames(cor_matrix)[pairs_idx[, 1]],
      isoform2 = rownames(cor_matrix)[pairs_idx[, 2]],
      correlation = upper_tri,
      stringsAsFactors = FALSE
    )
    stats$pairwise_details <- stats$pairwise_details[
      order(abs(stats$pairwise_details$correlation), decreasing = TRUE), 
    ]
  }
  
  return(stats)
}

#' Analyse co-expression conservation across cell types
#'
#' @description
#' Analyses how conserved co-expression patterns are across different cell types,
#' identifying isoform pairs with consistent relationships.
#'
#' @param integrated_scht IntegratedSCHT object
#' @param gene Gene to analyse
#' @param method Correlation method
#' @param min_cells Minimum cells per cell type
#' @param min_expression Minimum expression threshold
#' @param consistency_threshold Proportion of cell types needed for "conserved" classification
#' @param correlation_threshold Threshold for considering a correlation significant (default: 0.3)
#'
#' @return Conservation analysis results
#' @export
analyse_coexpression_conservation <- function(integrated_scht, 
                                            gene,
                                            method = "pearson",
                                            min_cells = 10,
                                            min_expression = 0,
                                            consistency_threshold = 0.7,
                                            correlation_threshold = 0.3) {
  
  if (!inherits(integrated_scht, "IntegratedSCHT")) {
    stop("Input must be an IntegratedSCHT object")
  }
  
  # Get correlations for all cell types
  coexpr_result <- calculate_gene_coexpression_all_celltypes(
    integrated_scht, gene, method, min_cells, min_expression
  )
  
  if (length(coexpr_result$cell_types) < 2) {
    stop("Need at least 2 cell types with sufficient data for conservation analysis")
  }
  
  # Extract all correlation matrices
  cor_matrices <- lapply(coexpr_result$cell_types, function(x) x$cor_matrix)
  
  # Get all unique isoform pairs
  all_isoforms <- unique(unlist(lapply(cor_matrices, rownames)))
  n_isoforms <- length(all_isoforms)
  
  if (n_isoforms < 2) {
    stop("Need at least 2 isoforms for conservation analysis")
  }
  
  # Initialize conservation matrix
  conservation_results <- list()
  
  # Analyze each pair
  for (i in 1:(n_isoforms-1)) {
    for (j in (i+1):n_isoforms) {
      iso1 <- all_isoforms[i]
      iso2 <- all_isoforms[j]
      pair_name <- paste(iso1, iso2, sep = "-")
      
      # Get correlations across cell types
      pair_correlations <- numeric()
      cell_types_with_pair <- character()
      
      for (ct_name in names(cor_matrices)) {
        cor_mat <- cor_matrices[[ct_name]]
        if (iso1 %in% rownames(cor_mat) && iso2 %in% rownames(cor_mat)) {
          cor_value <- cor_mat[iso1, iso2]
          # Skip NA correlations (from zero variance)
          if (!is.na(cor_value)) {
            pair_correlations <- c(pair_correlations, cor_value)
            cell_types_with_pair <- c(cell_types_with_pair, ct_name)
          }
        }
      }
      
      if (length(pair_correlations) >= 1) {
        # Classify conservation pattern
        n_positive <- sum(pair_correlations > correlation_threshold)
        n_negative <- sum(pair_correlations < -correlation_threshold)
        n_total <- length(pair_correlations)
        
        # Classification logic based on number of cell types
        if (n_total == 1) {
          # Pair only appears in one cell type
          pattern <- "Cell_Type_Specific"
        } else {
          # Pair appears in multiple cell types - check for mixed pattern first
          if (n_positive > 0 && n_negative > 0) {
            # Has both positive and negative correlations
            pattern <- "Mixed"
          } else if (n_positive / n_total >= consistency_threshold) {
            pattern <- "Conserved_Positive"
          } else if (n_negative / n_total >= consistency_threshold) {
            pattern <- "Conserved_Negative"
          } else {
            # Not conserved, not mixed - weak correlations across cell types
            pattern <- "Non_Conserved"
          }
        }
        
        conservation_results[[pair_name]] <- list(
          isoform1 = iso1,
          isoform2 = iso2,
          correlations = pair_correlations,
          cell_types = cell_types_with_pair,
          mean_correlation = mean(pair_correlations),
          sd_correlation = if(length(pair_correlations) > 1) sd(pair_correlations) else NA,
          pattern = pattern,
          n_cell_types = length(pair_correlations)
        )
      }
    }
  }
  
  # Summarize results
  summary_stats <- list(
    n_pairs_analysed = length(conservation_results),
    n_conserved_positive = sum(sapply(conservation_results, function(x) x$pattern == "Conserved_Positive")),
    n_conserved_negative = sum(sapply(conservation_results, function(x) x$pattern == "Conserved_Negative")),
    n_mixed = sum(sapply(conservation_results, function(x) x$pattern == "Mixed")),
    n_cell_type_specific = sum(sapply(conservation_results, function(x) x$pattern == "Cell_Type_Specific")),
    n_non_conserved = sum(sapply(conservation_results, function(x) x$pattern == "Non_Conserved"))
  )
  
  return(list(
    gene = gene,
    pairs = conservation_results,
    summary = summary_stats,
    cell_types = names(cor_matrices),
    consistency_threshold = consistency_threshold
  ))
}

#' Plot isoform co-expression heatmap
#'
#' @description
#' Creates a heatmap visualisation of isoform co-expression patterns.
#' This is the visualisation companion to calculate_isoform_coexpression.
#'
#' @param scht_obj SCHT or IntegratedSCHT object  
#' @param gene Gene name to visualize
#' @param method Correlation method: "pearson", "spearman", or "kendall"
#' @param display_numbers Whether to show correlation values
#' @param cluster_rows Whether to cluster rows
#' @param cluster_columns Whether to cluster columns
#' @param title Custom title for the heatmap
#'
#' @return A ComplexHeatmap object
#' @export
plot_isoform_coexpression <- function(scht_obj, 
                                    gene, 
                                    method = "pearson", 
                                    display_numbers = FALSE,
                                    cluster_rows = TRUE,
                                    cluster_columns = TRUE,
                                    title = NULL) {
  
  # Check if packages are available
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is needed for this function. Please install it with: BiocManager::install('ComplexHeatmap')")
  }
  
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is needed for this function. Please install it with: install.packages('circlize')")
  }
  
  # Calculate correlation
  cor_result <- calculate_isoform_coexpression(scht_obj, gene, method)
  cor_mat <- cor_result$cor_matrix
  
  # Create title if not provided
  if (is.null(title)) {
    title <- paste0(gene, " Isoform Co-expression (", 
                   toupper(substring(method, 1, 1)), substring(method, 2), ")")
  }
  
  # Create a diverging colour mapping
  heatmap_colours <- circlize::colorRamp2(
    c(-1, -0.5, 0, 0.5, 1),
    c("#6a0624", "#feab88", "#f7f7f7", "#6fafd2", "#053061")
  )
  
  # Format cell labels if display_numbers is TRUE
  if (display_numbers) {
    cell_fun <- function(j, i, x, y, width, height, fill) {
      if (i != j) {
        value <- sprintf("%.2f", cor_mat[i, j])
        grid::grid.text(value, x, y, gp = grid::gpar(
          fontsize = 10,
          col = ifelse(abs(cor_mat[i, j]) < 0.5, "black", "white")
        ))
      }
    }
  } else {
    cell_fun <- NULL
  }
  
  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    cor_mat,
    name = paste0(toupper(substring(method, 1, 1)), substring(method, 2), "\nCorrelation"),
    col = heatmap_colours,
    cell_fun = cell_fun,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    column_title = title,
    row_title = "Isoforms",
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    row_title_gp = grid::gpar(fontsize = 12),
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 10),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9)
    )
  )
  
  return(ht)
}

#' Plot co-expression patterns across cell types
#'
#' @description
#' Visualizes how co-expression patterns vary across different cell types
#' for a specific gene.
#'
#' @param coexpr_all_result Result from calculate_gene_coexpression_all_celltypes
#' @param pair_selection Which pairs to show: "all", "switching", "conserved"
#' @param threshold Threshold for pair selection
#'
#' @return A ggplot object
#' @export
plot_coexpression_across_celltypes <- function(coexpr_all_result,
                                             pair_selection = "all",
                                             threshold = 0.3) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }
  
  # Extract data
  ct_data <- coexpr_all_result$cell_types
  gene <- coexpr_all_result$gene
  
  # Collect all pairs and their correlations
  plot_data <- list()
  
  for (ct_name in names(ct_data)) {
    cor_mat <- ct_data[[ct_name]]$cor_matrix
    pairs_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
    
    for (k in 1:nrow(pairs_idx)) {
      i <- pairs_idx[k, 1]
      j <- pairs_idx[k, 2]
      pair_name <- paste(rownames(cor_mat)[i], rownames(cor_mat)[j], sep = " - ")
      
      plot_data[[length(plot_data) + 1]] <- data.frame(
        cell_type = ct_name,
        pair = pair_name,
        correlation = cor_mat[i, j],
        stringsAsFactors = FALSE
      )
    }
  }
  
  plot_df <- do.call(rbind, plot_data)
  
  # Filter based on selection
  if (pair_selection == "switching") {
    # Keep pairs that show switching (positive in some, negative in others)
    switching_pairs <- unique(plot_df$pair[plot_df$pair %in% 
      names(which(tapply(plot_df$correlation, plot_df$pair, 
                        function(x) any(x > threshold) && any(x < -threshold))))])
    plot_df <- plot_df[plot_df$pair %in% switching_pairs, ]
  } else if (pair_selection == "conserved") {
    # Keep pairs that are consistently positive or negative
    conserved_pairs <- unique(plot_df$pair[plot_df$pair %in% 
      names(which(tapply(plot_df$correlation, plot_df$pair, 
                        function(x) all(x > threshold) || all(x < -threshold))))])
    plot_df <- plot_df[plot_df$pair %in% conserved_pairs, ]
  }
  
  if (nrow(plot_df) == 0) {
    stop("No pairs match the selection criteria")
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = cell_type, y = correlation, 
                                             group = pair, color = pair)) +
    ggplot2::geom_line(size = 1, alpha = 0.7) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_hline(yintercept = c(-threshold, threshold), 
                       linetype = "dotted", color = "gray70") +
    ggplot2::scale_y_continuous(limits = c(-1, 1)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = paste("Co-expression Patterns Across Cell Types:", gene),
      subtitle = paste("Selection:", pair_selection),
      x = "Cell Type",
      y = "Correlation",
      color = "Isoform Pair"
    )
  
  return(p)
}

#' Calculate cell type co-expression similarity
#'
#' @description
#' Calculates similarity between cell types based on their co-expression patterns.
#'
#' @param correlation_list List of correlation matrices by cell type
#' @param method Distance method: "correlation", "euclidean", or "manhattan"
#'
#' @return Distance matrix between cell types
#' @export
calculate_celltype_coexpression_similarity <- function(correlation_list, 
                                                     method = "correlation") {
  
  if (length(correlation_list) < 2) {
    stop("Need at least 2 cell types for similarity analysis")
  }
  
  # Get all unique isoforms
  all_isoforms <- unique(unlist(lapply(correlation_list, rownames)))
  n_celltypes <- length(correlation_list)
  celltype_names <- names(correlation_list)
  
  # Create distance matrix
  dist_matrix <- matrix(NA, n_celltypes, n_celltypes)
  rownames(dist_matrix) <- colnames(dist_matrix) <- celltype_names
  
  # Calculate pairwise distances
  for (i in 1:(n_celltypes - 1)) {
    for (j in (i + 1):n_celltypes) {
      ct1 <- celltype_names[i]
      ct2 <- celltype_names[j]
      
      # Get common isoforms
      common_isos <- intersect(rownames(correlation_list[[ct1]]), 
                              rownames(correlation_list[[ct2]]))
      
      if (length(common_isos) >= 2) {
        # Extract common submatrices
        cor1 <- correlation_list[[ct1]][common_isos, common_isos]
        cor2 <- correlation_list[[ct2]][common_isos, common_isos]
        
        # Vectorize upper triangles
        vec1 <- cor1[upper.tri(cor1)]
        vec2 <- cor2[upper.tri(cor2)]
        
        # Calculate distance
        if (method == "correlation") {
          dist_val <- 1 - cor(vec1, vec2, use = "complete.obs")
        } else if (method == "euclidean") {
          dist_val <- sqrt(mean((vec1 - vec2)^2, na.rm = TRUE))
        } else if (method == "manhattan") {
          dist_val <- mean(abs(vec1 - vec2), na.rm = TRUE)
        }
        
        dist_matrix[i, j] <- dist_matrix[j, i] <- dist_val
      }
    }
  }
  
  # Set diagonal to 0
  diag(dist_matrix) <- 0
  
  return(as.dist(dist_matrix))
}

#' Export co-expression analysis results
#'
#' @description
#' Exports comprehensive co-expression analysis results to files.
#'
#' @param coexpr_result Co-expression analysis results
#' @param output_prefix Prefix for output files
#' @param formats Export formats: "csv", "xlsx", "rds"
#'
#' @export
export_coexpression_results <- function(coexpr_result, 
                                      output_prefix,
                                      formats = c("csv", "rds")) {
  
  # Create output directory if needed
  output_dir <- dirname(output_prefix)
  if (output_dir != "." && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Export based on result type
  if ("overall" %in% names(coexpr_result)) {
    # Results from calculate_gene_coexpression_all_celltypes
    
    # Export overall correlation
    if ("csv" %in% formats) {
      write.csv(coexpr_result$overall$cor_matrix, 
               paste0(output_prefix, "_overall_correlation.csv"))
    }
    
    # Export cell type specific
    for (ct in names(coexpr_result$cell_types)) {
      if ("csv" %in% formats) {
        write.csv(coexpr_result$cell_types[[ct]]$cor_matrix,
                 paste0(output_prefix, "_", ct, "_correlation.csv"))
      }
    }
  } else if ("pairs" %in% names(coexpr_result)) {
    # Results from conservation analysis
    
    # Convert to data frame
    conservation_df <- do.call(rbind, lapply(names(coexpr_result$pairs), function(pair) {
      x <- coexpr_result$pairs[[pair]]
      data.frame(
        pair = pair,
        isoform1 = x$isoform1,
        isoform2 = x$isoform2,
        mean_correlation = x$mean_correlation,
        sd_correlation = x$sd_correlation,
        pattern = x$pattern,
        n_cell_types = x$n_cell_types,
        stringsAsFactors = FALSE
      )
    }))
    
    if ("csv" %in% formats) {
      write.csv(conservation_df, 
               paste0(output_prefix, "_conservation_analysis.csv"),
               row.names = FALSE)
    }
  }
  
  # Export full object as RDS
  if ("rds" %in% formats) {
    saveRDS(coexpr_result, paste0(output_prefix, "_complete_results.rds"))
  }
  
  message("Results exported to: ", output_prefix, "_*")
}

#' Create a summary visualisation of co-expression conservation
#'
#' @description
#' Creates bar plots and visualisations summarising conservation patterns
#' across cell types.
#'
#' @param conservation_results Results from analyse_coexpression_conservation
#' @param output_file Output filename (default: "conservation_summary.pdf")
#' @param width Plot width in inches
#' @param height Plot height in inches
#'
#' @return Invisible NULL
#' @export
plot_conservation_summary <- function(conservation_results, 
                                    output_file = "conservation_summary.pdf",
                                    width = 8,
                                    height = 6) {
  
  # Handle both direct conservation results and formatted data frames
  if ("pairs" %in% names(conservation_results)) {
    # Results from analyse_coexpression_conservation function
    pattern_counts <- table(sapply(conservation_results$pairs, function(x) x$pattern))
    
    # Create data frame for plotting
    plot_data <- do.call(rbind, lapply(names(conservation_results$pairs), function(pair_name) {
      pair_data <- conservation_results$pairs[[pair_name]]
      data.frame(
        isoform_pair = pair_name,
        mean_correlation = pair_data$mean_correlation,
        sd_correlation = pair_data$sd_correlation,
        conservation_pattern = pair_data$pattern,
        stringsAsFactors = FALSE
      )
    }))
  } else if (is.data.frame(conservation_results)) {
    # Data frame format (e.g., from Shiny app)
    pattern_counts <- table(conservation_results$conservation_pattern)
    plot_data <- conservation_results
  } else {
    stop("Invalid conservation_results format")
  }
  
  # Close any existing devices first
  if (length(dev.list()) > 0) {
    dev.off()
  }
  
  pdf(output_file, width = width, height = height)
  
  tryCatch({
    par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))
    
    # Plot 1: Pattern distribution
    barplot(pattern_counts,
            main = "Co-expression Conservation Patterns",
            ylab = "Number of isoform pairs",
            col = c("Conserved_Positive" = "#6fafd2",
                    "Conserved_Negative" = "#feab88",
                    "Mixed" = "grey70",
                    "Cell_Type_Specific" = "#6a0624")[names(pattern_counts)],
            las = 2)
    
    # Plot 2: Top conserved pairs
    conserved_patterns <- c("Conserved_Positive", "Conserved_Negative")
    top_conserved <- plot_data[plot_data$conservation_pattern %in% conserved_patterns, ]
    
    if (nrow(top_conserved) > 0) {
      top_conserved <- head(top_conserved[order(abs(top_conserved$mean_correlation), 
                                              decreasing = TRUE), ], 10)
      
      plot(top_conserved$mean_correlation, 1:nrow(top_conserved),
           xlim = c(-1, 1),
           yaxt = "n",
           xlab = "Mean correlation",
           ylab = "",
           main = "Top Conserved Co-expression Pairs",
           pch = 19,
           col = ifelse(top_conserved$mean_correlation > 0, "#6fafd2", "#feab88"))
      
      axis(2, at = 1:nrow(top_conserved), 
           labels = top_conserved$isoform_pair,
           las = 2, cex.axis = 0.7)
      
      abline(v = 0, lty = 2, col = "grey50")
      
      # Add error bars (SD)
      segments(top_conserved$mean_correlation - top_conserved$sd_correlation,
               1:nrow(top_conserved),
               top_conserved$mean_correlation + top_conserved$sd_correlation,
               1:nrow(top_conserved))
    }
    
  }, error = function(e) {
    message("Error in plotting: ", e$message)
  }, finally = {
    dev.off()  # Always close the device
  })
  
  message(paste("Conservation summary saved to:", output_file))
  invisible(NULL)
}