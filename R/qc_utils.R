##############################################################
#  Author: [Siyuan Wu & Ulf Schmitz]                         #
#  Institution: [James Cook University]                      #
#  Date: Jul 29, 2025                                        #
#  Package: ScIsoX V1.1.0                                    #
##############################################################

#' Plot and analyse genes per cell distribution
#'
#' @description
#' Creates a detailed visualisation of genes per cell distribution and
#' calculates suggested quality control parameters. Optimised for both
#' dense and sparse matrices to accommodate large-scale single-cell datasets.
#'
#' @param gene_counts A matrix of gene counts (rows = genes, columns = cells).
#'        Supports both standard and sparse matrices from the Matrix package.
#' @param plot_type Character, one of "hist" (histogram) or "density" (density plot)
#' @param percentile_cutoffs Numeric vector of percentiles to show as cutoff lines
#' @param return_suggestions Logical, whether to return suggested QC parameters
#'
#' @return If return_suggestions=TRUE, returns a list containing:
#'   \itemize{
#'     \item min_genes_per_cell: Suggested minimum genes per cell
#'     \item max_genes_per_cell: Suggested maximum genes per cell
#'     \item median_genes: Median number of genes per cell
#'     \item summary_stats: Additional statistics
#'   }
#'
#' @examples
#' # Load example data
#' data(gene_counts_blood)
#' 
#' # Create histogram plot with default settings
#' qc_suggestions <- plot_genes_per_cell_distribution(
#'   gene_counts = gene_counts_blood,
#'   plot_type = "hist",
#'   percentile_cutoffs = c(0.05, 0.95),
#'   return_suggestions = TRUE
#' )
#' 
#' # View QC suggestions
#' print(qc_suggestions$min_genes_per_cell)
#' print(qc_suggestions$max_genes_per_cell)
#' 
#' # Create density plot
#' plot_genes_per_cell_distribution(
#'   gene_counts = gene_counts_blood,
#'   plot_type = "density",
#'   percentile_cutoffs = c(0.1, 0.9),
#'   return_suggestions = FALSE
#' )
#'
#' @export
plot_genes_per_cell_distribution <- function(gene_counts,
                                             plot_type = "hist",
                                             percentile_cutoffs = c(0.05, 0.95),
                                             return_suggestions = TRUE) {
  # Detect matrix type and process accordingly for memory efficiency
  is_sparse <- inherits(gene_counts, "sparseMatrix") ||
               inherits(gene_counts, "dgCMatrix") ||
               inherits(gene_counts, "dgTMatrix") ||
               inherits(gene_counts, "dgRMatrix")
  
  # Calculate genes per cell with optimised method based on matrix type
  if (is_sparse) {
    # For sparse matrices, use the Matrix package's specialised methods
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("The 'Matrix' package is required for sparse matrix support. Please install with: install.packages('Matrix')")
    }
    n_genes <- Matrix::colSums(gene_counts != 0)
  } else {
    # For standard matrices, use base R methods
    n_genes <- colSums(gene_counts > 0)
  }
  
  # Set up plotting parameters for optimal visualisation
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 4))
  
  # Create main plot based on user preference
  if (plot_type == "density") {
      plot(density(n_genes),
           main = "Density Distribution of Genes per Cell",
           xlab = "Number of Genes",
           ylab = "Density",
           col = "blue")
      polygon(density(n_genes), col = "lightblue", border = "blue")
  } else if (plot_type == "hist") {
      hist(n_genes,
           breaks = min(100, length(unique(n_genes))),
           main = "Distribution of Genes per Cell",
           xlab = "Number of Genes",
           ylab = "Number of Cells",
           col = "lightblue",
           border = "white")
  }
  
  # Add reference lines to highlight important distribution features
  median_genes <- median(n_genes)
  abline(v = median_genes, col = "red", lty = 2)
  text(median_genes, par("usr")[4],
       sprintf("Median: %d", round(median_genes)),
       pos = 3, col = "red")
  
  # Add percentile reference lines to indicate suggested filtering thresholds
  for (p in percentile_cutoffs) {
    q <- quantile(n_genes, p)
    abline(v = q, col = "blue", lty = 3)
  }
  
  # Add comprehensive summary statistics for dataset characterisation
  mtext(sprintf("5th: %.1f | Mean: %.1f | Median: %.1f | SD: %.1f | CV: %.2f | 95th: %.1f",
                quantile(n_genes, 0.05),
                mean(n_genes),
                median(n_genes),
                sd(n_genes),
                sd(n_genes)/mean(n_genes),
        quantile(n_genes, 0.95)), side = 3, line = 0, cex = 0.8)
  
  # Calculate suggested parameters for downstream quality control
  if (return_suggestions) {
    suggestions <- list(
      min_genes_per_cell = as.integer(quantile(n_genes, 0.05)),
      max_genes_per_cell = as.integer(quantile(n_genes, 0.95)),
      median_genes = as.integer(median_genes),
      summary_stats = list(
        mean = mean(n_genes),
        sd = sd(n_genes),
        cv = sd(n_genes)/mean(n_genes),
        quantiles = quantile(n_genes, probs = seq(0, 1, 0.1))
      )
    )
    return(suggestions)
  }
}

#' Generate QC parameter recommendations
#'
#' @description
#' Analyses the gene count distribution and provides detailed recommendations
#' for quality control parameter settings using different stringency levels.
#' Supports both dense and sparse matrices for efficient processing of
#' large-scale single-cell datasets.
#'
#' @param gene_counts A matrix of gene counts. Supports both standard and
#'        sparse matrices from the Matrix package.
#'
#' @return A list containing QC recommendations:
#'   \itemize{
#'     \item MAD_strategy: Conservative parameter settings based on median absolute deviation
#'     \item Interval_90: Moderate parameter settings using 5th and 95th percentiles
#'     \item Interval_80: Aggressive parameter settings using 10th and 90th percentiles
#'     \item explanation: Descriptions of each approach
#'   }
#'
#' @examples
#' # Load example data
#' data(gene_counts_blood)
#' 
#' # Get QC recommendations
#' qc_recommendations <- recommend_qc_parameters(gene_counts_blood)
#' 
#' # View different strategies
#' print("MAD Strategy (Conservative):")
#' print(qc_recommendations$MAD_strategy)
#' 
#' print("\n90% Interval Strategy (Moderate):")
#' print(qc_recommendations$Interval_90)
#' 
#' print("\n80% Interval Strategy (Aggressive):")
#' print(qc_recommendations$Interval_80)
#' 
#' # Use recommended parameters in create_scht
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' scht_obj <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   qc_params = c(
#'     qc_recommendations$MAD_strategy,
#'     list(
#'       min_cells_expressing = 0.02,
#'       min_expr = 1e-6
#'     )
#'   ),
#'   n_hvg = 3000,
#'   verbose = FALSE
#' )
#'
#' @export
recommend_qc_parameters <- function(gene_counts) {
  # Detect matrix type for memory-efficient processing
  is_sparse <- inherits(gene_counts, "sparseMatrix") ||
               inherits(gene_counts, "dgCMatrix") ||
               inherits(gene_counts, "dgTMatrix") ||
               inherits(gene_counts, "dgRMatrix")
  
  # Calculate number of genes per cell with optimised method based on matrix type
  if (is_sparse) {
    # For sparse matrices, use the Matrix package's specialised methods
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("The 'Matrix' package is required for sparse matrix support. Please install with: install.packages('Matrix')")
    }
    cell_gene_counts <- Matrix::colSums(gene_counts != 0)
  } else {
    # For standard matrices, use base R methods
    cell_gene_counts <- colSums(gene_counts > 0)
  }
  
  # Get distribution analysis using our dedicated function
  dist_analysis <- plot_genes_per_cell_distribution(gene_counts,
                                                    return_suggestions = TRUE)
  
  # Calculate MAD-based bounds for robust outlier detection
  median_genes <- dist_analysis$median_genes
  genes_mad <- mad(cell_gene_counts)
  MAD_strategy_min <- floor(median_genes - 3*genes_mad)
  MAD_strategy_max <- ceiling(median_genes + 3*genes_mad)
  
  # Generate multi-tier recommendations to accommodate different filtering strategies
  recommendations <- list(
    MAD_strategy = list(
      min_genes_per_cell = max(100, MAD_strategy_min),  # Lower bound to ensure minimal cell quality
      max_genes_per_cell = MAD_strategy_max
    ),
    Interval_90 = list(
      min_genes_per_cell = dist_analysis$min_genes_per_cell,
      max_genes_per_cell = dist_analysis$max_genes_per_cell
    ),
    Interval_80 = list(
      min_genes_per_cell = as.integer(quantile(cell_gene_counts, 0.1)),
      max_genes_per_cell = as.integer(quantile(cell_gene_counts, 0.9))
    ),
    explanation = c(
      "MAD_strategy: Uses median +/- 3 MAD, reduces risk of including poor quality cells whilst maintaining robustness to outliers",
      "Interval_90: Uses 5th and 95th percentiles, balances stringency with dataset preservation",
      "Interval_80: Uses 10th and 90th percentiles, provides more aggressive filtering for higher quality cell selection"
    )
  )
  
  return(recommendations)
}
