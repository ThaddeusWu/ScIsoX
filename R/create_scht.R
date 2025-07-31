##############################################################
#  Single-Cell Hierarchical Tensor (SCHT) Creation Pipeline  #
#  Enhanced version with normalization options and           #
#  comprehensive QC tracking                                 #                     
#                                                            #
#  Author: [Siyuan Wu & Ulf Schmitz]                         #
#  Institution: [James Cook University]                      #
#  Date: Jul 29, 2025                                        #
#  Package: ScIsoX V1.1.0                                    #
##############################################################

# This implementation supports efficient handling of large-scale
# single-cell long-read data with reduced memory footprint.
# Includes support for pre-normalized data and enhanced QC tracking.

########################
# Required Libraries   #
########################
#' @importFrom graphics abline hist mtext par polygon text
#' @importFrom stats density mad median na.omit quantile sd var setNames
#' @importFrom progress progress_bar
#' @importFrom methods is
#' @importFrom Matrix Matrix Diagonal
#' @importFrom utils packageVersion
NULL

##########################
# Core Utility Functions #
##########################

#' Process input data with matrix conversion
#'
#' @description
#' Validates and converts input data into an appropriate matrix format
#' for memory-efficient processing of large-scale count matrices.
#'
#' @param gene_counts Gene-level counts as matrix, data frame, or sparse matrix 
#' @param transcript_counts Transcript-level counts as matrix, data frame, or sparse matrix
#' @param transcript_info Data frame with transcript annotations
#' @param cell_info Optional data frame describing cells
#' @param n_hvg Number of highly variable genes to select
#' @param qc_params QC parameters list
#' @param require_cell_type Whether cell type information is required
#' @param sparsity_threshold Minimum sparsity to convert to sparse format (0-1)
#'
#' @return A list containing validated matrix inputs
#'
#' @keywords internal
.process_input_data <- function(gene_counts,
                                transcript_counts,
                                transcript_info,
                                cell_info = NULL,
                                n_hvg,
                                qc_params,
                                require_cell_type = TRUE,
                                sparsity_threshold = 0.4) {
  
  # Validate numeric parameters
  if (n_hvg < 1) stop("n_hvg must be positive")
  if (qc_params$min_genes_per_cell < 1) stop("min_genes_per_cell must be positive")
  if (qc_params$max_genes_per_cell < 1) stop("max_genes_per_cell must be positive")
  if (qc_params$min_cells_expressing <= 0 || qc_params$min_cells_expressing > 1) 
    stop("min_cells_expressing must be between 0 and 1")
  if (qc_params$min_expr < 0) stop("min_expr must be non-negative")
  
  # Efficient conversion to appropriate matrix format
  convert_to_matrix_format <- function(mat, name) {
    if (is(mat, "sparseMatrix")) {
      return(mat) # Already in sparse format
    } else if (is.data.frame(mat)) {
      # Check if numeric
      sample_cols <- min(5, ncol(mat))
      if (!all(sapply(mat[, 1:sample_cols, drop=FALSE], is.numeric))) {
        stop(paste("All columns in", name, "must be numeric"))
      }
      
      # Convert data frame to matrix first
      mat <- as.matrix(mat)
    } else if (!is.matrix(mat)) {
      stop(paste(name, "must be a matrix, data frame, or sparse matrix"))
    }
    
    # Count zeros to determine sparsity
    if (length(mat) > 1e6) { 
      # For very large matrices, estimate sparsity from a sample
      sample_size <- min(1e6, floor(length(mat) * 0.1))
      indices <- sample(length(mat), sample_size)
      zeros_ratio <- sum(mat[indices] == 0) / sample_size
    } else {
      zeros_ratio <- sum(mat == 0) / length(mat)
    }
    
    # Convert to sparse if sparsity exceeds threshold
    if (zeros_ratio >= sparsity_threshold) {
      # Use optimal sparse matrix format based on dimensions
      if (nrow(mat) > 1e4 || ncol(mat) > 1e4) {
        return(Matrix::Matrix(mat, sparse = TRUE))
      } else {
        return(Matrix::Matrix(mat, sparse = TRUE))
      }
    } else {
      # Keep as regular matrix if not sparse enough
      return(mat)
    }
  }
  
  # Convert inputs to appropriate matrix format
  gene_counts <- convert_to_matrix_format(gene_counts, "gene_counts")
  transcript_counts <- convert_to_matrix_format(transcript_counts, "transcript_counts")
  
  # Check transcript_info
  if (!is.data.frame(transcript_info)) {
    stop("'transcript_info' must be a data frame.")
  }
  
  # Check cell_info if provided
  if (!is.null(cell_info)) {
    if (!is.data.frame(cell_info)) {
      stop("When provided, 'cell_info' must be a data frame.")
    }
    
    # Check for cell type information
    if (require_cell_type && !"cell_type" %in% colnames(cell_info)) {
      stop("'cell_info' must contain a 'cell_type' column if require_cell_type = TRUE.")
    }
  } else if (require_cell_type) {
    stop("cell_info must be provided when require_cell_type = TRUE.")
  }
  
  # Check transcript_info for necessary columns
  required_cols <- c("transcript_id", "transcript_name", "gene_id", "gene_name")
  missing_cols <- setdiff(required_cols, colnames(transcript_info))
  if (length(missing_cols) > 0) {
    stop("transcript_info must contain columns: ",
         paste(required_cols, collapse = ", "))
  }
  
  # Return the processed inputs
  list(
    gene_counts = gene_counts,
    transcript_counts = transcript_counts,
    transcript_info = transcript_info,
    cell_info = cell_info
  )
}

#' Perform quality control on genes and transcripts
#'
#' @description
#' Performs initial quality control using matrix 
#' operations to filter low-quality genes and transcripts.
#'
#' @param gene_counts Gene-level counts (matrix)
#' @param transcript_counts Transcript-level counts (matrix)
#' @param transcript_info Data frame with transcript annotations
#' @param qc_params QC parameters list
#'
#' @return A list containing filtered data and QC metrics
#'
#' @keywords internal
.perform_qc <- function(gene_counts,
                        transcript_counts,
                        transcript_info,
                        qc_params) {
  
  # Calculate minimum cells threshold
  min_cells_threshold <- qc_params$min_cells_expressing * ncol(gene_counts)
  
  # For sparse matrices, use optimised methods
  if (is(gene_counts, "sparseMatrix")) {
    # Get non-zero counts per gene
    n_cells_per_gene <- Matrix::rowSums(gene_counts > 0)
    mean_expr_per_gene <- Matrix::rowMeans(gene_counts)
  } else {
    # For dense matrices use base functions
    n_cells_per_gene <- rowSums(gene_counts > 0)
    mean_expr_per_gene <- rowMeans(gene_counts)
  }
  
  # Combined filtering with single logical operation
  poor_genes <- (n_cells_per_gene < min_cells_threshold) | (mean_expr_per_gene < qc_params$min_expr)
  
  # Similar approach for transcripts
  if (is(transcript_counts, "sparseMatrix")) {
    n_cells_per_transcript <- Matrix::rowSums(transcript_counts > 0)
    mean_expr_per_transcript <- Matrix::rowMeans(transcript_counts)
  } else {
    n_cells_per_transcript <- rowSums(transcript_counts > 0)
    mean_expr_per_transcript <- rowMeans(transcript_counts)
  }
  
  # Direct NA checking
  poor_trans <- is.na(n_cells_per_transcript)
  
  # Filter matrices efficiently
  gene_counts_filtered <- gene_counts[!poor_genes, , drop = FALSE]
  transcript_counts_filtered <- transcript_counts[!poor_trans, , drop = FALSE]
  transcript_info_filtered <- transcript_info[!poor_trans, , drop = FALSE]
  
  # Return filtered data and QC metrics
  list(
    gene_counts_filtered = gene_counts_filtered,
    transcript_counts_filtered = transcript_counts_filtered,
    transcript_info_filtered = transcript_info_filtered,
    n_filtered_genes = sum(poor_genes),
    n_filtered_transcripts = sum(poor_trans),
    qc_metrics = list(
      genes = data.frame(
        n_cells = n_cells_per_gene,
        mean_expr = mean_expr_per_gene
      ),
      transcripts = data.frame(
        n_cells = n_cells_per_transcript,
        mean_expr = mean_expr_per_transcript
      )
    )
  )
}

#' Perform cell-level quality control
#'
#' @description
#' Identifies and removes low-quality cells using matrix operations
#' for efficient processing with large datasets.
#'
#' @param gene_counts Gene-level counts (matrix)
#' @param transcript_counts Transcript-level counts (matrix)
#' @param qc_params QC parameters list
#'
#' @return A list containing filtered data and QC metrics
#'
#' @keywords internal
.perform_cell_qc <- function(gene_counts,
                             transcript_counts,
                             qc_params) {
  
  # Vectorised calculation of genes per cell
  if (is(gene_counts, "sparseMatrix")) {
    n_genes_per_cell <- Matrix::colSums(gene_counts > 0)
  } else {
    n_genes_per_cell <- base::colSums(gene_counts > 0)
  }
  
  # Filter cells using vectorised operations
  poor_cells <- (n_genes_per_cell < qc_params$min_genes_per_cell) | 
    (n_genes_per_cell > qc_params$max_genes_per_cell)
  
  # Use logical indexing for matrices
  good_cells <- !poor_cells
  
  list(
    gene_counts_filtered = gene_counts[, good_cells, drop = FALSE],
    transcript_counts_filtered = transcript_counts[, good_cells, drop = FALSE],
    n_filtered_cells = sum(poor_cells),
    n_cells = sum(good_cells),
    qc_metrics = data.frame(
      n_genes = n_genes_per_cell
    )
  )
}

#' Perform normalisation
#'
#' @description
#' Normalises count data to counts per million (CPM) using
#' operations to maintain memory efficiency.
#'
#' @param gene_counts_filtered Filtered gene counts
#' @param transcript_counts_filtered Filtered transcript counts
#'
#' @return A list containing normalised data
#'
#' @keywords internal
.perform_normalisation <- function(gene_counts_filtered, 
                                   transcript_counts_filtered) {
  
  # Calculate scaling factors based on matrix type
  if (is(gene_counts_filtered, "sparseMatrix")) {
    gene_col_sums <- Matrix::colSums(gene_counts_filtered)
    transcript_col_sums <- Matrix::colSums(transcript_counts_filtered)
  } else {
    gene_col_sums <- base::colSums(gene_counts_filtered)
    transcript_col_sums <- base::colSums(transcript_counts_filtered)
  }
  
  gene_scaling_factors <- 1e6 / gene_col_sums
  transcript_scaling_factors <- 1e6 / transcript_col_sums
  
  # Store column names before operations
  gene_col_names <- colnames(gene_counts_filtered)
  transcript_col_names <- colnames(transcript_counts_filtered)
  
  
  # Use appropriate operations based on matrix type
  if (is(gene_counts_filtered, "sparseMatrix")) {
    # Column scaling for sparse matrices
    gene_counts_norm <- gene_counts_filtered %*% Diagonal(x = gene_scaling_factors)
    transcript_counts_norm <- transcript_counts_filtered %*% Diagonal(x = transcript_scaling_factors)
    
    # Restore column names
    colnames(gene_counts_norm) <- gene_col_names
    colnames(transcript_counts_norm) <- transcript_col_names
  } else {
    # For dense matrices
    gene_counts_norm <- sweep(gene_counts_filtered, 2, gene_scaling_factors, "*")
    transcript_counts_norm <- sweep(transcript_counts_filtered, 2, transcript_scaling_factors, "*")
  }
  
  list(
    gene_counts_norm = gene_counts_norm,
    transcript_counts_norm = transcript_counts_norm
  )
}

#' Perform log transformation
#'
#' @description
#' Applies log2(x+1) transformation to normalised count data with
#' special handling for different matrix types.
#'
#' @param gene_counts_norm Normalised gene counts
#' @param transcript_counts_norm Normalised transcript counts
#'
#' @return A list containing transformed data
#'
#' @keywords internal
.perform_log_transform <- function(gene_counts_norm,
                                   transcript_counts_norm,
                                   verbose) {
  
  gene_mat <- gene_counts_norm
  transcript_mat <- transcript_counts_norm
  
  # Handle different matrix types appropriately
  if (is(gene_counts_norm, "sparseMatrix")) {
    # log2(x+1) for sparse matrices efficiently preserves sparsity
    gene_mat@x <- log2(gene_mat@x + 1)
  } else {
    # For dense matrices
    gene_mat <- log2(gene_mat + 1)
  }
  
  if (is(transcript_mat, "sparseMatrix")) {
    transcript_mat@x <- log2(transcript_mat@x + 1)
  } else {
    transcript_mat <- log2(transcript_mat + 1)
  }
  
  if (verbose) {
    cat(
      "Matrix format after log transform:\n",
      "  class(gene_mat) =", class(gene_mat), "\n",
      "  dim(gene_mat)   =", paste(dim(gene_mat), collapse=", "), "\n",
      "  class(transcript_mat) =", class(transcript_mat), "\n",
      "  dim(transcript_mat)   =", paste(dim(transcript_mat), collapse=", "), "\n"
    )
  }
  
  list(
    gene_mat = gene_mat,
    transcript_mat = transcript_mat
  )
}

#' Select highly variable genes using second-moment approach
#'
#' @description
#' Identifies highly variable genes using matrix operations
#' for efficient analysis with large datasets. Variance is computed
#' by subtracting the square of the mean from the mean of squares.
#'
#' @param gene_mat Normalised gene expression matrix.
#' @param n_hvg Number of HVGs to select.
#' @param qc_params QC parameters list; must contain \code{min_cells_expressing}.
#' @param verbose Logical indicating whether to display progress messages.
#'
#' @return A list containing selected HVGs and metrics:
#'   \describe{
#'     \item{\code{var_genes}}{Character vector of selected HVG names.}
#'     \item{\code{dispersion}}{Named numeric vector of dispersion values (variance-to-mean ratio).}
#'     \item{\code{mean_expr}}{Named numeric vector of mean expression values for the selected genes.}
#'   }
#'
#' @keywords internal
.select_hvgs <- function(gene_mat,
                         n_hvg,
                         qc_params,
                         verbose = TRUE) {
  
  ## Validate the structure of the input matrix
  if (is.null(dim(gene_mat))) {
    stop("Input gene_mat must be a matrix with proper dimensions.")
  }
  
  ## Handle the empty matrix case
  if (nrow(gene_mat) == 0) {
    if (verbose) {
      message("Empty gene matrix provided. No highly variable genes to select.")
    }
    return(list(
      var_genes = character(0),
      dispersion = numeric(0),
      mean_expr = numeric(0)
    ))
  }
  
  ## Determine the minimum number of cells in which a gene must be expressed
  min_cells <- qc_params$min_cells_expressing * ncol(gene_mat)
  
  ## Compute the number of cells expressing each gene
  if (methods::is(gene_mat, "sparseMatrix")) {
    expr_cells_count <- Matrix::rowSums(gene_mat > 0)
  } else {
    expr_cells_count <- rowSums(gene_mat > 0)
  }
  
  ## Filter genes based on the minimum expression threshold
  keep_genes <- expr_cells_count >= min_cells
  
  ## If no genes pass the filter, select at least one fallback gene
  if (sum(keep_genes) == 0) {
    if (verbose) {
      message("Warning: No genes met the minimum cell expression threshold.")
      message("Using the gene with the highest number of expressing cells instead.")
    }
    top_idx <- which.max(expr_cells_count)
    keep_genes[top_idx] <- TRUE
  }
  
  ## Subset the matrix to only include the filtered genes
  gene_mat_filtered <- gene_mat[keep_genes, , drop = FALSE]
  
  if (verbose) {
    message("Selecting top highly variable genes (HVGs)...")
    message("Number of genes considered after min cell filter: ", nrow(gene_mat_filtered))
  }
  
  ## If only one gene is retained, return it directly with trivial dispersion
  if (nrow(gene_mat_filtered) == 1) {
    if (verbose) {
      message("Only one gene passed filtering. Returning it as the sole HVG.")
    }
    
    gene_name <- rownames(gene_mat_filtered)
    
    ## Compute the mean expression for this single gene
    mean_val <- if (methods::is(gene_mat_filtered, "sparseMatrix")) {
      as.numeric(Matrix::rowMeans(gene_mat_filtered))
    } else {
      as.numeric(rowMeans(gene_mat_filtered))
    }
    
    ## Return a structured result
    return(list(
      var_genes = gene_name,
      dispersion = stats::setNames(0, gene_name),  # Dispersion is zero when there's only one gene
      mean_expr = stats::setNames(mean_val, gene_name)
    ))
  }
  
  ## Compute row-wise means
  if (methods::is(gene_mat_filtered, "sparseMatrix")) {
    row_means <- Matrix::rowMeans(gene_mat_filtered)
  } else {
    row_means <- rowMeans(gene_mat_filtered)
  }
  
  ## Compute row-wise means of squares: E(X^2)
  if (methods::is(gene_mat_filtered, "sparseMatrix")) {
    # The '^2' operation is sparse-aware in Matrix package
    row_sq_means <- Matrix::rowMeans(gene_mat_filtered ^ 2)
  } else {
    row_sq_means <- rowMeans(gene_mat_filtered ^ 2)
  }
  
  ## Calculate variance from the second moment minus the square of the first moment
  row_vars <- row_sq_means - (row_means ^ 2)
  
  ## Compute the dispersion as the ratio of variance to mean
  dispersion <- row_vars / row_means
  dispersion[is.na(dispersion) | is.infinite(dispersion)] <- 0
  
  ## Name the dispersion values with the corresponding gene IDs
  names(dispersion) <- rownames(gene_mat_filtered)
  
  ## Rank genes by dispersion and select the top HVGs
  ordered_indices <- order(dispersion, decreasing = TRUE)
  top_n <- min(n_hvg, length(dispersion))
  top_indices <- ordered_indices[1:top_n]
  top_genes <- names(dispersion)[top_indices]
  
  if (verbose) {
    message("Selected ", length(top_genes), " highly variable genes.")
  }
  
  ## Return the results with gene names and associated metrics
  list(
    var_genes = top_genes,
    dispersion = dispersion[top_genes],
    mean_expr = row_means[top_indices]
  )
}


#' Create SCHT structure
#'
#' @description
#' Creates Single-Cell Hierarchical Tensor structure with efficient
#' data handling for large datasets.
#'
#' @param transcript_mat_final Normalised transcript expression
#' @param var_genes Selected highly variable genes
#' @param transcript_info_filtered Filtered transcript information
#' @param verbose Whether to show progress
#'
#' @return A list containing the SCHT structure
#'
#' @keywords internal
.create_SCHT <- function(transcript_mat_final,
                         var_genes,
                         transcript_info_filtered,
                         verbose = TRUE) {
  
  if (verbose) {
    message("Filtering transcripts to only those belonging to HVGs...")
    message(sprintf("Number of input var_genes: %d", length(var_genes)))
    message(sprintf("Number of unique gene names in transcript_info: %d", 
                    length(unique(transcript_info_filtered$gene_name))))
  }
  
  # Create lookup tables for performance
  gene_to_trans <- split(
    transcript_info_filtered$transcript_id, 
    transcript_info_filtered$gene_id
  )
  
  gene_id_to_name <- structure(
    transcript_info_filtered$gene_name,
    names = transcript_info_filtered$gene_id
  )
  gene_id_to_name <- gene_id_to_name[!duplicated(names(gene_id_to_name))]
  
  # Also create reverse lookup for gene names to IDs
  gene_name_to_id <- structure(
    names(gene_id_to_name),
    names = gene_id_to_name
  )
  
  # Create gene name to transcript lookup as well
  gene_name_to_trans <- split(
    transcript_info_filtered$transcript_id,
    transcript_info_filtered$gene_name
  )
  
  # Pre-allocate result list
  SCHT <- vector("list", length(var_genes))
  names_vec <- character(length(var_genes))
  
  # Track HVG isoform counts
  hvg_single_iso <- character()
  hvg_multi_iso <- character()
  
  # First pass: collect all gene names to identify conflicts
  gene_name_mapping <- list()
  for (i in seq_along(var_genes)) {
    current_gene <- var_genes[i]
    
    # Determine gene name
    gene_transcripts <- gene_to_trans[[current_gene]]
    gene_name_selected <- gene_id_to_name[current_gene]
    
    if (is.null(gene_transcripts)) {
      gene_transcripts <- gene_name_to_trans[[current_gene]]
      gene_name_selected <- current_gene
    }
    
    if (!is.null(gene_transcripts)) {
      # Store mapping
      if (gene_name_selected %in% names(gene_name_mapping)) {
        gene_name_mapping[[gene_name_selected]] <- c(gene_name_mapping[[gene_name_selected]], current_gene)
      } else {
        gene_name_mapping[[gene_name_selected]] <- current_gene
      }
    }
  }
  
  # Identify conflicts (gene names with multiple IDs)
  conflicted_gene_names <- names(gene_name_mapping)[sapply(gene_name_mapping, length) > 1]
  if (length(conflicted_gene_names) > 0 && verbose) {
    message(sprintf("Found %d gene names with multiple IDs, will use gene IDs for these", 
                   length(conflicted_gene_names)))
  }
  
  # Progress tracking
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "  Processing genes [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
      total = length(var_genes),
      clear = FALSE,
      width = 60
    )
  }
  
  # Get transcript rows once
  if (is(transcript_mat_final, "sparseMatrix")) {
    transcript_rows <- rownames(transcript_mat_final)
  } else {
    transcript_rows <- rownames(transcript_mat_final)
  }
  
  # Check if var_genes are gene IDs or gene names
  n_as_ids <- sum(var_genes %in% names(gene_to_trans))
  n_as_names <- sum(var_genes %in% names(gene_name_to_trans))
  
  if (verbose && n_as_names > n_as_ids) {
    message("Note: Input appears to use gene names rather than gene IDs. ",
            "Matching by gene names.")
  }
  
  # Process each gene
  for (i in seq_along(var_genes)) {
    current_gene <- var_genes[i]
    
    # Determine if current_gene is a gene ID or gene name
    # First try as gene ID
    gene_transcripts <- gene_to_trans[[current_gene]]
    gene_name_selected <- gene_id_to_name[current_gene]
    
    # If not found as gene ID, try as gene name
    if (is.null(gene_transcripts)) {
      gene_transcripts <- gene_name_to_trans[[current_gene]]
      gene_name_selected <- current_gene
    }
    
    if (!is.null(gene_transcripts)) {
      # Find transcripts present in data
      existing_transcripts <- intersect(gene_transcripts, transcript_rows)
      
      # Set the name: use gene ID if there's a conflict, otherwise use gene name
      if (gene_name_selected %in% conflicted_gene_names) {
        names_vec[i] <- current_gene  # Use gene ID for conflicted names
      } else {
        names_vec[i] <- gene_name_selected  # Use gene name for unique names
      }
      
      if (length(existing_transcripts) > 0) {
        # Extract transcript data
        result_data <- transcript_mat_final[existing_transcripts, , drop = FALSE]
        
        # Check for NAs and filter
        if (is(result_data, "sparseMatrix")) {
          # Sparse matrices handle NAs differently
          if (any(is.na(result_data@x))) {
            # Handle NAs in sparse matrix
            result_data <- na.omit(as.matrix(result_data))
            # Convert back to sparse if needed
            if (length(result_data) > 1000) {
              result_data <- Matrix(result_data, sparse = TRUE)
            }
          }
          
          # Calculate row sums
          row_sums <- Matrix::rowSums(result_data)
          non_zero_transcripts <- result_data[row_sums > 0, , drop = FALSE]
        } else {
          # For dense matrices
          result_data <- na.omit(result_data)
          row_sums <- rowSums(result_data)
          non_zero_transcripts <- result_data[row_sums > 0, , drop = FALSE]
        }
        
        # Track single vs multiple isoforms
        if (nrow(non_zero_transcripts) == 1) {
          hvg_single_iso <- c(hvg_single_iso, gene_name_selected)
        } else if (nrow(non_zero_transcripts) > 1) {
          SCHT[[i]] <- non_zero_transcripts
          hvg_multi_iso <- c(hvg_multi_iso, gene_name_selected)
        }
      }
    }
    
    # Update progress
    if (verbose) pb$tick()
  }
  
  # Remove NULL elements and add names
  not_null <- !sapply(SCHT, is.null)
  SCHT <- SCHT[not_null]
  names(SCHT) <- names_vec[not_null]
  
  # Filter out columns with all zeros
  for (j in seq_along(SCHT)) {
    gene_matrix <- SCHT[[j]]
    
    if (is(gene_matrix, "sparseMatrix")) {
      non_zero_cols <- Matrix::colSums(gene_matrix) > 0
    } else {
      non_zero_cols <- base::colSums(gene_matrix) > 0
    }
    
    if (any(non_zero_cols)) {
      SCHT[[j]] <- gene_matrix[, non_zero_cols, drop = FALSE]
    }
  }
  
  # Convert all sparse matrices to standard matrices
  if (verbose) {
    message("Converting to standard matrices...")
  }
  
  for (j in seq_along(SCHT)) {
    if (methods::is(SCHT[[j]], "sparseMatrix")) {
      SCHT[[j]] <- data.matrix(SCHT[[j]])
    }
  }
  
  # Final status
  if (verbose) {
    message(sprintf("Created isoform matrices for %d genes", length(SCHT)))
    message(sprintf("HVGs with single isoform: %d (removed)", length(hvg_single_iso)))
    message(sprintf("HVGs with multiple isoforms: %d (kept)", length(hvg_multi_iso)))
  }
  
  return(list(
    SCHT = SCHT,
    hvg_single_iso = hvg_single_iso,
    hvg_multi_iso = hvg_multi_iso
  ))
}

#' Build SCHT structure
#'
#' @description
#' Creates the final SCHT object structure with efficient data handling
#' and memory management.
#'
#' @param isoform_list List of gene-wise isoform expressions
#' @param transcript_info Transcript annotation data
#' @param n_cells Number of cells in dataset
#' @param cell_info Optional cell metadata
#' @param min_expr Minimum expression threshold
#' @param verbose Whether to show progress
#'
#' @return An SCHT object with attributes
#'
#' @keywords internal
.build_scht <- function(isoform_list,
                        transcript_info,
                        n_cells,
                        cell_info = NULL,
                        min_expr,
                        verbose = TRUE) {
  
  # Create basic structure
  scht <- isoform_list
  
  # Add metadata
  attr(scht, "creation_date") <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  attr(scht, "package_version") <- "1.3.0"
  attr(scht, "cell_info") <- cell_info
  
  # Compute statistics efficiently
  if (length(scht) > 0) {
    # Number of genes is straightforward
    n_genes <- length(scht)
    
    # Calculate total transcripts efficiently
    total_transcripts <- sum(vapply(scht, nrow, FUN.VALUE = integer(1)))
    
    # Calculate mean isoforms per gene
    mean_isoforms <- total_transcripts / n_genes
    
    # Optimised sparsity calculation
    sparsity_sum <- 0
    total_elements <- 0
    
    # Process in chunks for better memory management
    chunk_size <- min(100, length(scht))
    n_chunks <- ceiling(length(scht) / chunk_size)
    
    for (i in 1:n_chunks) {
      start_idx <- (i-1) * chunk_size + 1
      end_idx <- min(i * chunk_size, length(scht))
      chunk_range <- start_idx:end_idx
      
      for (j in chunk_range) {
        gene_matrix <- scht[[j]]
        # Use as.numeric to prevent integer overflow
        matrix_size <- as.numeric(nrow(gene_matrix)) * as.numeric(ncol(gene_matrix))
        
        if (is(gene_matrix, "sparseMatrix")) {
          # For sparse matrices, use nnzero for accurate non-zero count
          zeros_count <- matrix_size - Matrix::nnzero(gene_matrix)
        } else {
          # For dense matrices
          zeros_count <- sum(gene_matrix == 0)
        }
        
        sparsity_sum <- sparsity_sum + zeros_count
        total_elements <- total_elements + matrix_size
      }
    }
    
    # Calculate overall sparsity
    sparsity <- sparsity_sum / total_elements
  } else {
    n_genes <- 0
    total_transcripts <- 0
    mean_isoforms <- 0
    sparsity <- NA
  }
  
  # Set statistics
  attr(scht, "stats") <- list(
    n_genes = n_genes,
    n_cells = n_cells,
    total_transcripts = total_transcripts,
    mean_isoforms = mean_isoforms,
    sparsity = sparsity
  )
  
  # Set class
  class(scht) <- "SCHT"
  return(scht)
}

#' Generate cell type-specific SCHT (Optimized Version)
#'
#' @description
#' Creates cell type-specific data structures with optimised operations
#' for efficient handling of large datasets. This version inverts the loop
#' order to be more memory and stack efficient.
#'
#' @param scht Original SCHT object
#' @param cell_info Cell metadata with cell type information
#' @param qc_params QC parameters list
#'
#' @return A CellTypeSCHT object with cell type-specific data
#'
#' @keywords internal
.generate_cell_type_scht <- function(scht, cell_info, qc_params) {
  # Validate cell type information
  if (!"cell_type" %in% colnames(cell_info)) {
    stop("cell_info must contain 'cell_type' column")
  }
  
  # Get unique cell types and create a lookup for which cells belong to which type
  cell_types <- unique(as.character(cell_info$cell_type))
  
  # Ensure the first column of cell_info contains the cell IDs that match matrix colnames
  cell_ids <- as.character(cell_info[, 1])
  cells_by_cell_type <- split(cell_ids, cell_info$cell_type)
  
  # Pre-allocate a list structure for the results.
  cell_type_scht <- stats::setNames(vector("list", length(cell_types)), cell_types)
  for (ct in cell_types) {
    # Initialize each cell type with an empty list to hold gene matrices
    cell_type_scht[[ct]] <- list()
  }
  
  # --- OPTIMIZED LOOP: Iterate through GENES first, then assign to cell types ---
  for (gene_name in names(scht)) {
    original_mat <- scht[[gene_name]]
    
    # Get the column names (cell IDs) for the current gene's matrix once
    matrix_cells <- colnames(original_mat)
    
    # Now, for this single gene, iterate through the cell types
    for (current_cell_type in cell_types) {
      # Find which cells in this matrix belong to the current cell type
      cell_type_cells_for_this_gene <- intersect(matrix_cells, cells_by_cell_type[[current_cell_type]])
      
      if (length(cell_type_cells_for_this_gene) > 0) {
        # Extract the subset of the matrix for this cell type
        cell_type_mat <- original_mat[, cell_type_cells_for_this_gene, drop = FALSE]
        
        # Filter by expression (optional, but good practice)
        col_sums <- base::colSums(cell_type_mat)
        non_zero_cols <- which(col_sums > qc_params$min_expr)
        
        if (length(non_zero_cols) > 0) {
          # If there are cells left after filtering, add the matrix to our results
          cell_type_scht[[current_cell_type]][[gene_name]] <- cell_type_mat[, non_zero_cols, drop = FALSE]
        }
      }
    }
  }
  
  # --- Clean up the final structure ---
  
  # Remove cell types that ended up with zero genes
  non_empty_cell_types_mask <- sapply(cell_type_scht, function(x) length(x) > 0)
  cell_type_scht <- cell_type_scht[non_empty_cell_types_mask]
  
  # Create summary statistics
  summary_stats <- data.frame(
    cell_type = names(cell_type_scht),
    n_genes = sapply(cell_type_scht, length),
    n_cells = sapply(names(cell_type_scht), function(ct) length(cells_by_cell_type[[ct]])),
    row.names = NULL
  )
  
  # Create a simple data log
  data_log <- list()
  for(ct in names(cell_type_scht)){
    data_log[[ct]] <- list(
      initial_cells = summary_stats$n_cells[summary_stats$cell_type == ct],
      n_genes = summary_stats$n_genes[summary_stats$cell_type == ct]
    )
  }
  
  # Return the final result object
  result <- list(
    cell_type_matrices = cell_type_scht,
    summary = summary_stats,
    data_log = data_log
  )
  
  class(result) <- "CellTypeSCHT"
  return(result)
}

#' Create integrated SCHT
#'
#' @description
#' Creates an integrated data structure combining original SCHT and
#' cell type-specific analyses with efficient memory usage.
#'
#' @param scht_obj Original SCHT object
#' @param cell_type_result cell type-specific analysis results
#'
#' @return An IntegratedSCHT object
#'
#' @keywords internal
.create_integrated_scht <- function(scht_obj, cell_type_result) {
  # Create efficient structure preserving SCHT object
  final_result <- list(
    original_results = scht_obj,
    cell_type_matrices = cell_type_result$cell_type_matrices
  )
  
  # Set class and attributes efficiently
  class(final_result) <- "IntegratedSCHT"
  
  # Get current time once
  creation_date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # Set all attributes at once
  attributes(final_result) <- c(
    attributes(final_result),
    list(
      creation_date = creation_date,
      stats = list(
        original = attr(scht_obj, "stats"),
        cell_type_matrices = cell_type_result$summary
      ),
      data_log = cell_type_result$data_log
    )
  )
  
  return(final_result)
}

# Helper functions for enhanced QC tracking
.store_input_stats <- function(gene_counts, transcript_counts) {
  list(
    n_genes_original = nrow(gene_counts),
    n_transcripts_original = nrow(transcript_counts),
    n_cells_original = ncol(gene_counts),
    
    # Sparsity of original matrices
    gene_matrix_sparsity = if (inherits(gene_counts, "sparseMatrix")) {
      (1 - Matrix::nnzero(gene_counts) / length(gene_counts)) * 100
    } else {
      (sum(gene_counts == 0) / length(gene_counts)) * 100
    },
    
    transcript_matrix_sparsity = if (inherits(transcript_counts, "sparseMatrix")) {
      (1 - Matrix::nnzero(transcript_counts) / length(transcript_counts)) * 100
    } else {
      (sum(transcript_counts == 0) / length(transcript_counts)) * 100
    },
    
    # Basic statistics
    genes_per_cell = if(is(gene_counts, "sparseMatrix")) Matrix::colSums(gene_counts > 0) else base::colSums(gene_counts > 0),
    median_genes_per_cell = median(if(is(gene_counts, "sparseMatrix")) Matrix::colSums(gene_counts > 0) else base::colSums(gene_counts > 0)),
    median_transcripts_per_cell = median(if(is(transcript_counts, "sparseMatrix")) Matrix::colSums(transcript_counts > 0) else base::colSums(transcript_counts > 0)),
    median_counts_per_cell = median(if(is(gene_counts, "sparseMatrix")) Matrix::colSums(gene_counts) else base::colSums(gene_counts))
  )
}

.store_qc_parameters <- function(auto_params, used_params) {
  list(
    recommended = auto_params,
    applied = used_params
  )
}

.store_filtering_stats <- function(pre_qc_stats, existing_qc, input_stats, qc_params) {
  list(
    initial_qc = list(
      genes = list(
        original = pre_qc_stats$n_genes,
        removed = existing_qc$n_filtered_genes,
        remaining = pre_qc_stats$n_genes - existing_qc$n_filtered_genes,
        removal_rate = existing_qc$n_filtered_genes / pre_qc_stats$n_genes * 100
      ),
      transcripts = list(
        original = pre_qc_stats$n_transcripts,
        removed = existing_qc$n_filtered_transcripts,
        remaining = pre_qc_stats$n_transcripts - existing_qc$n_filtered_transcripts,
        removal_rate = existing_qc$n_filtered_transcripts / pre_qc_stats$n_transcripts * 100
      )
    ),
    cell_qc = list(
      cells = list(
        original = pre_qc_stats$n_cells,
        removed = existing_qc$n_filtered_cells,
        remaining = pre_qc_stats$n_cells - existing_qc$n_filtered_cells,
        removal_rate = existing_qc$n_filtered_cells / pre_qc_stats$n_cells * 100
      ),
      removal_reasons = list(
        too_few_genes = sum(input_stats$genes_per_cell < qc_params$min_genes_per_cell),
        too_many_genes = sum(input_stats$genes_per_cell > qc_params$max_genes_per_cell),
        percent_too_few = sum(input_stats$genes_per_cell < qc_params$min_genes_per_cell) / 
                         pre_qc_stats$n_cells * 100,
        percent_too_many = sum(input_stats$genes_per_cell > qc_params$max_genes_per_cell) / 
                          pre_qc_stats$n_cells * 100
      )
    )
  )
}

.store_hvg_stats <- function(pre_qc_stats, existing_qc, hvg_genes, scht_genes, 
                            single_iso_genes = NULL, multi_iso_genes = NULL) {
  list(
    total_genes_available = pre_qc_stats$n_genes - existing_qc$n_filtered_genes,
    hvg_requested = length(hvg_genes),
    hvg_selected = length(hvg_genes),
    hvg_in_scht = length(scht_genes),
    hvg_removed_single_isoform = length(hvg_genes) - length(scht_genes),
    selection_rate = length(hvg_genes) / (pre_qc_stats$n_genes - existing_qc$n_filtered_genes) * 100,
    # New detailed tracking
    hvg_with_single_isoform = if (!is.null(single_iso_genes)) length(single_iso_genes) else NA,
    hvg_with_multiple_isoforms = if (!is.null(multi_iso_genes)) length(multi_iso_genes) else NA,
    percent_multi_isoform = if (!is.null(multi_iso_genes) && length(hvg_genes) > 0) {
      length(multi_iso_genes) / length(hvg_genes) * 100
    } else NA
  )
}

.store_scht_structure_stats <- function(scht_list, n_cells_after_qc) {
  # Count isoforms per gene
  isoform_counts <- numeric()
  total_isoforms <- 0
  gene_names <- names(scht_list)
  
  for (gene in gene_names) {
    gene_mat <- scht_list[[gene]]
    if (!is.null(gene_mat) && length(gene_mat) > 0) {
      n_iso <- nrow(gene_mat)
      isoform_counts <- c(isoform_counts, n_iso)
      total_isoforms <- total_isoforms + n_iso
    }
  }
  
  # Calculate structure statistics
  list(
    n_genes_in_scht = length(isoform_counts),
    n_cells_after_qc = n_cells_after_qc,
    
    isoform_stats = list(
      total_isoforms = total_isoforms,
      max_isoforms_per_gene = max(isoform_counts),
      mean_isoforms_per_gene = mean(isoform_counts),
      median_isoforms_per_gene = median(isoform_counts),
      genes_with_single_isoform = sum(isoform_counts == 1),
      genes_with_multiple_isoforms = sum(isoform_counts > 1),
      percent_multi_isoform = sum(isoform_counts > 1) / length(isoform_counts) * 100
    ),
    
    isoform_distribution = table(isoform_counts)
  )
}

###################
# SCHT S3 Methods #
###################

#' Print method for SCHT objects
#'
#' @description
#' Provides a concise summary of an SCHT object, displaying key statistics
#' about genes, cells, and transcripts.
#'
#' @param x SCHT object to print
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
#' @examples
#' # Create example SCHT object
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
#'   require_cell_type = FALSE,  # Create basic SCHT
#'   verbose = FALSE
#' )
#' 
#' # Print the SCHT object
#' print(scht_obj)
#' @export
print.SCHT <- function(x, ...) {
  stats <- attr(x, "stats")
  
  # Efficient string concatenation
  cat(sprintf(
    "Single-Cell Hierarchical Tensor (SCHT) object\n%s\n%s\n%s\n%s\n%s\n",
    paste("Number of genes:", stats$n_genes),
    paste("Number of cells:", stats$n_cells),
    paste("Total number of transcripts:", stats$total_transcripts),
    paste("Mean isoforms per gene:", round(stats$mean_isoforms, 2)),
    paste("Overall sparsity:", round(stats$sparsity * 100, 1), "%")
  ))
  
  # Show matrix types if available
  # Note: SCHT objects are lists of gene matrices
  if (length(x) > 0 && !is.null(names(x))) {
    # Get first gene matrix
    first_gene <- names(x)[1]
    if (first_gene %in% names(x)) {
      sample_matrix <- x[[first_gene]]
      cat(sprintf("Matrix storage type: %s\n", class(sample_matrix)[1]))
    }
  }
}

#' Summary method for SCHT objects
#'
#' @description
#' Generates a detailed summary of an SCHT object, including preprocessing
#' information and data characteristics.
#'
#' @param object SCHT object to summarise
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
#' @examples
#' # Using the same SCHT object from print example
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
#'   require_cell_type = FALSE,
#'   verbose = FALSE
#' )
#' 
#' # Get summary of SCHT object
#' summary(scht_obj)
#' @export
summary.SCHT <- function(object, ...) {
  stats <- attr(object, "stats")
  preproc <- attr(object, "preprocessing")
  
  # Header
  cat("SCHT Object Summary:\n")
  cat("--------------------\n")
  
  # Basic statistics
  cat(sprintf(
    "  Cells: %d\n  Genes: %d\n  Total transcripts: %d\n  Mean isoforms: %.2f\n\n",
    stats$n_cells, stats$n_genes, stats$total_transcripts, stats$mean_isoforms
  ))
  
  # Preprocessing info if available
  if (!is.null(preproc)) {
    cat("Preprocessing Info:\n")
    cat(sprintf(
      "  HVGs selected: %d\n  QC-filtered genes: %d\n  QC-filtered transcripts: %d\n  QC-filtered cells: %d\n\n",
      length(preproc$hvg),
      preproc$qc_stats$n_filtered_genes,
      preproc$qc_stats$n_filtered_transcripts,
      preproc$qc_stats$n_filtered_cells
    ))
  }
  
  # Data characteristics
  cat("Data characteristics:\n")
  cat(sprintf("  Sparsity: %.1f%%\n", stats$sparsity * 100))
  
  cat(sprintf("  Created: %s\n", attr(object, "creation_date")))
  
  # Performance metrics if available
  perf <- attr(object, "performance")
  if (!is.null(perf)) {
    cat("\nPerformance metrics:\n")
    cat(sprintf(
      "  Processing time: %.2f seconds (%.2f minutes)\n  Memory utilised: %.2f MB\n",
      perf$total_time_sec, perf$total_time_sec/60, 
      ifelse(!is.null(perf$memory_used_mb), perf$memory_used_mb, perf$peak_memory_mb)
    ))
  }
}


############################
# CellTypeSCHT S3 Methods  #
############################

#' Print method for CellTypeSCHT objects
#'
#' @description
#' Displays a concise summary of cell type-specific SCHT analysis.
#'
#' @param x CellTypeSCHT object to print
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
#' @export
print.CellTypeSCHT <- function(x, ...) {
  cat("Cell Type-Specific SCHT Analysis\n")
  cat("===========================\n")
  cat(sprintf("Number of cell types: %d\n", length(x$cell_type_matrices)))
  
  # Print summary statistics
  cat("\nCell types summary:\n")
  print(x$summary)
}

#' Summary method for CellTypeSCHT objects
#'
#' @description
#' Generates a detailed summary of cell type-specific SCHT analysis.
#'
#' @param object CellTypeSCHT object to summarise
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
#' @export
summary.CellTypeSCHT <- function(object, ...) {
  cat("Cell Type-Specific SCHT Analysis Summary\n")
  cat("===================================\n")
  
  # Print summary statistics
  cat("\nCell Type statistics:\n")
  print(object$summary)
  
  # Print additional information about each cell type
  cat("\nDetailed cell type information:\n")
  for (cell_type_name in names(object$cell_type_matrices)) {
    cell_type_data <- object$cell_type_matrices[[cell_type_name]]
    cell_type_stats <- object$data_log[[cell_type_name]]
    
    # Calculate per-cell type metrics
    n_total_transcripts <- sum(vapply(cell_type_data, nrow, FUN.VALUE = integer(1)))
    mean_isoforms <- ifelse(length(cell_type_data) > 0, n_total_transcripts / length(cell_type_data), 0)
    
    cat(sprintf("\n  Cell Type: %s\n", cell_type_name))
    cat(sprintf("    Genes: %d\n", length(cell_type_data)))
    cat(sprintf("    Cells: %d\n", cell_type_stats$initial_cells))
    cat(sprintf("    Total transcripts: %d\n", n_total_transcripts))
    cat(sprintf("    Mean isoforms per gene: %.2f\n", mean_isoforms))
  }
}

##############################
# Integrated SCHT S3 Methods #
##############################

#' Print method for IntegratedSCHT objects
#'
#' @description
#' Displays a summary of an integrated SCHT object, including both original
#' and cell type-specific information.
#'
#' @param x IntegratedSCHT object to print
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
#' @examples
#' # Create IntegratedSCHT object (with cell types)
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
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
#'   require_cell_type = TRUE,  # This creates IntegratedSCHT
#'   verbose = FALSE
#' )
#' 
#' # Print the IntegratedSCHT object
#' print(integrated_scht)
#' @export
print.IntegratedSCHT <- function(x, ...) {
  cat("Integrated SCHT Object\n")
  cat("=====================\n")
  
  # Original SCHT summary
  cat("\nOriginal SCHT:\n")
  print(x$original_results)
  
  # Cell type-specific summary
  cat("\nCell type-specific analysis:\n")
  cat(sprintf("Number of cell types: %d\n", length(x$cell_type_matrices)))
  
  # Show cell type names
  if (length(x$cell_type_matrices) > 0) {
    cat("Cell Types:", paste(names(x$cell_type_matrices), collapse = ", "), "\n")
  }
  
  # Show creation time
  cat(sprintf("\nCreated: %s\n", attr(x, "creation_date")))
}

#' Summary method for IntegratedSCHT objects
#'
#' @description
#' Provides a comprehensive summary of an integrated SCHT object, including
#' original SCHT summary and cell type-specific analysis results.
#'
#' @param object IntegratedSCHT object to summarise
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
#' @examples
#' # Using the same IntegratedSCHT object
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
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
#'   require_cell_type = TRUE,
#'   verbose = FALSE
#' )
#' 
#' # Get summary of IntegratedSCHT object
#' summary(integrated_scht)
#' @export
summary.IntegratedSCHT <- function(object, ...) {
  cat("Integrated SCHT Summary\n")
  cat("======================\n")
  
  # Original SCHT summary
  cat("\nOriginal SCHT Summary:\n")
  cat("----------------------\n")
  summary(object$original_results)
  
  # Cell type-specific summary
  cat("\nCell type-specific Summary:\n")
  cat("----------------------\n")
  stats <- attr(object, "stats")$cell_type_matrices
  print(stats)
  
  
  # Performance metrics if available
  perf <- attr(object, "performance")
  if (!is.null(perf)) {
    cat("\nPerformance metrics:\n")
    cat(sprintf("  Total processing time: %.2f seconds (%.2f minutes)\n", 
                perf$total_time_sec, perf$total_time_sec/60))
    cat(sprintf("  Memory utilised: %.2f MB\n", 
                ifelse(!is.null(perf$memory_used_mb), perf$memory_used_mb, perf$peak_memory_mb)))
  }
  
  # Explanatory note
  cat("\nNote: The actual number of cells for each gene may vary from the shown n_cells,\n")
  cat("as cells with no expression for specific genes are removed from their respective matrices.\n")
  cat("This cell-wise filtering is performed independently for each gene to maintain data quality\n")
  cat("and avoid spurious zero expressions in the cell type-specific analyses.\n")
}

#' SCHT creation for large datasets with enhanced features
#'
#' @description
#' Creates a Single-Cell Hierarchical Tensor from raw count data
#' with efficient operations for large-scale datasets. Supports
#' pre-normalized data and comprehensive QC tracking.
#'
#' @param gene_counts Gene-level counts as matrix, data frame, or sparse matrix.
#'   Row names can be either gene IDs (e.g., ENSG00000000001) or gene names (e.g., GAPDH).
#'   The function automatically detects which format is used and handles both appropriately.
#' @param transcript_counts Transcript-level counts as matrix, data frame, or sparse matrix
#' @param transcript_info Data frame with transcript annotations
#' @param cell_info Optional data frame with cell metadata
#' @param n_hvg Number of highly variable genes to select
#' @param qc_params QC parameters list
#' @param require_cell_type Whether cell type information is required
#' @param verbose Whether to show progress
#' @param sparsity_threshold Minimum sparsity to use sparse representation (0-1)
#' @param input_type Type of input data: "raw_counts" or "normalised" (e.g., TPM, FPKM)
#'
#' @return An SCHT or IntegratedSCHT object with the following attributes:
#'   \itemize{
#'     \item \code{performance}: List containing performance metrics:
#'       \itemize{
#'         \item \code{total_time_sec}: Total processing time in seconds
#'         \item \code{memory_used_mb}: Memory utilised in megabytes
#'       }
#'     \item \code{preprocessing}: List containing QC statistics and enhanced QC report
#'     \item \code{creation_date}: Timestamp of object creation
#'     \item \code{package_version}: Version of ScIsoX used
#'   }
#' @export
#'
#' @examples
#' # Load example data
#' data(gene_counts_blood)
#' data(transcript_counts_blood)
#' data(transcript_info)
#' data(sample2stage)
#' 
#' # Example 1: Basic SCHT creation with default parameters
#' scht_basic <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   n_hvg = 2000,
#'   verbose = TRUE
#' )
#' 
#' # Examine the structure
#' print(scht_basic)
#' summary(scht_basic)
#' 
#' # Example 2: Using recommended QC parameters
#' # First, get QC recommendations
#' qc_recommendations <- recommend_qc_parameters(gene_counts_blood)
#' print(qc_recommendations$explanation)
#' 
#' # Use the moderate (90% interval) strategy
#' recommended_params <- qc_recommendations$Interval_90
#' # Add the missing parameters
#' recommended_params$min_cells_expressing <- 0.02
#' recommended_params$min_expr <- 1e-6
#' 
#' scht_recommended <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   n_hvg = 1500,
#'   qc_params = recommended_params,
#'   verbose = TRUE
#' )
#' 
#' # Example 3: Manual QC parameters for different experimental designs
#' # For high-depth sequencing or full-length protocols
#' scht_manual <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   n_hvg = 1000,
#'   qc_params = list(
#'     min_genes_per_cell = 4000,       
#'     max_genes_per_cell = 10000,      
#'     min_cells_expressing = 0.02,   
#'     min_expr = 1e-6
#'   ),
#'   verbose = TRUE
#' )
#' 
#' # Example 4: Working with normalised data
#' # If your data is already normalised (TPM/FPKM)
#' \donttest{
#' scht_normalised <- create_scht(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   cell_info = sample2stage,
#'   n_hvg = 1500,
#'   input_type = "normalised",
#'   verbose = TRUE
#' )
#' }
#' 
#' # Example 5: Accessing SCHT components
#' # Get list of cell types
#' cell_types <- names(scht_basic)
#' print(cell_types)
#' 
#' # Access specific cell type data
#' if ("AEC" %in% cell_types) {
#'   aec_data <- scht_basic[["AEC"]]
#'   print(paste("AEC dimensions:", paste(dim(aec_data), collapse=" x ")))
#' }
#' 
#' # Get highly variable genes used
#' hvgs <- attr(scht_basic, "hvg_genes")
#' print(head(hvgs))
#' print(paste("Total HVGs:", length(hvgs)))
#' 
#' # Check transcript usage per gene
#' transcript_usage <- attr(scht_basic, "transcript_gene_map")
#' print(head(transcript_usage))
create_scht <- function(gene_counts,
                        transcript_counts,
                        transcript_info,
                        cell_info = NULL,
                        n_hvg = 3000,
                        qc_params = list(
                          min_genes_per_cell = 200,       
                          max_genes_per_cell = 20000,      
                          min_cells_expressing = 0.02,   
                          min_expr = 1e-4),
                        require_cell_type = TRUE,
                        verbose = TRUE,
                        sparsity_threshold = 0.4,
                        input_type = c("raw_counts", "normalised")) {
  
  # Match input type
  input_type <- match.arg(input_type)
  
  # Start timing
  total_start_time <- Sys.time()
  mem_usage_start <- gc(reset = TRUE)
  
  if (verbose) {
    message("Starting SCHT creation with matrix support...")
    message(sprintf("Input type: %s", input_type))
    message("Checking package dependencies...")
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      warning("Package 'Matrix' is required for optimal performance with large datasets")
    }
  }
  
  # Step 0: Process and validate input data
  if (verbose) message("Step 0: Processing input data and converting to optimal format...")
  step0_time <- Sys.time()
  
  input_result <- .process_input_data(
    gene_counts,
    transcript_counts,
    transcript_info,
    cell_info,
    n_hvg,
    qc_params,
    require_cell_type,
    sparsity_threshold
  )
  
  gene_counts <- input_result$gene_counts
  transcript_counts <- input_result$transcript_counts
  transcript_info <- input_result$transcript_info
  cell_info <- input_result$cell_info
  
  # Store input statistics for enhanced QC
  input_stats <- .store_input_stats(gene_counts, transcript_counts)
  
  # Add cell type distribution if available
  if (!is.null(cell_info) && "cell_type" %in% colnames(cell_info)) {
    cell_type_dist <- table(cell_info$cell_type)
    input_stats$cell_type_distribution <- list(
      counts = as.list(cell_type_dist),
      percentages = as.list(prop.table(cell_type_dist) * 100),
      n_cell_types = length(unique(cell_info$cell_type))
    )
  }
  
  if (verbose) {
    mem_usage <- gc(reset = TRUE)
    is_sparse_gene <- is(gene_counts, "sparseMatrix")
    is_sparse_transcript <- is(transcript_counts, "sparseMatrix")
    
    message(sprintf("Input processing complete (%.2f sec). Memory usage: %.2f MB", 
                    as.numeric(difftime(Sys.time(), step0_time, units = "secs")),
                    sum(mem_usage[, 2] - mem_usage_start[, 2])))
    
    message(sprintf("Gene counts matrix: %d rows x %d cols (%s)", 
                    nrow(gene_counts), 
                    ncol(gene_counts),
                    ifelse(is_sparse_gene, "sparse", "dense")))
    
    message(sprintf("Transcript counts matrix: %d rows x %d cols (%s)", 
                    nrow(transcript_counts), 
                    ncol(transcript_counts),
                    ifelse(is_sparse_transcript, "sparse", "dense")))
  }
  
  # Always calculate QC recommendations for reporting purposes
  # This will provide the three strategies regardless of user input
  auto_qc_params <- NULL
  
  # Try to use recommend_qc_parameters if available
  if (exists("recommend_qc_parameters", where = asNamespace("ScIsoX")) || 
      exists("recommend_qc_parameters")) {
    tryCatch({
      # Try package namespace first
      if (exists("recommend_qc_parameters", where = asNamespace("ScIsoX"))) {
        auto_qc_params <- ScIsoX::recommend_qc_parameters(gene_counts)
      } else {
        auto_qc_params <- recommend_qc_parameters(gene_counts)
      }
        
        if (verbose) {
          message("Auto-detected QC parameters using three strategies:")
          message(sprintf("  MAD Strategy: min=%d, max=%d genes per cell",
                         auto_qc_params$MAD_strategy$min_genes_per_cell,
                         auto_qc_params$MAD_strategy$max_genes_per_cell))
          message(sprintf("  Interval 90: min=%d, max=%d genes per cell",
                         auto_qc_params$Interval_90$min_genes_per_cell,
                         auto_qc_params$Interval_90$max_genes_per_cell))
          message(sprintf("  Interval 80: min=%d, max=%d genes per cell",
                         auto_qc_params$Interval_80$min_genes_per_cell,
                         auto_qc_params$Interval_80$max_genes_per_cell))
        }
        
        # Use Interval_90 as default if not specified
        if (is.null(qc_params$min_genes_per_cell)) {
          qc_params$min_genes_per_cell <- auto_qc_params$Interval_90$min_genes_per_cell
        }
        if (is.null(qc_params$max_genes_per_cell)) {
          qc_params$max_genes_per_cell <- auto_qc_params$Interval_90$max_genes_per_cell
        }
      }, error = function(e) {
        # Fallback to manual calculation of all three strategies
        if (verbose) {
          message("Using fallback method for QC recommendations")
        }
        auto_qc_params <- list(
          MAD_strategy = list(
            min_genes_per_cell = as.integer(max(0, median(input_stats$genes_per_cell) - 3 * mad(input_stats$genes_per_cell))),
            max_genes_per_cell = as.integer(median(input_stats$genes_per_cell) + 3 * mad(input_stats$genes_per_cell))
          ),
          Interval_90 = list(
            min_genes_per_cell = as.integer(quantile(input_stats$genes_per_cell, 0.05)),
            max_genes_per_cell = as.integer(quantile(input_stats$genes_per_cell, 0.95))
          ),
          Interval_80 = list(
            min_genes_per_cell = as.integer(quantile(input_stats$genes_per_cell, 0.10)),
            max_genes_per_cell = as.integer(quantile(input_stats$genes_per_cell, 0.90))
          )
        )
        
        if (verbose) {
          message(sprintf("Auto-detected QC parameters (simple method): min=%d, max=%d genes per cell",
                         auto_qc_params$Interval_90$min_genes_per_cell,
                         auto_qc_params$Interval_90$max_genes_per_cell))
        }
        
        # Use auto-detected values if not provided
        if (is.null(qc_params$min_genes_per_cell)) {
          qc_params$min_genes_per_cell <- auto_qc_params$Interval_90$min_genes_per_cell
        }
        if (is.null(qc_params$max_genes_per_cell)) {
          qc_params$max_genes_per_cell <- auto_qc_params$Interval_90$max_genes_per_cell
        }
      })
  } else {
    # If function doesn't exist, calculate all three strategies manually
    auto_qc_params <- list(
      MAD_strategy = list(
        min_genes_per_cell = as.integer(max(0, median(input_stats$genes_per_cell) - 3 * mad(input_stats$genes_per_cell))),
        max_genes_per_cell = as.integer(median(input_stats$genes_per_cell) + 3 * mad(input_stats$genes_per_cell))
      ),
      Interval_90 = list(
        min_genes_per_cell = as.integer(quantile(input_stats$genes_per_cell, 0.05)),
        max_genes_per_cell = as.integer(quantile(input_stats$genes_per_cell, 0.95))
      ),
      Interval_80 = list(
        min_genes_per_cell = as.integer(quantile(input_stats$genes_per_cell, 0.10)),
        max_genes_per_cell = as.integer(quantile(input_stats$genes_per_cell, 0.90))
      )
    )
      
      if (verbose) {
        message(sprintf("Auto-detected QC parameters: min=%d, max=%d genes per cell",
                       auto_qc_params$Interval_90$min_genes_per_cell,
                       auto_qc_params$Interval_90$max_genes_per_cell))
      }
      
    # Use auto-detected values if not provided
    if (is.null(qc_params$min_genes_per_cell)) {
      qc_params$min_genes_per_cell <- auto_qc_params$Interval_90$min_genes_per_cell
    }
    if (is.null(qc_params$max_genes_per_cell)) {
      qc_params$max_genes_per_cell <- auto_qc_params$Interval_90$max_genes_per_cell
    }
  }
  
  # Store QC parameter info
  qc_params_info <- .store_qc_parameters(auto_qc_params, qc_params)
  
  # Track pre-QC statistics
  pre_qc_stats <- list(
    n_genes = nrow(gene_counts),
    n_transcripts = nrow(transcript_counts),
    n_cells = ncol(gene_counts)
  )
  
  # Step 1: Initial QC
  if (verbose) message("Step 1: Performing initial quality control...")
  step1_time <- Sys.time()
  
  qc_result <- .perform_qc(
    gene_counts, 
    transcript_counts, 
    transcript_info,
    qc_params
  )
  
  if (verbose) {
    mem_usage <- gc(reset = TRUE)
    message(sprintf("Initial QC complete (%.2f sec). Memory usage: %.2f MB", 
                    as.numeric(difftime(Sys.time(), step1_time, units = "secs")),
                    sum(mem_usage[, 2] - mem_usage_start[, 2])))
    
    message(sprintf("Filtered %d genes and %d transcripts",
                    qc_result$n_filtered_genes,
                    qc_result$n_filtered_transcripts))
  }
  
  # Step 2: Cell QC
  if (verbose) message("Step 2: Filtering poor-quality cells...")
  step2_time <- Sys.time()
  
  cell_qc_result <- .perform_cell_qc(
    qc_result$gene_counts_filtered,
    qc_result$transcript_counts_filtered,
    qc_params
  )
  
  if (verbose) {
    mem_usage <- gc(reset = TRUE)
    message(sprintf("Cell QC complete (%.2f sec). Memory usage: %.2f MB", 
                    as.numeric(difftime(Sys.time(), step2_time, units = "secs")),
                    sum(mem_usage[, 2] - mem_usage_start[, 2])))
    
    message(sprintf("Filtered %d cells, %d remaining",
                    cell_qc_result$n_filtered_cells,
                    cell_qc_result$n_cells))
  }
  
  # Step 3: CPM normalisation (conditional based on input_type)
  if (input_type == "raw_counts") {
    if (verbose) message("Step 3: Performing CPM normalisation...")
    step3_time <- Sys.time()
    
    norm_result <- .perform_normalisation(
      cell_qc_result$gene_counts_filtered,
      cell_qc_result$transcript_counts_filtered
    )
    
    gene_counts_norm <- norm_result$gene_counts_norm
    transcript_counts_norm <- norm_result$transcript_counts_norm
    
    if (verbose) {
      mem_usage <- gc(reset = TRUE)
      message(sprintf("CPM normalisation complete (%.2f sec). Memory usage: %.2f MB", 
                      as.numeric(difftime(Sys.time(), step3_time, units = "secs")),
                      sum(mem_usage[, 2] - mem_usage_start[, 2])))
    }
  } else {
    if (verbose) message("Step 3: Skipping CPM normalisation (data already normalised)...")
    gene_counts_norm <- cell_qc_result$gene_counts_filtered
    transcript_counts_norm <- cell_qc_result$transcript_counts_filtered
  }
  
  # Step 4: Log transformation
  if (verbose) message("Step 4: Performing log2 transformation...")
  step4_time <- Sys.time()
  
  log_result <- .perform_log_transform(
    gene_counts_norm,
    transcript_counts_norm,
    verbose
  )
  
  gene_mat_final <- log_result$gene_mat
  transcript_mat_final <- log_result$transcript_mat
  
  if (verbose) {
    mem_usage <- gc(reset = TRUE)
    message(sprintf("Log2 transformation complete (%.2f sec). Memory usage: %.2f MB", 
                    as.numeric(difftime(Sys.time(), step4_time, units = "secs")),
                    sum(mem_usage[, 2] - mem_usage_start[, 2])))
  }
  
  # Step 5: Select HVGs
  if (verbose) message("Step 5: Selecting highly variable genes...")
  step5_time <- Sys.time()
  
  hvg_result <- .select_hvgs(
    gene_mat_final,
    n_hvg = n_hvg,
    qc_params,
    verbose = verbose
  )
  
  if (verbose) {
    mem_usage <- gc(reset = TRUE)
    message(sprintf("HVG selection complete (%.2f sec). Memory usage: %.2f MB", 
                    as.numeric(difftime(Sys.time(), step5_time, units = "secs")),
                    sum(mem_usage[, 2] - mem_usage_start[, 2])))
    
    message(sprintf("Selected %d highly variable genes",
                    length(hvg_result$var_genes)))
  }
  
  # Step 6: Create SCHT
  if (verbose) message("Step 6: Creating isoform matrices for HVGs...")
  step6_time <- Sys.time()
  
  isoform_result <- .create_SCHT(
    transcript_mat_final,
    hvg_result$var_genes,
    qc_result$transcript_info_filtered,
    verbose = verbose
  )
  
  if (verbose) {
    mem_usage <- gc(reset = TRUE)
    message(sprintf("Isoform matrix creation complete (%.2f sec). Memory usage: %.2f MB", 
                    as.numeric(difftime(Sys.time(), step6_time, units = "secs")),
                    sum(mem_usage[, 2] - mem_usage_start[, 2])))
    
    message(sprintf("Created matrices for %d genes with multiple isoforms",
                    length(isoform_result$SCHT)))
  }
  
  # Step 7: Build SCHT structure
  if (verbose) message("Step 7: Building the final SCHT structure...")
  step7_time <- Sys.time()
  
  scht_obj <- .build_scht(
    isoform_result$SCHT,
    qc_result$transcript_info_filtered,
    cell_qc_result$n_cells,
    cell_info,
    qc_params$min_expr,
    verbose = verbose
  )
  
  if (verbose) {
    mem_usage <- gc(reset = TRUE)
    message(sprintf("SCHT structure building complete (%.2f sec). Memory usage: %.2f MB", 
                    as.numeric(difftime(Sys.time(), step7_time, units = "secs")),
                    sum(mem_usage[, 2] - mem_usage_start[, 2])))
  }
  
  # Store comprehensive filtering statistics
  filtering_stats <- .store_filtering_stats(
    pre_qc_stats, 
    list(
      n_filtered_genes = qc_result$n_filtered_genes,
      n_filtered_transcripts = qc_result$n_filtered_transcripts,
      n_filtered_cells = cell_qc_result$n_filtered_cells
    ),
    input_stats,
    qc_params
  )
  
  # Store HVG statistics with single/multi isoform tracking
  hvg_stats <- .store_hvg_stats(
    pre_qc_stats,
    list(
      n_filtered_genes = qc_result$n_filtered_genes,
      n_filtered_transcripts = qc_result$n_filtered_transcripts,
      n_filtered_cells = cell_qc_result$n_filtered_cells
    ),
    hvg_result$var_genes,
    names(isoform_result$SCHT),
    single_iso_genes = isoform_result$hvg_single_iso,
    multi_iso_genes = isoform_result$hvg_multi_iso
  )
  
  # Store SCHT structure statistics
  scht_structure_stats <- .store_scht_structure_stats(
    isoform_result$SCHT,
    cell_qc_result$n_cells
  )
  
  # Calculate sparsity if function is available
  if (exists("calculate_scht_sparsity")) {
    scht_sparsity <- tryCatch({
      calculate_scht_sparsity(scht_obj)
    }, error = function(e) {
      # If error, return NA values
      list(sparsity = NA, n_total = NA, n_nonzero = NA)
    })
    scht_structure_stats$sparsity <- scht_sparsity
  }
  
  # Perform comprehensive sparsity analysis if function is available
  comprehensive_sparsity <- NULL
  if (exists("analyse_sparsity_for_table")) {
    comprehensive_sparsity <- tryCatch({
      analyse_sparsity_for_table(
        gene_counts = gene_counts,  # Original
        transcript_counts = transcript_counts,  # Original
        transcript_info = transcript_info,
        scht_obj = scht_obj,
        filtered_gene_counts = cell_qc_result$gene_counts_filtered,
        filtered_transcript_counts = cell_qc_result$transcript_counts_filtered,
        dataset_name = "Current Dataset"
      )
    }, error = function(e) {
      NULL
    })
  }
  
  # Attach preprocessing metadata with enhanced QC report
  attr(scht_obj, "preprocessing") <- list(
    # Original QC stats (backward compatibility)
    qc_stats = list(
      n_filtered_genes = qc_result$n_filtered_genes,
      n_filtered_transcripts = qc_result$n_filtered_transcripts,
      n_filtered_cells = cell_qc_result$n_filtered_cells
    ),
    hvg = hvg_result$var_genes,
    hvg_stats = list(
      n_hvg = n_hvg,
      mean_expr = hvg_result$mean_expr,
      dispersion = hvg_result$dispersion
    ),
    params = list(
      n_hvg = n_hvg,
      min_genes_per_cell = qc_params$min_genes_per_cell,
      max_genes_per_cell = qc_params$max_genes_per_cell,
      min_cells_expressing = qc_params$min_cells_expressing,
      min_expr = qc_params$min_expr,
      require_cell_type = require_cell_type,
      verbose = verbose,
      sparsity_threshold = sparsity_threshold,
      input_type = input_type
    ),
    # NEW: Comprehensive QC report
    qc_report = list(
      # Input data characteristics
      input_data = input_stats,
      
      # QC parameters
      qc_parameters = qc_params_info,
      
      # Step-by-step filtering
      filtering = filtering_stats,
      
      # HVG selection
      hvg_selection = hvg_stats,
      
      # SCHT structure
      scht_structure = scht_structure_stats,
      
      # Metadata
      analysis_date = Sys.Date(),
      scisox_version = packageVersion("ScIsoX"),
      input_type = input_type,
      
      # Comprehensive sparsity analysis
      sparsity_analysis = comprehensive_sparsity
    )
  )
  
  # Perform cell type-specific analysis if required
  if (require_cell_type && !is.null(cell_info)) {
    # Step 8: Generate cell type-specific structures
    if (verbose) message("Step 8: Generating cell type-specific SCHT structures...")
    step8_time <- Sys.time()
    
    cell_type_result <- .generate_cell_type_scht(
      scht_obj, 
      cell_info, 
      qc_params
    )
    
    if (verbose) {
      mem_usage <- gc(reset = TRUE)
      message(sprintf("Cell type-specific structure generation complete (%.2f sec). Memory usage: %.2f MB",
                      as.numeric(difftime(Sys.time(), step8_time, units = "secs")),
                      sum(mem_usage[, 2] - mem_usage_start[, 2])))
    }
    
    # Step 9: Create integrated structure
    if (verbose) message("Step 9: Creating integrated SCHT structure...")
    step9_time <- Sys.time()
    
    scht_obj <- .create_integrated_scht(
      scht_obj, 
      cell_type_result
    )
    
    if (verbose) {
      mem_usage <- gc(reset = TRUE)
      message(sprintf("Integration complete (%.2f sec). Memory usage: %.2f MB", 
                      as.numeric(difftime(Sys.time(), step9_time, units = "secs")),
                      sum(mem_usage[, 2] - mem_usage_start[, 2])))
    }
    
    # Add cell type-specific stats to QC report
    if (exists("calculate_ct_scht_sparsity")) {
      ct_sparsity <- tryCatch({
        calculate_ct_scht_sparsity(scht_obj, verbose = FALSE)
      }, error = function(e) {
        NULL
      })
      
      if (!is.null(ct_sparsity)) {
        attr(scht_obj, "preprocessing")$qc_report$cell_type_sparsity <- ct_sparsity
      }
    }
  }
  
  # Final summary
  total_time <- as.numeric(difftime(Sys.time(), total_start_time, units = "secs"))
  final_mem <- gc(reset = TRUE)
  memory_increment <- sum(final_mem[, 2] - mem_usage_start[, 2])
  
  if (verbose) {
    message(sprintf("SCHT creation completed successfully in %.2f seconds (%.2f minutes)", 
                    total_time, total_time/60))
    message(sprintf("Memory utilised: %.2f MB", memory_increment))
  }
  
  # Add performance metrics to object
  attr(scht_obj, "performance") <- list(
    total_time_sec = total_time,
    memory_used_mb = memory_increment
  )
  
  # Update QC report with total time
  attr(scht_obj, "preprocessing")$qc_report$total_time <- total_time
  
  return(scht_obj)
}

