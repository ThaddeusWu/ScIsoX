##############################################################
#  Single-Cell Hierarchical Tensor (SCHT) Creation Pipeline  #                     
#                                                            #
#  Author: [Siyuan Wu & Ulf Schmitz]                         #
#  Institution: [James Cook University]                      #
#  Date: Apr 22, 2025                                        #
#  Package: ScIsoX V1.0.0                                    #
##############################################################

# This implementation supports efficient handling of large-scale
# single-cell long-read data with reduced memory footprint.

########################
# Required Libraries   #
########################
#' @importFrom graphics abline hist mtext par polygon text
#' @importFrom stats density mad median na.omit quantile sd var setNames
#' @importFrom progress progress_bar
#' @importFrom methods is
#' @importFrom Matrix Matrix Diagonal
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
                                require_cell_type = FALSE,
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
    n_genes_per_cell <- colSums(gene_counts > 0)
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
    gene_col_sums <- colSums(gene_counts_filtered)
    transcript_col_sums <- colSums(transcript_counts_filtered)
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
  
  # Pre-allocate result list
  SCHT <- vector("list", length(var_genes))
  names_vec <- character(length(var_genes))
  
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
  
  # Process each gene
  for (i in seq_along(var_genes)) {
    current_gene <- var_genes[i]
    
    # Get transcripts for current gene using lookup
    gene_transcripts <- gene_to_trans[[current_gene]]
    if (!is.null(gene_transcripts)) {
      # Find transcripts present in data
      existing_transcripts <- intersect(gene_transcripts, transcript_rows)
      
      # Get gene name efficiently
      gene_name_selected <- gene_id_to_name[current_gene]
      names_vec[i] <- gene_name_selected
      
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
        
        # Only keep genes with multiple isoforms
        if (nrow(non_zero_transcripts) > 1) {
          SCHT[[i]] <- non_zero_transcripts
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
      non_zero_cols <- colSums(gene_matrix) > 0
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
  }
  
  return(list(SCHT = SCHT))
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
        matrix_size <- nrow(gene_matrix) * ncol(gene_matrix)
        
        if (is(gene_matrix, "sparseMatrix")) {
          # For sparse matrices, we can directly get the number of non-zeros
          zeros_count <- matrix_size - length(gene_matrix@x)
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

#' Generate cell type-specific SCHT
#'
#' @description
#' Creates cell type-specific data structures with optimised operations
#' for efficient handling of large datasets.
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
  
  # Get unique cell type efficiently
  if (is.factor(cell_info$cell_type)) {
    cell_types <- levels(cell_info$cell_type)
  } else {
    cell_types <- unique(cell_info$cell_type)
  }
  
  # Pre-allocate result structures
  cell_type_scht <- vector("list", length(cell_types))
  names(cell_type_scht) <- cell_types
  data_log <- vector("list", length(cell_types))
  names(data_log) <- cell_types
  
  # Create cell type-to-cells lookup
  cells_by_cell_type <- split(cell_info[, 1], cell_info$cell_type)
  
  # Process each cell type efficiently
  for (i in seq_along(cell_types)) {
    current_cell_type <- cell_types[i]
    cell_type_cells <- cells_by_cell_type[[current_cell_type]]
    initial_cells <- length(cell_type_cells)
    
    # Pre-allocate for gene matrices
    cell_type_matrices <- vector("list", length(scht))
    names(cell_type_matrices) <- names(scht)
    genes_found <- 0
    
    # Process each gene
    for (j in seq_along(scht)) {
      gene_name <- names(scht)[j]
      original_mat <- scht[[gene_name]]
      
      # Find matching cells
      matching_cols <- which(colnames(original_mat) %in% cell_type_cells)
      
      if (length(matching_cols) > 0) {
        # Extract cell type-specific matrix
        cell_type_mat <- original_mat[, matching_cols, drop = FALSE]
        
        # Filter by expression
        col_sums <- colSums(cell_type_mat)
        non_zero_cols <- which(col_sums > qc_params$min_expr)
        
        if (length(non_zero_cols) > 0) {
          cell_type_mat <- cell_type_mat[, non_zero_cols, drop = FALSE]
          cell_type_matrices[[j]] <- cell_type_mat
          genes_found <- genes_found + 1
        }
      }
    }
    
    # Clean up empty entries
    cell_type_matrices <- cell_type_matrices[!sapply(cell_type_matrices, is.null)]
    
    # Store results
    cell_type_scht[[i]] <- cell_type_matrices
    data_log[[i]] <- list(
      initial_cells = initial_cells,
      n_genes = genes_found
    )
  }
  
  # Remove cell types with zero genes
  non_empty_cell_types <- !sapply(cell_type_scht, function(x) length(x) == 0)
  cell_type_scht <- cell_type_scht[non_empty_cell_types]
  data_log <- data_log[non_empty_cell_types]
  
  
  # Create summary statistics
  summary_stats <- data.frame(
    cell_type = names(cell_type_scht),
    n_genes = vapply(cell_type_scht, length, FUN.VALUE = integer(1)),
    n_cells = vapply(data_log, function(x) x$initial_cells, FUN.VALUE = integer(1)),
    row.names = NULL
  )
  
  # Return result
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
        cell_type_specific = cell_type_result$summary
      ),
      data_log = cell_type_result$data_log
    )
  )
  
  return(final_result)
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
  if (length(x$original_results) > 0) {
    sample_matrix <- x$original_results[[1]]
    cat(sprintf("Matrix storage type: %s\n", class(sample_matrix)[1]))
  }
}

#' Summary method for SCHT objects
#'
#' @description
#' Generates a detailed summary of an SCHT object, including preprocessing
#' information and data characteristics.
#'
#' @param object SCHT object to summarize
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
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
  cat(sprintf(
    "  Sparsity: %.1f%%\n  Created: %s\n",
    stats$sparsity * 100,
    attr(object, "creation_date")
  ))
  
  # Performance metrics if available
  perf <- attr(object, "performance")
  if (!is.null(perf)) {
    cat("\nPerformance metrics:\n")
    cat(sprintf(
      "  Processing time: %.2f seconds (%.2f minutes)\n  Peak memory usage: %.2f MB\n",
      perf$total_time_sec, perf$total_time_sec/60, perf$peak_memory_mb
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
#' @param object CellTypeSCHT object to summarize
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
#' @export
summary.CellTypeSCHTSCHT <- function(object, ...) {
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
#' @param object IntegratedSCHT object to summarize
#' @param ... Additional arguments (not used)
#'
#' @return None (prints to console)
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
  stats <- attr(object, "stats")$cell_type_specific
  print(stats)
  
  # Performance metrics if available
  perf <- attr(object, "performance")
  if (!is.null(perf)) {
    cat("\nPerformance metrics:\n")
    cat(sprintf("  Total processing time: %.2f seconds (%.2f minutes)\n", 
                perf$total_time_sec, perf$total_time_sec/60))
    cat(sprintf("  Peak memory usage: %.2f MB\n", perf$peak_memory_mb))
  }
  
  # Explanatory note
  cat("\nNote: The actual number of cells for each gene may vary from the shown n_cells,\n")
  cat("as cells with no expression for specific genes are removed from their respective matrices.\n")
  cat("This cell-wise filtering is performed independently for each gene to maintain data quality\n")
  cat("and avoid spurious zero expressions in the cell type-specific analyses.\n")
}

#' SCHT creation for large datasets
#'
#' @description
#' Creates a Single-Cell Hierarchical Tensor from raw count data
#' with efficient operations for large-scale datasets.
#'
#' @param gene_counts Gene-level counts as matrix, data frame, or sparse matrix
#' @param transcript_counts Transcript-level counts as matrix, data frame, or sparse matrix
#' @param transcript_info Data frame with transcript annotations
#' @param cell_info Optional data frame with cell metadata
#' @param n_hvg Number of highly variable genes to select
#' @param qc_params QC parameters list
#' @param require_cell_type Whether cell type information is required
#' @param verbose Whether to show progress
#' @param sparsity_threshold Minimum sparsity to use sparse representation (0-1)
#'
#' @return An SCHT or IntegratedSCHT object
#' @export
#'
#' @examples
#' # For large-scale single-cell data:
#' # scht_obj <- create_scht(
#' #   gene_counts = gene_counts,         # Large count matrix
#' #   transcript_counts = transcript_counts,
#' #   transcript_info = transcript_info,
#' #   cell_info = cell_info,
#' #   n_hvg = 3000,
#' #   qc_params = list(
#' #     min_genes_per_cell = 200,
#' #     max_genes_per_cell = 10000,
#' #     min_cells_expressing = 0.02,
#' #     min_expr = 1e-4
#' #   ),
#' #   require_cell_type = TRUE,
#' #   verbose = TRUE,
#' #   sparsity_threshold = 0.5  # Adjust based on data sparsity
#' # )
create_scht <- function(gene_counts,
                        transcript_counts,
                        transcript_info,
                        cell_info = NULL,
                        n_hvg = 3000,
                        qc_params = list(
                          min_genes_per_cell = 200,       
                          max_genes_per_cell = 10000,      
                          min_cells_expressing = 0.02,   
                          min_expr = 1e-4),
                        require_cell_type = FALSE,
                        verbose = TRUE,
                        sparsity_threshold = 0.4) {
  
  # Start timing
  total_start_time <- Sys.time()
  mem_usage_start <- gc(reset = TRUE)
  
  if (verbose) {
    message("Starting SCHT creation with matrix support...")
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
  
  # Step 3: CPM normalisation
  if (verbose) message("Step 3: Performing CPM normalisation...")
  step3_time <- Sys.time()
  
  norm_result <- .perform_normalisation(
    cell_qc_result$gene_counts_filtered,
    cell_qc_result$transcript_counts_filtered
  )
  
  if (verbose) {
    mem_usage <- gc(reset = TRUE)
    message(sprintf("CPM normalisation complete (%.2f sec). Memory usage: %.2f MB", 
                    as.numeric(difftime(Sys.time(), step3_time, units = "secs")),
                    sum(mem_usage[, 2] - mem_usage_start[, 2])))
  }
  
  # Step 4: Log transformation
  if (verbose) message("Step 4: Performing log2 transformation...")
  step4_time <- Sys.time()
  
  log_result <- .perform_log_transform(
    norm_result$gene_counts_norm,
    norm_result$transcript_counts_norm,
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
  
  # Attach preprocessing metadata
  attr(scht_obj, "preprocessing") <- list(
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
      min_cells_expressing = qc_params$min_cells_expressing,
      min_expr = qc_params$min_expr,
      require_cell_type = require_cell_type,
      verbose = verbose,
      sparsity_threshold = sparsity_threshold
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
  }
  
  # Final summary
  total_time <- as.numeric(difftime(Sys.time(), total_start_time, units = "secs"))
  final_mem <- gc(reset = TRUE)
  peak_mem_usage <- sum(final_mem[, 2] - mem_usage_start[, 2])
  
  if (verbose) {
    message(sprintf("SCHT creation completed successfully in %.2f seconds (%.2f minutes)", 
                    total_time, total_time/60))
    message(sprintf("Peak memory usage: %.2f MB", peak_mem_usage))
  }
  
  # Add performance metrics to object
  attr(scht_obj, "performance") <- list(
    total_time_sec = total_time,
    peak_memory_mb = peak_mem_usage
  )
  
  return(scht_obj)
}
