##############################################################
#  SCHT Sparsity Calculation Functions                       #
#                                                            #
#  Functions to calculate sparsity statistics for SCHT       #
#  data structures                                           #
#                                                            #
#  Author: Siyuan Wu & Ulf Schmitz                           #
#  Date: Jul 29, 2025                                        #
#  Package: ScIsoX V1.1.0                                    #
##############################################################

#' Calculate sparsity statistics for SCHT structure
#'
#' @description
#' Calculates the total number of elements, non-zero elements, and sparsity
#' percentage for an SCHT object. This demonstrates the memory efficiency
#' of the SCHT structure compared to alternative representations.
#'
#' @param scht_obj An SCHT object created by create_scht
#' @param verbose Print detailed statistics
#'
#' @return A list containing:
#'   \itemize{
#'     \item n_total: Total number of elements in SCHT
#'     \item n_nonzero: Number of non-zero elements
#'     \item n_zeros: Number of zero elements
#'     \item sparsity: Sparsity percentage
#'     \item n_genes: Number of genes in SCHT
#'     \item gene_stats: Per-gene statistics (if verbose)
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
#' # Calculate sparsity statistics
#' sparsity_stats <- calculate_scht_sparsity(scht_obj)
#' print(paste("SCHT sparsity:", round(sparsity_stats$sparsity, 2), "%"))
#' 
#' # With verbose output
#' sparsity_verbose <- calculate_scht_sparsity(scht_obj, verbose = TRUE)
calculate_scht_sparsity <- function(scht_obj, verbose = FALSE) {
  
  # Validate input and extract SCHT data based on object type
  if (inherits(scht_obj, "IntegratedSCHT")) {
    # For IntegratedSCHT, use the original results which is the SCHT
    scht_obj <- scht_obj$original_results
  } else if (!inherits(scht_obj, "SCHT")) {
    stop("Input must be an SCHT or IntegratedSCHT object")
  }
  
  # Extract SCHT data
  if (!is.null(scht_obj$scht)) {
    scht_data <- scht_obj$scht
  } else {
    # SCHT object is itself the gene list
    scht_data <- scht_obj
  }
  n_genes <- length(scht_data)
  
  # Initialize counters
  total_elements <- 0
  total_nonzero <- 0
  gene_stats <- list()
  
  # Process each gene
  for (i in seq_along(scht_data)) {
    gene_name <- names(scht_data)[i]
    gene_matrix <- scht_data[[i]]
    
    if (!is.null(gene_matrix) && length(gene_matrix) > 0) {
      # Calculate elements for this gene
      n_elem <- length(gene_matrix)
      total_elements <- total_elements + n_elem
      
      # Count non-zeros
      if (inherits(gene_matrix, "sparseMatrix")) {
        n_nz <- Matrix::nnzero(gene_matrix)
      } else {
        n_nz <- sum(gene_matrix != 0)
      }
      total_nonzero <- total_nonzero + n_nz
      
      # Store per-gene stats if verbose
      if (verbose) {
        gene_stats[[gene_name]] <- list(
          n_isoforms = nrow(gene_matrix),
          n_cells = ncol(gene_matrix),
          n_total = n_elem,
          n_nonzero = n_nz,
          n_zeros = n_elem - n_nz,
          sparsity = ((n_elem - n_nz) / n_elem) * 100
        )
      }
    }
  }
  
  # Calculate overall statistics
  n_zeros <- total_elements - total_nonzero
  sparsity <- (n_zeros / total_elements) * 100
  
  # Print verbose output
  if (verbose) {
    cat("SCHT Sparsity Statistics:\n")
    cat("========================\n")
    cat(sprintf("Number of genes: %d\n", n_genes))
    cat(sprintf("Total elements: %s\n", format(total_elements, big.mark = ",")))
    cat(sprintf("Non-zero elements: %s\n", format(total_nonzero, big.mark = ",")))
    cat(sprintf("Zero elements: %s\n", format(n_zeros, big.mark = ",")))
    cat(sprintf("Sparsity: %.2f%%\n", sparsity))
    
    # Per-gene summary
    if (length(gene_stats) > 0) {
      sparsity_values <- sapply(gene_stats, function(x) x$sparsity)
      cat("\nPer-gene sparsity:\n")
      cat(sprintf("  Min: %.2f%%\n", min(sparsity_values)))
      cat(sprintf("  Mean: %.2f%%\n", mean(sparsity_values)))
      cat(sprintf("  Median: %.2f%%\n", median(sparsity_values)))
      cat(sprintf("  Max: %.2f%%\n", max(sparsity_values)))
    }
  }
  
  # Return results
  result <- list(
    n_total = total_elements,
    n_nonzero = total_nonzero,
    n_zeros = n_zeros,
    sparsity = sparsity,
    n_genes = n_genes
  )
  
  if (verbose) {
    result$gene_stats <- gene_stats
  }
  
  return(result)
}

#' Calculate sparsity statistics for cell type-specific SCHT
#'
#' @description
#' Calculates sparsity statistics for an IntegratedSCHT object,
#' including both overall and per-cell-type statistics.
#'
#' @param integrated_scht An IntegratedSCHT object
#' @param verbose Print detailed statistics
#'
#' @return A list containing overall and per-cell-type sparsity statistics
#'
#' @export
calculate_ct_scht_sparsity <- function(integrated_scht, verbose = FALSE) {
  
  # Validate input
  if (!inherits(integrated_scht, "IntegratedSCHT")) {
    stop("Input must be an IntegratedSCHT object")
  }
  
  # Extract cell type data
  ct_data <- integrated_scht$cell_type_specific
  n_cell_types <- length(ct_data)
  
  if (verbose) {
    cat("Cell Type-Specific SCHT Sparsity:\n")
    cat("=================================\n")
    cat(sprintf("Number of cell types: %d\n\n", n_cell_types))
  }
  
  # Calculate statistics for each cell type
  ct_stats <- list()
  total_elements_all <- 0
  total_nonzero_all <- 0
  
  for (ct_name in names(ct_data)) {
    if (verbose) {
      cat(sprintf("\nCell type: %s\n", ct_name))
    }
    
    # Calculate sparsity for this cell type
    ct_scht <- ct_data[[ct_name]]
    ct_sparsity <- calculate_scht_sparsity(ct_scht)
    
    ct_stats[[ct_name]] <- ct_sparsity
    total_elements_all <- total_elements_all + ct_sparsity$n_total
    total_nonzero_all <- total_nonzero_all + ct_sparsity$n_nonzero
    
    if (verbose) {
      cat(sprintf("  Elements: %s\n", format(ct_sparsity$n_total, big.mark = ",")))
      cat(sprintf("  Sparsity: %.2f%%\n", ct_sparsity$sparsity))
    }
  }
  
  # Calculate overall statistics
  n_zeros_all <- total_elements_all - total_nonzero_all
  sparsity_all <- (n_zeros_all / total_elements_all) * 100
  
  if (verbose) {
    cat("\nOverall Statistics:\n")
    cat(sprintf("Total elements: %s\n", format(total_elements_all, big.mark = ",")))
    cat(sprintf("Total non-zero: %s\n", format(total_nonzero_all, big.mark = ",")))
    cat(sprintf("Total zeros: %s\n", format(n_zeros_all, big.mark = ",")))
    cat(sprintf("Overall sparsity: %.2f%%\n", sparsity_all))
  }
  
  # Return results
  list(
    n_cell_types = n_cell_types,
    overall = list(
      n_total = total_elements_all,
      n_nonzero = total_nonzero_all,
      n_zeros = n_zeros_all,
      sparsity = sparsity_all
    ),
    cell_type_stats = ct_stats
  )
}

#' Compare sparsity across different data representations
#'
#' @description
#' Compares the sparsity and memory efficiency of SCHT structure
#' against original matrices and theoretical tensor representations.
#'
#' @param scht_obj SCHT object
#' @param original_transcript_counts Original transcript count matrix
#' @param filtered_transcript_counts Filtered transcript count matrix (optional)
#'
#' @return Comparison statistics
#'
#' @export
compare_sparsity <- function(scht_obj, 
                           original_transcript_counts,
                           filtered_transcript_counts = NULL) {
  
  cat("Sparsity Comparison Analysis\n")
  cat("============================\n\n")
  
  # 1. Original matrix sparsity
  if (inherits(original_transcript_counts, "sparseMatrix")) {
    orig_nonzero <- Matrix::nnzero(original_transcript_counts)
  } else {
    orig_nonzero <- sum(original_transcript_counts != 0)
  }
  
  orig_total <- length(original_transcript_counts)
  orig_zeros <- orig_total - orig_nonzero
  orig_sparsity <- (orig_zeros / orig_total) * 100
  
  cat("Original Transcript Matrix:\n")
  cat(sprintf("  Total elements: %s\n", format(orig_total, big.mark = ",")))
  cat(sprintf("  Sparsity: %.2f%%\n", orig_sparsity))
  
  # 2. SCHT sparsity
  scht_stats <- calculate_scht_sparsity(scht_obj)
  
  cat("\nSCHT Structure:\n")
  cat(sprintf("  Total elements: %s\n", format(scht_stats$n_total, big.mark = ",")))
  cat(sprintf("  Sparsity: %.2f%%\n", scht_stats$sparsity))
  
  # 3. Memory savings
  memory_reduction <- (1 - scht_stats$n_total / orig_total) * 100
  zeros_avoided <- orig_zeros - scht_stats$n_zeros
  
  cat("\nMemory Efficiency:\n")
  cat(sprintf("  Size reduction: %.1f%%\n", memory_reduction))
  cat(sprintf("  Zero elements avoided: %s\n", format(zeros_avoided, big.mark = ",")))
  
  # Return comparison
  list(
    original = list(
      total = orig_total,
      nonzero = orig_nonzero,
      zeros = orig_zeros,
      sparsity = orig_sparsity
    ),
    scht = scht_stats,
    efficiency = list(
      size_reduction = memory_reduction,
      zeros_avoided = zeros_avoided
    )
  )
}

#' Comprehensive sparsity analysis for manuscript table
#'
#' @description
#' Generates comprehensive sparsity statistics matching the manuscript table format,
#' comparing original matrices, filtered matrices, theoretical tensor, and SCHT structure.
#'
#' @param gene_counts Original gene count matrix
#' @param transcript_counts Original transcript count matrix
#' @param transcript_info Transcript information data frame
#' @param scht_obj SCHT object created by create_scht
#' @param filtered_gene_counts Gene counts after QC (optional)
#' @param filtered_transcript_counts Transcript counts after QC (optional)
#' @param dataset_name Name of the dataset for reporting
#'
#' @return List containing comprehensive sparsity statistics at all levels
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
#' # Analyse sparsity across representations
#' results <- analyse_sparsity_for_table(
#'   gene_counts = gene_counts_blood,
#'   transcript_counts = transcript_counts_blood,
#'   transcript_info = transcript_info,
#'   scht_obj = scht_obj,
#'   dataset_name = "Blood Cell Data"
#' )
#' 
#' # Access specific statistics
#' print(paste("SCHT uses only", 
#'             round(results$scht$total / results$original$total * 100, 2),
#'             "% of original matrix size"))
analyse_sparsity_for_table <- function(gene_counts, 
                                     transcript_counts,
                                     transcript_info,
                                     scht_obj,
                                     filtered_gene_counts = NULL,
                                     filtered_transcript_counts = NULL,
                                     dataset_name = "Dataset") {
  
  cat(sprintf("\n=== Sparsity Analysis for %s ===\n", dataset_name))
  
  results <- list(dataset = dataset_name)
  
  # 1. Original Transcript Matrix
  cat("\n1. Original Transcript Matrix:\n")
  
  # Get gene and isoform counts
  n_genes_original <- nrow(gene_counts)  # Number of genes from gene count matrix
  n_isoforms_original <- nrow(transcript_counts)
  n_cells_original <- ncol(transcript_counts)
  
  # Calculate sparsity
  if (inherits(transcript_counts, "sparseMatrix")) {
    nonzero_original <- Matrix::nnzero(transcript_counts)
  } else {
    nonzero_original <- sum(transcript_counts != 0)
  }
  
  total_elements_original <- as.numeric(n_isoforms_original) * as.numeric(n_cells_original)
  zero_elements_original <- total_elements_original - nonzero_original
  sparsity_original <- (zero_elements_original / total_elements_original) * 100
  
  results$original <- list(
    n_genes = n_genes_original,
    n_isoforms = n_isoforms_original,
    n_cells = n_cells_original,
    nonzero = nonzero_original,
    zeros = zero_elements_original,
    total = total_elements_original,
    sparsity = sparsity_original
  )
  
  cat(sprintf("  Number of genes: %s\n", format(n_genes_original, big.mark = ",")))
  cat(sprintf("  Number of isoforms: %s\n", format(n_isoforms_original, big.mark = ",")))
  cat(sprintf("  Number of cells: %s\n", format(n_cells_original, big.mark = ",")))
  cat(sprintf("  Non-zero elements: %s\n", format(nonzero_original, big.mark = ",")))
  cat(sprintf("  Zero elements: %s\n", format(zero_elements_original, big.mark = ",")))
  cat(sprintf("  Total elements: %s\n", format(total_elements_original, big.mark = ",")))
  cat(sprintf("  Sparsity: %.2f%%\n", sparsity_original))
  
  # 2. Filtered Transcript Matrix (Post-QC HVG)
  cat("\n2. Filtered Transcript Matrix (Post-QC HVG):\n")
  
  # Get preprocessing info
  preproc <- attr(scht_obj, "preprocessing")
  
  # Check SCHT object structure and get the data
  if (inherits(scht_obj, "IntegratedSCHT")) {
    # For IntegratedSCHT, original_results IS the SCHT object
    scht_list <- scht_obj$original_results
    n_cells_filtered <- attr(scht_obj$original_results, "stats")$n_cells
  } else if (inherits(scht_obj, "SCHT")) {
    # SCHT object is itself the gene list
    scht_list <- scht_obj
    # Get n_cells from attributes
    n_cells_filtered <- attr(scht_obj, "stats")$n_cells
  } else {
    stop("Unable to determine SCHT structure")
  }
  
  # Get actual genes in SCHT (it's a named list)
  genes_in_scht <- names(scht_list)
  n_genes_filtered <- length(genes_in_scht)
  
  # Get SCHT gene names/IDs (which now may be a mix)
  scht_gene_identifiers <- names(scht_list)
  
  # For each identifier in SCHT, find corresponding transcripts
  all_hvg_transcripts <- character()
  
  for (identifier in scht_gene_identifiers) {
    # Check if it's a gene ID (contains ':') or gene name
    if (grepl(":", identifier)) {
      # It's a gene ID
      matching_transcripts <- transcript_info$transcript_id[transcript_info$gene_id == identifier]
    } else {
      # It's a gene name
      matching_transcripts <- transcript_info$transcript_id[transcript_info$gene_name == identifier]
    }
    all_hvg_transcripts <- c(all_hvg_transcripts, matching_transcripts)
  }
  
  unique_scht_transcripts <- unique(all_hvg_transcripts)
  
  # Collect what's actually in SCHT
  all_scht_transcripts <- character()
  for (gene in genes_in_scht) {
    gene_matrix <- scht_list[[gene]]
    if (!is.null(gene_matrix) && length(gene_matrix) > 0) {
      all_scht_transcripts <- c(all_scht_transcripts, rownames(gene_matrix))
    }
  }
  actual_scht_unique <- unique(all_scht_transcripts)
  
  # Use actual SCHT transcripts for the count
  n_isoforms_filtered <- length(actual_scht_unique)
  
  # Calculate non-zeros from the original transcript matrix for these specific isoforms and cells
  if (!is.null(transcript_counts)) {
    # Get all unique isoform IDs from SCHT for filtering
    scht_rows <- which(rownames(transcript_counts) %in% actual_scht_unique)
    
    if (length(scht_rows) > 0) {
      # Get cells that are actually in SCHT (post-QC cells)
      scht_cells <- unique(unlist(lapply(genes_in_scht, function(gene) {
        gene_mat <- scht_list[[gene]]
        if (!is.null(gene_mat)) colnames(gene_mat) else NULL
      })))
      
      # Find these cells in the original transcript matrix
      cell_cols <- which(colnames(transcript_counts) %in% scht_cells)
      
      if (length(cell_cols) > 0) {
        # Extract subset with only SCHT isoforms and cells from ORIGINAL matrix
        scht_subset <- transcript_counts[scht_rows, cell_cols, drop = FALSE]
        if (inherits(scht_subset, "sparseMatrix")) {
          nonzero_filtered <- Matrix::nnzero(scht_subset)
        } else {
          nonzero_filtered <- sum(scht_subset != 0)
        }
      } else {
        # Fallback - sum from SCHT matrices
        nonzero_filtered <- 0
        for (gene in genes_in_scht) {
          gene_matrix <- scht_list[[gene]]
          if (!is.null(gene_matrix)) {
            if (inherits(gene_matrix, "sparseMatrix")) {
              nonzero_filtered <- nonzero_filtered + Matrix::nnzero(gene_matrix)
            } else {
              nonzero_filtered <- nonzero_filtered + sum(gene_matrix != 0)
            }
          }
        }
      }
    } else {
      # Fallback - sum from SCHT matrices
      nonzero_filtered <- 0
      for (gene in genes_in_scht) {
        gene_matrix <- scht_list[[gene]]
        if (!is.null(gene_matrix)) {
          if (inherits(gene_matrix, "sparseMatrix")) {
            nonzero_filtered <- nonzero_filtered + Matrix::nnzero(gene_matrix)
          } else {
            nonzero_filtered <- nonzero_filtered + sum(gene_matrix != 0)
          }
        }
      }
    }
  } else {
    # Fallback if transcript_counts not provided - sum from SCHT matrices
    nonzero_filtered <- 0
    for (gene in genes_in_scht) {
      gene_matrix <- scht_list[[gene]]
      if (!is.null(gene_matrix)) {
        if (inherits(gene_matrix, "sparseMatrix")) {
          nonzero_filtered <- nonzero_filtered + Matrix::nnzero(gene_matrix)
        } else {
          nonzero_filtered <- nonzero_filtered + sum(gene_matrix != 0)
        }
      }
    }
  }
  
  total_elements_filtered <- as.numeric(n_isoforms_filtered) * as.numeric(n_cells_filtered)
  zero_elements_filtered <- total_elements_filtered - nonzero_filtered
  sparsity_filtered <- (zero_elements_filtered / total_elements_filtered) * 100
  
  results$filtered <- list(
    n_genes = n_genes_filtered,
    n_isoforms = n_isoforms_filtered,
    n_cells = n_cells_filtered,
    nonzero = nonzero_filtered,
    zeros = zero_elements_filtered,
    total = total_elements_filtered,
    sparsity = sparsity_filtered
  )
  
  cat(sprintf("  Number of genes: %s\n", format(n_genes_filtered, big.mark = ",")))
  cat(sprintf("  Number of isoforms: %s\n", format(n_isoforms_filtered, big.mark = ",")))
  cat(sprintf("  Number of cells: %s\n", format(n_cells_filtered, big.mark = ",")))
  cat(sprintf("  Non-zero elements: %s\n", format(nonzero_filtered, big.mark = ",")))
  cat(sprintf("  Zero elements: %s\n", format(zero_elements_filtered, big.mark = ",")))
  cat(sprintf("  Total elements: %s\n", format(total_elements_filtered, big.mark = ",")))
  cat(sprintf("  Sparsity: %.2f%%\n", sparsity_filtered))
  
  # 3. SCHT Structure (Post-QC HVG) - Calculate this first
  cat("\n3. SCHT Structure (Post-QC HVG):\n")
  
  # Calculate actual SCHT statistics manually since we have the data
  total_scht_elements <- 0
  total_scht_nonzero <- 0
  
  for (gene in genes_in_scht) {
    gene_matrix <- scht_list[[gene]]
    if (!is.null(gene_matrix) && length(gene_matrix) > 0) {
      n_elem <- length(gene_matrix)
      total_scht_elements <- total_scht_elements + n_elem
      
      if (inherits(gene_matrix, "sparseMatrix")) {
        total_scht_nonzero <- total_scht_nonzero + Matrix::nnzero(gene_matrix)
      } else {
        total_scht_nonzero <- total_scht_nonzero + sum(gene_matrix != 0)
      }
    }
  }
  
  scht_stats <- list(
    n_total = total_scht_elements,
    n_nonzero = total_scht_nonzero,
    n_zeros = total_scht_elements - total_scht_nonzero,
    sparsity = ((total_scht_elements - total_scht_nonzero) / total_scht_elements) * 100
  )
  
  results$scht <- list(
    nonzero = scht_stats$n_nonzero,
    zeros = scht_stats$n_zeros,
    total = scht_stats$n_total,
    sparsity = scht_stats$sparsity
  )
  
  cat(sprintf("  Non-zero elements: %s\n", format(scht_stats$n_nonzero, big.mark = ",")))
  cat(sprintf("  Zero elements: %s\n", format(scht_stats$n_zeros, big.mark = ",")))
  cat(sprintf("  Total elements: %s\n", format(scht_stats$n_total, big.mark = ",")))
  cat(sprintf("  Sparsity: %.2f%%\n", scht_stats$sparsity))
  
  # 4. Naive 3D Tensor (Post-QC HVG)
  cat("\n4. Naive 3D Tensor (Post-QC HVG):\n")
  
  # For naive 3D tensor based on SCHT dimensions:
  # - Number of genes in SCHT
  # - Maximum isoforms per gene in SCHT 
  # - Number of cells after QC
  
  # Use SCHT genes (not all HVGs)
  n_tensor_genes <- n_genes_filtered  # Same as SCHT
  
  # Find max isoforms in SCHT
  max_isoforms_in_scht <- 0
  for (gene in genes_in_scht) {
    gene_matrix <- scht_list[[gene]]
    if (!is.null(gene_matrix) && length(gene_matrix) > 0) {
      max_isoforms_in_scht <- max(max_isoforms_in_scht, nrow(gene_matrix))
    }
  }
  
  # Calculate theoretical tensor dimensions
  tensor_dims <- c(n_tensor_genes, max_isoforms_in_scht, n_cells_filtered)
  total_elements_tensor <- prod(as.numeric(tensor_dims))
  
  # For naive 3D tensor, non-zeros should equal SCHT non-zeros
  # This is because the tensor would contain the same data as SCHT, 
  # just arranged in a 3D structure with padding
  nonzero_tensor <- total_scht_nonzero  # Same as SCHT
  
  zero_elements_tensor <- total_elements_tensor - nonzero_tensor
  sparsity_tensor <- (zero_elements_tensor / total_elements_tensor) * 100
  
  results$tensor <- list(
    dimensions = tensor_dims,
    nonzero = nonzero_tensor,
    zeros = zero_elements_tensor,
    total = total_elements_tensor,
    sparsity = sparsity_tensor
  )
  
  cat(sprintf("  Required dimensions: %s x %s x %s\n", 
              format(tensor_dims[1], big.mark = ","),
              tensor_dims[2],
              format(tensor_dims[3], big.mark = ",")))
  cat(sprintf("  Non-zero elements: %s\n", format(nonzero_tensor, big.mark = ",")))
  cat(sprintf("  Zero elements: %s\n", format(zero_elements_tensor, big.mark = ",")))
  cat(sprintf("  Total elements: %s\n", format(total_elements_tensor, big.mark = ",")))
  cat(sprintf("  Sparsity: %.2f%%\n", sparsity_tensor))
  
  # 5. Zero Elements Avoided by SCHT
  cat("\n5. Zero Elements Avoided by SCHT:\n")
  
  zeros_avoided_original <- zero_elements_original - scht_stats$n_zeros
  zeros_avoided_filtered <- zero_elements_filtered - scht_stats$n_zeros
  zeros_avoided_tensor <- zero_elements_tensor - scht_stats$n_zeros
  
  results$zeros_avoided <- list(
    vs_original = zeros_avoided_original,
    vs_filtered = zeros_avoided_filtered,
    vs_tensor = zeros_avoided_tensor
  )
  
  cat(sprintf("  vs Original Matrix: %s\n", format(zeros_avoided_original, big.mark = ",")))
  cat(sprintf("  vs Filtered Matrix: %s\n", format(zeros_avoided_filtered, big.mark = ",")))
  cat(sprintf("  vs Naive 3D Tensor: %s\n", format(zeros_avoided_tensor, big.mark = ",")))
  
  # Summary
  cat("\n6. Memory Efficiency Summary:\n")
  cat(sprintf("  SCHT uses only %.2f%% of original matrix size\n", 
              (scht_stats$n_total / total_elements_original) * 100))
  cat(sprintf("  SCHT uses only %.2f%% of filtered matrix size\n", 
              (scht_stats$n_total / total_elements_filtered) * 100))
  cat(sprintf("  SCHT uses only %.2f%% of naive tensor size\n", 
              (scht_stats$n_total / total_elements_tensor) * 100))
  
  return(results)
}