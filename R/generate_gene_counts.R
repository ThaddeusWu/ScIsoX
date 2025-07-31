#################################################################
#  Generate Gene Counts from Isoform Expression Data            #
#                                                               #
#  Author: [Siyuan Wu & Ulf Schmitz]                            #
#  Institution: [James Cook University]                         #
#  Date: Jul 29, 2025                                           #
#  Package: ScIsoX V1.1.0                                       #
#################################################################

#' Generate gene counts from annotated isoform counts
#'
#' @name generate_gene_counts
#'
#' @description
#' Generates gene-level counts and transcript information from a fully annotated
#' isoform count matrix. This function converts an isoform-level expression matrix
#' into a gene-level expression matrix by summing up counts of isoforms belonging
#' to the same gene. The function assumes isoform names are in the format 
#' "GENE-NUMBER" (e.g., "Fasl-201").
#'
#' @param isoform_counts Matrix or data frame of isoform-level counts with
#'   isoform names as rownames. Each row represents an isoform, and each column
#'   represents a cell or sample.
#' @param show_progress Logical indicating whether to show progress messages
#'   (default: TRUE).
#'
#' @return A list containing two elements:
#'   \itemize{
#'     \item gene_counts: Data frame of gene-level counts where rows are genes
#'           and columns are cells/samples
#'     \item transcript_info: Data frame with columns:
#'           \itemize{
#'             \item transcript_id: Original isoform identifier
#'             \item transcript_name: Original isoform name
#'             \item gene_id: Associated gene identifier
#'             \item gene_name: Associated gene name
#'           }
#'   }
#'
#' @examples
#' # IMPORTANT: This function requires transcript names in "GeneName-TranscriptID" format
#' # For example: "Sox2-201", "Nanog-001", "Oct4-202"
#' 
#' # Example 1: Create synthetic data with correct naming format
#' # Simulate an isoform count matrix
#' set.seed(123)
#' n_genes <- 50
#' n_isoforms <- 120  # Some genes have multiple isoforms
#' n_cells <- 100
#' 
#' # Generate transcript names in required format
#' gene_names <- paste0("Gene", 1:n_genes)
#' transcript_names <- character(n_isoforms)
#' idx <- 1
#' 
#' for (i in 1:n_genes) {
#'   n_iso <- sample(1:4, 1)  # 1-4 isoforms per gene
#'   for (j in 1:n_iso) {
#'     transcript_names[idx] <- paste0(gene_names[i], "-", sprintf("%03d", j))
#'     idx <- idx + 1
#'     if (idx > n_isoforms) break
#'   }
#'   if (idx > n_isoforms) break
#' }
#' 
#' # Remove empty entries and create count matrix
#' transcript_names <- transcript_names[transcript_names != ""]
#' isoform_counts <- matrix(
#'   rpois(length(transcript_names) * n_cells, lambda = 5),
#'   nrow = length(transcript_names),
#'   dimnames = list(transcript_names,
#'                  paste0("Cell", 1:n_cells))
#' )
#' 
#' # Generate gene counts
#' result <- generate_gene_counts(isoform_counts, show_progress = FALSE)
#' 
#' # Check results
#' print(paste("Input transcripts:", nrow(isoform_counts)))
#' print(paste("Output genes:", nrow(result$gene_counts)))
#' print(head(result$transcript_info))
#' 
#' # Example 2: Working with real data that needs reformatting
#' # If your transcript names are not in the required format:
#' \donttest{
#' # Current format: ENSMUST00000193812
#' # Need to convert to: GeneName-TranscriptID format
#' # This would require mapping transcript IDs to gene names first
#' 
#' # For demonstration, create a small subset with correct format
#' demo_transcripts <- matrix(
#'   rpois(500, lambda = 3),
#'   nrow = 10,
#'   dimnames = list(
#'     c("Actb-201", "Actb-202", "Gapdh-201", "Gapdh-202", "Gapdh-203",
#'       "Tubb5-201", "Ppia-201", "Ppia-202", "B2m-201", "Hprt-201"),
#'     paste0("Cell", 1:50)
#'   )
#' )
#' 
#' demo_result <- generate_gene_counts(demo_transcripts)
#' print(demo_result$gene_counts[, 1:5])
#' }
#' 
#' # Example 3: Understanding the pseudo-transcript_info
#' # Even without create_transcript_info(), this function creates
#' # a basic transcript_info for gene-transcript mapping
#' print("Pseudo-transcript_info structure:")
#' print(str(result$transcript_info))
#' 
#' # This ensures compatibility with downstream SCHT analysis
#' # The mapping can be used directly in create_scht()
#'
#' @details
#' The function performs the following steps:
#' 1. Validates input format and extracts gene names from isoform identifiers
#' 2. Converts data to efficient format for processing
#' 3. Aggregates isoform counts to gene level
#' 4. Creates necessary transcript information
#' 5. Returns results in format compatible with SCHT analysis
#'
#' @importFrom data.table .SD
#' @import data.table
#' 
#' @export
utils::globalVariables(c(".SD", "gene"))
generate_gene_counts <- function(isoform_counts, show_progress = TRUE) {
  # Input validation
  if (!is.matrix(isoform_counts) && !is.data.frame(isoform_counts)) {
    stop("'isoform_counts' must be a matrix or data frame")
  }
  
  if (is.null(rownames(isoform_counts))) {
    stop("'isoform_counts' must have rownames containing isoform identifiers")
  }
  
  if (!is.logical(show_progress)) {
    stop("'show_progress' must be TRUE or FALSE")
  }
  
  # Progress message
  if (show_progress) {
    message("Starting gene count generation...")
  }
  
  # Extract gene names from isoform names
  isoform_names <- rownames(isoform_counts)
  gene_names <- sapply(strsplit(isoform_names, "-"), `[`, 1)
  
  if (show_progress) {
    message("Converting data format...")
  }
  
  # Convert to data.frame and add necessary columns
  dt <- as.data.frame(isoform_counts)
  dt$isoform_name <- isoform_names
  dt$gene <- gene_names
  
  # Convert to data.table for faster processing
  dt <- data.table::as.data.table(dt)
  
  if (show_progress) {
    message("Calculating gene-level counts...")
    pb <- progress::progress_bar$new(
      format = "Processing [:bar] :percent eta: :eta",
      total = 1,
      clear = FALSE,
      width = 60
    )
  }
  
  # Calculate gene counts
  cols_to_sum <- setdiff(names(dt), c("isoform_name", "gene"))
  gene_counts <- dt[, lapply(.SD, sum), 
                    by = gene,
                    .SDcols = cols_to_sum]
  
  if (show_progress) {
    pb$tick()
    message("Creating output format...")
  }
  
  # Format gene counts as data frame with proper row names
  gene_counts_df <- as.data.frame(gene_counts)
  rownames(gene_counts_df) <- gene_counts_df$gene
  gene_counts_df$gene <- NULL
  
  # Create transcript information data frame
  transcript_info <- data.frame(
    transcript_id = isoform_names,
    transcript_name = isoform_names,
    gene_id = gene_names,
    gene_name = gene_names,
    stringsAsFactors = FALSE,
    row.names = isoform_names
  )
  
  if (show_progress) {
    message("Process completed successfully.")
  }
  
  # Return results
  return(list(
    gene_counts = gene_counts_df,
    transcript_info = transcript_info
  ))
}
