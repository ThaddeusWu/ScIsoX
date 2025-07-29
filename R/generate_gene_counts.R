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
#' \dontrun{
#' # Create example data
#' isoform_mat <- matrix(
#'   rnorm(20), nrow = 4,
#'   dimnames = list(
#'     c("Gene1-201", "Gene1-202", "Gene2-201", "Gene2-202"),
#'     c("cell1", "cell2", "cell3", "cell4", "cell5")
#'   )
#' )
#'
#' # Generate gene counts
#' result <- generate_gene_counts(isoform_mat)
#'
#' # Access results
#' gene_counts <- result$gene_counts
#' transcript_info <- result$transcript_info
#' }
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
