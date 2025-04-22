##############################################################
#  Author: [Siyuan Wu & Ulf Schmitz]                         #
#  Institution: [James Cook University]                      #
#  Date: Apr 22, 2025                                        #
#  Package: ScIsoX V1.0.0                                    #
##############################################################

#' Create transcript information data frame from GTF file
#' 
#' @name create_transcript_info
#'
#' @description
#' Imports a GTF file and creates a standardised transcript information data frame
#' for use in the SCHT pipeline. Can optionally remove version numbers from
#' transcript and gene IDs.
#'
#' @param gtf_path Character string specifying the path to the GTF file
#' @param remove_version Logical indicating whether to remove version numbers from
#'   transcript and gene IDs (default: TRUE)
#' @param progress Logical indicating whether to show progress messages (default: TRUE)
#'
#' @return A data frame containing transcript information with columns:
#'   \itemize{
#'     \item transcript_id: Transcript identifier
#'     \item transcript_name: Transcript name
#'     \item gene_id: Gene identifier
#'     \item gene_name: Gene name
#'     \item transcript_type: Type of transcript (e.g., protein_coding)
#'     \item gene_type: Type of gene
#'   }
#'
#' @examples
#' \dontrun{
#' # With version number removal (default)
#' transcript_info <- create_transcript_info("path/to/gencode.gtf")
#'
#' # Without version number removal
#' transcript_info <- create_transcript_info("path/to/gencode.gtf", 
#'                                         remove_version = FALSE)
#' }
#'
#' @importFrom rtracklayer import
#' @importFrom tools file_ext file_path_sans_ext
#' @export
create_transcript_info <- function(gtf_path, 
                                   remove_version = TRUE,
                                   progress = TRUE) {
  # Input validation
  if (!file.exists(gtf_path)) {
    stop("GTF file not found at specified path: ", gtf_path)
  }
  
  # Check file extension
  file_ext <- tolower(tools::file_ext(gtf_path))
  if (!file_ext %in% c("gtf", "gz")) {
    stop("File must be a GTF file (*.gtf or *.gtf.gz)")
  }
  
  # For .gz files, check if it's actually a GTF file
  if (file_ext == "gz") {
    base_name <- tools::file_path_sans_ext(gtf_path)
    if (tolower(tools::file_ext(base_name)) != "gtf") {
      stop("Compressed file must be a GTF file (*.gtf.gz)")
    }
  }
  
  if (!is.logical(remove_version)) {
    stop("remove_version must be TRUE or FALSE")
  }
  
  # Show progress message
  if (progress) message("Importing GTF file...")
  
  # Import GTF file
  tryCatch({
    gtf_data <- rtracklayer::import(gtf_path, format = "gtf")
  }, error = function(e) {
    stop("Error importing GTF file: ", e$message)
  })
  
  if (progress) message("Creating transcript information data frame...")
  
  # Extract transcript information
  gtf_transcripts <- gtf_data[gtf_data$type == "transcript", ]
  transcript_info <- data.frame(
    transcript_id = rtracklayer::mcols(gtf_transcripts)$transcript_id,
    transcript_name = rtracklayer::mcols(gtf_transcripts)$transcript_name,
    gene_id = rtracklayer::mcols(gtf_transcripts)$gene_id,
    gene_name = rtracklayer::mcols(gtf_transcripts)$gene_name,
    transcript_type = rtracklayer::mcols(gtf_transcripts)$transcript_type,
    gene_type = rtracklayer::mcols(gtf_transcripts)$gene_type
  )
  
  # Remove version numbers if requested
  if (remove_version) {
    if (progress) message("Removing version numbers from identifiers...")
    transcript_info$transcript_id <- gsub("\\.\\d+$", "", transcript_info$transcript_id)
    transcript_info$gene_id <- gsub("\\.\\d+$", "", transcript_info$gene_id)
  } else {
    # Check if version numbers exist
    has_version_transcript <- any(grepl("\\.\\d+$", transcript_info$transcript_id))
    has_version_gene <- any(grepl("\\.\\d+$", transcript_info$gene_id))
    
    if (has_version_transcript || has_version_gene) {
      warning(
        "Version numbers detected in identifiers. This might cause issues with ID matching.\n",
        "Consider setting remove_version = TRUE to ensure compatibility."
      )
    }
  }
  
  # Set row names
  rownames(transcript_info) <- transcript_info$transcript_id
  
  # Final check for version numbers if they should have been removed
  if (remove_version) {
    has_version_after <- any(grepl("\\.\\d+$", transcript_info$transcript_id)) ||
      any(grepl("\\.\\d+$", transcript_info$gene_id))
    if (has_version_after) {
      warning(
        "Some version numbers may still be present in the identifiers ",
        "despite removal attempt. Please check your data."
      )
    }
  }
  
  if (progress) message("Transcript information processing complete.")
  
  return(transcript_info)
}
