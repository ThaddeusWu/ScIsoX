#' Blood dataset gene expression counts
#'
#' Gene-level expression count matrix from blood single-cell RNA-seq data.
#' This data is derived from Wang et al. (2022) study on alternative splicing 
#' during hematopoietic stem cell formation.
#'
#' @name gene_counts_blood
#' @docType data
#' @format A data frame
#' @source Wang et al. (2022). Science Advances 8(1):eabg5369. Data from https://zenodo.org/records/5706781. Licensed under CC BY 4.0.
#' @references Wang F, et al. (2022). Single-cell architecture and functional requirement of alternative splicing during hematopoietic stem cell formation. Science Advances 8(1):eabg5369.
#' @keywords datasets
NULL

#' Blood dataset transcript expression counts
#'
#' Transcript-level expression count matrix from blood single-cell RNA-seq data.
#' This data is derived from Wang et al. (2022) study on alternative splicing 
#' during hematopoietic stem cell formation.
#'
#' @name transcript_counts_blood
#' @docType data
#' @format A data frame
#' @source Wang et al. (2022). Science Advances 8(1):eabg5369. Data from https://zenodo.org/records/5706781. Licensed under CC BY 4.0.
#' @references Wang F, et al. (2022). Single-cell architecture and functional requirement of alternative splicing during hematopoietic stem cell formation. Science Advances 8(1):eabg5369.
#' @keywords datasets
NULL

#' Blood dataset cell metadata
#'
#' Cell metadata including cell type annotations for the blood dataset.
#'
#' @name sample2stage
#' @docType data
#' @format A data frame with columns:
#' \describe{
#'   \item{sample}{Sample identifier}
#'   \item{cell_type}{Cell type annotation}
#' }
#' @source Wang et al. (2022). Science Advances 8(1):eabg5369. Data from https://zenodo.org/records/5706781. Licensed under CC BY 4.0.
#' @references Wang F, et al. (2022). Single-cell architecture and functional requirement of alternative splicing during hematopoietic stem cell formation. Science Advances 8(1):eabg5369.
#' @keywords datasets
NULL

#' Transcript information
#'
#' Transcript and gene annotation information from mouse genome.
#' This data was processed using the create_transcript_info() function from ScIsoX package.
#'
#' @name transcript_info
#' @docType data
#' @format A data frame with columns:
#' \describe{
#'   \item{transcript_id}{Transcript identifier}
#'   \item{transcript_name}{Transcript name}
#'   \item{gene_id}{Gene identifier}
#'   \item{gene_name}{Gene symbol}
#'   \item{transcript_type}{Transcript biotype}
#'   \item{gene_type}{Gene biotype}
#' }
#' @source GENCODE Mouse Release M36 from https://www.gencodegenes.org/mouse/release_M36.html. Processed using create_transcript_info() function from ScIsoX package.
#' @seealso \code{\link{create_transcript_info}} for the function used to generate this data
#' @keywords datasets
NULL