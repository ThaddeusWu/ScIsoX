##################################################################
#  Comprehensive QC Report Generation for ScIsoX                 #
#                                                                #
#  Generate detailed quality control reports from SCHT objects   #
#                                                                #
#  Author: Siyuan Wu & Ulf Schmitz                               #
#  Date: Jul 29, 2025                                            #
#  Package: ScIsoX V1.1.0                                        #
##################################################################

#' Generate comprehensive QC report from SCHT object
#'
#' @description
#' Creates a detailed quality control report including data characteristics,
#' filtering statistics, HVG selection, sparsity analysis, and performance metrics.
#'
#' @param scht_obj SCHT object with enhanced QC information
#' @param output_dir Directory to save report files
#' @param format Output format: "html" or "markdown"
#' @param dataset_name Optional dataset name to include in the output filename (e.g., "blood_data")
#'
#' @return Path to generated report
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
#' # Generate report in temporary directory
#' temp_dir <- tempdir()
#' report_path <- generate_qc_report(
#'   scht_obj = scht_obj,
#'   output_dir = temp_dir,
#'   format = "html",
#'   dataset_name = "blood_cells_example"
#' )
#' 
#' # Check if report was created
#' print(paste("Report created at:", report_path))
#' print(file.exists(report_path))
#' 
#' # List files created
#' list.files(temp_dir, pattern = "blood_cells_example")
#' @export
generate_qc_report <- function(scht_obj, 
                              output_dir = "qc_report",
                              format = c("html", "markdown"),
                              dataset_name = NULL) {
  
  format <- match.arg(format)
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract QC information based on object type
  if (inherits(scht_obj, "IntegratedSCHT")) {
    qc_info <- attr(scht_obj$original_results, "preprocessing")
  } else {
    qc_info <- attr(scht_obj, "preprocessing")
  }
  
  # Check if enhanced QC information is available
  has_enhanced_qc <- !is.null(qc_info$qc_report)
  
  if (!has_enhanced_qc) {
    warning("SCHT object does not contain enhanced QC information. ",
            "Report will be limited. Consider recreating SCHT with latest version.")
  }
  
  # Generate report sections
  report_sections <- list()
  
  # 1. Header and summary
  report_sections$header <- generate_header_section(scht_obj, qc_info)
  
  # 2. Input data characteristics
  report_sections$input_data <- generate_input_data_section(scht_obj, qc_info)
  
  # 3. QC parameters
  report_sections$qc_params <- generate_qc_parameters_section(qc_info)
  
  # 4. Filtering summary
  report_sections$filtering <- generate_filtering_section(qc_info)
  
  # 5. HVG selection
  report_sections$hvg <- generate_hvg_section(qc_info)
  
  # 6. SCHT structure
  report_sections$scht_structure <- generate_scht_structure_section(scht_obj, qc_info)
  
  # 7. Sparsity analysis
  report_sections$sparsity <- generate_sparsity_section(scht_obj)
  
  # 8. Performance metrics
  report_sections$performance <- generate_performance_section(scht_obj)
  
  # Compile report
  if (format == "html") {
    report_path <- compile_html_report(report_sections, output_dir, scht_obj, dataset_name)
  } else if (format == "markdown") {
    report_path <- compile_markdown_report(report_sections, output_dir, dataset_name)
  }
  
  message("QC report generated: ", report_path)
  return(report_path)
}

# Section generators

generate_header_section <- function(scht_obj, qc_info) {
  
  # Get basic information - prefer QC report data for accuracy
  # Get n_genes from QC report first (this accounts for duplicates correctly)
  if (!is.null(qc_info$qc_report$scht_structure$n_genes_in_scht)) {
    n_genes <- qc_info$qc_report$scht_structure$n_genes_in_scht
  } else {
    # Fallback: count from SCHT object
    if (inherits(scht_obj, "IntegratedSCHT")) {
      original_scht <- scht_obj$original_results
      gene_names <- setdiff(names(original_scht), c("n_genes", "n_cells", "n_isoforms", 
                                                     "median_isoforms", "mean_isoforms", 
                                                     "sparsity", "class"))
      n_genes <- length(gene_names)
    } else {
      gene_names <- setdiff(names(scht_obj), c("n_genes", "n_cells", "n_isoforms", 
                                                "median_isoforms", "mean_isoforms", 
                                                "sparsity", "class"))
      n_genes <- length(gene_names)
    }
  }
  
  # Get n_cells from QC report
  if (!is.null(qc_info$qc_report$scht_structure$n_cells_after_qc)) {
    n_cells <- qc_info$qc_report$scht_structure$n_cells_after_qc
  } else if (!is.null(qc_info$qc_stats$n_filtered_cells) && !is.null(qc_info$qc_report$input_data$n_cells_original)) {
    n_cells <- qc_info$qc_report$input_data$n_cells_original - qc_info$qc_stats$n_filtered_cells
  } else {
    n_cells <- "N/A"
  }
  
  header <- c(
    "# ScIsoX Quality Control Report",
    "",
    paste("**Generated on:**", Sys.Date()),
    paste("**ScIsoX version:**", packageVersion("ScIsoX")),
    "",
    "## Summary",
    "",
    paste("- **Genes in SCHT:**", format(n_genes, big.mark = ",")),
    paste("- **Cells after QC:**", format(n_cells, big.mark = ",")),
    ""
  )
  
  return(header)
}

generate_input_data_section <- function(scht_obj, qc_info) {
  
  section <- c(
    "## Input Data Characteristics",
    ""
  )
  
  # Add input data type at the beginning - check multiple locations
  input_type <- NULL
  
  # First check qc_report
  if (!is.null(qc_info$qc_report$input_type)) {
    input_type <- qc_info$qc_report$input_type
  }
  # Then check params
  else if (!is.null(qc_info$params$input_type)) {
    input_type <- qc_info$params$input_type
  }
  # Then check root level
  else if (!is.null(qc_info$input_type)) {
    input_type <- qc_info$input_type
  }
  
  # Format for display
  input_type_display <- ifelse(!is.null(input_type),
                              ifelse(input_type == "raw_counts", 
                                     "Raw Count Matrices", 
                                     "Normalised Count Matrices"),
                              "Unknown")
  
  section <- c(section,
    sprintf("**Input Data Type:** %s", input_type_display),
    ""
  )
  
  if (!is.null(qc_info$qc_report$input_data)) {
    data <- qc_info$qc_report$input_data
    
    section <- c(section,
      "| Metric | Value |",
      "|--------|-------|",
      sprintf("| Original genes | %s |", format(data$n_genes_original, big.mark = ",")),
      sprintf("| Original transcripts | %s |", format(data$n_transcripts_original, big.mark = ",")),
      sprintf("| Original cells | %s |", format(data$n_cells_original, big.mark = ",")),
      sprintf("| Gene matrix sparsity | %.2f%% |", data$gene_matrix_sparsity),
      sprintf("| Transcript matrix sparsity | %.2f%% |", data$transcript_matrix_sparsity),
      sprintf("| Median genes/cell | %s |", format(data$median_genes_per_cell, big.mark = ",")),
      sprintf("| Median transcripts/cell | %s |", format(data$median_transcripts_per_cell, big.mark = ",")),
      ""
    )
    
    # Add cell type distribution if available
    if (!is.null(data$cell_type_distribution)) {
      section <- c(section,
        "",
        "### Cell Type Distribution",
        sprintf("- Number of cell types: %d", data$cell_type_distribution$n_cell_types),
        "",
        "| Cell Type | Count | Percentage |",
        "|-----------|-------|------------|"
      )
      
      for (ct in names(data$cell_type_distribution$counts)) {
        section <- c(section,
          sprintf("| %s | %s | %.2f%% |", 
                  ct, 
                  format(data$cell_type_distribution$counts[[ct]], big.mark = ","), 
                  data$cell_type_distribution$percentages[[ct]]))
      }
      section <- c(section, "")
    }
  } else {
    # Fallback for older SCHT objects
    section <- c(section,
      "*Enhanced QC information not available for input data.*",
      ""
    )
  }
  
  return(section)
}

generate_qc_parameters_section <- function(qc_info) {
  
  section <- c(
    "## QC Parameters",
    ""
  )
  
  if (!is.null(qc_info$qc_report$qc_parameters)) {
    params <- qc_info$qc_report$qc_parameters
    
    # Create table for parameters with all three strategies
    section <- c(section,
      "| Parameter | Applied Value | MAD Strategy | Interval 90 | Interval 80 |",
      "|-----------|---------------|--------------|-------------|-------------|"
    )
    
    # Min genes per cell
    mad_min <- ifelse(!is.null(params$recommended) && !is.null(params$recommended$MAD_strategy), 
                      as.character(params$recommended$MAD_strategy$min_genes_per_cell), "N/A")
    int90_min <- ifelse(!is.null(params$recommended) && !is.null(params$recommended$Interval_90), 
                        as.character(params$recommended$Interval_90$min_genes_per_cell), "N/A")
    int80_min <- ifelse(!is.null(params$recommended) && !is.null(params$recommended$Interval_80), 
                        as.character(params$recommended$Interval_80$min_genes_per_cell), "N/A")
    
    section <- c(section,
      sprintf("| Min genes per cell | %d | %s | %s | %s |", 
              params$applied$min_genes_per_cell, 
              mad_min, int90_min, int80_min)
    )
    
    # Max genes per cell
    mad_max <- ifelse(!is.null(params$recommended) && !is.null(params$recommended$MAD_strategy), 
                      as.character(params$recommended$MAD_strategy$max_genes_per_cell), "N/A")
    int90_max <- ifelse(!is.null(params$recommended) && !is.null(params$recommended$Interval_90), 
                        as.character(params$recommended$Interval_90$max_genes_per_cell), "N/A")
    int80_max <- ifelse(!is.null(params$recommended) && !is.null(params$recommended$Interval_80), 
                        as.character(params$recommended$Interval_80$max_genes_per_cell), "N/A")
    
    section <- c(section,
      sprintf("| Max genes per cell | %d | %s | %s | %s |", 
              params$applied$max_genes_per_cell, 
              mad_max, int90_max, int80_max)
    )
    
    # Min cells expressing
    section <- c(section,
      sprintf("| Min cells expressing | %.2f%% | - | - | - |", 
              params$applied$min_cells_expressing * 100)
    )
    
    # Min expression
    section <- c(section,
      sprintf("| Min expression | %.1e | - | - | - |", params$applied$min_expr)
    )
    
    section <- c(section, "", 
      "### Strategy Explanations:",
      "",
      "- MAD Strategy: Uses median +/- 3 MAD, reduces risk of including poor quality cells whilst maintaining robustness to outliers",
      "- Interval 90: Uses 5th and 95th percentiles, balances stringency with dataset preservation", 
      "- Interval 80: Uses 10th and 90th percentiles, provides more aggressive filtering for higher quality cell selection",
      ""
    )
    
  } else {
    # Use basic info in table format
    section <- c(section,
      "| Parameter | Value |",
      "|-----------|-------|",
      sprintf("| Min genes per cell | %d |", qc_info$params$min_genes_per_cell),
      sprintf("| Max genes per cell | %d |", qc_info$params$max_genes_per_cell),
      ""
    )
  }
  
  return(section)
}

generate_filtering_section <- function(qc_info) {
  
  section <- c(
    "## Filtering Summary",
    ""
  )
  
  # Basic filtering info in table
  section <- c(section,
    "| Category | Count Removed |",
    "|----------|---------------|",
    sprintf("| Genes | %s |", format(qc_info$qc_stats$n_filtered_genes, big.mark = ",")),
    sprintf("| Transcripts | %s |", format(qc_info$qc_stats$n_filtered_transcripts, big.mark = ",")),
    sprintf("| Cells | %s |", format(qc_info$qc_stats$n_filtered_cells, big.mark = ",")),
    ""
  )
  
  # Detailed reasons if available
  if (!is.null(qc_info$qc_report$filtering$cell_qc$removal_reasons)) {
    reasons <- qc_info$qc_report$filtering$cell_qc$removal_reasons
    
    section <- c(section,
      "### Cell Removal Reasons",
      "",
      "| Reason | Cell Count | Percentage |",
      "|--------|------------|------------|",
      sprintf("| Too few genes | %s | %.2f%% |", 
              format(reasons$too_few_genes, big.mark = ","), reasons$percent_too_few),
      sprintf("| Too many genes | %s | %.2f%% |", 
              format(reasons$too_many_genes, big.mark = ","), reasons$percent_too_many),
      ""
    )
  }
  
  return(section)
}

generate_hvg_section <- function(qc_info) {
  
  section <- c(
    "## Highly Variable Gene Selection",
    ""
  )
  
  # Main HVG stats table
  section <- c(section,
    "| Metric | Value |",
    "|--------|-------|",
    sprintf("| HVGs requested | %s |", format(qc_info$params$n_hvg, big.mark = ",")),
    sprintf("| HVGs selected | %s |", format(length(qc_info$hvg), big.mark = ",")),
    ""
  )
  
  if (!is.null(qc_info$qc_report$hvg_selection)) {
    hvg_data <- qc_info$qc_report$hvg_selection
    
    section <- c(section,
      "### HVG Filtering Details",
      "",
      "| Description | Count | Status |",
      "|-------------|-------|--------|",
      sprintf("| Total genes available after QC | %s | - |", format(hvg_data$total_genes_available, big.mark = ",")),
      sprintf("| HVGs with single isoform | %s | Removed |", 
              format(ifelse(is.na(hvg_data$hvg_with_single_isoform), 0, hvg_data$hvg_with_single_isoform), big.mark = ",")),
      sprintf("| HVGs with multiple isoforms | %s | Kept |", 
              format(ifelse(is.na(hvg_data$hvg_with_multiple_isoforms), 0, hvg_data$hvg_with_multiple_isoforms), big.mark = ",")),
      sprintf("| Percentage multi-isoform HVGs | %.2f%% | - |", 
              ifelse(is.na(hvg_data$percent_multi_isoform), 0, hvg_data$percent_multi_isoform)),
      sprintf("| Final genes in SCHT | %s | - |", format(hvg_data$hvg_in_scht, big.mark = ",")),
      ""
    )
  }
  
  return(section)
}

generate_scht_structure_section <- function(scht_obj, qc_info) {
  
  section <- c(
    "## SCHT Structure",
    ""
  )
  
  if (!is.null(qc_info$qc_report$scht_structure)) {
    struct <- qc_info$qc_report$scht_structure
    iso_stats <- struct$isoform_stats
    
    section <- c(section,
      "| Metric | Value |",
      "|--------|-------|",
      sprintf("| Genes in SCHT | %s |", format(struct$n_genes_in_scht, big.mark = ",")),
      sprintf("| Cells after QC | %s |", format(struct$n_cells_after_qc, big.mark = ",")),
      sprintf("| Total isoforms | %s |", format(iso_stats$total_isoforms, big.mark = ",")),
      sprintf("| Max isoforms per gene | %d |", iso_stats$max_isoforms_per_gene),
      sprintf("| Mean isoforms per gene | %.2f |", iso_stats$mean_isoforms_per_gene),
      ""
    )
  } else {
    # Calculate from SCHT based on object type
    if (inherits(scht_obj, "IntegratedSCHT")) {
      original_scht <- scht_obj$original_results
      gene_names <- setdiff(names(original_scht), c("n_genes", "n_cells", "n_isoforms", 
                                                     "median_isoforms", "mean_isoforms", 
                                                     "sparsity", "class"))
      n_genes <- length(gene_names)
    } else {
      gene_names <- setdiff(names(scht_obj), c("n_genes", "n_cells", "n_isoforms", 
                                                "median_isoforms", "mean_isoforms", 
                                                "sparsity", "class"))
      n_genes <- length(gene_names)
    }
    section <- c(section,
      "| Metric | Value |",
      "|--------|-------|",
      sprintf("| Genes in SCHT | %s |", format(n_genes, big.mark = ",")),
      "| Enhanced info | Not available |",
      ""
    )
  }
  
  return(section)
}

generate_sparsity_section <- function(scht_obj) {
  
  section <- c(
    "## Sparsity Analysis",
    ""
  )
  
  # Check if comprehensive sparsity analysis is available
  if (inherits(scht_obj, "IntegratedSCHT")) {
    qc_info <- attr(scht_obj$original_results, "preprocessing")
  } else {
    qc_info <- attr(scht_obj, "preprocessing")
  }
  if (!is.null(qc_info$qc_report$sparsity_analysis)) {
    sparsity_data <- qc_info$qc_report$sparsity_analysis
    
    section <- c(section,
      "### Comprehensive Sparsity Comparison",
      "",
      "| Matrix Type | Elements | Non-zeros | Zeros | Sparsity % |",
      "|-------------|----------|-----------|-------|------------|"
    )
    
    # Original matrix
    if (!is.null(sparsity_data$original)) {
      orig <- sparsity_data$original
      section <- c(section,
        sprintf("| Original Transcript Matrix | %s | %s | %s | %.2f%% |",
                format(orig$total, big.mark = ","),
                format(orig$nonzero, big.mark = ","),
                format(orig$zeros, big.mark = ","),
                orig$sparsity))
    }
    
    # Filtered matrix
    if (!is.null(sparsity_data$filtered)) {
      filt <- sparsity_data$filtered
      section <- c(section,
        sprintf("| Filtered Matrix (Post-QC) | %s | %s | %s | %.2f%% |",
                format(filt$total, big.mark = ","),
                format(filt$nonzero, big.mark = ","),
                format(filt$zeros, big.mark = ","),
                filt$sparsity))
    }
    
    # Naive 3D tensor
    if (!is.null(sparsity_data$tensor)) {
      tens <- sparsity_data$tensor
      section <- c(section,
        sprintf("| Naive 3D Tensor | %s | %s | %s | %.2f%% |",
                format(tens$total, big.mark = ","),
                format(tens$nonzero, big.mark = ","),
                format(tens$zeros, big.mark = ","),
                tens$sparsity))
    }
    
    # SCHT structure
    if (!is.null(sparsity_data$scht)) {
      scht <- sparsity_data$scht
      section <- c(section,
        sprintf("| SCHT Structure | %s | %s | %s | %.2f%% |",
                format(scht$total, big.mark = ","),
                format(scht$nonzero, big.mark = ","),
                format(scht$zeros, big.mark = ","),
                scht$sparsity))
    }
    
    section <- c(section, "")
    
    # Memory efficiency
    if (!is.null(sparsity_data$zeros_avoided)) {
      section <- c(section,
        "### Zero Padding Reduction",
        "",
        "| Comparison | Zero Elements Avoided |",
        "|------------|----------------------|",
        sprintf("| vs Original Matrix | %s |", 
                format(sparsity_data$zeros_avoided$vs_original, big.mark = ",")),
        sprintf("| vs Filtered Matrix | %s |", 
                format(sparsity_data$zeros_avoided$vs_filtered, big.mark = ",")),
        sprintf("| vs Naive 3D Tensor | %s |", 
                format(sparsity_data$zeros_avoided$vs_tensor, big.mark = ",")),
        ""
      )
    }
    
  } else {
    # Fallback to basic sparsity calculation
    if (exists("calculate_scht_sparsity")) {
      sparsity <- tryCatch({
        calculate_scht_sparsity(scht_obj)
      }, error = function(e) {
        list(n_total = NA, n_nonzero = NA, sparsity = NA)
      })
      
      section <- c(section,
        sprintf("- Total elements in SCHT: %s", format(sparsity$n_total, big.mark = ",")),
        sprintf("- Non-zero elements: %s", format(sparsity$n_nonzero, big.mark = ",")),
        sprintf("- Sparsity: %.2f%%", sparsity$sparsity),
        ""
      )
    } else {
      section <- c(section,
        "*Sparsity analysis function not available.*",
        ""
      )
    }
  }
  
  return(section)
}

generate_performance_section <- function(scht_obj) {
  
  section <- c(
    "## Performance Metrics",
    ""
  )
  
  # Performance metrics are stored on the main object
  perf <- attr(scht_obj, "performance")
  
  if (!is.null(perf)) {
    section <- c(section,
      "| Metric | Value |",
      "|--------|-------|",
      sprintf("| Total processing time | %.2f seconds (%.2f minutes) |", 
              perf$total_time_sec, perf$total_time_sec / 60),
      sprintf("| Memory used | %.2f MB |", perf$memory_used_mb),
      ""
    )
  } else {
    section <- c(section,
      "*Performance metrics not available.*",
      ""
    )
  }
  
  return(section)
}

# Compilation functions

compile_markdown_report <- function(sections, output_dir, dataset_name = NULL) {
  
  # Combine all sections
  full_report <- c()
  
  for (section_name in names(sections)) {
    full_report <- c(full_report, sections[[section_name]], "")
  }
  
  # Write to file
  # Generate filename based on dataset_name
  if (!is.null(dataset_name) && dataset_name != "") {
    output_filename <- paste0("qc_report_", dataset_name, ".md")
  } else {
    output_filename <- "qc_report.md"
  }
  output_file <- file.path(output_dir, output_filename)
  writeLines(full_report, output_file)
  
  return(output_file)
}

compile_html_report <- function(sections, output_dir, scht_obj = NULL, dataset_name = NULL) {
  
  # Get HTML styles
  style_content <- get_html_style()
  
  # Create HTML structure
  html_content <- c(
    "<!DOCTYPE html>",
    "<html lang='en'>",
    "<head>",
    "<meta charset='UTF-8'>",
    "<meta name='viewport' content='width=device-width, initial-scale=1.0'>",
    "<title>ScIsoX Quality Control Report</title>",
    style_content,
    "</head>",
    "<body>",
    "<div class='container'>",
    "<div class='report-header'>",
    "<h1>ScIsoX Quality Control Report</h1>",
    "<div class='report-meta'>",
    paste0("Generated on ", Sys.Date(), " | ScIsoX v", packageVersion("ScIsoX")),
    "</div>",
    "</div>"
  )
  
  # Add summary cards if we have the data
  if (!is.null(scht_obj)) {
    # Extract key metrics from SCHT object
    if (inherits(scht_obj, "IntegratedSCHT")) {
      qc_info <- attr(scht_obj$original_results, "preprocessing")
    } else {
      qc_info <- attr(scht_obj, "preprocessing")
    }
    
    # Get basic stats
    # Get n_genes from QC report first (this accounts for duplicates correctly)
    if (!is.null(qc_info$qc_report$scht_structure$n_genes_in_scht)) {
      n_genes <- qc_info$qc_report$scht_structure$n_genes_in_scht
    } else {
      # Fallback: count from SCHT object
      if (inherits(scht_obj, "IntegratedSCHT")) {
        original_scht <- scht_obj$original_results
        gene_names <- setdiff(names(original_scht), c("n_genes", "n_cells", "n_isoforms", 
                                                       "median_isoforms", "mean_isoforms", 
                                                       "sparsity", "class"))
        n_genes <- length(gene_names)
      } else {
        gene_names <- setdiff(names(scht_obj), c("n_genes", "n_cells", "n_isoforms", 
                                                  "median_isoforms", "mean_isoforms", 
                                                  "sparsity", "class"))
        n_genes <- length(gene_names)
      }
    }
    
    # Get n_cells from QC report
    if (!is.null(qc_info$qc_report$scht_structure$n_cells_after_qc)) {
      n_cells <- qc_info$qc_report$scht_structure$n_cells_after_qc
    } else if (!is.null(qc_info$qc_stats$n_filtered_cells) && !is.null(qc_info$qc_report$input_data$n_cells_original)) {
      n_cells <- qc_info$qc_report$input_data$n_cells_original - qc_info$qc_stats$n_filtered_cells
    } else {
      n_cells <- "N/A"
    }
    
    # Get SCHT list for isoform counting
    if (inherits(scht_obj, "IntegratedSCHT")) {
      original_scht <- scht_obj$original_results
      gene_names <- setdiff(names(original_scht), c("n_genes", "n_cells", "n_isoforms", 
                                                     "median_isoforms", "mean_isoforms", 
                                                     "sparsity", "class"))
      scht_list <- original_scht[gene_names]
    } else {
      gene_names <- setdiff(names(scht_obj), c("n_genes", "n_cells", "n_isoforms", 
                                                "median_isoforms", "mean_isoforms", 
                                                "sparsity", "class"))
      scht_list <- scht_obj[gene_names]
    }
    
    # Calculate total isoforms
    total_isoforms <- 0
    if (!is.null(scht_list)) {
      for (gene in names(scht_list)) {
        gene_matrix <- scht_list[[gene]]
        if (!is.null(gene_matrix) && length(gene_matrix) > 0) {
          total_isoforms <- total_isoforms + nrow(gene_matrix)
        }
      }
    }
    
    # Get sparsity from QC report or calculate it
    scht_sparsity <- "N/A"
    if (!is.null(qc_info$qc_report$sparsity_analysis$scht$sparsity)) {
      scht_sparsity <- sprintf("%.2f%%", qc_info$qc_report$sparsity_analysis$scht$sparsity)
    } else if (!is.null(qc_info$qc_report$scht_structure$sparsity$sparsity)) {
      scht_sparsity <- sprintf("%.2f%%", qc_info$qc_report$scht_structure$sparsity$sparsity)
    } else if (!is.null(qc_info$qc_report$scht_structure$sparsity)) {
      # Sometimes it's stored directly as a number
      scht_sparsity <- sprintf("%.2f%%", qc_info$qc_report$scht_structure$sparsity)
    }
    
    html_content <- c(html_content,
      "<div class='summary-grid'>",
      "<div class='summary-card'>",
      "<div class='label'>Genes in SCHT</div>",
      "<div class='value'>", format(n_genes, big.mark = ","), "</div>",
      "</div>",
      "<div class='summary-card'>",
      "<div class='label'>Cells after QC</div>",
      "<div class='value'>", format(n_cells, big.mark = ","), "</div>",
      "</div>",
      "<div class='summary-card'>",
      "<div class='label'>Total Isoforms</div>",
      "<div class='value'>", format(total_isoforms, big.mark = ","), "</div>",
      "</div>",
      "<div class='summary-card'>",
      "<div class='label'>SCHT Sparsity</div>",
      "<div class='value'>", scht_sparsity, "</div>",
      "</div>",
      "</div>"
    )
  }
  
  # Process each section
  for (section_name in names(sections)) {
    if (section_name != "header") {
      # Convert markdown to HTML
      in_list <- FALSE
      in_table <- FALSE
      
      for (line in sections[[section_name]]) {
        # Handle lists
        if (startsWith(line, "- ") && !in_list) {
          html_content <- c(html_content, "<ul>")
          in_list <- TRUE
        } else if (!startsWith(line, "- ") && in_list && line != "") {
          html_content <- c(html_content, "</ul>")
          in_list <- FALSE
        }
        
        # Process line
        if (line == "") {
          if (in_list) {
            html_content <- c(html_content, "</ul>")
            in_list <- FALSE
          }
          if (in_table) {
            html_content <- c(html_content, "</table>")
            in_table <- FALSE
          }
          html_content <- c(html_content, "")
        } else if (startsWith(line, "# ")) {
          html_content <- c(html_content, paste0("<h1>", substring(line, 3), "</h1>"))
        } else if (startsWith(line, "## ")) {
          html_content <- c(html_content, paste0("<h2>", substring(line, 4), "</h2>"))
        } else if (startsWith(line, "### ")) {
          html_content <- c(html_content, paste0("<h3>", substring(line, 5), "</h3>"))
        } else if (startsWith(line, "- ")) {
          html_content <- c(html_content, paste0("<li>", substring(line, 3), "</li>"))
        } else if (startsWith(line, "|")) {
          # Handle table rows
          if (grepl("\\|[-:]+\\|", line)) {
            # Header separator - mark table start
            if (!in_table) {
              html_content <- c(html_content, "<table>")
              in_table <- TRUE
            }
          } else {
            # Table row
            cells <- strsplit(line, "\\|")[[1]]
            cells <- cells[cells != ""]  # Remove empty cells
            cells <- trimws(cells)
            
            if (!in_table) {
              html_content <- c(html_content, "<table>")
              in_table <- TRUE
            }
            
            # Check if this is header row (comes before separator)
            if (length(cells) > 0) {
              next_line_idx <- which(sections[[section_name]] == line) + 1
              is_header <- FALSE
              if (next_line_idx <= length(sections[[section_name]])) {
                next_line <- sections[[section_name]][next_line_idx]
                if (grepl("\\|[-:]+\\|", next_line)) {
                  is_header <- TRUE
                }
              }
              
              html_content <- c(html_content, "<tr>")
              for (cell in cells) {
                if (is_header) {
                  html_content <- c(html_content, paste0("<th>", cell, "</th>"))
                } else {
                  html_content <- c(html_content, paste0("<td>", cell, "</td>"))
                }
              }
              html_content <- c(html_content, "</tr>")
            }
          }
        } else {
          # Convert markdown bold to HTML
          processed_line <- gsub("\\*\\*(.+?)\\*\\*", "<strong>\\1</strong>", line)
          html_content <- c(html_content, paste0("<p>", processed_line, "</p>"))
        }
      }
      
      # Close any open lists or tables
      if (in_list) {
        html_content <- c(html_content, "</ul>")
      }
      if (in_table) {
        html_content <- c(html_content, "</table>")
      }
    }
  }
  
  # Add footer
  html_content <- c(html_content,
    "<div class='footer'>",
    paste0("<p>Generated on ", Sys.Date(), " using <a href='https://github.com/ThaddeusWu/ScIsoX'>ScIsoX</a> v", 
           packageVersion("ScIsoX"), "</p>"),
    "<p>Single-cell Transcriptomic Complexity Analysis</p>",
    "</div>",
    "</div>", # Close container
    "</body>",
    "</html>"
  )
  
  # Write HTML file
  # Generate filename based on dataset_name
  if (!is.null(dataset_name) && dataset_name != "") {
    html_filename <- paste0("qc_report_", dataset_name, ".html")
  } else {
    html_filename <- "qc_report.html"
  }
  html_file <- file.path(output_dir, html_filename)
  writeLines(html_content, html_file)
  
  return(html_file)
}

# Enhanced HTML styling for QC reports
get_html_style <- function() {
  return(c(
    "<style>",
    "/* Reset and base styles */",
    "* { margin: 0; padding: 0; box-sizing: border-box; }",
    "body {",
    "  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;",
    "  line-height: 1.6;",
    "  color: #2c3e50;",
    "  background-color: #f5f7fa;",
    "}",
    "",
    "/* Container */",
    ".container {",
    "  max-width: 1200px;",
    "  margin: 0 auto;",
    "  padding: 20px;",
    "  background-color: white;",
    "  min-height: 100vh;",
    "}",
    "",
    "/* Header styles */",
    ".report-header {",
    "  background: linear-gradient(135deg, #6B73FF 0%, #000DFF 100%);",
    "  color: white;",
    "  padding: 60px 40px;",
    "  margin: -20px -20px 40px -20px;",
    "  text-align: center;",
    "  box-shadow: 0 4px 20px rgba(0,0,0,0.1);",
    "}",
    "",
    ".report-header h1 {",
    "  font-size: 3em;",
    "  font-weight: 300;",
    "  margin-bottom: 15px;",
    "  color: white !important;",
    "  text-shadow: 2px 2px 4px rgba(0,0,0,0.2);",
    "}",
    "",
    ".report-meta {",
    "  font-size: 1.1em;",
    "  opacity: 0.95;",
    "  margin-top: 10px;",
    "}",
    "",
    "/* Section styles */",
    "h1 {",
    "  color: #2c3e50;",
    "  font-size: 2.2em;",
    "  font-weight: 400;",
    "  margin: 50px 0 25px 0;",
    "  padding-bottom: 15px;",
    "  border-bottom: 3px solid #6B73FF;",
    "  position: relative;",
    "}",
    "",
    "h2 {",
    "  color: #34495e;",
    "  font-size: 1.6em;",
    "  font-weight: 500;",
    "  margin: 35px 0 20px 0;",
    "  padding-left: 15px;",
    "  border-left: 4px solid #6B73FF;",
    "}",
    "",
    "h3 {",
    "  color: #555;",
    "  font-size: 1.3em;",
    "  font-weight: 600;",
    "  margin: 25px 0 15px 0;",
    "}",
    "",
    "/* Summary cards */",
    ".summary-grid {",
    "  display: grid;",
    "  grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));",
    "  gap: 20px;",
    "  margin: 30px 0;",
    "}",
    "",
    ".summary-card {",
    "  background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);",
    "  padding: 25px;",
    "  border-radius: 12px;",
    "  text-align: center;",
    "  box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
    "  transition: transform 0.3s ease;",
    "}",
    "",
    ".summary-card:hover {",
    "  transform: translateY(-5px);",
    "  box-shadow: 0 6px 20px rgba(0,0,0,0.12);",
    "}",
    "",
    ".summary-card .value {",
    "  font-size: 2.5em;",
    "  font-weight: 700;",
    "  color: #6B73FF;",
    "  margin: 10px 0;",
    "}",
    "",
    ".summary-card .label {",
    "  color: #666;",
    "  font-size: 0.95em;",
    "  font-weight: 500;",
    "  text-transform: uppercase;",
    "  letter-spacing: 0.5px;",
    "}",
    "",
    "/* Tables */",
    "table {",
    "  width: 100%;",
    "  border-collapse: collapse;",
    "  margin: 25px 0;",
    "  background: white;",
    "  box-shadow: 0 2px 10px rgba(0,0,0,0.08);",
    "  border-radius: 10px;",
    "  overflow: hidden;",
    "}",
    "",
    "th {",
    "  background: linear-gradient(135deg, #6B73FF 0%, #5a63e6 100%);",
    "  color: white;",
    "  font-weight: 600;",
    "  padding: 15px 20px;",
    "  text-align: left;",
    "  font-size: 0.95em;",
    "  text-transform: uppercase;",
    "  letter-spacing: 0.5px;",
    "}",
    "",
    "td {",
    "  padding: 12px 20px;",
    "  border-bottom: 1px solid #e8ecf1;",
    "  color: #555;",
    "}",
    "",
    "tr:last-child td {",
    "  border-bottom: none;",
    "}",
    "",
    "tr:nth-child(even) {",
    "  background-color: #f8f9fc;",
    "}",
    "",
    "tr:hover {",
    "  background-color: #eef2f7;",
    "  transition: background-color 0.2s;",
    "}",
    "",
    "/* Number formatting in tables */",
    "td:nth-child(n+2) {",
    "  text-align: right;",
    "  font-family: 'SF Mono', Monaco, 'Courier New', monospace;",
    "}",
    "",
    "/* Lists */",
    "ul, ol {",
    "  margin: 20px 0 20px 40px;",
    "}",
    "",
    "li {",
    "  margin: 10px 0;",
    "  color: #555;",
    "}",
    "",
    "li strong {",
    "  color: #2c3e50;",
    "}",
    "",
    "/* Info sections */",
    ".info-section {",
    "  background: #f8f9fc;",
    "  padding: 25px;",
    "  border-radius: 10px;",
    "  margin: 25px 0;",
    "  border-left: 4px solid #6B73FF;",
    "}",
    "",
    ".info-section h3 {",
    "  margin-top: 0;",
    "  color: #2c3e50;",
    "}",
    "",
    "/* Code blocks */",
    "code {",
    "  background-color: #f0f3f7;",
    "  padding: 3px 8px;",
    "  border-radius: 4px;",
    "  font-family: 'SF Mono', Monaco, 'Courier New', monospace;",
    "  font-size: 0.9em;",
    "  color: #e74c3c;",
    "}",
    "",
    "/* Paragraphs */",
    "p {",
    "  margin: 15px 0;",
    "  color: #555;",
    "  line-height: 1.7;",
    "}",
    "",
    "/* Footer */",
    ".footer {",
    "  margin-top: 80px;",
    "  padding: 30px 0;",
    "  border-top: 2px solid #e8ecf1;",
    "  text-align: center;",
    "  color: #999;",
    "  font-size: 0.9em;",
    "}",
    "",
    ".footer a {",
    "  color: #6B73FF;",
    "  text-decoration: none;",
    "}",
    "",
    ".footer a:hover {",
    "  text-decoration: underline;",
    "}",
    "",
    "/* Responsive design */",
    "@media (max-width: 768px) {",
    "  .container { padding: 10px; }",
    "  .report-header { padding: 30px 20px; margin: -10px -10px 30px -10px; }",
    "  .report-header h1 { font-size: 2em; }",
    "  .summary-grid { grid-template-columns: 1fr; }",
    "  table { font-size: 0.85em; }",
    "  th, td { padding: 10px; }",
    "}",
    "",
    "/* Print styles */",
    "@media print {",
    "  body { background: white; }",
    "  .container { box-shadow: none; max-width: 100%; }",
    "  .report-header { ",
    "    background: none !important; ",
    "    color: #2c3e50 !important; ",
    "    border-bottom: 3px solid #2c3e50;",
    "    print-color-adjust: exact;",
    "    -webkit-print-color-adjust: exact;",
    "  }",
    "  .summary-card { break-inside: avoid; }",
    "  table { break-inside: avoid; }",
    "}",
    "</style>"
  ))
}

