# ScIsoX 1.1.1 (2025-08-01)

## Major Improvements

* **Runnable Examples**: Replaced all `\dontrun{}` blocks with fully executable examples using built-in datasets
* **Example Datasets**: Added four example datasets from Wang et al. (2022) Science Advances 8(1):eabg5369 (CC BY 4.0)
  - `gene_counts_blood`: Gene expression matrix
  - `transcript_counts_blood`: Transcript expression matrix
  - `transcript_info`: Transcript-to-gene mapping
  - `sample2stage`: Cell type annotations
* **Comprehensive Vignettes**: Created 5 detailed vignettes covering all aspects of the package
  - Getting Started with ScIsoX
  - Data Import and Quality Control
  - Understanding Transcriptomic Complexity Metrics
  - Visualisation Gallery and Interpretation
  - Co-expression Analysis and Isoform Switching
* **pkgdown Website**: Built comprehensive documentation website with improved design and navigation

## Bug Fixes

* Fixed `plot_isoform_profile()` error when no multi-isoform genes exist in IntegratedSCHT
* Fixed duplicate `plot_isoform_coexpression()` function definition
* Fixed `plot_complexity_radar()` parameter name mismatch (`gene_colors` vs `gene_colours`)
* Fixed incorrect classification names in complexity metrics (e.g., "high"/"low" → proper category names)
* Fixed NAMESPACE typo: `CellTypeSCHTSCHT` → `CellTypeSCHT`
* Replaced deprecated `class()` usage with `inherits()` for S3 class checking
* Fixed vignette errors with incorrect function parameters

## Documentation

* Added comprehensive Input Data Recommendations section to README
* Added Visualisation Gallery showcasing all plot types
* Added 5-Minute Quick Start guide
* All functions now have runnable examples
* Improved function documentation with clearer parameter descriptions
* Added proper attribution for example data throughout documentation

# ScIsoX 1.1.0 (2025-07-29)

## Bug Fixes

* Fixed issue in `create_scht()` where the function would fail when gene count matrix used gene names instead of gene IDs (#2). The function now automatically detects whether gene IDs or gene names are being used and handles both cases appropriately.
* Fixed issue in `plot_compare_multiple_genes_radar_cell_type()` where the function would fail when encountering NA values (#1). The function now removes genes/cell types with NA values and provides appropriate warnings.
* Fixed complexity metrics showing 100% NA for cell type-related metrics due to field name mismatch between `cell_type_specific` and `cell_type_matrices`
* Fixed conservation analysis Mixed pattern detection logic to check Mixed patterns before Conserved patterns
* Fixed sparsity analysis to correctly count isoforms in filtered matrix when using mixed gene names/IDs in SCHT
* Resolved gene name conflicts in SCHT by using gene IDs for conflicted names while preserving gene names for unique entries

## Enhancements

* Added graceful handling of missing suggested packages with informative messages
* Added `install_scisox_suggests()` function for easy installation of all optional packages
* Improved package startup messages to inform users about missing optional dependencies
* Added British English spelling consistency throughout the package
* **Performance tracking**: All core functions now automatically track and report processing time and memory utilisation
* **Enhanced S3 methods**: `summary()` methods now display detailed performance metrics and data structure efficiency statistics
* **Sparsity reporting**: SCHT objects now report both original and hierarchical structure sparsity, demonstrating memory efficiency gains
* **Enhanced visualisation flexibility**: 
  - `plot_tc_landscape()` now supports highlighting both top and bottom genes via the new `n_label` and `label_direction` parameters
  - `plot_isoform_transitions()` now includes a `selected_isoforms` parameter to plot specific isoforms of interest

## New Features

* Added `simple_benchmark.R` script in `inst/scripts/` for performance evaluation
* Added `generate_qc_report.R` script in `inst/scripts/` for quality control reporting
* Added `sparsity_analysis.R` script in `inst/scripts/` for tensor efficiency analysis
* **Quality Control Report Generation**: Added `generate_qc_report()` function for comprehensive QC reporting
  - Customisable dataset names for output files
  - Accurate gene counting from SCHT structure
  - Integrated HTML styling within the function
* **Sparsity Analysis Functions**: Added comprehensive sparsity analysis tools
  - `calculate_scht_sparsity()`: Calculates sparsity statistics for SCHT structures
  - `calculate_ct_scht_sparsity()`: Cell type-specific sparsity analysis
  - `analyse_sparsity_for_table()`: Comprehensive comparison across data representations
  - Demonstrates memory efficiency of SCHT vs naive tensor approaches
* **Co-expression Analysis Suite**: Added comprehensive isoform co-expression analysis functionality
  - `calculate_isoform_coexpression()`: Core correlation calculation for single genes
  - `calculate_gene_coexpression_all_celltypes()`: Multi-cell type correlation analysis
  - `analyse_coexpression_conservation()`: Identifies conserved vs cell-type-specific patterns with bootstrap and FDR
  - `detect_isoform_switching()`: Identifies antagonistic isoform relationships
  - `plot_isoform_coexpression()`: Creates correlation heatmaps using ComplexHeatmap
  - `plot_coexpression_across_celltypes()`: Line plots showing correlation dynamics across cell types
  - `calculate_celltype_coexpression_similarity()`: Cell type similarity based on co-expression
  - `plot_conservation_summary()`: Bar charts of conservation pattern distributions
  - `export_coexpression_results()`: Export functionality for all results
* **Interactive Shiny Application**: Added `launch_coexpression_app()` for interactive co-expression exploration
  - Overall and cell type-specific co-expression heatmaps
  - Isoform switching detection through negative correlations
  - Conservation analysis with mixed pattern detection
  - Statistical analysis: confidence intervals, bootstrap stability (100 iterations), FDR correction
  - Interactive plotly visualisations with detailed hover information
  - Comprehensive export functionality for all analysis results

## Internal Changes

* Created `package_utils.R` with helper functions for package management
* Created `zzz.R` for package startup hooks
* Added fallback colour palettes when RColorBrewer is not available
* Moved `viridis` and `RColorBrewer` from Suggests to Imports for core visualisation functionality
* Renamed internal performance tracking from "peak memory" to "memory utilised" for clarity
* Added performance attributes to SCHT and transcriptomic_complexity objects
* Added validation scripts in `inst/scripts/` for robustness testing
* Standardised code formatting and documentation