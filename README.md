# ScIsoX

## Single-Cell Hierarchical Tensor (SCHT) Creation Pipeline and Transcriptomic Complexity Analysis in R

*A comprehensive toolkit for analysing isoform expression patterns at single-cell resolution*

```
Wu S and Schmitz U (2025). ScIsoX: A Multidimensional Framework for Measuring Transcriptomic Complexity in Single-Cell Long-Read Sequencing Data.
```

---

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Package Architecture](#package-architecture)
4. [Basic Tutorial](#basic-tutorial)
5. [Advanced Usage](#advanced-usage)
6. [Visualisation Gallery](#visualisation-gallery)
7. [Troubleshooting](#troubleshooting)
8. [Citation](#citation)
9. [Contact & Support](#contact--support)

---

## Introduction

ScIsoX provides a robust computational framework for investigating transcriptomic complexity at single-cell resolution. The package enables researchers to analyse alternative splicing patterns across cell types, revealing cell type-specific isoform usage and co-expression dynamics.

### Key Features

- **Single-Cell Hierarchical Tensor**: Generate multi-dimensional representation of isoform expression
- **Complexity Analysis**: Calculate seven core complexity metrics that capture different aspects of transcriptomic diversity
- **Advanced Visualisations**: Create advanced figures for complexity analysis
- **Cell Type Comparison**: Analyse cell type-specific isoform usage patterns

### Theoretical Framework

ScIsoX implements a novel analytical framework based on hierarchical tensor decomposition to quantify transcriptomic complexity across multiple dimensions. The seven core metrics capture:

1. **Intra-cellular isoform diversity**: Co-expression of multiple isoforms within individual cells
2. **Inter-cellular isoform  diversity**: Distribution of different isoforms across the cellular population
3. **Intra-cell-type heterogeneity**: Cell-to-cell variation in isoform usage within a cell type
4. **Inter-cell-type dpecificity**: Cell-type-specific patterns of isoform usage
5. **Intra-cell-type heterogeneity variability**: Variation in cellular heterogeneity across cell types
6. **Inter-cell-type hifference variability**: Identification of cell type pairs with distinctive isoform usage
7. **Cell-type-specific co-expression variability**: Variation in cellular co-expression mechanisms across cell types

---

## Installation

### Prerequisites

ScIsoX requires R version 4.0.0 or higher and depends on several packages for efficient data processing and visualisation.

### From GitHub

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install ScIsoX
devtools::install_github("ThaddeusWu/ScIsoX")
```

### Dependencies

ScIsoX requires several R packages to function properly. Below is a comprehensive list of all dependencies:

#### Required Dependencies (Imports)

These packages will be automatically installed when you install ScIsoX:

```r
# Core dependencies
required_packages <- c(
  "progress",   # For progress bars
  "stats",      # For statistical functions
  "graphics",   # For base plotting
  "rtracklayer", # For handling GTF files
  "tools",      # For file utilities
  "data.table", # For efficient data manipulation
  "Matrix",     # For sparse matrix support
  "methods",    # For S3/S4 method handling
  "dplyr",      # For data manipulation
  "mclust",     # For clustering and mixture models
  "ggplot2",    # For plotting
  "moments",    # For statistical moments calculation
  "utils",      # For utility functions
  "knitr",      # For documentation
  "car",        # For statistical transformations
  "diptest",    # For multimodality tests
  "scales",     # For scale transformations
  "gridExtra",  # For arranging multiple plots
  "magrittr",    # For pipe operations
  "MASS"             # For statistical functions
)

# Check and install missing required packages
missing_required <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_required) > 0) {
  install.packages(missing_required)
}
```

#### Optional Dependencies (Suggests)

These packages provide enhanced functionality:

```r
# Optional packages for advanced visualisations
suggested_packages <- c(
  "ggridges",         # For ridge plots
  "ggrepel",          # For repelling text labels
  "ggExtra",          # For marginal plots
  "viridis",          # For colour palettes
  "ComplexHeatmap",   # For heatmaps
  "cowplot",          # For plot composition
  "grid",             # For advanced graphics
  "tidyr",            # For data reshaping
  "RColorBrewer",     # For colour palettes
  "circlize",   # For colour palettes
  "patchwork"         # For combining plots
)

# Check and install missing suggested packages
missing_suggested <- suggested_packages[!sapply(suggested_packages, requireNamespace, quietly = TRUE)]
if (length(missing_suggested) > 0) {
  install.packages(missing_suggested)
}
```

#### Other Packages

Some advanced features require packages:

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install ggradar for radar charts
if (!requireNamespace("ggradar", quietly = TRUE)) {
  devtools::install_github("ricardo-bion/ggradar")
}

# Install ComplexHeatmap for heat maps
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

To ensure maximum compatibility and functionality, you can install all dependencies with:

```r
# Install all dependencies at once
all_packages <- c(
  required_packages, 
  suggested_packages
)

# Check and install missing packages
missing_all <- all_packages[!sapply(all_packages, requireNamespace, quietly = TRUE)]
if (length(missing_all) > 0) {
  install.packages(missing_all)
}

# Install GitHub packages
if (!requireNamespace("ggradar", quietly = TRUE)) {
  devtools::install_github("ricardo-bion/ggradar")
}
```

Note: Installing all suggested packages is recommended for full functionality of ScIsoX.

---

## Package Architecture

ScIsoX is organised into several functional modules:

1. **Data Import & Preprocessing**
   - `create_transcript_info()`: Generate standardised transcript annotations from GTF files
   - `generate_gene_counts()`: Create gene-level counts from isoform-level expression data
   - `plot_genes_per_cell_distribution()`: Visualise and analyse genes per cell distribution
   - `recommend_qc_parameters()`: Generate data-driven QC parameter recommendations

2. **SCHT Creation Pipeline**
   - `create_scht()`: Process raw count data into Single-Cell Hierarchical Tensor objects
     
3. **Complexity Analysis**
   - `calculate_isoform_complexity_metrics()`: Calculate all seven core complexity metrics
   - `select_genes_of_interest()`: Filter genes based on complexity classifications
   - `find_complexity_pattern()`: Identify genes with specific complexity patterns
   - `compare_gene_metrics()`: Extract and compare metrics across multiple genes

4. **Visualisation & Exploration (13 specialised functions)**
   - `plot_threshold_visualisations()`: Create bar charts for comparing distributions and threshold choices across different complexity metrics
   - `plot_tc_landscape()`: Create complexity landscape plots with marginal distributions
   - `plot_tc_density()`: Generate contour plots for complexity landscapes
   - `plot_diversity_comparison()`: Compare intra-cellular and inter-cellular diversity
   - `plot_complexity_ridges()`: Create ridge plots for complexity metrics
   - `plot_isoform_coexpression()`: Visualise co-expression patterns between isoforms
   - `plot_isoform_profile()`: Create stacked bar charts of isoform usage across cell types
   - `plot_isoform_transitions()`: Visualise isoform usage transitions across cell types
   - `plot_complexity_radar()`: Create radar charts for multi-dimensional complexity comparison
   - `plot_single_gene_radar_cell_type()`: Create radar chart for a single gene across cell types
   - `plot_compare_multiple_genes_radar_cell_type()`: Compare multiple genes across cell types with radar charts
   - `plot_compare_tc_density_difference()`: Compare complexity density differences between groups
   - `plot_compare_tc_complexity_heatmap()`: Create comparative heatmaps of complexity differences

---

## Basic Tutorial

### 1. Data Preparation

ScIsoX works with single-cell gene and transcript expression matrices. Let's start by setting up the required data formats:

```r
library(ScIsoX)

# If working with a GTF file, create transcript information
transcript_info <- create_transcript_info(
  gtf_path = "path/to/gencode.gtf",
  remove_version = TRUE  # Remove version numbers from transcript and gene IDs for better compatibility
)

# If starting with isoform counts, generate gene counts
count_results <- generate_gene_counts(
  isoform_counts = your_isoform_matrix,
  show_progress = TRUE
)

# Extracted output
gene_counts <- count_results$gene_counts
transcript_counts <- your_isoform_matrix  # Original isoform counts
```

### 2. Quality Control Analysis

Before creating the SCHT object, analyse the distribution of genes per cell to determine appropriate QC parameters:

```r
# Visualise genes per cell distribution
qc_suggestions <- plot_genes_per_cell_distribution(
  gene_counts = gene_counts,
  plot_type = "density",  # Options: "density" or "hist"
  percentile_cutoffs = c(0.05, 0.95),  # Customise percentile lines to show
  return_suggestions = TRUE
)

# Get data-driven QC recommendations with multiple strategies
qc_recommendations <- recommend_qc_parameters(gene_counts)

# View recommendations
print(qc_recommendations)
```

The `recommend_qc_parameters()` function provides three different QC strategies:

1. **MAD_strategy**: A conservative approach using median ± 3 MAD (Median Absolute Deviation). This reduces the risk of including poor quality cells while remaining robust to outliers.

2. **Interval_90**: A balanced approach using 5th and 95th percentiles. This strategy balances stringency with dataset preservation.

3. **Interval_80**: An aggressive approach using 10th and 90th percentiles. This provides more stringent filtering for higher quality cell selection.

You can choose the strategy that best suits your dataset or use your own expertise to set custom thresholds.

```r
# Select a QC strategy
selected_strategy <- qc_recommendations$MAD_strategy  # Conservative approach

# Or use custom parameters based on experience
custom_params <- list(
  min_genes_per_cell = 300,  # Based on prior knowledge of your dataset
  max_genes_per_cell = 5000, # Based on prior knowledge of your dataset
  min_cells_expressing = 0.02,
  min_expr = 1e-6
)
```

### 3. Create SCHT Object

Use the QC parameters to create the Single-Cell Hierarchical Tensor object:

```r
# Create cell info data frame with cell type information
cell_info <- data.frame(
  cell_id = colnames(gene_counts),
  cell_type = your_cell_type_vector,
  stringsAsFactors = FALSE
)

# Create SCHT object using selected QC strategy
scht_obj <- create_scht(
  gene_counts = gene_counts,
  transcript_counts = transcript_counts,
  transcript_info = transcript_info,
  cell_info = cell_info,
  n_hvg = 3000,  # Number of highly variable genes to select (can be adjusted based on dataset size)
  qc_params = selected_strategy,  # Or use custom_params if preferred
  require_cell_type = TRUE,  # Set to FALSE if cell type information is not available
  verbose = TRUE,
  sparsity_threshold = 0.4  # Adjust to control when sparse matrices are used (lower = more memory efficient)
)

# Examine the SCHT object
print(scht_obj)
summary(scht_obj)
```

### 4. Calculate Complexity Metrics

Calculate the seven core complexity metrics with flexibility in threshold determination:

```r
# Option 1: Calculate transcriptomic complexity metrics with data-driven thresholds
tc_results <- calculate_isoform_complexity_metrics(
  scht_obj = scht_obj,
  data_driven_thresholds = TRUE,  # Automatically determine optimal thresholds from the data
  visualise = TRUE,  # Generate threshold determination visualisations
  verbose = TRUE
)

# Option 2: Calculate with custom thresholds based on prior knowledge
custom_thresholds <- list(
  intra_cellular_isoform_diversity = 0.5,
  inter_cellular_isoform_diversity = 0.5,
  intra_cell_type_heterogeneity = 0.3,
  inter_cell_type_specificity = 0.7,
  intra_cell_type_heterogeneity_variability = 0.4,
  inter_cell_type_difference_variability = 0.3,
  cell_type_coexpression_variability = 0.4
)

tc_results_custom <- calculate_isoform_complexity_metrics(
  scht_obj = scht_obj,
  default_thresholds = custom_thresholds,
  data_driven_thresholds = FALSE,  # Use provided thresholds instead of data-driven ones
  visualise = TRUE,
  verbose = TRUE
)

# View summary of complexity metrics
summary(tc_results)

# Access complexity metrics data frame
complexity_df <- tc_results$metrics
head(complexity_df)
```

The function supports both data-driven threshold determination, which automatically calculates optimal thresholds based on your dataset's distribution, and user-defined thresholds based on prior knowledge or specific research requirements. Data-driven thresholds are particularly useful for exploratory analysis, while custom thresholds may be preferred for comparative studies or when specific biological knowledge guides threshold selection.

### 5. Basic Visualisation

Create standard visualisations for exploring transcriptomic complexity:

```r
# Create complexity landscape plot with default settings
landscape_plot <- plot_tc_landscape(
  tc_results = tc_results,
  x_metric = "inter_cellular_isoform_diversity",
  y_metric = "inter_cell_type_isoform_specificity",
  highlight_genes = NULL,  # Optionally specify genes to highlight
  label_annotation = "intra_cell_type_heterogeneity",  # Metric used for colour intensity
  label_top = 10,  # Number of genes to label
  use_thresholds = TRUE,  # Use thresholds from tc_results
  point_transparency = 0.85
)

# Display the plot
print(landscape_plot)

# Create density plot for complexity landscape
density_plot <- plot_tc_density(
  tc_results = tc_results,
  x_metric = "inter_cellular_isoform_diversity",
  y_metric = "inter_cell_type_specificity",
  use_thresholds = TRUE,
  show_threshold_lines = TRUE
)

# Display the plot
print(density_plot)

# Visualise diversity comparison
diversity_plot <- plot_diversity_comparison(
  tc_results = tc_results,
  label_top = 10,  # Label top genes below the diagonal
  use_thresholds = TRUE
)

# Display the plot
print(diversity_plot)
```

---

## Advanced Usage

### Exploring Genes with Specific Complexity Patterns

Identify genes with interesting complexity patterns:

```r
# Find genes with high intra-cellular diversity and cell type specificity
complex_genes <- find_complexity_pattern(
  tc_results$metrics,
  pattern = list(
    intra_cellular_isoform_diversity_class = "Strong Isoform Co-expression",
    inter_cell_type_specificity_class = "Cell-Type-Specific Expression"
  ),
  top_n = 20,
  sort_by = "inter_cell_type_specificity"  # Sort by this metric
)

# Select genes with high cellular heterogeneity
heterogeneous_genes <- select_genes_of_interest(
  tc_results$metrics,
  category = "High Cellular Heterogeneity",
  column = "intra_cell_type_heterogeneity_class",
  top_n = 15,
  sort_by = "intra_cell_type_heterogeneity"  # Sort by most heterogeneous
)

# Compare metrics for selected genes
gene_comparison <- compare_gene_metrics(
  tc_metrics = tc_results,
  gene_names = heterogeneous_genes,
  include_mean = TRUE  # Include mean values for easier comparison
)
```

The `find_complexity_pattern()` function enables you to identify genes matching specific combinations of complexity classifications across multiple dimensions. You can specify as many criteria as needed and sort the results by your metric of interest.

The `select_genes_of_interest()` function provides a simpler interface for filtering genes based on a single classification dimension, useful for quickly isolating genes with specific characteristics.

### Multi-dimensional Complexity Visualisation

Create advanced visualisations for multi-dimensional complexity analysis:

```r
# Create dual diversity comparison plot
diversity_plot <- plot_diversity_comparison(
  tc_results = tc_results,
  label_top = 15,  # Number of genes to label
  point_transparency = 0.85,
  use_thresholds = TRUE
)

# Create radar chart for selected genes
radar_plot <- plot_complexity_radar(
  tc_metrics = tc_results,
  genes = complex_genes[1:5],  # Compare up to 5 genes
  scale_type = "per_metric"  # Options: "global", "per_metric", or "none"
)

# Create ridge plots for visualising complexity distributions
ridge_plot <- plot_complexity_ridges(
  tc_results = tc_results,
  type = "global",  # Options: "global" or "cell_type"
  metrics = NULL,  # NULL will use all metrics, or specify a subset
  scale_values = TRUE
)

# Create cell type-specific ridge plots
ct_ridges <- plot_complexity_ridges(
  tc_results = tc_results,
  type = "cell_type",
  n_celltypes = 5,  # Maximum number of cell types to show
  label_y_axis = FALSE
)

# Create radar charts for a single gene across cell types
single_gene_radar <- plot_single_gene_radar_cell_type(
  tc_results = tc_results,
  gene_name = complex_genes[1],
  metrics = NULL,  # NULL will use default metrics, or specify custom metrics
  scale_values = TRUE
)

# Compare multiple genes across cell types with radar charts
multi_gene_radar <- plot_compare_multiple_genes_radar_cell_type(
  tc_results = tc_results,
  gene_names = complex_genes[1:3],
  cell_types = NULL,  # NULL will auto-detect, or specify cell types
  metrics = NULL,  # NULL will use default metrics
  scale_type = "per_cell_type",  # Options: "global" or "per_cell_type"
  ncol = 3  # Number of columns in the grid layout
)
```

These visualisations provide different perspectives on gene complexity. The radar charts are particularly useful for comparing multiple metrics simultaneously, while ridge plots reveal the distribution of each metric across your dataset or cell types.

### Isoform Usage Analysis

Analyse isoform usage patterns for genes of interest:

```r
# Create heatmap for isoform co-expression
coexp_heatmap <- plot_isoform_coexpression(
  scht_obj = scht_obj,
  gene = complex_genes[1],
  display_numbers = TRUE  # Show correlation values in cells
)

# Create stacked bar chart for isoform usage across cell types
isoform_profile <- plot_isoform_profile(
  scht_obj = scht_obj,
  gene = complex_genes[1],
  cell_type_order = NULL,  # Optional ordering of cell types
  min_prop = 0.05,  # Minimum proportion to display (minor isoforms grouped as "Other")
  color_palette = NULL  # Optional custom colour palette
)

# Create line plot for isoform transitions across cell types
# For developmental trajectory analysis
cell_type_order <- c("Progenitor", "Intermediate", "Mature")
transition_plot <- plot_isoform_transitions(
  scht_obj = scht_obj,
  gene = complex_genes[1],
  cell_type_order = cell_type_order,  # Required for meaningful transitions
  min_prop = 0.05,  # Minimum proportion to display
  smooth = TRUE,  # Apply smoothing to lines
  color_palette = NULL  # Optional custom colour palette
)
```

These functions help you understand isoform usage dynamics:

- `plot_isoform_coexpression()` reveals correlations between isoforms, indicating coordinated or mutually exclusive expression
- `plot_isoform_profile()` shows the relative abundance of each isoform across cell types
- `plot_isoform_transitions()` is particularly useful for developmental trajectories or other ordered progressions, showing how isoform usage changes across stages

### Comparing Complexity Across Conditions

For experiments comparing different conditions or treatments:

```r
# Assuming you have SCHT objects for different conditions
scht_obj_A <- scht_obj  # Control
scht_obj_B <- your_treatment_scht_obj

# Calculate complexity metrics for each condition
tc_results_A <- calculate_isoform_complexity_metrics(scht_obj_A)
tc_results_B <- calculate_isoform_complexity_metrics(scht_obj_B)

# Create a list of results
tc_results_list <- list(
  Control = tc_results_A,
  Treatment = tc_results_B
)

# Compare complexity density differences
diff_plot <- plot_compare_tc_density_difference(
  tc_results_list = tc_results_list,
  group_names = c("Control", "Treatment"),
  x_metric = "inter_cellular_diversity",
  y_metric = "inter_cell_type_specificity",
  grid_size = 100,  # Resolution of density estimation
  show_threshold_lines = TRUE
)

# Create heatmaps for metrics across conditions
heatmap_results <- plot_compare_tc_complexity_heatmap(
  tc_results_list = tc_results_list,
  groups = c("Control", "Treatment"),
  metrics = NULL,  # NULL will use all metrics, or specify a subset
  n_top_genes = 50,  # Number of genes to include
  selection_method = "variance",  # Options: "variance", "magnitude", or "custom"
  custom_genes = NULL,  # Used when selection_method = "custom"
  cluster_genes = FALSE,  # Whether to cluster genes
  show_changes = TRUE  # Whether to show changes between groups
)

# Display a heatmap
print(heatmap_results$heatmaps$intra_cellular_diversity)
```

These comparative functions are designed for experimental designs with multiple conditions:

- `plot_compare_tc_density_difference()` highlights regions in the complexity landscape where genes shift between conditions
- `plot_compare_tc_complexity_heatmap()` creates heatmaps for visualising multiple complexity metrics and their changes across different groups or conditions.

Each function offers different ways to select genes of interest, including variance-based selection (genes with highest variation across conditions), magnitude-based selection (genes with largest absolute changes), or custom selection (specific genes of interest).

---

## Troubleshooting

### Common Issues and Solutions

1. **Memory Issues with Large Datasets**
   
   **Problem**: R crashes or shows high memory usage when creating SCHT objects from large datasets.
   
   **Solution**: 
   - Set `sparsity_threshold` to a lower value (e.g., 0.2) to force sparse matrix usage
   - Process data in batches using smaller values
   - For very large datasets, pre-filter genes with low expression before creating the SCHT object

2. **Missing Cell Type Information**
   
   **Problem**: Error about missing cell types when calculating complexity metrics.
   
   **Solution**:
   - Ensure `cell_info` data frame contains a column named "cell_type"
   - Make sure cell IDs in `cell_info` match column names in count matrices
   - Set `require_cell_type = FALSE` if cell type information is not available
   - Check if cell type names contain special characters or spaces that might cause issues

3. **Missing Dependencies for Visualisations**
   
   **Problem**: Errors when trying to create visualisations.
   
   **Solution**:
   - Install all suggested packages for full visualisation support
   - Check for specific missing packages mentioned in error messages
   - For radar charts, ensure `ggradar` is installed from GitHub
   - Use `requireNamespace("package_name", quietly = TRUE)` to check if needed packages are installed

4. **Performance Optimisation**
   
   **Problem**: Processing is slow for large datasets.
   
   **Solution**:
   - Use smaller `n_hvg` values (e.g., 1000-2000)
   - Increase `min_cells_expressing` threshold to focus on well-expressed genes
   - Process cell types separately if working with many cell types
   - Consider pre-filtering the dataset to include only genes of interest

### Getting Help

If you encounter problems not addressed here, please:

1. Submit issues via the GitHub repository
2. Contact the package maintainers with detailed information about your issue, including:
   - R session information (`sessionInfo()`)
   - Error messages (full traceback if possible)
   - A minimal reproducible example if possible

---

## Citation

If you use ScIsoX in your research, please cite:

```
Scisox: a multidimensional framework for measuring transcriptomic complexity in single-cell long-read sequencing data
Wu, S and Schmitz, U
bioRxiv, 2025
DOI: 10.1101/2025.04.28.650897
url: https://www.biorxiv.org/content/10.1101/2025.04.28.650897v1
```

BibTeX entry:

```bibtex

@article{Wu2025.04.28.650897,
	author = {Wu, Siyuan and Schmitz, Ulf},
	doi = {10.1101/2025.04.28.650897},
	eprint = {https://www.biorxiv.org/content/early/2025/05/01/2025.04.28.650897.full.pdf},
	journal = {bioRxiv},
	title = {{ScIsoX}: A Multidimensional Framework for Measuring Transcriptomic Complexity in Single-Cell Long-Read Sequencing Data},
	url = {https://www.biorxiv.org/content/early/2025/05/01/2025.04.28.650897},
	year = {2025}
}
```
---

## Contact & Support

### Maintainer

Siyuan Wu (thaddeus.wu@jcu.edu.au)

### Contributors

- Siyuan Wu (James Cook University)
- Ulf Schmitz (James Cook University)

### Reporting Issues

Please report bugs and feature requests via the [GitHub issue tracker](https://github.com/siyuanwu/ScIsoX/issues).

### License

This package is licensed under the MIT License.

---

*ScIsoX: Unlocking the complexity of the transcriptome, one cell at a time.*
