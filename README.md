# ScIsoX

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/Version-1.1.2-blue.svg)](https://github.com/ThaddeusWu/ScIsoX)
[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue.svg)](https://www.r-project.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16569859.svg)](https://doi.org/10.5281/zenodo.16569859)

## Single-Cell Hierarchical Tensor (SCHT) Creation Pipeline and Transcriptomic Complexity Analysis in R

*A comprehensive toolkit for analysing isoform expression patterns at single-cell resolution*

```
Wu S and Schmitz U (2025). ScIsoX: A Multidimensional Framework for Measuring Isoform-Level Transcriptomic Complexity in Single Cells.
```

---

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Package Architecture](#package-architecture)
4. [5-Minute Quick Start](#5-minute-quick-start)
5. [Basic Tutorial](#basic-tutorial)
6. [Advanced Usage](#advanced-usage)
7. [Visualisation Gallery](#visualisation-gallery)
8. [Troubleshooting](#troubleshooting)
9. [Citation](#citation)
10. [Contact & Support](#contact--support)

---

## Introduction

ScIsoX provides a robust computational framework for investigating transcriptomic complexity at single-cell resolution. The package enables researchers to analyse alternative splicing patterns across cell types, revealing cell type-specific isoform usage and co-expression dynamics.

### Key Features

- **Single-Cell Hierarchical Tensor**: Generate multi-dimensional representation of isoform expression
- **Complexity Analysis**: Calculate seven core complexity metrics that capture different aspects of transcriptomic diversity
- **Advanced Visualisations**: Create advanced figures for complexity analysis with 13+ specialised plotting functions
- **Cell Type Comparison**: Analyse cell type-specific isoform usage patterns
- **Co-expression Analysis**: Interactive Shiny application for exploring isoform relationships
- **Performance Tracking**: Built-in memory and runtime monitoring with sparsity analysis

### Theoretical Framework

ScIsoX implements a novel analytical framework based on hierarchical tensor decomposition to quantify transcriptomic complexity across multiple dimensions. The seven core metrics capture:

1. **Intra-cellular isoform diversity**: Co-expression of multiple isoforms within individual cells
2. **Inter-cellular isoform  diversity**: Distribution of different isoforms across the cellular population
3. **Intra-cell-type heterogeneity**: Cell-to-cell variation in isoform usage within a cell type
4. **Inter-cell-type specificity**: Cell-type-specific patterns of isoform usage
5. **Intra-cell-type heterogeneity variability**: Variation in cellular heterogeneity across cell types
6. **Inter-cell-type difference variability**: Identification of cell type pairs with distinctive isoform usage
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

### First Time Setup

When you first load ScIsoX, it will check for optional packages and provide helpful information:

```r
library(ScIsoX)

# If some optional packages are missing, you'll see:
# ===============================================
# ScIsoX: Single-Cell Isoform Complexity Analysis
# ===============================================
# 
# Some optional packages are not installed:
#   - For visualisation: ggridges, ggrepel, ggExtra
#   - For advanced_plots: ComplexHeatmap, ggradar, cowplot, patchwork
#   - For data_manipulation: tidyr
# 
# Some visualisation features may not be available.
# To install all optional packages, run:
#   install_scisox_suggests()

# Install all optional packages with one command:
install_scisox_suggests()

# Or install only specific types:
install_scisox_suggests(include_bioc = FALSE)  # Skip Bioconductor packages
install_scisox_suggests(include_github = FALSE)  # Skip GitHub packages
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
  "grDevices",  # For graphics devices and colours
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
  "magrittr",   # For pipe operations
  "MASS",       # For statistical functions
  "viridis",    # For colour palettes (frequently used)
  "RColorBrewer", # For colour palettes (frequently used)
  "shiny",      # For interactive applications
  "shinydashboard", # For Shiny dashboard layouts
  "plotly",     # For interactive plots
  "DT"          # For interactive data tables
)

# These are automatically installed with ScIsoX
```

#### Optional Dependencies (Suggests)

These packages provide enhanced functionality. **ScIsoX is designed to work without these packages**, but some visualisations and features will have reduced functionality:

```r
# Optional packages for advanced visualisations
suggested_packages <- c(
  "ggridges",         # For ridge plots (required for plot_complexity_ridges)
  "ggrepel",          # For repelling text labels (falls back to regular text)
  "ggExtra",          # For marginal plots (skipped if not available)
  "ComplexHeatmap",   # For heatmaps (required for plot_isoform_coexpression)
  "cowplot",          # For plot composition (required for multi-panel layouts)
  "grid",             # For advanced graphics
  "tidyr",            # For data reshaping (falls back to base R reshape)
  "circlize",         # For colour palettes in heatmaps
  "patchwork"         # For combining plots
)
```

**Quick Installation:** ScIsoX provides a convenient function to install all optional packages:

```r
# After installing ScIsoX, run:
library(ScIsoX)
install_scisox_suggests()
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

#### Package Behaviour with Missing Dependencies

ScIsoX is designed to be robust when optional packages are not installed:

- **Core functionality always works**: All complexity calculations and basic analyses
- **Graceful degradation**: When optional packages are missing:
  - `ggrepel`: Labels use standard text positioning
  - `ggExtra`: Marginal distributions are skipped
  - `tidyr`: Base R reshape functions are used
  - Colour palettes: Built-in alternatives are provided
- **Clear error messages**: Functions requiring specific packages will provide installation instructions
- **Startup messages**: When loading ScIsoX, you'll see which optional packages are missing

```r
# Example: Using ScIsoX without all optional packages
library(ScIsoX)
# You'll see a message about missing packages and how to install them

# The core analysis still works perfectly:
# First create SCHT object
data(gene_counts_blood)
data(transcript_counts_blood)
data(transcript_info)
data(sample2stage)

scht_obj <- create_scht(
  gene_counts = gene_counts_blood,
  transcript_counts = transcript_counts_blood,
  transcript_info = transcript_info,
  cell_info = sample2stage,
  n_hvg = 3000,
  qc_params = list(
    min_genes_per_cell = 4000,       
    max_genes_per_cell = 10000,      
    min_cells_expressing = 0.02,   
    min_expr = 1e-6
  ),
  verbose = FALSE
)

tc_results <- calculate_isoform_complexity_metrics(scht_obj)

# Some visualisations may have reduced features:
plot_tc_landscape(tc_results)  # Works, but labels might overlap without ggrepel
```

Note: For the best experience, we recommend installing all suggested packages using `install_scisox_suggests()`.

---

## Package Architecture

ScIsoX is organised into several functional modules:

1. **Data Import & Preprocessing**
   - `create_transcript_info()`: Generate standardised transcript annotations from GTF files
   - `generate_gene_counts()`: Create gene-level counts from isoform-level expression data
   - `plot_genes_per_cell_distribution()`: Visualise and analyse genes per cell distribution
   - `recommend_qc_parameters()`: Generate data-driven QC parameter recommendations

2. **SCHT Creation Pipeline**
   - `create_scht()`: Process raw count data into Single-Cell Hierarchical Tensor objects with automatic sparsity detection
     
3. **Complexity Analysis**
   - `calculate_isoform_complexity_metrics()`: Calculate all seven core complexity metrics
   - `select_genes_of_interest()`: Filter genes based on complexity classifications
   - `find_complexity_pattern()`: Identify genes with specific complexity patterns
   - `compare_gene_metrics()`: Extract and compare metrics across multiple genes

4. **Visualisation & Exploration (13+ specialised functions)**
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

5. **Co-expression Analysis Suite**
   - `calculate_isoform_coexpression()`: Core correlation calculation for single genes
   - `calculate_gene_coexpression_all_celltypes()`: Multi-cell type correlation analysis
   - `analyse_coexpression_conservation()`: Identifies conserved vs cell-type-specific patterns
   - `detect_isoform_switching()`: Identifies antagonistic isoform relationships
   - `launch_coexpression_app()`: Interactive Shiny application for exploration

6. **Performance & Quality Control**
   - `calculate_scht_sparsity()`: Analyse memory efficiency of SCHT structures
   - `generate_qc_report()`: Generate comprehensive QC reports

---

## 5-Minute Quick Start

Want to get started quickly? Here's a minimal example using the included example data:

```r
# Load ScIsoX
library(ScIsoX)

# Load example data
data(gene_counts_blood)
data(transcript_counts_blood)
data(transcript_info)
data(sample2stage)

# Create SCHT object with recommended parameters
scht_obj <- create_scht(
  gene_counts = gene_counts_blood,
  transcript_counts = transcript_counts_blood,
  transcript_info = transcript_info,
  cell_info = sample2stage,
  qc_params = list(
    min_genes_per_cell = 4000,
    max_genes_per_cell = 10000,
    min_cells_expressing = 0.02,
    min_expr = 1e-6
  ),
  n_hvg = 3000,
  verbose = FALSE
)

# Calculate complexity metrics
tc_results <- calculate_isoform_complexity_metrics(scht_obj, verbose = FALSE)

# Visualise complexity landscape
plot_tc_landscape(tc_results, n_label = 5)

# Find genes with high complexity
complex_genes <- tc_results$metrics[
  tc_results$metrics$intra_cellular_isoform_diversity_class == "high" &
  tc_results$metrics$inter_cellular_isoform_diversity_class == "high",
]
print(paste("Found", nrow(complex_genes), "highly complex genes"))

# Launch interactive co-expression app (optional)
# launch_coexpression_app(scht_obj)
```

That's it! You've just:
- Created a Single-Cell Hierarchical Tensor (SCHT)
- Calculated transcriptomic complexity metrics
- Visualised the complexity landscape
- Identified highly complex genes

For more details, continue to the tutorials below.

---

## Basic Tutorial

### Input data recommendations
The accuracy of _ScIsoX_'s complexity metrics is contingent upon the quality of the input isoform-by-cell count matrix. Users should be aware that artefacts from library preparation, sequencing errors, or read misalignment can lead to the spurious identification of transcripts, potentially inflating complexity metrics. Rigorous upstream quality control and filtering are therefore essential. For instance, for long-read data, several tools can be employed to collapse redundant transcripts and filter artefacts such as TALON or IsoQuant. Similarly, for isoform-level quantification from short-read data, established tools include Salmon, Kallisto-Bustools, and RSEM. The performance of these and other tools can vary depending on the dataset, and we encourage readers to consult comprehensive benchmarking studies, such as Dong et al. for long-read RNA-seq and Westoby et al. for short-read RNA-seq, to select the most appropriate method for their specific needs. For experimental designs with multiple batches, we recommend applying batch correction using established methods (e.g., Harmony, Seurat's integration workflow) before analysis with _ScIsoX_. For all data types, starting with a high-quality reference annotation and applying appropriate expression thresholds to filter out extremely low-abundance isoforms is critical. While _ScIsoX_'s internal filtering of low-expression genes and its focus on highly variable genes help to mitigate the impact of residual noise, the validity of novel isoform analysis ultimately depends on robust upstream quantification.

### 1. Data Preparation

ScIsoX works with single-cell gene and transcript expression matrices. Let's start by setting up the required data formats:

```r
library(ScIsoX)

# Example data is included with the package
data(gene_counts_blood)
data(transcript_counts_blood)
data(transcript_info)
data(sample2stage)

# View the structure of example data
str(gene_counts_blood)    # Gene-level count matrix
str(transcript_counts_blood)  # Transcript-level count matrix
str(transcript_info)      # Transcript annotation information
str(sample2stage)         # Cell metadata
```

#### Working with Your Own Data

```r
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

### Data Attribution

The example dataset included in ScIsoX is derived from:

> Wang Q, et al. (2022). Single-cell transcriptomic atlas of the human endometrium during the menstrual cycle. *Science Advances* 8(1):eabg5369. DOI: 10.1126/sciadv.abg5369

This data is provided under CC BY 4.0 license and has been processed for use as example data in this package.

### 2. Quality Control Analysis

Before creating the SCHT object, analyse the distribution of genes per cell to determine appropriate QC parameters:

```r
# Visualise genes per cell distribution using example data
qc_suggestions <- plot_genes_per_cell_distribution(
  gene_counts = gene_counts_blood,
  plot_type = "density",  # Options: "density" or "hist"
  percentile_cutoffs = c(0.05, 0.95),  # Customise percentile lines to show
  return_suggestions = TRUE
)

# Get data-driven QC recommendations with multiple strategies
qc_recommendations <- recommend_qc_parameters(gene_counts_blood)

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
  min_expr = 1e-4
)
```

### 3. Create SCHT Object

Use the QC parameters to create the Single-Cell Hierarchical Tensor object:

```r
# Using the example data
data(gene_counts_blood)
data(transcript_counts_blood)
data(transcript_info)
data(sample2stage)  # This contains cell type information

# Create SCHT object using recommended parameters
scht_obj <- create_scht(
  gene_counts = gene_counts_blood,
  transcript_counts = transcript_counts_blood,
  transcript_info = transcript_info,
  cell_info = sample2stage,
  n_hvg = 3000,
  qc_params = list(
    min_genes_per_cell = 4000,       
    max_genes_per_cell = 10000,      
    min_cells_expressing = 0.02,   
    min_expr = 1e-6
  ),
  require_cell_type = TRUE,
  verbose = TRUE,
  sparsity_threshold = 0.4
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
  y_metric = "inter_cell_type_specificity",
  highlight_genes = NULL,  # Optionally specify genes to highlight
  label_annotation = "intra_cell_type_heterogeneity",  # Metric used for colour intensity
  n_label = 10,  # Number of genes to label
  label_direction = "top",  # "top" or "bottom"
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
  use_thresholds = TRUE
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
    inter_cell_type_specificity_class = "Cell-Type-Specific Isoform Expression"
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
# Find a gene with multiple isoforms for demonstration
multi_iso_genes <- names(which(table(transcript_info$gene_name) > 2))
if(length(multi_iso_genes) > 0) {
  example_gene <- multi_iso_genes[1]
  
  # Create heatmap for isoform co-expression
  coexp_heatmap <- plot_isoform_coexpression(
    scht_obj = scht_obj,
    gene = example_gene,
    display_numbers = TRUE  # Show correlation values in cells
  )
  
  # Create stacked bar chart for isoform usage across cell types
  isoform_profile <- plot_isoform_profile(
    scht_obj = scht_obj,
    gene = example_gene,
    cell_type_order = NULL,  # Optional ordering of cell types
    min_prop = 0.05,  # Minimum proportion to display
    color_palette = NULL  # Optional custom colour palette
  )
  
  # Create line plot for isoform transitions across cell types
  # Get actual cell types from the data
  cell_types <- sort(unique(scht_obj$cell_info$cell_type))
  if(length(cell_types) >= 2) {
    transition_plot <- plot_isoform_transitions(
      scht_obj = scht_obj,
      gene = example_gene,
      cell_type_order = cell_types,  # Use actual cell types
      min_prop = 0.05,  # Minimum proportion to display
      smooth = TRUE,  # Apply smoothing to lines
      color_palette = NULL  # Optional custom colour palette
    )
  }
}
```

These functions help you understand isoform usage dynamics:

- `plot_isoform_coexpression()` reveals correlations between isoforms, indicating coordinated or mutually exclusive expression
- `plot_isoform_profile()` shows the relative abundance of each isoform across cell types
- `plot_isoform_transitions()` is particularly useful for developmental trajectories or other ordered progressions, showing how isoform usage changes across stages

### Comparing Complexity Across Conditions

For experiments comparing different conditions:

```r
# For demonstration, we'll simulate different conditions by subsetting our data
# In practice, you would have separate datasets for each condition

# Create subsets to simulate conditions (e.g., different developmental stages)
stages <- unique(sample2stage$cell_type)
if(length(stages) >= 2) {
  # Create SCHT objects for different stages/conditions
  stage1_cells <- rownames(sample2stage)[sample2stage$cell_type == stages[1]]
  stage2_cells <- rownames(sample2stage)[sample2stage$cell_type == stages[2]]
  
  # Create subset SCHT objects (in practice, these would be your different conditions)
  # Note: This is for demonstration only
  scht_stage1 <- create_scht(
    gene_counts = gene_counts_blood[, stage1_cells],
    transcript_counts = transcript_counts_blood[, stage1_cells],
    transcript_info = transcript_info,
    cell_info = sample2stage[stage1_cells, , drop = FALSE],
    n_hvg = 1000,  # Fewer HVGs for smaller dataset
    qc_params = list(
      min_genes_per_cell = 2000,  # Adjusted for subset
      max_genes_per_cell = 8000,
      min_cells_expressing = 0.05,
      min_expr = 1e-6
    ),
    verbose = FALSE
  )
  
  # Calculate complexity metrics
  tc_results_stage1 <- calculate_isoform_complexity_metrics(scht_stage1, verbose = FALSE)
  tc_results_stage2 <- tc_results_stage1  # For demo purposes
  
  # Create a list of results
  tc_results_list <- list(
    Stage1 = tc_results_stage1,
    Stage2 = tc_results_stage2  # In practice, calculate separately
  )
  
  # Compare complexity (demonstration with same data)
  # In practice, these would show real differences between conditions
  diff_plot <- plot_compare_tc_density_difference(
    tc_results_list = tc_results_list,
    group_names = c("Stage 1", "Stage 2"),
    x_metric = "inter_cellular_isoform_diversity",
    y_metric = "inter_cell_type_specificity",
    grid_size = 50  # Lower resolution for demo
  )
}
```

These comparative functions are designed for experimental designs with multiple conditions:

- `plot_compare_tc_density_difference()` highlights regions in the complexity landscape where genes shift between conditions
- `plot_compare_tc_complexity_heatmap()` creates heatmaps for visualising multiple complexity metrics and their changes across different groups or conditions.

Each function offers different ways to select genes of interest, including variance-based selection (genes with highest variation across conditions), magnitude-based selection (genes with largest absolute changes), or custom selection (specific genes of interest).

### Co-expression Analysis

ScIsoX provides comprehensive tools for analysing isoform co-expression patterns:

```r
# Find genes with multiple isoforms for co-expression analysis
multi_iso_genes <- names(which(table(transcript_info$gene_name) > 2))

if(length(multi_iso_genes) > 0) {
  # Calculate co-expression for a single gene
  coexp_result <- calculate_isoform_coexpression(
    scht_obj = scht_obj,
    gene = multi_iso_genes[1],
    method = "pearson",  # Options: "pearson", "spearman"
    min_cells = 10
  )
  
  # Analyse co-expression across all cell types
  multi_coexp <- calculate_gene_coexpression_all_celltypes(
    scht_obj = scht_obj,
    gene = multi_iso_genes[1],
    method = "pearson"
  )
}

# Identify conservation patterns
if(exists("multi_coexp") && length(multi_iso_genes) > 0) {
  conservation_results <- analyse_coexpression_conservation(
    integrated_scht = scht_obj,
    gene = multi_iso_genes[1],
    method = "pearson",
    min_cells = 10,
    consistency_threshold = 0.7,
    correlation_threshold = 0.3
  )
  
  # View conservation summary
  print(conservation_results$summary)
  
  # Plot conservation summary
  plot_conservation_summary(
    conservation_results,
    output_file = "conservation_summary.pdf",
    width = 8,
    height = 6
  )
}

# Detect isoform switching (antagonistic relationships)
if(exists("coexp_result")) {
  switching_results <- detect_isoform_switching(
    cor_result = coexp_result,
    threshold = -0.3,
    strong_threshold = -0.5
  )
  
  # View switching pairs
  if(switching_results$n_switching_pairs > 0) {
    print(switching_results$switching_pairs)
  }
}

# Calculate co-expression statistics
if(exists("coexp_result")) {
  coexp_stats <- calculate_coexpression_stats(
    cor_result = coexp_result,
    include_pairwise = TRUE
  )
  print(paste("Mean correlation:", round(coexp_stats$mean_correlation, 3)))
  print(paste("Positive correlations:", coexp_stats$pct_positive, "%"))
}

# Export co-expression results
if(exists("multi_coexp")) {
  export_coexpression_results(
    coexpr_result = multi_coexp,
    output_prefix = "coexpression_analysis",
    formats = c("csv", "rds")
  )
}

# Launch interactive Shiny application
launch_coexpression_app(scht_obj)
```

The co-expression analysis suite provides:
- **Correlation analysis**: Calculate pairwise correlations between isoforms
- **Conservation patterns**: Identify conserved, variable, and mixed patterns across cell types
- **Isoform switching**: Detect antagonistic isoform relationships
- **Interactive exploration**: Shiny app with heatmaps, statistics, and export functionality

### Performance Analysis & Quality Control

Monitor performance and generate comprehensive QC reports:

```r
# Analyse SCHT sparsity and memory efficiency
sparsity_stats <- calculate_scht_sparsity(scht_obj)
print(sparsity_stats)

# Compare with cell type-specific analysis
ct_sparsity <- calculate_ct_scht_sparsity(scht_obj)

# Generate comprehensive QC report
generate_qc_report(
  scht_obj = scht_obj,
  output_dir = "qc_reports",
  dataset_name = "my_dataset",
  include_plots = TRUE,
  include_summary_stats = TRUE
)

# Analyse sparsity across different data representations
sparsity_comparison <- analyse_sparsity_for_table(
  gene_counts = gene_counts_blood,
  transcript_counts = transcript_counts_blood,
  transcript_info = transcript_info,
  scht_obj = scht_obj,
  dataset_name = "Blood Cell Example Data"
)
```

These functions help you:
- **Understand memory usage**: See how SCHT's hierarchical structure saves memory
- **Monitor performance**: Track processing time and memory utilisation
- **Quality control**: Generate detailed reports on data quality and filtering results
- **Compare efficiency**: Demonstrate advantages over naive tensor approaches

---

## Visualisation Gallery

ScIsoX provides over 13 specialised visualisation functions for comprehensive exploration of transcriptomic complexity. Here are detailed examples:

### Threshold Visualisations

```r
# View threshold determination plots from complexity calculation
if(tc_results$visualise) {
  plot_threshold_visualisations(
    threshold_plots = tc_results$threshold_plots,
    ncol = 3,
    title = "Threshold Determination for Complexity Metrics"
  )
}
```

### Complexity Landscape Plots

```r
# Create landscape plot with marginal distributions
landscape <- plot_tc_landscape(
  tc_results = tc_results,
  x_metric = "inter_cellular_isoform_diversity",
  y_metric = "inter_cell_type_specificity",
  highlight_genes = names(sort(tc_results$metrics$inter_cellular_isoform_diversity, 
                              decreasing = TRUE))[1:5],
  label_annotation = "intra_cell_type_heterogeneity",
  n_label = 10,
  label_direction = "top",  # "top" or "bottom"
  use_thresholds = TRUE,
  x_threshold = NULL,  # Use threshold from tc_results
  y_threshold = NULL,  # Use threshold from tc_results
  point_transparency = 0.85
)

# Create density contour plot
density <- plot_tc_density(
  tc_results = tc_results,
  x_metric = "inter_cellular_isoform_diversity",
  y_metric = "inter_cell_type_specificity",
  use_thresholds = TRUE,
  x_threshold = 0.6,
  y_threshold = 0.6
)
```

### Diversity Analysis Plots

```r
# Compare intra- vs inter-cellular diversity
diversity_comp <- plot_diversity_comparison(
  tc_results = tc_results,
  label_top = 20,
  point_transparency = 0.85,
  use_thresholds = TRUE,
  x_threshold = 0.6,
  y_threshold = 0.6
)

# Create ridge plots for metric distributions
ridge_global <- plot_complexity_ridges(
  tc_results = tc_results,
  type = "global",
  metrics = NULL,  # Use all metrics
  scale_values = TRUE
)

# Cell type-specific ridge plots
ridge_ct <- plot_complexity_ridges(
  tc_results = tc_results,
  type = "cell_type",
  cell_types = NULL,  # Auto-detect
  metrics = c("intra_cellular_isoform_diversity", 
              "inter_cellular_isoform_diversity"),
  n_celltypes = 5,
  label_y_axis = TRUE
)
```

### Isoform Co-expression Visualisations

```r
# Find genes with multiple isoforms
multi_iso_genes <- names(which(table(transcript_info$gene_name) > 2))

if(length(multi_iso_genes) > 0) {
  # Co-expression heatmap
  coexp_heatmap <- plot_isoform_coexpression(
    scht_obj = scht_obj,
    gene = multi_iso_genes[1],
    method = "pearson",
    display_numbers = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    title = NULL  # Auto-generated
  )
  
  # Isoform usage profile across cell types
  usage_profile <- plot_isoform_profile(
    scht_obj = scht_obj,
    gene = multi_iso_genes[1],
    cell_type_order = NULL,  # Auto-order
    min_prop = 0.05,
    colour_palette = NULL  # Default palette
  )
  
  # Isoform transitions (if ordered cell types)
  cell_types <- sort(unique(scht_obj$cell_info$cell_type))
  if(length(cell_types) >= 2) {
    transitions <- plot_isoform_transitions(
      scht_obj = scht_obj,
      gene = multi_iso_genes[1],
      cell_type_order = cell_types,
      selected_isoforms = NULL,  # All isoforms
      min_prop = 0.05,
      smooth = TRUE,
      colour_palette = NULL
    )
  }
}
```

### Radar Charts for Multi-dimensional Comparison

```r
# Compare genes from the results
# Use genes from the middle of the results to ensure they exist
available_genes <- rownames(tc_results$metrics)
complex_genes <- available_genes[c(100, 200, 300, 400, 500)]

# Multi-gene radar chart
if(requireNamespace("ggradar", quietly = TRUE)) {
  radar_multi <- plot_complexity_radar(
    tc_metrics = tc_results,
    genes = complex_genes,
    scale_type = "global"  # "global", "per_metric", or "none"
  )
  
  # Single gene across cell types
  radar_single <- plot_single_gene_radar_cell_type(
    tc_results = tc_results,
    gene_name = complex_genes[1],
    metrics = NULL,  # Default metrics
    scale_values = TRUE
  )
  
  # Multiple genes across cell types
  radar_compare <- plot_compare_multiple_genes_radar_cell_type(
    tc_results = tc_results,
    gene_names = complex_genes[1:3],
    cell_types = NULL,  # All cell types
    metrics = NULL,
    scale_type = "per_cell_type",
    ncol = 3
  )
}
```

### Comparative Analysis Visualisations

```r
# For comparing conditions (requires multiple tc_results)
# Example with simulated conditions
if(length(unique(sample2stage$cell_type)) >= 2) {
  # Create comparison list (in practice, use different conditions)
  tc_results_list <- list(
    Baseline = tc_results,
    Treatment = tc_results  # Would be different condition
  )
  
  # Density difference plots
  density_diff <- plot_compare_tc_density_difference(
    tc_results_list = tc_results_list,
    group_names = c("Baseline", "Treatment"),
    pair_indices = list(c(1, 2)),
    x_metric = "inter_cellular_isoform_diversity",
    y_metric = "inter_cell_type_specificity",
    grid_size = 50
  )
  
  # Complexity heatmaps
  if(requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    heatmap_comp <- plot_compare_tc_complexity_heatmap(
      tc_results_list = tc_results_list,
      groups = c("Baseline", "Treatment"),
      metrics = c("intra_cellular_isoform_diversity", 
                  "inter_cellular_isoform_diversity"),
      n_top_genes = 30,
      selection_method = "variance",
      cluster_genes = TRUE,
      show_changes = TRUE
    )
    
    # Display specific metric heatmap
    print(heatmap_comp$heatmaps$intra_cellular_isoform_diversity)
  }
}
```

### Co-expression Pattern Visualisations

```r
# Plot co-expression across cell types
if(exists("multi_coexp") && length(multi_coexp$cell_types) > 0) {
  coexp_pattern <- plot_coexpression_across_celltypes(
    coexpr_all_result = multi_coexp,
    pair_selection = "all",  # "all", "switching", "conserved"
    threshold = 0.3
  )
  print(coexp_pattern)
}

# Cell type similarity based on co-expression
if(exists("multi_coexp") && length(multi_coexp$cell_types) >= 2) {
  # Extract correlation matrices
  cor_list <- lapply(multi_coexp$cell_types, function(x) x$cor_matrix)
  
  # Calculate similarity
  ct_similarity <- calculate_celltype_coexpression_similarity(
    correlation_list = cor_list,
    method = "correlation"  # "correlation", "euclidean", "manhattan"
  )
  
  # Visualise with hierarchical clustering
  hc <- hclust(ct_similarity)
  plot(hc, main = "Cell Type Similarity Based on Co-expression Patterns")
}
```

### Export Visualisations

```r
# Save plots to PDF
pdf("scisox_visualisations.pdf", width = 10, height = 8)

# Add your plots here
print(landscape)
print(diversity_comp)
print(ridge_global)

dev.off()

# Save individual plots
ggsave("complexity_landscape.pdf", landscape, width = 10, height = 8)
ggsave("diversity_comparison.pdf", diversity_comp, width = 8, height = 8)
```

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

**Manuscript:**
```
ScIsoX: A Multidimensional Framework for Measuring Isoform-Level Transcriptomic Complexity in Single Cells
Wu, S and Schmitz, U
bioRxiv, 2025
DOI: 10.1101/2025.04.28.650897
url: https://www.biorxiv.org/content/10.1101/2025.04.28.650897v2
```

**Software:**
```
Wu, S and Schmitz, U (2025). ScIsoX: Single-Cell Hierarchical Tensor (SCHT) Creation Pipeline 
and Transcriptomic Complexity Analysis in R. R package version 1.1.1.
DOI: 10.5281/zenodo.16569860
```

BibTeX entry:

```bibtex

@article{Wu2025.04.28.650897,
	author = {Wu, Siyuan and Schmitz, Ulf},
	doi = {10.1101/2025.04.28.650897},
	eprint = {https://www.biorxiv.org/content/early/2025/05/01/2025.04.28.650897.full.pdf},
	journal = {bioRxiv},
	title = {{ScIsoX}: A Multidimensional Framework for Measuring Isoform-Level Transcriptomic Complexity in Single Cells},
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

Please report bugs and feature requests via the [GitHub issue tracker](https://github.com/ThaddeusWu/ScIsoX/issues).


### License

This package is licensed under the MIT License.

### Acknowledgments

We gratefully acknowledge Wang et al. (2022) for making their single-cell multiomics data publicly available. The example datasets included in ScIsoX are derived from their publication:

> Wang W, et al. (2022). Single-cell multiomics defines tolerogenic extrathymic Aire-expressing populations. *Science Advances* 8(1):eabg5369. DOI: 10.1126/sciadv.abg5369

These datasets are used under the Creative Commons Attribution 4.0 International License (CC BY 4.0).

---

*ScIsoX: Unlocking the complexity of the transcriptome, one cell at a time.*
