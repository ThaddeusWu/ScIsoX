% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualisation_analysis.R
\name{plot_compare_tc_difference_heatmap}
\alias{plot_compare_tc_difference_heatmap}
\title{Create heatmaps for transcriptome complexity metrics across groups}
\usage{
plot_compare_tc_difference_heatmap(
  tc_results_list,
  groups = NULL,
  metrics = c("intra_cellular_diversity", "inter_cellular_diversity",
    "intra_cell_type_heterogeneity", "inter_cell_type_specificity",
    "intra_cell_type_heterogeneity_variability",
    "inter_cell_type_difference_variability",
    "cell_type_diversity_mechanism_variability"),
  n_top_genes = 50,
  selection_method = "variance",
  custom_genes = NULL,
  cluster_genes = FALSE,
  show_changes = TRUE
)
}
\arguments{
\item{tc_results_list}{List of transcriptomic complexity results objects}

\item{groups}{Vector of names for each group (if NULL, generic names are used)}

\item{metrics}{Vector of metric names to include in heatmaps}

\item{n_top_genes}{Number of top genes to include in heatmaps}

\item{selection_method}{Method for selecting genes ("variance", "magnitude", or "custom")}

\item{custom_genes}{Vector of custom gene names to use when selection_method is "custom"}

\item{cluster_genes}{Whether to cluster genes in heatmaps}

\item{show_changes}{Whether to create separate heatmaps showing changes between groups}
}
\value{
A list containing:
\item{heatmaps}{List of pheatmap objects for each metric}
\item{change_heatmaps}{List of pheatmap objects showing changes between groups}
\item{metric_matrices}{List of matrices with metric values}
\item{change_matrices}{List of matrices with change values}
\item{top_genes}{Vector of selected gene names}
\item{selection_method}{Method used for gene selection}
}
\description{
This function creates heatmaps for visualising multiple complexity metrics
across different groups or conditions. It can also show changes in metrics
between consecutive groups. The function supports different methods for
selecting genes to display.
}
\examples{
\dontrun{
# Assuming you have transcriptomic complexity results for 3 conditions
tc_results_list <- list(condition_A = tc_results_A, 
                        condition_B = tc_results_B,
                        condition_C = tc_results_C)
                        
# Create heatmaps with default settings
heatmap_results <- plot_compare_tc_difference_heatmap(tc_results_list, 
                  groups = c("Condition A", "Condition B", "Condition C"))
                  
# Display the first heatmap
print(heatmap_results$heatmaps$intra_cellular_diversity)

# Use custom gene selection
custom_heatmaps <- plot_compare_tc_difference_heatmap(tc_results_list,
                  groups = c("Condition A", "Condition B", "Condition C"),
                  selection_method = "custom",
                  custom_genes = c("GENE1", "GENE2", "GENE3"))
}

}
