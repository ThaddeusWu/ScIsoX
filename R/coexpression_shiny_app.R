##################################################################
#  Interactive Co-expression Analysis Shiny Application          #
#  for Single-cell Hierarchical Tensor (SCHT) structure          #
#                                                                #
#  Author: [Siyuan Wu & Ulf Schmitz]                             #
#  Institution: [James Cook University]                          #
#  Date: Jul 29, 2025                                            #
#  Package: ScIsoX V1.1.0                                        #
##################################################################

########################
# Required Libraries   #
########################
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny selectInput numericInput actionButton downloadButton 
#' @importFrom shiny renderPlot renderTable renderUI observeEvent reactive req
#' @importFrom shiny updateSelectInput updateNumericInput downloadHandler
#' @importFrom shiny tabsetPanel tabPanel fluidRow column wellPanel br hr h4 p
#' @importFrom shiny selectizeInput updateSelectizeInput renderText textOutput
#' @importFrom shiny renderDataTable dataTableOutput onStop
#' @importFrom stats cor p.adjust sd quantile hclust pt complete.cases
#' @importFrom utils write.csv zip
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom grid unit gpar
#' @importFrom graphics par plot.new
#' @importFrom grDevices png dev.off
#' @importFrom plotly plot_ly layout add_trace colorbar renderPlotly plotlyOutput
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridis viridis
#' @importFrom DT datatable renderDT DTOutput
#' @importFrom plyr rbind.fill

# Define global variables to avoid R CMD check notes
utils::globalVariables(c(".scht_obj_temp"))

## Launch interactive co-expression analysis app
#'
#' @title Interactive Co-expression Analysis Application
#' @description Launch a Shiny application for interactive exploration of isoform co-expression patterns
#' @param scht_obj An IntegratedSCHT object created by create_scht()
#' @param port Port number for the app (default: NULL, uses random port)
#' @param launch.browser Whether to launch browser automatically (default: TRUE)
#' @return Launches Shiny application
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
#' # Launch the Shiny app
#' \dontrun{
#' launch_coexpression_app(scht_obj)
#' }
#' 
#' # Launch on specific port without opening browser
#' \dontrun{
#' launch_coexpression_app(scht_obj, port = 8080, launch.browser = FALSE)
#' }
#' @export
launch_coexpression_app <- function(scht_obj, port = NULL, launch.browser = TRUE) {
  
  # Check required libraries
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Please install it.")
  }
  if (!requireNamespace("shinydashboard", quietly = TRUE)) {
    stop("Package 'shinydashboard' is required. Please install it.")
  }
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("Package 'DT' is required. Please install it.")
  }
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Please install it.")
  }
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required. Please install it.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required. Please install it.")
  }
  
  # Validate input
  if (!inherits(scht_obj, "IntegratedSCHT")) {
    stop("scht_obj must be an IntegratedSCHT object")
  }
  
  # Get available genes
  available_genes <- names(scht_obj$original_results)
  
  # Store objects in global environment
  .GlobalEnv$.scht_obj_temp <- scht_obj
  
  # Define UI
  ui <- shinydashboard::dashboardPage(
    
    # Header
    shinydashboard::dashboardHeader(
      title = "ScIsoX Co-expression Analysis",
      titleWidth = 350
    ),
    
    # Sidebar
    shinydashboard::dashboardSidebar(
      width = 300,
      
      # Gene selection with server-side selectize
      shiny::selectizeInput(
        "selected_gene",
        "Select Gene:",
        choices = NULL,  # Will be populated server-side
        options = list(
          placeholder = 'Type to search genes...',
          maxOptions = 50,
          searchField = "label"
        )
      ),
      
      # Gene info
      shiny::tags$div(
        style = "padding: 10px; font-size: 12px; color: #666; background-color: #f4f4f4; border-radius: 5px; margin: 10px;",
        shiny::textOutput("gene_info")
      ),
      
      # Analysis parameters
      shiny::tags$div(
        style = "padding: 0 10px;",
        shiny::h4("Analysis Parameters"),
        
        shiny::sliderInput(
          "min_cells",
          "Minimum cells per cell type:",
          min = 1,
          max = 50,
          value = 5,
          step = 1
        ),
        
        shiny::sliderInput(
          "min_expression",
          "Minimum expression threshold:",
          min = 0,
          max = 1,
          value = 0.01,
          step = 0.01
        ),
        
        shiny::selectInput(
          "cor_method",
          "Correlation method:",
          choices = c("pearson", "spearman", "kendall"),
          selected = "pearson"
        )
      ),
      
      # Conservation parameters
      shiny::tags$div(
        style = "padding: 0 10px;",
        shiny::h4("Conservation Analysis"),
        
        shiny::sliderInput(
          "cor_threshold",
          "Correlation threshold:",
          min = 0,
          max = 1,
          value = 0.3,
          step = 0.05
        ),
        
        shiny::sliderInput(
          "consistency_threshold",
          "Consistency (%):",
          min = 0,
          max = 100,
          value = 60,
          step = 5
        ),
        
        shiny::helpText(
          style = "font-size: 11px; color: #666; margin-top: -10px;",
          "Consistency: % of cell types where the isoform pair shows correlation above threshold"
        ),
        
        shiny::actionButton(
          "run_conservation",
          "Run Analysis",
          class = "btn-success",
          style = "width: 80%; margin-left: 10%; margin-right: 10%;"
        )
      ),
      
      # Visualisation options
      shiny::tags$div(
        style = "padding: 0 10px; margin-top: 20px;",
        shiny::h4("Visualisation Options"),
        
        shiny::checkboxInput(
          "show_values",
          "Show correlation values",
          value = FALSE
        ),
        
        shiny::checkboxInput(
          "cluster_isoforms",
          "Cluster isoforms",
          value = TRUE
        )
      )
    ),
    
    # Main body
    shinydashboard::dashboardBody(
      
      # Custom CSS
      shiny::tags$head(
        shiny::tags$style(shiny::HTML("
          .content-wrapper, .right-side {
            background-color: #f4f4f4;
          }
          .nav-tabs > li.active > a {
            border-top: 3px solid #3c8dbc;
          }
          .info-box {
            min-height: 60px;
          }
          .conservation-legend {
            background-color: #fff;
            padding: 15px;
            border-radius: 5px;
            margin: 10px 0;
            border: 1px solid #ddd;
          }
          .switching-box {
            background-color: #f9f9f9;
            padding: 15px;
            border-radius: 5px;
            border: 1px solid #e0e0e0;
            margin: 10px 0;
          }
          .stat-box {
            background-color: #f0f8ff;
            padding: 10px;
            border-radius: 5px;
            margin: 5px 0;
          }
        "))
      ),
      
      # Main tabs
      shiny::tabsetPanel(
        
        # Overview tab
        shiny::tabPanel(
          "Overview",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Gene Summary",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              shiny::htmlOutput("gene_summary")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Cell Type Information",
              width = 12,
              DT::DTOutput("celltype_table")
            )
          )
        ),
        
        # Overall co-expression tab
        shiny::tabPanel(
          "Overall Co-expression",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Gene-wide Co-expression (All Cell Types Combined)",
              width = 12,
              status = "info",
              solidHeader = TRUE,
              shiny::plotOutput("overall_coexpression_heatmap", height = "700px")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Overall Statistics",
              width = 6,
              shiny::tableOutput("overall_stats")
            ),
            shinydashboard::box(
              title = "Expression Distribution",
              width = 6,
              plotly::plotlyOutput("overall_expression_plot", height = "400px")
            )
          )
        ),
        
        # Individual heatmaps tab
        shiny::tabPanel(
          "Cell Type Heatmaps",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Select Cell Type",
              width = 12,
              shiny::selectInput(
                "selected_celltype",
                "Cell Type:",
                choices = NULL,
                width = "100%"
              )
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Co-expression Heatmap",
              width = 8,
              shiny::plotOutput("coexpression_heatmap", height = "600px")
            ),
            shinydashboard::box(
              title = "Heatmap Statistics",
              width = 4,
              shiny::tableOutput("heatmap_stats")
            )
          )
        ),
        
        # Isoform Switching Analysis tab
        shiny::tabPanel(
          "Isoform Switching",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Isoform Switching Analysis",
              width = 12,
              status = "warning",
              solidHeader = TRUE,
              shiny::div(
                class = "switching-box",
                shiny::h4("What is Isoform Switching?"),
                shiny::p("Isoform switching occurs when cells change which isoform of a gene they predominantly express. 
                         This is often indicated by negative correlations between isoforms of the same gene."),
                shiny::p("Key indicators:"),
                shiny::tags$ul(
                  shiny::tags$li(shiny::strong("Strong negative correlation"), " between isoforms (< -0.3)"),
                  shiny::tags$li(shiny::strong("Mutually exclusive"), " expression patterns"),
                  shiny::tags$li(shiny::strong("Different dominant isoforms"), " in different cell types")
                ),
                shiny::br(),
                shiny::div(
                  class = "alert alert-warning",
                  shiny::p(shiny::strong("How to interpret the results:")),
                  shiny::tags$ul(
                    shiny::tags$li("Red cells in the heatmap indicate negative correlations"),
                    shiny::tags$li("The darker the red, the stronger the switching relationship"),
                    shiny::tags$li("Check the table for specific correlation values and cell type patterns")
                  )
                )
              )
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Isoform Relationships",
              width = 8,
              shiny::plotOutput("isoform_switching_plot", height = "600px")
            ),
            shinydashboard::box(
              title = "Switching Statistics",
              width = 4,
              shiny::tableOutput("switching_stats")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Negative Correlation Pairs",
              width = 12,
              DT::DTOutput("switching_pairs_table")
            )
          )
        ),
        
        # Cell Type Similarity tab
        shiny::tabPanel(
          "Cell Type Similarity",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Cell Type Similarity Analysis",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              shiny::div(
                style = "padding: 10px;",
                shiny::h4("How Cell Type Similarity Works"),
                shiny::p("This analysis compares cell types based on their co-expression patterns:"),
                shiny::tags$ul(
                  shiny::tags$li("Similar cell types show similar isoform co-expression patterns"),
                  shiny::tags$li("Can reveal functional relationships between cell types"),
                  shiny::tags$li("Helps identify cell type-specific regulatory programs")
                ),
                shiny::actionButton(
                  "run_celltype_similarity",
                  "Calculate Cell Type Similarity",
                  class = "btn-primary"
                )
              )
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Cell Type Dendrogram",
              width = 6,
              shiny::plotOutput("celltype_dendrogram", height = "500px")
            ),
            shinydashboard::box(
              title = "Similarity Matrix",
              width = 6,
              shiny::plotOutput("celltype_similarity_heatmap", height = "500px")
            )
          )
        ),
        
        # Advanced Statistics tab with detailed results
        shiny::tabPanel(
          "Advanced Statistics",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Statistical Analysis Options",
              width = 12,
              status = "info",
              solidHeader = TRUE,
              shiny::div(
                class = "stat-box",
                shiny::h4("Available Statistical Analyses"),
                shiny::radioButtons(
                  "stat_analysis_type",
                  "Select analysis:",
                  choices = list(
                    "Correlation confidence intervals" = "conf_int",
                    "Bootstrap stability analysis" = "bootstrap",
                    "FDR correction for multiple testing" = "fdr"
                  ),
                  selected = "conf_int"
                ),
                shiny::br(),
                shiny::conditionalPanel(
                  condition = "input.stat_analysis_type == 'conf_int'",
                  shiny::div(
                    class = "alert alert-info",
                    shiny::h5("About Confidence Intervals"),
                    shiny::p("Confidence intervals show the reliability range of correlation coefficients. 
                            A 95% CI means we're 95% confident the true correlation lies within this range."),
                    shiny::tags$ul(
                      shiny::tags$li(shiny::strong("Significant:"), " CI doesn't include 0 (blue points)"),
                      shiny::tags$li(shiny::strong("Not significant:"), " CI includes 0 (grey points)"),
                      shiny::tags$li(shiny::strong("Narrower CI:"), " More precise estimate")
                    )
                  )
                ),
                shiny::conditionalPanel(
                  condition = "input.stat_analysis_type == 'bootstrap'",
                  shiny::div(
                    class = "alert alert-info",
                    shiny::h5("About Bootstrap Stability"),
                    shiny::p("Bootstrap analysis tests how stable correlations are by resampling cells. 
                            Lower SD means more stable relationships."),
                    shiny::tags$ul(
                      shiny::tags$li(shiny::strong("Very stable (SD < 0.1):"), " Robust correlation"),
                      shiny::tags$li(shiny::strong("Stable (SD 0.1-0.2):"), " Moderately robust"),
                      shiny::tags$li(shiny::strong("Unstable (SD > 0.2):"), " May be driven by outliers")
                    ),
                    shiny::p("Lighter colours in the heatmap indicate more stable relationships.")
                  )
                ),
                shiny::conditionalPanel(
                  condition = "input.stat_analysis_type == 'fdr'",
                  shiny::div(
                    class = "alert alert-info",
                    shiny::h5("About FDR Correction"),
                    shiny::p("When testing many correlations, some may appear significant by chance. 
                            FDR correction controls the expected proportion of false discoveries."),
                    shiny::tags$ul(
                      shiny::tags$li(shiny::strong("***"), " q < 0.001 (very high confidence)"),
                      shiny::tags$li(shiny::strong("**"), " q < 0.01 (high confidence)"),
                      shiny::tags$li(shiny::strong("*"), " q < 0.05 (moderate confidence)")
                    ),
                    shiny::p("Darker colours in the heatmap indicate stronger significance after correction.")
                  )
                ),
                shiny::actionButton(
                  "run_stats",
                  "Run Statistical Analysis",
                  class = "btn-info"
                )
              )
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Statistical Summary",
              width = 12,
              shiny::uiOutput("stats_results")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Detailed Results by Isoform Pair",
              width = 12,
              DT::DTOutput("stats_detailed_table")
            )
          )
        ),
        
        # Conservation analysis tab
        shiny::tabPanel(
          "Conservation Analysis",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Conservation Analysis Overview",
              width = 12,
              status = "success",
              solidHeader = TRUE,
              shiny::uiOutput("conservation_status"),
              shiny::div(
                class = "alert alert-success",
                shiny::h5("What is Conservation Analysis?"),
                shiny::p("Conservation analysis identifies co-expression relationships that are maintained across different cell types. 
                        These conserved patterns often indicate fundamental regulatory relationships."),
                shiny::p("The consistency threshold determines how strict the conservation criteria are:")
              ),
              shiny::div(
                class = "conservation-legend",
                shiny::h4("Conservation Pattern Definitions:"),
                shiny::div(
                  class = "legend-item",
                  shiny::span(class = "legend-color", 
                             style = "background-color: #6fafd2;"),
                  shiny::strong("Conserved Positive:"),
                  " Isoform pairs showing consistent positive correlation across cell types"
                ),
                shiny::div(
                  class = "legend-item",
                  shiny::span(class = "legend-color", 
                             style = "background-color: #feab88;"),
                  shiny::strong("Conserved Negative:"),
                  " Isoform pairs showing consistent negative correlation across cell types"
                ),
                shiny::div(
                  class = "legend-item",
                  shiny::span(class = "legend-color", 
                             style = "background-color: #808080;"),
                  shiny::strong("Mixed:"),
                  " Isoform pairs with both positive and negative correlations"
                ),
                shiny::div(
                  class = "legend-item",
                  shiny::span(class = "legend-color", 
                             style = "background-color: #6a0624;"),
                  shiny::strong("Cell Type Specific:"),
                  " Isoform pairs showing correlation in only few cell types"
                )
              )
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Conservation Heatmap",
              width = 8,
              shiny::plotOutput("conservation_heatmap", height = "600px")
            ),
            shinydashboard::box(
              title = "Conservation Statistics",
              width = 4,
              plotly::plotlyOutput("conservation_plot", height = "400px")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Conserved Isoform Pairs",
              width = 12,
              DT::DTOutput("conservation_table")
            )
          )
        ),
        
        # Export tab
        shiny::tabPanel(
          "Export",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Export Options",
              width = 12,
              shiny::h4("Download Heatmaps"),
              shiny::p("Download all co-expression heatmaps as a ZIP file."),
              shiny::downloadButton(
                "download_heatmaps",
                "Download All Heatmaps (PNG)",
                class = "btn-primary"
              ),
              shiny::br(),
              shiny::br(),
              shiny::h4("Download Analysis Results"),
              shiny::p("Download all analysis results including conservation and statistical data."),
              shiny::downloadButton(
                "download_results",
                "Download Complete Results (CSV)",
                class = "btn-info"
              )
            )
          )
        )
      )
    )
  )
  
  # Define server
  server <- function(input, output, session) {
    
    # Get objects from global environment
    scht_obj <- .GlobalEnv$.scht_obj_temp
    
    # Reactive values
    values <- shiny::reactiveValues(
      correlation_list = NULL,
      overall_correlation = NULL,
      conservation_results = NULL,
      valid_celltypes = NULL,
      celltype_similarity = NULL,
      stats_results = NULL,
      stats_detailed = NULL
    )
    
    # Initialize gene selector with server-side processing
    shiny::updateSelectizeInput(
      session, 
      "selected_gene",
      choices = available_genes,
      selected = available_genes[1],
      server = TRUE  # Enable server-side selectize
    )
    
    # Gene info
    output$gene_info <- shiny::renderText({
      paste("Total genes available:", length(available_genes))
    })
    
    # Update cell type choices when gene changes
    shiny::observe({
      gene <- input$selected_gene
      if (!is.null(gene) && gene != "") {
        # Get cell types with this gene
        ct_with_gene <- names(scht_obj$cell_type_matrices)[
          sapply(scht_obj$cell_type_matrices, function(x) gene %in% names(x))
        ]
        
        shiny::updateSelectInput(
          session,
          "selected_celltype",
          choices = ct_with_gene,
          selected = ct_with_gene[1]
        )
      }
    })
    
    # Process gene data
    shiny::observe({
      gene <- input$selected_gene
      if (is.null(gene) || gene == "") return()
      
      # Get overall matrix from original results
      if (gene %in% names(scht_obj$original_results)) {
        overall_matrix <- scht_obj$original_results[[gene]]
        
        # Filter based on expression
        mean_expr <- rowMeans(overall_matrix)
        keep_isoforms <- mean_expr >= input$min_expression
        
        if (sum(keep_isoforms) >= 2) {
          overall_matrix_filtered <- overall_matrix[keep_isoforms, , drop = FALSE]
          overall_cor <- cor(t(overall_matrix_filtered), 
                           method = input$cor_method, 
                           use = "pairwise.complete.obs")
          values$overall_correlation <- overall_cor
        } else {
          values$overall_correlation <- NULL
        }
      }
      
      # Get correlation data for all cell types using the analysis function
      all_result <- calculate_gene_coexpression_all_celltypes(
        scht_obj,
        gene,
        method = input$cor_method,
        min_cells = input$min_cells,
        min_expression = input$min_expression
      )
      
      # Extract cell type data
      correlation_data <- list()
      valid_cts <- character()
      
      if (!is.null(all_result$cell_types)) {
        for (ct_name in names(all_result$cell_types)) {
          ct_data <- all_result$cell_types[[ct_name]]
          
          # Get expression matrix for this cell type
          expr_matrix <- scht_obj$cell_type_matrices[[ct_name]][[gene]]
          mean_expr <- rowMeans(expr_matrix)
          keep_isoforms <- mean_expr >= input$min_expression
          expr_matrix_filtered <- expr_matrix[keep_isoforms, , drop = FALSE]
          
          correlation_data[[ct_name]] <- list(
            cor_matrix = ct_data$cor_matrix,
            expr_matrix = expr_matrix_filtered,
            n_cells = ct_data$n_cells,
            n_isoforms = ct_data$n_isoforms
          )
          valid_cts <- c(valid_cts, ct_name)
        }
      }
      
      values$correlation_list <- correlation_data
      values$valid_celltypes <- valid_cts
    })
    
    # Gene summary
    output$gene_summary <- shiny::renderUI({
      gene <- input$selected_gene
      ct_data <- values$correlation_list
      
      if (is.null(ct_data) || length(ct_data) == 0) {
        return(shiny::HTML("<p>No data available for current settings</p>"))
      }
      
      n_celltypes <- length(ct_data)
      total_cells <- sum(sapply(ct_data, function(x) x$n_cells))
      avg_isoforms <- mean(sapply(ct_data, function(x) x$n_isoforms))
      
      # Get overall stats if available
      n_overall_isoforms <- if (!is.null(values$overall_correlation)) {
        nrow(values$overall_correlation)
      } else {
        0
      }
      
      shiny::fluidRow(
        shinydashboard::infoBox(
          "Selected Gene",
          gene,
          icon = shiny::icon("dna"),
          color = "blue",
          width = 3
        ),
        shinydashboard::infoBox(
          "Cell Types",
          n_celltypes,
          icon = shiny::icon("layer-group"),
          color = "green",
          width = 3
        ),
        shinydashboard::infoBox(
          "Total Cells",
          format(total_cells, big.mark = ","),
          icon = shiny::icon("table-cells"),
          color = "yellow",
          width = 3
        ),
        shinydashboard::infoBox(
          "Total Isoforms",
          n_overall_isoforms,
          icon = shiny::icon("code-branch"),
          color = "red",
          width = 3
        )
      )
    })
    
    # Cell type table
    output$celltype_table <- DT::renderDT({
      ct_data <- values$correlation_list
      
      if (is.null(ct_data) || length(ct_data) == 0) {
        return(data.frame(Message = "No data available"))
      }
      
      summary_df <- data.frame(
        CellType = names(ct_data),
        Cells = sapply(ct_data, function(x) x$n_cells),
        Isoforms = sapply(ct_data, function(x) x$n_isoforms),
        stringsAsFactors = FALSE
      )
      
      DT::datatable(summary_df, 
                    options = list(pageLength = 10),
                    rownames = FALSE)
    })
    
    # Overall co-expression heatmap
    output$overall_coexpression_heatmap <- shiny::renderPlot({
      if (is.null(values$overall_correlation)) {
        plot.new()
        text(0.5, 0.5, "Insufficient data for overall co-expression analysis", 
             cex = 1.5, col = "gray")
        return()
      }
      
      cor_matrix <- values$overall_correlation
      
      # Create heatmap
      colours <- c("#6a0624", "#feab88", "#f7f7f7", "#6fafd2", "#053061")
      
      ht <- ComplexHeatmap::Heatmap(
        cor_matrix,
        name = "Correlation",
        column_title = paste(input$selected_gene, "- All Cell Types"),
        col = circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), colours),
        cell_fun = if (input$show_values) {
          function(j, i, x, y, width, height, fill) {
            grid::grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y,
                           gp = grid::gpar(fontsize = 10))
          }
        } else {
          NULL
        },
        cluster_rows = input$cluster_isoforms,
        cluster_columns = input$cluster_isoforms,
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10)
      )
      
      ComplexHeatmap::draw(ht)
    })
    
    # Overall statistics
    output$overall_stats <- shiny::renderTable({
      if (is.null(values$overall_correlation)) {
        return(data.frame(Metric = "No data", Value = NA))
      }
      
      cor_matrix <- values$overall_correlation
      upper_tri <- cor_matrix[upper.tri(cor_matrix)]
      
      stats_df <- data.frame(
        Metric = c("Number of isoforms", 
                   "Number of pairs", 
                   "Mean correlation", 
                   "Median correlation",
                   "Strong positive (> 0.5)", 
                   "Strong negative (< -0.5)"),
        Value = c(
          nrow(cor_matrix),
          length(upper_tri),
          round(mean(upper_tri), 3),
          round(median(upper_tri), 3),
          sum(upper_tri > 0.5),
          sum(upper_tri < -0.5)
        )
      )
      
      stats_df
    }, digits = 3)
    
    # Overall expression plot
    output$overall_expression_plot <- plotly::renderPlotly({
      gene <- input$selected_gene
      if (is.null(gene) || !gene %in% names(scht_obj$original_results)) return(NULL)
      
      overall_matrix <- scht_obj$original_results[[gene]]
      
      # Create expression summary
      expr_summary <- data.frame(
        isoform = rownames(overall_matrix),
        mean_expression = rowMeans(overall_matrix),
        cells_expressed = rowSums(overall_matrix > 0),
        stringsAsFactors = FALSE
      )
      
      p <- plotly::plot_ly(expr_summary, 
                          x = ~mean_expression, 
                          y = ~cells_expressed,
                          text = ~isoform,
                          type = 'scatter',
                          mode = 'markers',
                          marker = list(
                            size = 10,
                            color = ~mean_expression,
                            colorscale = 'Viridis',
                            showscale = TRUE
                          )) %>%
        plotly::layout(
          title = "Isoform Expression Overview",
          xaxis = list(title = "Mean Expression"),
          yaxis = list(title = "Number of Cells Expressing"),
          hovermode = 'closest'
        )
      
      p
    })
    
    # Individual cell type heatmap
    output$coexpression_heatmap <- shiny::renderPlot({
      ct <- input$selected_celltype
      
      if (is.null(ct) || is.null(values$correlation_list[[ct]])) return()
      
      cor_matrix <- values$correlation_list[[ct]]$cor_matrix
      
      # Create heatmap
      colours <- c("#6a0624", "#feab88", "#f7f7f7", "#6fafd2", "#053061")
      
      ht <- ComplexHeatmap::Heatmap(
        cor_matrix,
        name = "Correlation",
        column_title = paste(input$selected_gene, "-", ct),
        col = circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), colours),
        cell_fun = if (input$show_values) {
          function(j, i, x, y, width, height, fill) {
            grid::grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y,
                           gp = grid::gpar(fontsize = 10))
          }
        } else {
          NULL
        },
        cluster_rows = input$cluster_isoforms,
        cluster_columns = input$cluster_isoforms
      )
      
      ComplexHeatmap::draw(ht)
    })
    
    # Heatmap statistics
    output$heatmap_stats <- shiny::renderTable({
      ct <- input$selected_celltype
      
      if (is.null(ct) || is.null(values$correlation_list[[ct]])) {
        return(data.frame(Metric = "No data", Value = NA))
      }
      
      cor_matrix <- values$correlation_list[[ct]]$cor_matrix
      upper_tri <- cor_matrix[upper.tri(cor_matrix)]
      
      stats_df <- data.frame(
        Metric = c("Number of isoforms", "Number of pairs", 
                   "Mean correlation", "Median correlation",
                   "Min correlation", "Max correlation"),
        Value = c(
          nrow(cor_matrix),
          length(upper_tri),
          round(mean(upper_tri), 3),
          round(median(upper_tri), 3),
          round(min(upper_tri), 3),
          round(max(upper_tri), 3)
        )
      )
      
      stats_df
    }, digits = 3)
    
    # Isoform Switching Analysis
    output$isoform_switching_plot <- shiny::renderPlot({
      if (is.null(values$overall_correlation)) {
        plot.new()
        text(0.5, 0.5, "No data available", cex = 1.5, col = "gray")
        return()
      }
      
      cor_matrix <- values$overall_correlation
      
      # Find negative correlations (potential switching)
      negative_pairs <- which(cor_matrix < -0.3 & upper.tri(cor_matrix), arr.ind = TRUE)
      
      if (nrow(negative_pairs) == 0) {
        plot.new()
        text(0.5, 0.5, "No significant negative correlations found", 
             cex = 1.5, col = "gray")
        return()
      }
      
      # Create a heatmap focusing on negative correlations
      colours <- c("#8B0000", "#CD5C5C", "#FFA07A", "#FFFFFF", "#87CEEB", "#4169E1", "#00008B")
      
      ht <- ComplexHeatmap::Heatmap(
        cor_matrix,
        name = "Correlation",
        column_title = "Isoform Switching Analysis",
        col = circlize::colorRamp2(c(-1, -0.7, -0.3, 0, 0.3, 0.7, 1), colours),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (cor_matrix[i, j] < -0.3 && i != j) {
            grid::grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y,
                           gp = grid::gpar(fontsize = 10, fontface = "bold"))
          }
        }
      )
      
      ComplexHeatmap::draw(ht)
    })
    
    # Switching statistics
    output$switching_stats <- shiny::renderTable({
      if (is.null(values$overall_correlation)) {
        return(data.frame(Metric = "No data", Value = NA))
      }
      
      cor_matrix <- values$overall_correlation
      upper_tri <- cor_matrix[upper.tri(cor_matrix)]
      
      # Count negative correlations
      strong_negative <- sum(upper_tri < -0.5)
      moderate_negative <- sum(upper_tri < -0.3 & upper_tri >= -0.5)
      
      stats_df <- data.frame(
        Metric = c(
          "Total isoform pairs",
          "Strong negative (< -0.5)",
          "Moderate negative (-0.3 to -0.5)",
          "Potential switching pairs",
          "Mean negative correlation"
        ),
        Value = c(
          length(upper_tri),
          strong_negative,
          moderate_negative,
          strong_negative + moderate_negative,
          if (any(upper_tri < 0)) round(mean(upper_tri[upper_tri < 0]), 3) else 0
        )
      )
      
      stats_df
    })
    
    # Negative correlation pairs table
    output$switching_pairs_table <- DT::renderDT({
      if (is.null(values$overall_correlation)) return(NULL)
      
      cor_matrix <- values$overall_correlation
      
      # Get all negative pairs
      negative_pairs <- which(cor_matrix < -0.3 & upper.tri(cor_matrix), arr.ind = TRUE)
      
      if (nrow(negative_pairs) == 0) {
        return(data.frame(Message = "No negative correlations found"))
      }
      
      # Create data frame
      pairs_df <- data.frame(
        Isoform1 = rownames(cor_matrix)[negative_pairs[,1]],
        Isoform2 = colnames(cor_matrix)[negative_pairs[,2]],
        Correlation = round(cor_matrix[negative_pairs], 3),
        stringsAsFactors = FALSE
      )
      
      # Sort by correlation strength
      pairs_df <- pairs_df[order(pairs_df$Correlation), ]
      
      # Add interpretation
      pairs_df$Interpretation <- ifelse(
        pairs_df$Correlation < -0.5,
        "Strong switching",
        "Moderate switching"
      )
      
      DT::datatable(pairs_df,
                    options = list(pageLength = 10),
                    rownames = FALSE) %>%
        DT::formatStyle(
          "Correlation",
          backgroundColor = DT::styleInterval(c(-0.5), c("#ffcccc", "#ff9999"))
        )
    })
    
    # Cell Type Similarity Analysis
    shiny::observeEvent(input$run_celltype_similarity, {
      shiny::withProgress(message = 'Calculating cell type similarity...', value = 0, {
        
        if (is.null(values$correlation_list) || length(values$correlation_list) < 2) {
          shiny::showNotification("Need at least 2 cell types for similarity analysis", 
                                type = "warning")
          return()
        }
        
        # Calculate similarity between cell types based on correlation patterns
        ct_names <- names(values$correlation_list)
        n_ct <- length(ct_names)
        similarity_matrix <- matrix(NA, n_ct, n_ct)
        rownames(similarity_matrix) <- colnames(similarity_matrix) <- ct_names
        
        shiny::incProgress(0.3, detail = "Comparing cell types...")
        
        for (i in 1:n_ct) {
          for (j in i:n_ct) {
            if (i == j) {
              similarity_matrix[i, j] <- 1
            } else {
              # Get common isoforms
              iso1 <- rownames(values$correlation_list[[i]]$cor_matrix)
              iso2 <- rownames(values$correlation_list[[j]]$cor_matrix)
              common_iso <- intersect(iso1, iso2)
              
              if (length(common_iso) >= 2) {
                # Compare correlation patterns
                cor1 <- values$correlation_list[[i]]$cor_matrix[common_iso, common_iso]
                cor2 <- values$correlation_list[[j]]$cor_matrix[common_iso, common_iso]
                
                # Calculate correlation of correlations
                similarity <- cor(as.vector(cor1), as.vector(cor2), 
                                use = "complete.obs")
                similarity_matrix[i, j] <- similarity_matrix[j, i] <- similarity
              } else {
                similarity_matrix[i, j] <- similarity_matrix[j, i] <- NA
              }
            }
          }
        }
        
        values$celltype_similarity <- similarity_matrix
        shiny::incProgress(1, detail = "Complete!")
      })
    })
    
    # Cell type dendrogram
    output$celltype_dendrogram <- shiny::renderPlot({
      if (is.null(values$celltype_similarity)) {
        plot.new()
        text(0.5, 0.5, "Run cell type similarity analysis first", 
             cex = 1.5, col = "gray")
        return()
      }
      
      # Create distance matrix and hierarchical clustering
      dist_matrix <- as.dist(1 - values$celltype_similarity)
      hc <- hclust(dist_matrix, method = "complete")
      
      plot(hc, main = "Cell Type Clustering Based on Co-expression Patterns",
           xlab = "", sub = "", ylab = "Distance (1 - Correlation)")
    })
    
    # Cell type similarity heatmap
    output$celltype_similarity_heatmap <- shiny::renderPlot({
      if (is.null(values$celltype_similarity)) {
        plot.new()
        text(0.5, 0.5, "Run cell type similarity analysis first", 
             cex = 1.5, col = "gray")
        return()
      }
      
      # Create heatmap
      colours <- c("#053061", "#6fafd2", "#f7f7f7", "#feab88", "#6a0624")
      
      ht <- ComplexHeatmap::Heatmap(
        values$celltype_similarity,
        name = "Similarity",
        column_title = "Cell Type Similarity Matrix",
        col = circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), colours),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid::grid.text(sprintf("%.2f", values$celltype_similarity[i, j]), x, y,
                         gp = grid::gpar(fontsize = 10))
        },
        cluster_rows = TRUE,
        cluster_columns = TRUE
      )
      
      ComplexHeatmap::draw(ht)
    })
    
    # Advanced Statistics with detailed results
    shiny::observeEvent(input$run_stats, {
      shiny::withProgress(message = 'Running statistical analysis...', value = 0, {
        
        if (is.null(values$overall_correlation)) {
          shiny::showNotification("No correlation data available", type = "warning")
          return()
        }
        
        stat_type <- input$stat_analysis_type
        cor_matrix <- values$overall_correlation
        n_cells <- ncol(scht_obj$original_results[[input$selected_gene]])
        
        # Get all pairs
        pairs_idx <- which(upper.tri(cor_matrix), arr.ind = TRUE)
        isoform_names <- rownames(cor_matrix)
        
        if (stat_type == "conf_int") {
          # Calculate confidence intervals for correlations
          
          # Initialize results
          ci_results <- data.frame(
            Isoform1 = isoform_names[pairs_idx[,1]],
            Isoform2 = isoform_names[pairs_idx[,2]],
            Correlation = numeric(nrow(pairs_idx)),
            CI_Lower = numeric(nrow(pairs_idx)),
            CI_Upper = numeric(nrow(pairs_idx)),
            CI_Width = numeric(nrow(pairs_idx)),
            Significant = logical(nrow(pairs_idx)),
            stringsAsFactors = FALSE
          )
          
          for (i in 1:nrow(pairs_idx)) {
            row_idx <- pairs_idx[i, 1]
            col_idx <- pairs_idx[i, 2]
            cor_val <- cor_matrix[row_idx, col_idx]
            
            # Fisher transformation
            z_score <- 0.5 * log((1 + cor_val) / (1 - cor_val))
            se_z <- 1 / sqrt(n_cells - 3)
            
            # 95% confidence intervals
            z_lower <- z_score - 1.96 * se_z
            z_upper <- z_score + 1.96 * se_z
            
            # Back transform
            ci_lower <- (exp(2 * z_lower) - 1) / (exp(2 * z_lower) + 1)
            ci_upper <- (exp(2 * z_upper) - 1) / (exp(2 * z_upper) + 1)
            
            ci_results[i, "Correlation"] <- round(cor_val, 3)
            ci_results[i, "CI_Lower"] <- round(ci_lower, 3)
            ci_results[i, "CI_Upper"] <- round(ci_upper, 3)
            ci_results[i, "CI_Width"] <- round(ci_upper - ci_lower, 3)
            ci_results[i, "Significant"] <- sign(ci_lower) == sign(ci_upper)
          }
          
          # Sort by absolute correlation
          ci_results <- ci_results[order(abs(ci_results$Correlation), decreasing = TRUE), ]
          
          values$stats_results <- list(type = "conf_int")
          values$stats_detailed <- ci_results
          
        } else if (stat_type == "bootstrap") {
          # Bootstrap analysis
          gene_matrix <- scht_obj$original_results[[input$selected_gene]]
          n_boot <- 100
          
          # Initialize results
          boot_results <- data.frame(
            Isoform1 = isoform_names[pairs_idx[,1]],
            Isoform2 = isoform_names[pairs_idx[,2]],
            Original_Cor = numeric(nrow(pairs_idx)),
            Bootstrap_Mean = numeric(nrow(pairs_idx)),
            Bootstrap_SD = numeric(nrow(pairs_idx)),
            Stability = character(nrow(pairs_idx)),
            stringsAsFactors = FALSE
          )
          
          # Store bootstrap correlations
          boot_cors <- matrix(NA, nrow(pairs_idx), n_boot)
          
          for (b in 1:n_boot) {
            # Resample cells
            boot_idx <- sample(ncol(gene_matrix), replace = TRUE)
            boot_matrix <- gene_matrix[, boot_idx]
            
            # Filter
            mean_expr <- rowMeans(boot_matrix)
            keep <- mean_expr >= input$min_expression
            
            if (sum(keep) >= 2) {
              boot_cor <- cor(t(boot_matrix[keep, ]), method = input$cor_method)
              
              # Extract correlations for pairs that exist in bootstrap
              for (i in 1:nrow(pairs_idx)) {
                iso1 <- isoform_names[pairs_idx[i,1]]
                iso2 <- isoform_names[pairs_idx[i,2]]
                
                if (iso1 %in% rownames(boot_cor) && iso2 %in% colnames(boot_cor)) {
                  boot_cors[i, b] <- boot_cor[iso1, iso2]
                }
              }
            }
            
            shiny::incProgress(b/n_boot, detail = paste("Bootstrap", b, "of", n_boot))
          }
          
          # Calculate statistics
          for (i in 1:nrow(pairs_idx)) {
            boot_results[i, "Original_Cor"] <- round(cor_matrix[pairs_idx[i,1], pairs_idx[i,2]], 3)
            boot_results[i, "Bootstrap_Mean"] <- round(mean(boot_cors[i,], na.rm = TRUE), 3)
            boot_results[i, "Bootstrap_SD"] <- round(sd(boot_cors[i,], na.rm = TRUE), 3)
            
            # Classify stability
            sd_val <- as.numeric(boot_results[i, "Bootstrap_SD"])
            boot_results[i, "Stability"] <- ifelse(
              sd_val < 0.1, "Very stable",
              ifelse(sd_val < 0.2, "Stable", "Unstable")
            )
          }
          
          # Sort by stability (SD)
          boot_results <- boot_results[order(boot_results$Bootstrap_SD), ]
          
          values$stats_results <- list(type = "bootstrap")
          values$stats_detailed <- boot_results
          
        } else if (stat_type == "fdr") {
          # FDR correction for correlation significance
          
          # Calculate p-values for correlations
          fdr_results <- data.frame(
            Isoform1 = isoform_names[pairs_idx[,1]],
            Isoform2 = isoform_names[pairs_idx[,2]],
            Correlation = numeric(nrow(pairs_idx)),
            P_Value = numeric(nrow(pairs_idx)),
            Q_Value = numeric(nrow(pairs_idx)),
            Significant = logical(nrow(pairs_idx)),
            stringsAsFactors = FALSE
          )
          
          for (i in 1:nrow(pairs_idx)) {
            row_idx <- pairs_idx[i, 1]
            col_idx <- pairs_idx[i, 2]
            cor_val <- cor_matrix[row_idx, col_idx]
            
            # Calculate t-statistic and p-value
            t_stat <- cor_val * sqrt((n_cells - 2) / (1 - cor_val^2))
            p_val <- 2 * pt(abs(t_stat), df = n_cells - 2, lower.tail = FALSE)
            
            fdr_results[i, "Correlation"] <- round(cor_val, 3)
            fdr_results[i, "P_Value"] <- p_val
          }
          
          # FDR correction
          fdr_results$Q_Value <- p.adjust(fdr_results$P_Value, method = "BH")
          fdr_results$Significant <- fdr_results$Q_Value < 0.05
          
          # Add significance stars
          fdr_results$Significance <- ifelse(fdr_results$Q_Value < 0.001, "***",
                                           ifelse(fdr_results$Q_Value < 0.01, "**",
                                                ifelse(fdr_results$Q_Value < 0.05, "*", "")))
          
          # Format p-values
          fdr_results$P_Value <- format(fdr_results$P_Value, scientific = TRUE, digits = 3)
          fdr_results$Q_Value <- round(fdr_results$Q_Value, 4)
          
          # Sort by q-value
          fdr_results <- fdr_results[order(fdr_results$Q_Value), ]
          
          values$stats_results <- list(type = "fdr")
          values$stats_detailed <- fdr_results
        }
        
        shiny::incProgress(1, detail = "Complete!")
      })
    })
    
    # Display statistics results summary
    output$stats_results <- shiny::renderUI({
      if (is.null(values$stats_results)) {
        return(shiny::HTML("<p>Select an analysis and click 'Run Statistical Analysis'</p>"))
      }
      
      results <- values$stats_results
      detailed <- values$stats_detailed
      
      if (results$type == "conf_int") {
        n_sig <- sum(detailed$Significant)
        n_total <- nrow(detailed)
        
        shiny::tagList(
          shiny::h4("Confidence Interval Analysis Results"),
          shiny::p(sprintf("Total correlations: %d", n_total)),
          shiny::p(sprintf("Significant correlations (CI doesn't include 0): %d (%.1f%%)", 
                         n_sig, 100 * n_sig / n_total)),
          shiny::p(sprintf("Mean CI width: %.3f", mean(detailed$CI_Width))),
          shiny::br(),
          plotly::plotlyOutput("ci_plot", height = "500px"),
          shiny::br(),
          shiny::p("See the table below for detailed results by isoform pair.")
        )
        
      } else if (results$type == "bootstrap") {
        stable_counts <- table(detailed$Stability)
        
        shiny::tagList(
          shiny::h4("Bootstrap Stability Analysis Results"),
          shiny::p(sprintf("Total isoform pairs: %d", nrow(detailed))),
          shiny::p("Stability distribution:"),
          shiny::tags$ul(
            shiny::tags$li(sprintf("Very stable (SD < 0.1): %d", 
                                 stable_counts["Very stable"] %||% 0)),
            shiny::tags$li(sprintf("Stable (SD 0.1-0.2): %d", 
                                 stable_counts["Stable"] %||% 0)),
            shiny::tags$li(sprintf("Unstable (SD > 0.2): %d", 
                                 stable_counts["Unstable"] %||% 0))
          ),
          shiny::p(sprintf("Mean stability (SD): %.3f", mean(detailed$Bootstrap_SD))),
          shiny::br(),
          plotly::plotlyOutput("bootstrap_heatmap", height = "500px")
        )
        
      } else if (results$type == "fdr") {
        n_sig <- sum(detailed$Significant)
        n_total <- nrow(detailed)
        
        shiny::tagList(
          shiny::h4("FDR Correction Results"),
          shiny::p(sprintf("Total correlations tested: %d", n_total)),
          shiny::p(sprintf("Significant after FDR correction: %d (%.1f%%)", 
                         n_sig, 100 * n_sig / n_total)),
          shiny::p("FDR-corrected p-values (q-values) control the expected proportion of false positives."),
          shiny::p("Significance levels: *** (q < 0.001), ** (q < 0.01), * (q < 0.05)"),
          shiny::br(),
          plotly::plotlyOutput("fdr_heatmap", height = "500px")
        )
      }
    })
    
    # Detailed statistics table
    output$stats_detailed_table <- DT::renderDT({
      if (is.null(values$stats_detailed)) return(NULL)
      
      detailed <- values$stats_detailed
      
      if (values$stats_results$type == "conf_int") {
        DT::datatable(detailed,
                      options = list(pageLength = 10),
                      rownames = FALSE) %>%
          DT::formatStyle(
            "Significant",
            backgroundColor = DT::styleEqual(TRUE, "lightgreen")
          )
      } else if (values$stats_results$type == "bootstrap") {
        DT::datatable(detailed,
                      options = list(pageLength = 10),
                      rownames = FALSE) %>%
          DT::formatStyle(
            "Stability",
            backgroundColor = DT::styleEqual(
              c("Very stable", "Stable", "Unstable"),
              c("lightgreen", "lightyellow", "lightcoral")
            )
          )
      } else if (values$stats_results$type == "fdr") {
        DT::datatable(detailed,
                      options = list(pageLength = 10),
                      rownames = FALSE) %>%
          DT::formatStyle(
            "Significant",
            backgroundColor = DT::styleEqual(TRUE, "lightgreen")
          ) %>%
          DT::formatStyle(
            "Q_Value",
            color = DT::styleInterval(c(0.01, 0.05), 
                                    c("darkgreen", "orange", "red"))
          )
      }
    })
    
    # Run conservation analysis
    shiny::observeEvent(input$run_conservation, {
      shiny::withProgress(message = 'Running conservation analysis...', value = 0, {
        
        gene <- input$selected_gene
        if (is.null(values$correlation_list) || length(values$correlation_list) < 2) {
          shiny::showNotification("Need at least 2 cell types for conservation analysis", 
                                type = "warning")
          return()
        }
        
        shiny::incProgress(0.3, detail = "Analyzing conservation patterns...")
        
        # Use the analyse_coexpression_conservation function
        tryCatch({
          conservation_result <- analyse_coexpression_conservation(
            scht_obj,
            gene,
            method = input$cor_method,
            min_cells = input$min_cells,
            min_expression = input$min_expression,
            consistency_threshold = input$consistency_threshold / 100,
            correlation_threshold = input$cor_threshold
          )
          
          shiny::incProgress(0.7, detail = "Formatting results...")
          
          # Convert to data frame format expected by the UI
          if (length(conservation_result$pairs) > 0) {
            conservation_df <- do.call(rbind, lapply(names(conservation_result$pairs), function(pair_name) {
              pair_data <- conservation_result$pairs[[pair_name]]
              data.frame(
                isoform_pair = pair_name,
                isoform1 = pair_data$isoform1,
                isoform2 = pair_data$isoform2,
                mean_correlation = pair_data$mean_correlation,
                sd_correlation = pair_data$sd_correlation,
                conservation_pattern = pair_data$pattern,
                n_cell_types = pair_data$n_cell_types,
                consistency = length(which(abs(pair_data$correlations) > input$cor_threshold)) / length(pair_data$correlations),
                stringsAsFactors = FALSE
              )
            }))
            
            values$conservation_results <- conservation_df
          } else {
            values$conservation_results <- NULL
          }
          
        }, error = function(e) {
          shiny::showNotification(
            paste("Error in conservation analysis:", e$message),
            type = "error",
            duration = 10
          )
          values$conservation_results <- NULL
        })
        
        shiny::incProgress(1, detail = "Complete!")
      })
    })
    
    # Conservation status
    output$conservation_status <- shiny::renderUI({
      if (is.null(values$conservation_results)) {
        return(shiny::HTML("<p>Click 'Run Analysis' to analyse conservation patterns across cell types.</p>"))
      }
      
      results <- values$conservation_results
      n_conserved <- sum(results$conservation_pattern %in% c("Conserved_Positive", "Conserved_Negative"))
      n_total <- nrow(results)
      
      shiny::HTML(sprintf(
        "<div class='alert alert-success'>
          <strong>Conservation Analysis Complete!</strong><br>
          Analysed %d isoform pairs across %d cell types.<br>
          Found %d conserved pairs (%.1f%% of total).
        </div>",
        n_total, length(values$correlation_list), n_conserved, 
        100 * n_conserved / n_total
      ))
    })
    
    # Conservation heatmap
    output$conservation_heatmap <- shiny::renderPlot({
      if (is.null(values$conservation_results)) {
        plot.new()
        text(0.5, 0.5, "Run conservation analysis first", cex = 1.5, col = "gray")
        return()
      }
      
      # Create matrix of mean correlations
      results <- values$conservation_results
      all_isoforms <- unique(c(results$isoform1, results$isoform2))
      
      cor_matrix <- matrix(NA, length(all_isoforms), length(all_isoforms))
      rownames(cor_matrix) <- colnames(cor_matrix) <- all_isoforms
      
      for (i in 1:nrow(results)) {
        iso1 <- results$isoform1[i]
        iso2 <- results$isoform2[i]
        val <- results$mean_correlation[i]
        cor_matrix[iso1, iso2] <- val
        cor_matrix[iso2, iso1] <- val
      }
      
      diag(cor_matrix) <- 1
      
      # Create heatmap
      colours <- c("#6a0624", "#feab88", "#f7f7f7", "#6fafd2", "#053061")
      
      ht <- ComplexHeatmap::Heatmap(
        cor_matrix,
        name = "Mean\nCorrelation",
        column_title = paste(input$selected_gene, "- Conservation Across Cell Types"),
        col = circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), colours),
        na_col = "gray90",
        cluster_rows = TRUE,
        cluster_columns = TRUE
      )
      
      ComplexHeatmap::draw(ht)
    })
    
    # Conservation plot
    output$conservation_plot <- plotly::renderPlotly({
      if (is.null(values$conservation_results)) return(NULL)
      
      pattern_counts <- table(values$conservation_results$conservation_pattern)
      
      df <- data.frame(
        pattern = factor(names(pattern_counts), 
                        levels = c("Conserved_Positive", "Conserved_Negative", 
                                  "Mixed", "Cell_Type_Specific")),
        count = as.numeric(pattern_counts)
      )
      
      colours <- c(
        "Conserved_Positive" = "#6fafd2",
        "Conserved_Negative" = "#feab88",
        "Mixed" = "#808080",
        "Cell_Type_Specific" = "#6a0624"
      )
      
      p <- plotly::plot_ly(df, x = ~pattern, y = ~count, type = 'bar',
                          marker = list(color = colours[as.character(df$pattern)])) %>%
        plotly::layout(
          title = "Conservation Patterns",
          xaxis = list(title = ""),
          yaxis = list(title = "Number of Isoform Pairs"),
          showlegend = FALSE
        )
      
      p
    })
    
    # Conservation table
    output$conservation_table <- DT::renderDT({
      if (is.null(values$conservation_results)) return(NULL)
      
      df <- values$conservation_results
      
      df <- df[order(abs(df$mean_correlation), decreasing = TRUE), ]
      
      DT::datatable(
        df[, c("isoform_pair", "mean_correlation", "conservation_pattern", 
               "n_cell_types", "consistency")],
        options = list(pageLength = 10),
        rownames = FALSE
      ) %>%
        DT::formatRound(c("mean_correlation", "consistency"), 3) %>%
        DT::formatStyle(
          "conservation_pattern",
          backgroundColor = DT::styleEqual(
            c("Conserved_Positive", "Conserved_Negative", "Conserved_Neutral", "Cell_Type_Specific"),
            c("#e6f2ff", "#fff5e6", "#f0f0f0", "#ffe6e6")
          )
        )
    })
    
    # Download handlers
    output$download_heatmaps <- shiny::downloadHandler(
      filename = function() {
        paste0("coexpression_", input$selected_gene, "_", Sys.Date(), ".zip")
      },
      content = function(file) {
        # Create temporary directory
        temp_dir <- tempfile()
        dir.create(temp_dir)
        
        # Save overall heatmap
        if (!is.null(values$overall_correlation)) {
          png_file <- file.path(temp_dir, paste0(input$selected_gene, "_Overall.png"))
          png(png_file, width = 800, height = 800)
          
          cor_matrix <- values$overall_correlation
          colours <- c("#6a0624", "#feab88", "#f7f7f7", "#6fafd2", "#053061")
          
          ht <- ComplexHeatmap::Heatmap(
            cor_matrix, 
            name = "Correlation",
            column_title = paste(input$selected_gene, "- All Cell Types"),
            col = circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), colours)
          )
          ComplexHeatmap::draw(ht)
          
          dev.off()
        }
        
        # Save individual heatmaps
        for (ct in names(values$correlation_list)) {
          png_file <- file.path(temp_dir, paste0(input$selected_gene, "_", ct, ".png"))
          png(png_file, width = 800, height = 800)
          
          cor_matrix <- values$correlation_list[[ct]]$cor_matrix
          colours <- c("#6a0624", "#feab88", "#f7f7f7", "#6fafd2", "#053061")
          
          ht <- ComplexHeatmap::Heatmap(
            cor_matrix, 
            name = "Correlation",
            column_title = paste(input$selected_gene, "-", ct),
            col = circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1), colours)
          )
          ComplexHeatmap::draw(ht)
          
          dev.off()
        }
        
        # Zip files
        old_wd <- getwd()
        setwd(temp_dir)
        zip(file, list.files(pattern = "\\.png$"))
        setwd(old_wd)
      }
    )
    
    output$download_results <- shiny::downloadHandler(
      filename = function() {
        paste0("coexpression_results_", input$selected_gene, "_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        
        # Check if openxlsx is available
        if (!requireNamespace("openxlsx", quietly = TRUE)) {
          # Fall back to CSV
          all_results <- list()
          
          if (!is.null(values$conservation_results)) {
            all_results$conservation <- values$conservation_results
          }
          
          if (!is.null(values$stats_detailed)) {
            all_results$statistics <- values$stats_detailed
          }
          
          # Combine results
          combined <- do.call(plyr::rbind.fill, all_results)
          write.csv(combined, file, row.names = FALSE)
          
        } else {
          # Use Excel format
          if (!requireNamespace("openxlsx", quietly = TRUE)) {
            stop("Package 'openxlsx' is required for Excel export. Please install it.")
          }
          wb <- openxlsx::createWorkbook()
          
          # Add conservation results
          if (!is.null(values$conservation_results)) {
            openxlsx::addWorksheet(wb, "Conservation")
            openxlsx::writeData(wb, "Conservation", values$conservation_results)
          }
          
          # Add statistical results
          if (!is.null(values$stats_detailed)) {
            sheet_name <- switch(values$stats_results$type,
                               "conf_int" = "Confidence_Intervals",
                               "bootstrap" = "Bootstrap_Stability",
                               "fdr" = "FDR_Correction")
            openxlsx::addWorksheet(wb, sheet_name)
            openxlsx::writeData(wb, sheet_name, values$stats_detailed)
          }
          
          # Add switching pairs if available
          if (!is.null(values$overall_correlation)) {
            cor_matrix <- values$overall_correlation
            negative_pairs <- which(cor_matrix < -0.3 & upper.tri(cor_matrix), arr.ind = TRUE)
            
            if (nrow(negative_pairs) > 0) {
              switching_df <- data.frame(
                Isoform1 = rownames(cor_matrix)[negative_pairs[,1]],
                Isoform2 = colnames(cor_matrix)[negative_pairs[,2]],
                Correlation = round(cor_matrix[negative_pairs], 3),
                stringsAsFactors = FALSE
              )
              
              openxlsx::addWorksheet(wb, "Isoform_Switching")
              openxlsx::writeData(wb, "Isoform_Switching", switching_df)
            }
          }
          
          openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
        }
      }
    )
    
    # CI plot with isoform pairs
    output$ci_plot <- plotly::renderPlotly({
      if (is.null(values$stats_results) || values$stats_results$type != "conf_int") return(NULL)
      
      ci_data <- values$stats_detailed
      
      # Select top correlations for visualisation
      n_show <- min(20, nrow(ci_data))
      plot_data <- ci_data[1:n_show, ]
      
      # Create hover text
      plot_data$hover_text <- paste(
        "Pair: ", plot_data$Isoform1, " - ", plot_data$Isoform2, "<br>",
        "Correlation: ", round(plot_data$Correlation, 3), "<br>",
        "CI: [", round(plot_data$CI_Lower, 3), ", ", round(plot_data$CI_Upper, 3), "]<br>",
        "Width: ", round(plot_data$CI_Width, 3), "<br>",
        ifelse(plot_data$Significant, "Significant", "Not significant"),
        sep = ""
      )
      
      # Create plotly figure
      p <- plot_ly() %>%
        # Add error bars
        add_trace(
          x = plot_data$Correlation,
          y = 1:n_show,
          type = 'scatter',
          mode = 'markers',
          error_x = list(
            type = 'data',
            symmetric = FALSE,
            array = plot_data$CI_Upper - plot_data$Correlation,
            arrayminus = plot_data$Correlation - plot_data$CI_Lower,
            color = ifelse(plot_data$Significant, "#2196F3", "#9E9E9E"),
            thickness = 2
          ),
          marker = list(
            size = 10,
            color = ifelse(plot_data$Significant, "#2196F3", "#9E9E9E")
          ),
          text = plot_data$hover_text,
          hoverinfo = 'text',
          showlegend = FALSE
        ) %>%
        layout(
          title = "Top Correlations with 95% Confidence Intervals",
          xaxis = list(
            title = "Correlation Coefficient",
            range = c(-1, 1),
            zeroline = FALSE
          ),
          yaxis = list(
            title = "",
            tickmode = "array",
            tickvals = 1:n_show,
            ticktext = paste(plot_data$Isoform1, "-", plot_data$Isoform2),
            automargin = TRUE
          ),
          margin = list(l = 200),
          shapes = list(
            list(
              type = "line",
              x0 = 0, x1 = 0,
              y0 = 0, y1 = n_show + 1,
              line = list(color = "gray50", dash = "dash")
            )
          )
        )
      
      p
    })
    
    # Bootstrap stability heatmap
    output$bootstrap_heatmap <- plotly::renderPlotly({
      if (is.null(values$stats_results) || values$stats_results$type != "bootstrap") return(NULL)
      
      detailed <- values$stats_detailed
      cor_matrix <- values$overall_correlation
      
      # Create stability matrix
      stability_matrix <- matrix(NA, nrow(cor_matrix), ncol(cor_matrix))
      rownames(stability_matrix) <- rownames(cor_matrix)
      colnames(stability_matrix) <- colnames(cor_matrix)
      
      # Create hover text matrix
      hover_matrix <- matrix("", nrow(cor_matrix), ncol(cor_matrix))
      
      for (i in 1:nrow(detailed)) {
        iso1_idx <- which(rownames(cor_matrix) == detailed$Isoform1[i])
        iso2_idx <- which(colnames(cor_matrix) == detailed$Isoform2[i])
        
        if (length(iso1_idx) > 0 && length(iso2_idx) > 0) {
          sd_val <- as.numeric(detailed$Bootstrap_SD[i])
          if (!is.na(sd_val)) {
            stability_matrix[iso1_idx, iso2_idx] <- sd_val
            stability_matrix[iso2_idx, iso1_idx] <- sd_val
          }
        }
        
        if (length(iso1_idx) > 0 && length(iso2_idx) > 0) {
          hover_text <- paste(
            "Pair: ", detailed$Isoform1[i], " - ", detailed$Isoform2[i], "<br>",
            "Correlation: ", round(as.numeric(detailed$Original_Cor[i]), 3), "<br>",
            "Bootstrap SD: ", round(as.numeric(detailed$Bootstrap_SD[i]), 3), "<br>",
            "Stability: ", detailed$Stability[i],
            sep = ""
          )
          
          hover_matrix[iso1_idx, iso2_idx] <- hover_text
          hover_matrix[iso2_idx, iso1_idx] <- hover_text
        }
      }
      
      # Set diagonal
      diag(stability_matrix) <- NA
      for (i in 1:nrow(cor_matrix)) {
        hover_matrix[i, i] <- rownames(cor_matrix)[i]
      }
      
      # Check if we have any valid data
      if (all(is.na(stability_matrix))) {
        return(plotly::plot_ly() %>% 
          plotly::layout(
            title = "No data available for bootstrap stability heatmap",
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE)
          ))
      }
      
      # Create plotly heatmap
      p <- plot_ly(
        x = colnames(stability_matrix),
        y = rownames(stability_matrix),
        z = stability_matrix,
        text = hover_matrix,
        hovertemplate = "%{text}<extra></extra>",
        type = "heatmap",
        colorscale = list(
          c(0, "#F9DDDAFF"),
          c(0.125, "#F3BEC7FF"),
          c(0.25, "#E8A0BCFF"),
          c(0.375, "#D785B5FF"),
          c(0.5, "#BF6DB0FF"),
          c(0.625, "#A159A9FF"),
          c(0.75, "#7C489CFF"),
          c(1, "#573B88FF")
        ),
        zmin = 0,
        zmax = 0.5,
        colorbar = list(title = "Bootstrap SD")
      ) %>%
        layout(
          title = "Bootstrap Stability (SD of correlations)",
          xaxis = list(title = "", tickangle = -45),
          yaxis = list(title = ""),
          plot_bgcolor = 'rgba(0,0,0,0)',
          paper_bgcolor = 'rgba(0,0,0,0)'
        )
      
      p
    })
    
    # FDR significance heatmap
    output$fdr_heatmap <- plotly::renderPlotly({
      if (is.null(values$stats_results) || values$stats_results$type != "fdr") return(NULL)
      
      detailed <- values$stats_detailed
      cor_matrix <- values$overall_correlation
      
      # Create significance matrix
      sig_matrix <- matrix(0, nrow(cor_matrix), ncol(cor_matrix))
      rownames(sig_matrix) <- rownames(cor_matrix)
      colnames(sig_matrix) <- colnames(cor_matrix)
      
      # Create hover text matrix
      hover_matrix <- matrix("", nrow(cor_matrix), ncol(cor_matrix))
      
      for (i in 1:nrow(detailed)) {
        iso1_idx <- which(rownames(cor_matrix) == detailed$Isoform1[i])
        iso2_idx <- which(colnames(cor_matrix) == detailed$Isoform2[i])
        
        # Convert significance to numeric levels
        sig_level <- switch(detailed$Significance[i],
                           "***" = 3,
                           "**" = 2,
                           "*" = 1,
                           0)
        
        sig_matrix[iso1_idx, iso2_idx] <- sig_level
        sig_matrix[iso2_idx, iso1_idx] <- sig_level
        
        hover_text <- paste(
          "Pair: ", detailed$Isoform1[i], " - ", detailed$Isoform2[i], "<br>",
          "Correlation: ", round(detailed$Correlation[i], 3), "<br>",
          "P-value: ", detailed$P_Value[i], "<br>",
          "Q-value: ", detailed$Q_Value[i], "<br>",
          "Significance: ", ifelse(detailed$Significance[i] == "", "Not significant", detailed$Significance[i]),
          sep = ""
        )
        
        hover_matrix[iso1_idx, iso2_idx] <- hover_text
        hover_matrix[iso2_idx, iso1_idx] <- hover_text
      }
      
      # Set diagonal
      diag(sig_matrix) <- NA
      for (i in 1:nrow(cor_matrix)) {
        hover_matrix[i, i] <- rownames(cor_matrix)[i]
      }
      
      # Create plotly heatmap
      p <- plot_ly(
        x = colnames(sig_matrix),
        y = rownames(sig_matrix),
        z = sig_matrix,
        text = hover_matrix,
        hovertemplate = "%{text}<extra></extra>",
        type = "heatmap",
        colorscale = list(
          c(0, "#F9DDDAFF"),
          c(0.33, "#E8A0BCFF"),
          c(0.66, "#BF6DB0FF"),
          c(1, "#7C489CFF")
        ),
        zmin = 0,
        zmax = 3,
        colorbar = list(
          title = "Significance",
          tickmode = "array",
          tickvals = c(0, 1, 2, 3),
          ticktext = c("NS", "*", "**", "***")
        )
      ) %>%
        layout(
          title = "FDR Significance Levels",
          xaxis = list(title = "", tickangle = -45),
          yaxis = list(title = ""),
          plot_bgcolor = 'rgba(0,0,0,0)',
          paper_bgcolor = 'rgba(0,0,0,0)'
        )
      
      p
    })
  }
  
  # Clean up on exit
  onStop(function() {
    if (exists(".scht_obj_temp", envir = .GlobalEnv)) {
      rm(.scht_obj_temp, envir = .GlobalEnv)
    }
  })
  
  # Run the app
  shiny::runApp(
    shiny::shinyApp(ui = ui, server = server),
    port = port,
    launch.browser = launch.browser
  )
}

# Null default operator
`%||%` <- function(x, y) if (is.null(x)) y else x