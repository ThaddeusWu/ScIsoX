#################################################################
#  Package Utility Functions for ScIsoX                         #
#                                                               #
#  Author: [Siyuan Wu & Ulf Schmitz]                            #
#  Institution: [James Cook University]                         #
#  Date: Jul 29, 2025                                           #
#  Package: ScIsoX V1.1.0                                       #
#################################################################

#' @importFrom utils install.packages
NULL

#' Check if a package is available
#' 
#' Internal function to check package availability with informative messages
#' 
#' @param pkg_name Name of the package to check
#' @param stop_if_missing Logical, whether to stop execution if package is missing
#' @param install_message Optional custom installation message
#' 
#' @return Logical indicating whether package is available
#' 
#' @keywords internal
.check_package <- function(pkg_name, stop_if_missing = FALSE, install_message = NULL) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    msg <- sprintf("Package '%s' is required but not installed.", pkg_name)
    
    if (!is.null(install_message)) {
      msg <- paste0(msg, " ", install_message)
    } else {
      # Special messages for Bioconductor packages
      if (pkg_name %in% c("ComplexHeatmap", "circlize")) {
        msg <- paste0(msg, sprintf(" Please install it using: BiocManager::install('%s')", pkg_name))
      } else if (pkg_name == "ggradar") {
        msg <- paste0(msg, " Please install it using: devtools::install_github('ricardo-bion/ggradar')")
      } else {
        msg <- paste0(msg, sprintf(" Please install it using: install.packages('%s')", pkg_name))
      }
    }
    
    if (stop_if_missing) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg)
      return(FALSE)
    }
  }
  return(TRUE)
}

#' Install all suggested packages for ScIsoX
#' 
#' This function helps users install all optional packages that enhance ScIsoX functionality.
#' It provides a convenient way to ensure all visualisation and analysis features are available.
#' 
#' @param include_bioc Logical, whether to include Bioconductor packages (default: TRUE)
#' @param include_github Logical, whether to include GitHub packages (default: TRUE)
#' 
#' @return NULL (invisibly). Messages indicate installation progress.
#' 
#' @examples
#' \dontrun{
#' # Install all suggested packages
#' install_scisox_suggests()
#' 
#' # Install only CRAN packages
#' install_scisox_suggests(include_bioc = FALSE, include_github = FALSE)
#' }
#' 
#' @export
install_scisox_suggests <- function(include_bioc = TRUE, include_github = TRUE) {
  # CRAN packages
  cran_pkgs <- c("ggridges", "ggrepel", "ggExtra", "viridis", 
                 "RColorBrewer", "cowplot", "tidyr", "patchwork")
  
  # Check which are not installed
  to_install <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
  
  if (length(to_install) > 0) {
    message("Installing CRAN packages: ", paste(to_install, collapse = ", "))
    install.packages(to_install)
  } else {
    message("All CRAN suggested packages are already installed.")
  }
  
  # Bioconductor packages
  if (include_bioc) {
    bioc_pkgs <- c("ComplexHeatmap", "circlize")
    bioc_missing <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
    
    if (length(bioc_missing) > 0) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        message("Installing BiocManager...")
        install.packages("BiocManager")
      }
      message("Installing Bioconductor packages: ", paste(bioc_missing, collapse = ", "))
      BiocManager::install(bioc_missing, update = FALSE, ask = FALSE)
    } else {
      message("All Bioconductor suggested packages are already installed.")
    }
  }
  
  # GitHub packages
  if (include_github) {
    if (!requireNamespace("ggradar", quietly = TRUE)) {
      if (!requireNamespace("devtools", quietly = TRUE)) {
        message("Installing devtools...")
        install.packages("devtools")
      }
      message("Installing ggradar from GitHub...")
      devtools::install_github("ricardo-bion/ggradar", upgrade = "never")
    } else {
      message("ggradar is already installed.")
    }
  }
  
  message("\nInstallation complete! All suggested packages have been processed.")
  invisible(NULL)
}

#' Get alternative colour palette
#' 
#' Internal function to provide fallback colour palettes when packages are not available
#' 
#' @param n Number of colours needed
#' @param type Type of palette ("sequential", "diverging", "qualitative")
#' 
#' @return Vector of colour codes
#' 
#' @keywords internal
.get_fallback_colours <- function(n, type = "qualitative") {
  if (type == "sequential") {
    # Fallback for viridis-like sequential palette
    colours <- c("#440154", "#482878", "#3e4989", "#31688e", "#26828e", 
                "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725")
    return(colorRampPalette(colours)(n))
  } else if (type == "diverging") {
    # Fallback for diverging palette
    colours <- c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
                "#fddbc7", "#f4a582", "#d6604d", "#b2182b")
    return(colorRampPalette(colours)(n))
  } else {
    # Fallback for qualitative palette (like RColorBrewer Set3)
    if (n <= 12) {
      colours <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", 
                  "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
                  "#ccebc5", "#ffed6f")
      return(colours[1:n])
    } else {
      # For larger n, use a combination of palettes
      base_colours <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", 
                       "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
                       "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4", "#b2df8a",
                       "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00")
      if (n <= 20) {
        return(base_colours[1:n])
      } else {
        # Generate additional colours
        return(colorRampPalette(base_colours)(n))
      }
    }
  }
}