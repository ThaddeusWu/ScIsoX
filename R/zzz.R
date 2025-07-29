#################################################################
#  Package Startup Functions for ScIsoX                         #
#                                                               #
#  Author: [Siyuan Wu & Ulf Schmitz]                            #
#  Institution: [James Cook University]                         #
#  Date: Jul 29, 2025                                           #
#  Package: ScIsoX V1.1.0                                       #
#################################################################

#' Package attach hook
#' 
#' This function runs when the package is attached, checking for suggested packages
#' and providing helpful messages to users.
#' 
#' @param libname Library name
#' @param pkgname Package name
#' 
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  # Define suggested packages and their uses
  suggested_pkgs <- list(
    visualisation = c("ggridges", "ggrepel", "ggExtra"),
    advanced_plots = c("ggradar", "cowplot", "patchwork"),
    data_manipulation = c("tidyr")
  )
  
  # Check which packages are not installed
  missing_pkgs <- list()
  for (category in names(suggested_pkgs)) {
    pkgs <- suggested_pkgs[[category]]
    not_installed <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
    if (length(not_installed) > 0) {
      missing_pkgs[[category]] <- not_installed
    }
  }
  
  # Only show message if there are missing packages
  if (length(missing_pkgs) > 0) {
    packageStartupMessage(
      "\n",
      "===============================================\n",
      "ScIsoX: Single-Cell Isoform Complexity Analysis\n",
      "===============================================\n",
      "\n",
      "Some optional packages are not installed:\n"
    )
    
    for (category in names(missing_pkgs)) {
      packageStartupMessage(
        sprintf("  - For %s: %s\n", 
                category, 
                paste(missing_pkgs[[category]], collapse = ", "))
      )
    }
    
    packageStartupMessage(
      "\nSome visualisation features may not be available.\n",
      "To install all optional packages, run:\n",
      "  install_scisox_suggests()\n",
      "\n",
      "For more information, visit:\n",
      "  https://github.com/ThaddeusWu/ScIsoX\n"
    )
  } else {
    # All packages installed - show brief welcome message
    packageStartupMessage(
      "\n",
      "ScIsoX v", utils::packageVersion(pkgname), " successfully loaded.\n",
      "All optional packages are installed.\n"
    )
  }
}