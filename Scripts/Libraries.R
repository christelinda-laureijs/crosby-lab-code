# Load Libraries

# Data Management
suppressPackageStartupMessages({
  library(here) # Easier file management
  library(abftools) # Required to read .ABF files
  library(tidyr) # Reformat data frames from wide to long and vice versa
  library(dplyr) # Data manipulation
  library(stringr) # Manipulating text objects
  library(glue) # Required for easy mixing of variables and text when naming objects
  library(purrr) # Required to map over lists and reduce repetitive code
  library(reactable) # Required for display tables
  library(htmltools) # Required for reactable
  library(reactablefmtr) # required for sparklines; also requires dataui to be installed
  # Plotting
  
  library(ggplot2) # Make plots
  library(viridis) # Accessible and beautiful colour palettes
  library(extrafont) # Required for custom fonts in plots
  library(grid) # Required for rasterGrob function to have a background image on a plot
  library(gridExtra) # Required for more rasterGrob functions on cell coordinates plot
  library(ggtext) # Required for formatting plot text (e.g. coloured text in title)
  library(ggforce) # Required for sina plots
  library(patchwork) # Required for multi-plot layout
  library(ggsignif) # Required for significance stars and brackets on plots
  library(ggh4x) # Required to create nested x-axis labels
  library(ggpubr) # Required for significance brackets
  library(ggblend) # Required for blending colours on overlapping points
  
  # Statistical packages
  
  library(ggfortify) # Autoplot multiple diagnostic plots at once
  library(nlme) # Required for linear mixed models and other statistical tests
  library(MASS) # Get studentized residuals for linear regression
  library(rcompanion) # Plot residuals against a normalized histogram
  library(car) # Provides ncvTest and other functions
  library(heplots) # Required for Box's M test for MANOVA
  library(RVAideMemoire) # Required for multivariate normality tests
  library(broom) # Required to format statistical test output into tidy dataframes
  library(rstatix) # Required to use pipes with statistical tests
  library(knitr) # Required to access better table formatting through kable
  library(kableExtra) # Required for full width tables
  library(lazyWeave) # Provides the pvalString() function to format publication-ready p-values in a table
  
})


# Redefine the function select() to prioritize dplyr select()
# This is because the MASS select() function masks dplyr's select()

select <- dplyr::select
filter <- dplyr::filter



