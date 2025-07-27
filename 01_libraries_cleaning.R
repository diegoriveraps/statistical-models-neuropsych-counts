# =============================================================================
# Project: Linear and Generalized Linear Models with Poisson Distribution 
#          Applied to Count Data in Neuropsychology: A Comparative Analysis
# Author: Diego Rivera
# Date: July 27, 2025
# Script: Package loading, auxiliary functions, data import and preprocessing
# =============================================================================

# ---- Package installation and loading (automatic) ----

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  BAS, BlandAltmanLeh, caret, coda, dplyr, foreign, future.apply,
  knitr, kableExtra, patchwork, PerformanceAnalytics, pROC,
  psych, RCurl, rjags, runjags, tidyverse, tidyr
)


# ---- Define folder paths ----

folders <- c(
  "data",        # Raw and processed data
  "models",      # Fitted model objects (.Rdata, etc.)
  "results",     # Tables, summaries, exportable figures
  "scripts",     # Modular R scripts (01, 02, 03...)
  "utils",       # Utility functions (e.g., DBDA2E-utilities.R)
  "figures"     # Diagnostic and final plots
)

# ---- Create folders if they don't exist ----

for (f in folders) {
  if (!dir.exists(f)) dir.create(f, recursive = TRUE)
}


# ---- Auxiliary functions ----
# The following script is sourced from:
# Kruschke, J. K. (2015). *Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan (2nd ed.)*.
# File: DBDA2E-utilities.R
# It includes custom functions for MCMC diagnostics, plotting, and summarizing posterior samples.
# Download from: https://sites.google.com/site/doingbayesiandataanalysis/
# (see "Programs & Data" section for the Second Edition)
source("utils/DBDA2E-utilities.R")  # Path must point to the downloaded script




# ---- Load and subset data ----

data <- read.spss("data/Special Issue trabajo15072015_4.sav",
                  use.value.labels = TRUE,
                  to.data.frame = TRUE)

# Select relevant variables and subset to Mexican participants
data <- data[, c(1, 3, 4, 6, 7, 10, 11, 31:33, 70, 46:52)]
data <- data[data$Pais == "Mexico", ]
data <- data %>%
  filter(Stroop_palabras != 100)
data <- na.exclude(data)


# ---- Split datasets ----

set.seed(4567)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
df <- data[ind == 1, ]
df2 <- data[ind == 2, ]


# ---- Create polynomial age terms ----

age <- poly(df$Edad, 2)
df$age1 <- age[, 1]
df$age2 <- age[, 2]

