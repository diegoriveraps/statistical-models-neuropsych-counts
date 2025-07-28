# =============================================================================
# Script: 03_score_to_percentile.R
# Project: Linear and Generalized Linear Models with Poisson Distribution
#          Applied to Count Data in Neuropsychology: A Comparative Analysis
# Author: Diego Rivera
# Date: July 27, 2025
# Purpose: Generate normative percentiles for LM and Poisson GLM models
# Dependencies: data, df, ind, age basis, model_specs_lm, model_specs_glm loaded from previous modules
# =============================================================================

# ---- 0. Load required package for Excel import ----
library(readxl)

# ---- 1. Prepare combined dataset ----
# Subset test data
df2 <- data[ind == 2, ]
df2 <- df2[, 3:18]

# Load clinical sample data
data_dx_al <- read_excel(file.path("data", "df_alz.xlsx")) %>%
  na.exclude()
data_dx <- data_dx_al[, c(1, 3:18)] %>%
  na.exclude()

# Align columns and tag healthy controls
df2 <- df2 %>%
  mutate(Grupo = "HC") %>%
  select(Grupo, everything())

# Combine clinical and healthy samples
combined_df <- bind_rows(
  data_dx,
  df2
)

# Add polynomial age terms (using basis 'age' from module 01)
ed_poly <- predict(age, newdata = combined_df$Edad)
combined_df <- combined_df %>%
  mutate(
    age1 = ed_poly[, 1],
    age2 = ed_poly[, 2]
  )

# Rename for consistency
combined_df <- combined_df %>%
  rename(
    edu = Escolaridad,
    sex = Sexo
  )

# ---- 2. Load model specifications ----
# Expects model_specs_lm and model_specs_glm in model_specs.RData
load(file.path("models", "model_specs.RData"))  # or adjust path

# ---- 3. Function: predict percentile using Gaussian LM posterior ----
predict_percentile_lm <- function(samples, data_df, score_var, predictors) {
  # Combine chains into matrix
  samp_mat <- as.matrix(do.call(rbind, samples))
  beta <- samp_mat[, grep("^beta\\[", colnames(samp_mat))]
  sigma <- samp_mat[, "sig"]
  
  # Initialize output
  percentiles <- numeric(nrow(data_df))
  
  for (i in seq_len(nrow(data_df))) {
    row <- data_df[i, ]
    # Build design vector: intercept + predictors
    X <- c(1)
    for (p in predictors) {
      col_val <- case_when(
        p == "edu" ~ log(row$edu),
        p == "sex" ~ as.numeric(row$sex == "Female"),
        str_detect(p, "_") ~ {
          parts <- str_split(p, "_", simplify = TRUE)
          v1 <- if (parts[1] == "edu") log(row$edu) else if (parts[1] == "sex") as.numeric(row$sex == "Female") else row[[parts[1]]]
          v2 <- if (parts[2] == "edu") log(row$edu) else if (parts[2] == "sex") as.numeric(row$sex == "Female") else row[[parts[2]]]
          v1 * v2
        },
        TRUE ~ row[[p]]
      )
      X <- c(X, col_val)
    }
    # Check dimensions
    if (length(X) != ncol(beta)) {
      stop("Dimension mismatch in row ", i)
    }
    # Compute linear predictor distribution
    mu_i <- beta %*% X
    # Observed score
    x <- row[[score_var]]
    # CDF under normal posterior for each sample
    ps <- pnorm(x, mean = mu_i, sd = sigma)
    percentiles[i] <- mean(ps)
  }
  return(percentiles)
}

# ---- 4. Compute LM percentiles for each score ----
scores <- names(model_specs_lm)
percentiles_lm <- tibble(id = seq_len(nrow(combined_df)))
for (sc in scores) {
  message("Computing LM percentiles for ", sc)
  samples <- readRDS(file.path("models", "modelos_lineales_jags", paste0(sc, "_linear_mcmc.rds")))
  preds   <- predict_percentile_lm(samples, combined_df, sc, model_specs_lm[[sc]])
  percentiles_lm[[paste0(sc, "_percentile_lm")]] <- preds
}

# Save LM percentiles
write_csv(percentiles_lm, file.path("results", "percentiles_lm.csv"))

# ---- 5. Function: predict percentile using Poisson GLM posterior ----
predict_percentile_glm <- function(samples, data_df, score_var, predictors) {
  samp_mat <- as.matrix(do.call(rbind, samples))
  beta <- samp_mat[, grep("^beta\\[", colnames(samp_mat))]
  
  percentiles <- numeric(nrow(data_df))
  for (i in seq_len(nrow(data_df))) {
    row <- data_df[i, ]
    X <- c(1)
    for (p in predictors) {
      col_val <- case_when(
        p == "edu" ~ log(row$edu),
        p == "sex" ~ as.numeric(row$sex == "Female"),
        str_detect(p, "_") ~ {
          parts <- str_split(p, "_", simplify = TRUE)
          v1 <- if (parts[1] == "edu") log(row$edu) else if (parts[1] == "sex") as.numeric(row$sex == "Female") else row[[parts[1]]]
          v2 <- if (parts[2] == "edu") log(row$edu) else if (parts[2] == "sex") as.numeric(row$sex == "Female") else row[[parts[2]]]
          v1 * v2
        },
        TRUE ~ row[[p]]
      )
      X <- c(X, col_val)
    }
    if (length(X) != ncol(beta)) stop("Dimension mismatch in row ", i)
    lambda <- exp(beta %*% X)
    x <- row[[score_var]]
    ps <- ppois(x, lambda = lambda)
    percentiles[i] <- mean(ps)
  }
  return(percentiles)
}

# ---- 6. Compute GLM percentiles ----
percentiles_glm <- tibble(id = seq_len(nrow(combined_df)))
for (sc in scores) {
  message("Computing GLM percentiles for ", sc)
  samples <- readRDS(file.path("models", "modelos_glm_jags", paste0(sc, "_glm_mcmc.rds")))
  preds   <- predict_percentile_glm(samples, combined_df, sc, model_specs_glm[[sc]])
  percentiles_glm[[paste0(sc, "_percentile_glm")]] <- preds
}

# Save GLM percentiles
write_csv(percentiles_glm, file.path("results", "percentiles_glm.csv"))

# ---- 7. Add group variable to output tables ----
percentiles_lm  <- percentiles_lm  %>% mutate(Grupo = combined_df$Grupo)
percentiles_glm <- percentiles_glm %>% mutate(Grupo = combined_df$Grupo)

# Optionally inspect
print(head(percentiles_lm))
print(head(percentiles_glm))
