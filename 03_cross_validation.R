# =============================================================================
# Script: 03_cross_validation.R
# Project: Linear and Generalized Linear Models with Poisson Distribution
#          Applied to Count Data in Neuropsychology: A Comparative Analysis
# Author: Diego Rivera
# Date: July 27, 2025
# Purpose: Perform cross-validation for Bayesian Gaussian linear models and Poisson GLMs
# =============================================================================



# ————————————————————
# Reproducibility: set a global seed
set.seed(202506)

# Note: All required packages (rjags, coda, tidyverse, caret, future.apply, parallel) are loaded in module 01.

# ---- Prepare validation dataset ----
# Subset test data and generate polynomial terms for df2

df2 <- data[ind == 2, ]
df2 <- df2[, 3:18]
# Reuse polynomial basis from df
age_basis <- poly(df$Edad, 2)
df2_poly  <- predict(age_basis, newdata = df2$Edad)
df2 <- df2 %>%
  mutate(
    age1 = df2_poly[, 1],
    age2 = df2_poly[, 2]
  )

# ---- 1. Cross-validation for Gaussian linear models ----
# Configure parallel backend
n_cores   <- parallel::detectCores()
n_workers <- max(n_cores - 1, 1)
plan(multisession, workers = n_workers)
message("Parallel plan: multisession with ", n_workers, " workers")

# Initialization for LM
init_values_lm <- function(predictors) {
  list(
    beta  = rnorm(length(predictors) + 1, 0, 5),
    sig   = runif(1, 0, 1000)
  )
}

# JAGS model generator for LM
generate_jags_model_lm <- function(varname, predictors) {
  pred_exprs <- map_chr(predictors, ~ if (str_detect(.x, "_")) {
    parts <- str_split(.x, "_", simplify = TRUE)
    paste0(parts[1], "[i] * ", parts[2], "[i]")
  } else {
    paste0(.x, "[i]")
  })
  terms <- map2_chr(seq_along(pred_exprs) + 1, pred_exprs,
                    ~ paste0("beta[", .x, "] * ", .y))
  pred_line <- if (length(terms)) {
    paste("mu[i] <- beta[1] +", paste(terms, collapse = " + "))
  } else {
    "mu[i] <- beta[1]"
  }
  
  model_lines <- c(
    "model {",
    "  for (i in 1:n) {",
    "    y[i] ~ dnorm(mu[i], tau_y)",
    paste0("    ", pred_line),
    "  }",
    sprintf("  for (j in 1:%d) {", length(predictors) + 1),
    "    beta[j] ~ dnorm(0.0, 0.0001)",
    "  }",
    "  tau_y <- 1 / (sig * sig)",
    "  sig   ~ dunif(0, 1000)",
    "}"
  )
  
  file <- paste0("model_", varname, "_lm.txt")
  write_lines(model_lines, file)
  file
}

# Cross-validation function for LM
cross_validate_jags_lm_parallel <- function(varname, predictors, k = 10) {
  folds <- caret::createFolds(df2[[varname]], k = k)
  future_lapply(seq_along(folds), function(i) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(df2)), test_idx)
    df_train  <- df2[train_idx, ]
    df_test   <- df2[test_idx, ]
    
    # Fit JAGS model
    model_file <- generate_jags_model_lm(varname, predictors)
    jags_data  <- list(y = df_train[[varname]], n = nrow(df_train))
    base_vars  <- unique(unlist(str_split(predictors, "_")))
    for (v in base_vars) {
      jags_data[[v]] <- case_when(
        v == "sex" ~ as.numeric(df_train$Sexo == "Female"),
        v == "edu" ~ log(df_train$Escolaridad),
        TRUE        ~ df_train[[v]]
      )
    }
    jm <- jags.model(file = model_file,
                     data = jags_data,
                     inits = function() init_values_lm(predictors),
                     n.chains = 3,
                     n.adapt  = 5000)
    update(jm, n.iter = 1000)
    samp <- coda.samples(jm,
                         variable.names = c("beta", "sig"),
                         n.iter = 10000,
                         thin   = 50)
    betas <- summary(samp)$statistics %>%
      .[grep("beta", rownames(.)), "Mean"]
    
    # Build X_test
    X_test <- matrix(1, nrow = nrow(df_test), ncol = 1)
    for (pred in predictors) {
      newcol <- if (str_detect(pred, "_")) {
        parts <- str_split(pred, "_", simplify = TRUE)
        v1 <- if (parts[1] == "edu") log(df_test$Escolaridad) else if (parts[1] == "sex") as.numeric(df_test$Sexo == "Female") else df_test[[parts[1]]]
        v2 <- if (parts[2] == "edu") log(df_test$Escolaridad) else if (parts[2] == "sex") as.numeric(df_test$Sexo == "Female") else df_test[[parts[2]]]
        v1 * v2
      } else if (pred == "sex") {
        as.numeric(df_test$Sexo == "Female")
      } else if (pred == "edu") {
        log(df_test$Escolaridad)
      } else {
        df_test[[pred]]
      }
      X_test <- cbind(X_test, newcol)
    }
    
    preds <- X_test %*% betas
    rmse  <- sqrt(mean((df_test[[varname]] - preds)^2))
    mae   <- mean(abs(df_test[[varname]] - preds))
    tibble(model = varname, type = "LM", fold = i, RMSE = rmse, MAE = mae)
  }, future.seed = TRUE) %>%
    bind_rows()
}

# Loop over LM specs
cv_all_results <- map(model_specs_lm, ~ cross_validate_jags_lm_parallel(.y, .x, k = 10))
cv_lm_df     <- bind_rows(cv_all_results, .id = "model")
cv_summary_lm <- cv_lm_df %>%
  group_by(model, type) %>%
  summarise(mean_RMSE = mean(RMSE), mean_MAE = mean(MAE), .groups = "drop")
readr::write_csv(cv_summary_lm, file.path("results", "cross_validation_LM_summary_parallel.csv"))

# ---- 2. Cross-validation for Poisson GLMs ----
# Initialization for GLM
init_values_glm <- function(predictors) {
  list(beta = rnorm(length(predictors) + 1, 0, 0.01))
}

# JAGS model generator for GLM
generate_jags_model_glm <- function(varname, predictors) {
  pred_exprs <- map_chr(predictors, ~ if (str_detect(.x, "_")) {
    parts <- str_split(.x, "_", simplify = TRUE)
    paste0(parts[1], "[i] * ", parts[2], "[i]")
  } else {
    paste0(.x, "[i]")
  })
  terms <- map2_chr(seq_along(pred_exprs) + 1, pred_exprs,
                    ~ paste0("beta[", .x, "] * ", .y))
  pred_line <- if (length(terms)) {
    paste("log(lambda[i]) <- beta[1] +", paste(terms, collapse = " + "))
  } else {
    "log(lambda[i]) <- beta[1]"
  }
  model_lines <- c(
    "model {",
    "  for (i in 1:n) {",
    "    y[i] ~ dpois(lambda[i])",
    paste0("    ", pred_line),
    "  }",
    sprintf("  for (j in 1:%d) {", length(predictors) + 1),
    "    beta[j] ~ dnorm(0.0, 0.0001)",
    "  }",
    "}"
  )
  file <- paste0("model_", varname, "_glm.txt")
  write_lines(model_lines, file)
  file
}

# Cross-validation for GLM
cross_validate_jags_glm_parallel <- function(varname, predictors, k = 10) {
  folds <- caret::createFolds(df2[[varname]], k = k)
  future_lapply(seq_along(folds), function(i) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(df2)), test_idx)
    df_train  <- df2[train_idx, ]
    df_test   <- df2[test_idx, ]
    
    model_file <- generate_jags_model_glm(varname, predictors)
    jags_data  <- list(y = df_train[[varname]], n = nrow(df_train))
    base_vars  <- unique(unlist(str_split(predictors, "_")))
    for (v in base_vars) {
      jags_data[[v]] <- case_when(
        v == "sex" ~ as.numeric(df_train$Sexo == "Female"),
        v == "edu" ~ log(df_train$Escolaridad),
        TRUE        ~ df_train[[v]]
      )
    }
    jm <- jags.model(file = model_file,
                     data = jags_data,
                     inits = function() init_values_glm(predictors),
                     n.chains = 3,
                     n.adapt  = 5000)
    update(jm, n.iter = 1000)
    samp <- coda.samples(jm,
                         variable.names = "beta",
                         n.iter = 10000,
                         thin   = 50)
    betas <- summary(samp)$statistics %>%
      .[grep("beta", rownames(.)), "Mean"]
    
    X_test <- matrix(1, nrow = nrow(df_test), ncol = 1)
    for (pred in predictors) {
      newcol <- if (str_detect(pred, "_")) {
        parts <- str_split(pred, "_", simplify = TRUE)
        v1 <- if (parts[1] == "edu") log(df_test$Escolaridad) else if (parts[1] == "sex") as.numeric(df_test$Sexo == "Female") else df_test[[parts[1]]]
        v2 <- if (parts[2] == "edu") log(df_test$Escolaridad) else if (parts[2] == "sex") as.numeric(df_test$Sexo == "Female") else df_test[[parts[2]]]
        v1 * v2
      } else if (pred == "sex") {
        as.numeric(df_test$Sexo == "Female")
      } else if (pred == "edu") {
        log(df_test$Escolaridad)
      } else {
        df_test[[pred]]
      }
      X_test <- cbind(X_test, newcol)
    }
    
    log_lambda <- X_test %*% betas
    preds      <- exp(log_lambda)
    rmse       <- sqrt(mean((df_test[[varname]] - preds)^2))
    mae        <- mean(abs(df_test[[varname]] - preds))
    tibble(model = varname, type = "GLM", fold = i, RMSE = rmse, MAE = mae)
  }, future.seed = TRUE) %>%
    bind_rows()
}

# Loop over GLM specs
cv_all_results_glm <- map(model_specs_glm, ~ cross_validate_jags_glm_parallel(.y, .x, k = 10))
cv_glm_df        <- bind_rows(cv_all_results_glm, .id = "model")
cv_summary_glm   <- cv_glm_df %>%
  group_by(model, type) %>%
  summarise(mean_RMSE = mean(RMSE), mean_MAE = mean(MAE), .groups = "drop")
readr::write_csv(cv_summary_glm, file.path("results", "cross_validation_GLM_summary_parallel.csv"))

