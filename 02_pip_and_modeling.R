# =============================================================================
# Script: 02_pip_modeling.R
# Project: Linear and Generalized Linear Models with Poisson Distribution 
#          Applied to Count Data in Neuropsychology: A Comparative Analysis
# Author: Diego Rivera
# Date: July 27, 2025
# Script: Bayesian variable selection and model fitting 
# =============================================================================


# ---- Bayesian variable selection and model fitting using BAS  ----

compute <- TRUE

if (compute) {
  
  # ---- Model fitting function for LM and Poisson GLM ----
  
  fit_model <- function(IDscore = NULL, data) {
    dd <- data.frame(
      YY   = df[[IDscore]],
      age1 = df$age1,
      age2 = df$age2,
      edu  = df$Escolaridad,
      sex  = df$Sexo
    )
    
    formula <- "YY ~ (age1 + age2 + log(edu) + sex)^2 - age1:age2"
    
    model.lm  <- bas.lm(formula, data = dd)
    model.glm <- bas.glm(formula, family = poisson(), data = dd)
    
    return(list(model.lm = model.lm, model.glm = model.glm))
  }
  
  # ---- Target neuropsychological variables (count scores) ----
  
  aux <- c(
    "Stroop_palabras", "Stroop_colores", "Stroop_PC", "SDMT",
    "Fonologica_F", "Fonologica_A", "Fonologica_S", "Fonologica_M",
    "Semantica_animales", "Semantica_frutas", "Semantica_profesiones"
  )
  
  # ---- Parallel model fitting ----
  
  plan("multisession")
  models <- future_lapply(aux, function(x) fit_model(IDscore = x, data = df), future.seed = TRUE)
  names(models) <- aux
  
  # ---- Save fitted models ----
  
  save(models, file = "models/PIP_models.Rdata")
  
} else {
  
  # ---- Load precomputed models ----
  
  load("models/PIP_models.Rdata")
}


# ---- Extract predictors with PIP > 0.5 ----

linear_models <- lapply(models, function(x) sort(summary(x[[1]])[, 1], decreasing = TRUE))
linear_models <- lapply(linear_models, function(x) x[x > 0.5])

genlinear_models <- lapply(models, function(x) sort(summary(x[[2]])[, 1], decreasing = TRUE))
genlinear_models <- lapply(genlinear_models, function(x) x[x > 0.5])



# =============================================================================
# Fit Bayesian Gaussian linear models and Poisson GLMs in JAGS
# =============================================================================


# ---- Note ----
# Dependencies (rjags, coda, tidyverse, etc.) are loaded in module 01.

# ---- 1. Model specifications for Gaussian linear regression ----
model_specs_lm <- list(
  Stroop_palabras       = c("age1", "age2", "edu", "age2_edu"),
  Stroop_colores        = c("age1", "age2", "edu", "sex", "age1_edu", "age2_edu", "age1_sex"),
  Stroop_PC             = c("age1", "edu",  "age2", "age1_edu", "age2_edu"),
  SDMT                  = c("age1", "age2", "edu", "age2_edu"),
  Fonologica_F          = c("age1", "age2", "edu", "age2_edu"),
  Fonologica_A          = c("age1", "age2", "edu", "age1_edu", "age2_edu"),
  Fonologica_S          = c("age1", "age2", "edu", "age2_edu"),
  Fonologica_M          = c("age1", "age2", "edu", "age1_edu"),
  Semantica_animales    = c("age1", "age2", "edu", "age2_edu"),
  Semantica_frutas      = c("age1", "age2", "edu", "sex", "age1_edu", "age2_edu"),
  Semantica_profesiones = c("age1", "age2", "edu", "sex", "age2_edu", "age1_sex")
)

# ---- 2. Converter for predictor names to JAGS syntax ----
jags_safe_var <- function(var) {
  if (str_detect(var, "_")) {
    parts <- str_split(var, "_", simplify = TRUE)
    paste0(parts[1], "[i] * ", parts[2], "[i]")
  } else {
    paste0(var, "[i]")
  }
}

# ---- 3. Generate JAGS model file for Gaussian linear regression ----
generate_jags_model <- function(varname, predictors) {
  pred_exprs <- map_chr(predictors, jags_safe_var)
  terms      <- map2_chr(seq_along(pred_exprs) + 1, pred_exprs,
                         ~ paste0("beta[", .x, "] * ", .y))
  
  pred_line <- if (length(terms) > 0) {
    paste("mu[i] <- beta[1] +", paste(terms, collapse = " + "))
  } else {
    "mu[i] <- beta[1]"
  }
  
  model_string <- c(
    "model {",
    "  for (i in 1:n) {",
    "    y[i] ~ dnorm(mu[i], tau_y)",
    paste0("    ", pred_line),
    "  }",
    paste0("  for (j in 1:", length(predictors) + 1, ") {"),
    "    beta[j] ~ dunif(-1e6, 1e6)",
    "  }",
    "  tau_y <- 1 / (sigma * sigma)",
    "  sigma   ~ dunif(0, 10000)",
    "}"
  )
  
  model_file <- paste0("model_", varname, "_linear.txt")
  write_lines(model_string, model_file)
  model_file
}

# ---- 4. Function to run a Gaussian linear model ----
fit_variable_model <- function(varname, predictors) {
  message("Fitting model: ", varname)
  model_file <- generate_jags_model(varname, predictors)
  
  # Prepare data for JAGS
  jags_data <- list(
    y = df[[varname]],
    n = nrow(df)
  )
  base_vars <- unique(unlist(str_split(predictors, "_")))
  base_vars <- base_vars[base_vars != ""]
  
  for (v in base_vars) {
    if (v == "sex") {
      jags_data$sex <- as.numeric(df$Sexo == "Female")
    } else if (v == "edu") {
      jags_data$edu <- log(df$Escolaridad)
    } else {
      jags_data[[v]] <- df[[v]]
    }
  }
  
  # Initial values
  init_values <- function() {
    list(
      beta = rnorm(length(predictors) + 1, 0, 5),
      sigma = runif(1, 0, 1000)
    )
  }
  
  jm <- jags.model(
    file     = model_file,
    data     = jags_data,
    inits    = init_values,
    n.chains = 3,
    n.adapt  = 10000
  )
  
  update(jm, n.iter = 2000)
  
  coda.samples(
    jm,
    variable.names = c("beta", "sigma"),
    n.iter = 50000,
    thin   = 100
  )
}

# ---- 5. Run all Gaussian linear models and save outputs ----
results <- map(model_specs_lm, ~ fit_variable_model(.y, .x))
names(results) <- names(model_specs_lm)

# ---- Create JAGS output folder under models/ if missing ----
jags_out_lm <- file.path("models", "modelos_lineales_jags")
if (!dir.exists(jags_out_lm)) {
  dir.create(jags_out_lm, recursive = TRUE)
}

# Save each MCMC result as RDS in models/...
iwalk(results, ~ saveRDS(.x,
                         file = file.path(jags_out_lm,
                                          paste0(.y, "_linear_mcmc.rds"))))

# ---- 6. Model specifications for Poisson GLMs ----
model_specs_glm <- list(
  Stroop_palabras       = c("age1", "age2", "edu"),
  Stroop_colores        = c("age1", "age2", "edu", "sex", "edu_sex", "age2_sex", "age1_edu"),
  Stroop_PC             = c("age1", "age2", "edu", "age2_edu"),
  SDMT                  = c("age1", "age2", "edu", "age1_edu"),
  Fonologica_F          = c("age1", "age2", "edu", "age2_edu"),
  Fonologica_A          = c("age1", "age2", "edu", "age2_edu"),
  Fonologica_S          = c("age1", "age2", "edu", "age2_edu"),
  Fonologica_M          = c("age1", "age2", "edu", "age2_edu"),
  Semantica_animales    = c("age1", "age2", "edu", "age2_edu"),
  Semantica_frutas      = c("age1", "age2", "edu", "sex", "age1_edu", "age2_edu"),
  Semantica_profesiones = c("age1", "age2", "edu", "sex", "age1_sex")
)

# ---- 7. Generate JAGS model file for Poisson GLMs ----
generate_jags_model_glm <- function(varname, predictors) {
  pred_exprs <- map_chr(predictors, jags_safe_var)
  terms      <- map2_chr(seq_along(pred_exprs) + 1, pred_exprs,
                         ~ paste0("beta[", .x, "] * ", .y))
  
  pred_line <- if (length(terms) > 0) {
    paste("log(lambda[i]) <- beta[1] +", paste(terms, collapse = " + "))
  } else {
    "log(lambda[i]) <- beta[1]"
  }
  
  model_string <- c(
    "model {",
    "  for (i in 1:n) {",
    "    y[i] ~ dpois(lambda[i])",
    paste0("    ", pred_line),
    "  }",
    paste0("  for (j in 1:", length(predictors) + 1, ") {"),
    "    beta[j] ~ dnorm(0.0, 0.0001)",
    "  }",
    "}"
  )
  
  model_file <- paste0("model_", varname, "_glm.txt")
  write_lines(model_string, model_file)
  model_file
}

# ---- 8. Function to fit Poisson GLMs ----
fit_variable_model_glm <- function(varname, predictors) {
  message("Fitting GLM: ", varname)
  model_file <- generate_jags_model_glm(varname, predictors)
  
  jags_data <- list(
    y = df[[varname]],
    n = nrow(df)
  )
  
  base_vars <- unique(unlist(str_split(predictors, "_")))
  base_vars <- base_vars[base_vars != ""]
  for (v in base_vars) {
    if (v == "sex") {
      jags_data$sex <- as.numeric(df$Sexo == "Female")
    } else if (v == "edu") {
      jags_data$edu <- log(df$Escolaridad)
    } else {
      jags_data[[v]] <- df[[v]]
    }
  }
  
  init_values <- function() {
    list(beta = rnorm(length(predictors) + 1, 0, 0.01))
  }
  
  jm <- jags.model(
    file     = model_file,
    data     = jags_data,
    inits    = init_values,
    n.chains = 3,
    n.adapt  = 10000
  )
  
  update(jm, n.iter = 2000)
  
  coda.samples(
    jm,
    variable.names = "beta",
    n.iter         = 50000,
    thin           = 100
  )
}

# ---- 9. Run all Poisson GLMs and save outputs ----
glm_results <- map(model_specs_glm, ~ fit_variable_model_glm(.y, .x))
names(glm_results) <- names(model_specs_glm)

# ---- Create JAGS output folder under models/ if missing ----
jags_out_glm <- file.path("models", "modelos_glm_jags")
if (!dir.exists(jags_out_glm)) {
  dir.create(jags_out_glm, recursive = TRUE)
}

iwalk(glm_results, ~ saveRDS(.x,
                             file = file.path(jags_out_glm,
                                              paste0(.y, "_glm_mcmc.rds"))))
