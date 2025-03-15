
# part3-HMMTrainingandTesting.R
#
# CMPT 318 Term Project - Part 3
# HMM Training & Testing (Multivariate)
#

r <- getOption("repos")
if (is.null(r) || r["CRAN"] == "@CRAN@") {
  r["CRAN"] <- "https://cran.rstudio.com/"
  options(repos = r)
}

# 1) Directory structure and package loading
OUTPUT_DIR <- "output"
PART3_DIR  <- file.path(OUTPUT_DIR, "part3")   # Outputs from Part 3
MODELS_DIR <- file.path(OUTPUT_DIR, "models")  # Where we store trained models

# Create directories if they don't exist
dir_list <- c(OUTPUT_DIR, PART3_DIR, MODELS_DIR)
for (d in dir_list) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Packages required
packages <- c("depmixS4", "dplyr", "lubridate", "ggplot2")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)  # uses the mirror we set above
  }
  library(pkg, character.only = TRUE)
}

# 2) Load your preprocessed, scaled dataset (from Task 1)
scaled_data_path <- "data/processed/TermProjectData_Standardized.csv"
if (!file.exists(scaled_data_path)) {
  stop("Scaled data not found at ", scaled_data_path, 
       ". Please run Task 1 first to generate the standardized CSV.")
}

sdf <- read.csv(scaled_data_path, stringsAsFactors = FALSE)
cat("Loaded scaled dataset with", nrow(sdf), "rows and", ncol(sdf), "columns.\n")

# 2a) Try parsing 'DateTime' in two stages:


if (!"DateTime" %in% colnames(sdf)) {
  stop("No 'DateTime' column found in the CSV. Check your Task 1 output.")
}

parsed1 <- as.POSIXct(sdf$DateTime, format="%Y-%m-%d %H:%M:%S", tz="UTC")
needs_fallback <- is.na(parsed1)

# For those that were NA, try the second format
parsed2 <- as.POSIXct(sdf$DateTime[needs_fallback], format="%d/%m/%Y %H:%M:%S", tz="UTC")

# Combine
sdf$DateTime <- parsed1
sdf$DateTime[needs_fallback] <- parsed2

# Check how many remain NA after both attempts
bad_count <- sum(is.na(sdf$DateTime))
if (bad_count > 0) {
  warning("Could not parse DateTime for ", bad_count, " rows. Removing those rows...\n")
  sdf <- sdf[!is.na(sdf$DateTime), ]
  cat("New data size after removing unparsed rows:", nrow(sdf), "\n")
}

# Now create a Weekday column
sdf$Weekday <- weekdays(sdf$DateTime)

# 3) Choose a weekday/weekend and a time window
SELECTED_WEEKDAY <- "Monday"   
START_TIME       <- "18:00:00"
END_TIME         <- "22:00:00"
# Adjust these as needed.

# 4) Partition function:

partition_data <- function(df, weekday, start_t, end_t) {
  # Filter for chosen weekday
  df_filtered <- df[df$Weekday == weekday, ]
  
  # Filter by time range
  df_filtered$TimeHM <- format(df_filtered$DateTime, "%H:%M:%S")
  
  df_filtered <- df_filtered %>%
    dplyr::filter(TimeHM >= start_t & TimeHM <= end_t)
  
  # Extract numeric years from DateTime
  df_filtered$Year <- lubridate::year(df_filtered$DateTime)
  all_years <- sort(unique(df_filtered$Year))
  cat("Data covers years:", paste(all_years, collapse=", "), "\n")
  
  # Need at least 4 distinct years for 3-year train + 1-year test
  if (length(all_years) < 4) {
    stop("Not enough distinct years in the dataset for a 3-year train + 1-year test partition.")
  }
  
  train_years <- all_years[1:3]
  test_year   <- all_years[4]
  train_data <- df_filtered[df_filtered$Year %in% train_years, ]
  test_data  <- df_filtered[df_filtered$Year == test_year, ]
  
  cat("Training data: years", paste(train_years, collapse=", "), 
      "→", nrow(train_data), "rows\n")
  cat("Test data: year", test_year, "→", nrow(test_data), "rows\n")
  
  list(train_data = train_data, test_data = test_data)
}

partitioned <- partition_data(sdf, SELECTED_WEEKDAY, START_TIME, END_TIME)
train_data  <- partitioned$train_data
test_data   <- partitioned$test_data

# 5) Choose variables for a MULTIVARIATE HMM
#    Typically, you'd pick these based on PCA from Task 2.
selected_vars <- c("Global_active_power", "Voltage", "Sub_metering_3")

# Sanity check
if (any(!selected_vars %in% names(train_data))) {
  stop("One or more of the selected variables are missing from the data.\n",
       "Check 'selected_vars' and columns in your CSV.")
}

# Build train/test data frames
train_df <- train_data[, selected_vars, drop=FALSE]
test_df  <- test_data[, selected_vars, drop=FALSE]

cat("\nSelected variables:", paste(selected_vars, collapse=", "), "\n")
cat("Train set dimension:", dim(train_df)[1], "rows x", dim(train_df)[2], "cols\n")
cat("Test set dimension:", dim(test_df)[1], "rows x", dim(test_df)[2], "cols\n")

# 6) Define a helper function to train a multivariate HMM using depmixS4
library(depmixS4)

train_hmm <- function(train_df, nstates) {
  # Each column is a separate response, ~ 1 means intercept-only
  # Adjust 'family' if your columns are not all continuous (gaussian).
  hmm_spec <- depmix(
    list(
      train_df[[1]] ~ 1,
      train_df[[2]] ~ 1,
      train_df[[3]] ~ 1
      # Add more lines if you have more selected_vars
    ),
    data    = train_df,
    nstates = nstates,
    family  = list(gaussian(), gaussian(), gaussian())
  )
  
  set.seed(123)  # reproducibility
  hmm_fit <- fit(hmm_spec, verbose=FALSE, emcontrol=em.control(maxit=200))
  return(hmm_fit)
}

# 7) Train multiple models for different numbers of states
state_candidates <- c(4, 6, 8, 10, 12)
results_df <- data.frame(
  states         = integer(),
  logLik_train   = numeric(),
  train_norm_ll  = numeric(),  # NEW
  BIC            = numeric(),
  logLik_test    = numeric(),
  test_norm_ll   = numeric(),  # NEW
  stringsAsFactors = FALSE
)

hmm_list <- list()

for (nstates in state_candidates) {
  cat("\n--------------------------------------------------------\n")
  cat("Training HMM with", nstates, "states...\n")
  
  # (A) Train on the training set
  hmm_fit <- train_hmm(train_df, nstates)
  ll_train <- logLik(hmm_fit)
  bic_val  <- -2 * ll_train + npar(hmm_fit) * log(nrow(train_df))
  
  # Normalized log-likelihood on train
  train_norm_ll <- as.numeric(ll_train) / nrow(train_df)
  
  # (B) Evaluate on the test set
  test_spec <- depmix(
    list(
      test_df[[1]] ~ 1,
      test_df[[2]] ~ 1,
      test_df[[3]] ~ 1
    ),
    data    = test_df,
    nstates = nstates,
    family  = list(gaussian(), gaussian(), gaussian())
  )
  test_spec <- setpars(test_spec, getpars(hmm_fit))
  fb_test   <- forwardbackward(test_spec)
  ll_test   <- fb_test$logLike
  
  # Normalized log-likelihood on test
  test_norm_ll <- as.numeric(ll_test) / nrow(test_df)
  
  # (C) Save results
  results_df <- rbind(
    results_df,
    data.frame(
      states        = nstates,
      logLik_train  = as.numeric(ll_train),
      train_norm_ll = train_norm_ll,
      BIC           = as.numeric(bic_val),
      logLik_test   = as.numeric(ll_test),
      test_norm_ll  = test_norm_ll,
      stringsAsFactors = FALSE
    )
  )
  hmm_list[[as.character(nstates)]] <- hmm_fit
  
  cat("Train logLik     =", as.numeric(ll_train), "\n")
  cat("Train normalized =", train_norm_ll, "\n")
  cat("Train BIC        =", as.numeric(bic_val), "\n")
  cat("Test logLik      =", as.numeric(ll_test), "\n")
  cat("Test normalized  =", test_norm_ll, "\n")
}

# 8) Model selection:

valid_rows <- results_df[results_df$logLik_train < 0, ]
if (nrow(valid_rows) > 0) {
  best_idx <- which.min(valid_rows$BIC)
  best_states <- valid_rows$states[best_idx]
  cat("\nBest model among negative-train-LL models by BIC has", best_states, "states.\n")
} else {
  cat("\nNo model had negative train log-likelihood. Selecting minimal BIC overall.\n")
  best_idx <- which.min(results_df$BIC)
  best_states <- results_df$states[best_idx]
}

best_model <- hmm_list[[as.character(best_states)]]

# Print final info
cat("\n=== BEST MODEL SELECTION ===\n")
cat("Best model:", best_states, "states\n")
cat("Train LL:  ", results_df$logLik_train[results_df$states == best_states], "\n")
cat("Train BIC: ", results_df$BIC[results_df$states == best_states], "\n")
cat("Test LL:   ", results_df$logLik_test[results_df$states == best_states], "\n")
cat("Train Norm LL:", results_df$train_norm_ll[results_df$states == best_states], "\n")
cat("Test Norm LL: ", results_df$test_norm_ll[results_df$states == best_states], "\n")

# Save outputs
csv_path <- file.path(PART3_DIR, "hmm_training_results.csv")
write.csv(results_df, csv_path, row.names=FALSE)
cat("Saved model comparison to:", csv_path, "\n")

rds_best_model <- file.path(MODELS_DIR, "best_multivariate_hmm.rds")
saveRDS(best_model, rds_best_model)
cat("Saved best model as:", rds_best_model, "\n")

rds_best_states <- file.path(MODELS_DIR, "best_multivariate_states.rds")
saveRDS(best_states, file.path(MODELS_DIR, "best_states.rds"))
cat("Saved best model's state count as:", rds_best_states, "\n")

cat("\n=== Part 3: HMM Training & Testing Complete ===\n")
