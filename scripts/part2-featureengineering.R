# Part2 code

# Create necessary directories
if (!dir.exists("data/processed")) {
  dir.create("data/processed", recursive = TRUE)
}
if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}

# Load required packages
packages <- c("depmixS4", "dplyr", "lubridate", "ggplot2", "stats", "reshape2")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(depmixS4)  # For HMM models
library(dplyr)     # For data manipulation
library(lubridate) # For date/time handling
library(ggplot2)   # For plotting
library(stats)     # For PCA
library(reshape2)  # For data reshaping


if (!requireNamespace("ggbiplot", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("vqv/ggbiplot")
}
suppressWarnings(library(ggbiplot))  # For PCA visualization

# Check if standardized data exists, if not run Part I feature scaling
if (!file.exists("data/processed/TermProjectData_Standardized.csv")) {
  cat("Standardized data not found. Running feature scaling first...\n")
  
  # Load raw data
  data_raw <- read.csv("data/TermProjectData.txt", header = TRUE, stringsAsFactors = FALSE)
  
  # Convert date and time to DateTime format
  data_raw$DateTime <- dmy_hms(paste(data_raw$Date, data_raw$Time))
  if (all(is.na(data_raw$DateTime))) {
    # Try alternate format
    data_raw$DateTime <- as.POSIXct(paste(data_raw$Date, data_raw$Time), format="%d/%m/%Y %H:%M:%S")
  }
  
  # Convert to epoch time
  datetimes <- as.integer(data_raw$DateTime)
  data_raw$Epoch <- datetimes
  
  # Save original Date and Time
  data_raw$DateOriginal <- data_raw$Date
  data_raw$TimeOriginal <- data_raw$Time
  
  # Function to handle missing values using linear interpolation
  replace_na_with_lerped <- function(ind, dep) {
    if (!any(is.na(dep))) {
      return(NULL)
    }
    new_values <- approx(ind, dep, ind[which(is.na(dep))])
    return(new_values)
  }
  
  # Create a copy of the data for processing
  cdf <- data_raw
  
  # Interpolate missing values for each variable
  variables <- c("Global_active_power", "Global_reactive_power", "Voltage", 
                "Global_intensity", "Sub_metering_1", "Sub_metering_2", "Sub_metering_3")
  
  for (var in variables) {
    nv <- replace_na_with_lerped(data_raw$Epoch, data_raw[[var]])
    if (!is.null(nv)) {
      cdf[[var]][which(is.na(data_raw[[var]]))] <- nv$y
    }
  }
  
  # Standardize the data (z-score)
  sdf <- cdf
  for (var in variables) {
    sdf[[var]] <- scale(cdf[[var]], center = TRUE, scale = TRUE)
  }
  
  # Check for any remaining NA or infinite values
  any_na <- sapply(sdf[, variables], function(x) any(is.na(x) | is.infinite(x)))
  
  # Replace any remaining problematic values with means
  if(any(any_na)) {
    cat("Found NA or infinite values after standardization. Replacing them with mean values.\n")
    for(var in names(any_na)[any_na]) {
      var_mean <- mean(sdf[[var]], na.rm = TRUE)
      sdf[[var]][is.na(sdf[[var]]) | is.infinite(sdf[[var]])] <- var_mean
      cat(paste("Replaced values in", var, "with mean:", var_mean, "\n"))
    }
  }
  
  # Add back Date and Time columns
  sdf$Date <- data_raw$DateOriginal
  sdf$Time <- data_raw$TimeOriginal
  
  # Save standardized data
  write.csv(sdf, "data/processed/TermProjectData_Standardized.csv", row.names = FALSE)
  cat("Standardized data saved to data/processed/TermProjectData_Standardized.csv\n")
} else {
  cat("Using existing standardized data...\n")
}

# Step 1: Load the standardized data
cat("\nStep 1: Loading standardized data...\n")
sdf <- read.csv("data/processed/TermProjectData_Standardized.csv", stringsAsFactors = FALSE)

# Add DateTime and Weekday columns
sdf$DateTime <- as.POSIXct(paste(sdf$Date, sdf$Time), format="%d/%m/%Y %H:%M:%S")
if (all(is.na(sdf$DateTime))) {
  # Try alternate format
  sdf$DateTime <- as.POSIXct(paste(sdf$Date, sdf$Time), format="%Y-%m-%d %H:%M:%S")
}
sdf$Weekday <- weekdays(sdf$DateTime)

# Step 2: Perform PCA for feature selection
cat("Step 2: Performing PCA analysis...\n")
variables <- c("Global_active_power", "Global_reactive_power", "Voltage", 
              "Global_intensity", "Sub_metering_1", "Sub_metering_2", "Sub_metering_3")
features <- sdf[, variables]

# Check for NA values before PCA
na_count <- colSums(is.na(features))
if(sum(na_count) > 0) {
  cat("Found NA values in the data. Imputing with column means...\n")
  for(col in names(features)) {
    if(na_count[col] > 0) {
      features[[col]][is.na(features[[col]])] <- mean(features[[col]], na.rm = TRUE)
    }
  }
}

# Run PCA
pca_result <- prcomp(features, center = TRUE, scale. = TRUE)
summary_pca <- summary(pca_result)

# Print PCA summary
cat("PCA Summary:\n")
print(summary_pca)

# Save PCA results to file
capture.output(print(summary_pca), file = "output/pca_summary.txt")

# Visualize PCA results - Scree plot
var_explained <- summary_pca$importance[2,] * 100
cum_var_explained <- summary_pca$importance[3,] * 100

scree_data <- data.frame(
  Component = factor(1:length(var_explained), levels = 1:length(var_explained)),
  Variance = var_explained,
  Cumulative = cum_var_explained
)

scree_plot <- ggplot(scree_data, aes(x = Component, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(y = Cumulative, group = 1), color = "red", linewidth = 1) +
  geom_point(aes(y = Cumulative), color = "red", size = 3) +
  labs(title = "PCA Scree Plot", 
       x = "Principal Component", 
       y = "Percent of Variance Explained") +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Cumulative Percent")) +
  theme_minimal()

print(scree_plot)
ggsave("output/pca_scree_plot.png", scree_plot, width = 10, height = 6)

# Try biplot with tryCatch to handle potential errors
tryCatch({
  biplot <- ggbiplot(pca_result, labels = rownames(features), groups = NULL) +
    labs(title = "PCA Biplot of Electricity Consumption Variables")
  
  print(biplot)
  ggsave("output/pca_biplot.png", biplot, width = 10, height = 8)
}, error = function(e) {
  cat("Error generating biplot:", e$message, "\n")
  cat("Skipping biplot generation...\n")
})

# Loading plot
loadings <- as.data.frame(pca_result$rotation)
loadings$Variable <- rownames(loadings)
loadings_long <- reshape2::melt(loadings, id.vars = "Variable", 
                               variable.name = "Component", 
                               value.name = "Loading")

# Focus on first 3 components
loadings_long_subset <- loadings_long[loadings_long$Component %in% c("PC1", "PC2", "PC3"), ]

loading_plot <- ggplot(loadings_long_subset, aes(x = Variable, y = Loading, fill = Component)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "PCA Loading Plot", 
       x = "Variable", 
       y = "Loading") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(loading_plot)
ggsave("output/pca_loading_plot.png", loading_plot, width = 10, height = 6)

# Based on PCA results, select variables for HMM
# Select the top 3 variables that contribute most to PC1
loadings_pc1 <- abs(pca_result$rotation[, "PC1"])
selected_variables <- names(loadings_pc1)[order(loadings_pc1, decreasing = TRUE)[1:3]]

cat("Selected variables for HMM based on PCA:", paste(selected_variables, collapse=", "), "\n")

# Save variable selection information
cat(paste("Selected variables for HMM based on PCA:", paste(selected_variables, collapse=", ")), 
    file = "output/selected_variables.txt")

# Step 3: Select time window and partition data
cat("Step 3: Selecting time window and partitioning data...\n")
# Choose a weekday and time window
selected_weekday <- "Monday"  # You can change this to any day
start_time <- "18:00:00"
end_time <- "22:00:00"

# Function to check if a time falls within the selected window
is_within_timewindow <- function(time_str, start_time, end_time) {
  # 日本語の曜日名を含む文字列をParseするために、単純な比較に切り替える
  hour <- as.numeric(substr(time_str, 1, 2))
  start_hour <- as.numeric(substr(start_time, 1, 2))
  end_hour <- as.numeric(substr(end_time, 1, 2))
  
  return(hour >= start_hour && hour <= end_hour)
}

# Filter data for the selected weekday and time window
selected_data <- sdf[sdf$Weekday == selected_weekday & 
                     sapply(sdf$Time, is_within_timewindow, start_time, end_time), ]

# Extract year from DateTime for partitioning
selected_data$Year <- year(selected_data$DateTime)

# Check the number of data points for each year
year_count <- table(selected_data$Year)
cat("Data points by year:\n")
print(year_count)

# Partition data into training (first 3 years) and testing (4th year)
all_years <- sort(unique(selected_data$Year))
if(length(all_years) >= 4) {
  training_years <- all_years[1:3]
  testing_year <- all_years[4]
  
  train_data <- selected_data[selected_data$Year %in% training_years, ]
  test_data <- selected_data[selected_data$Year == testing_year, ]
  
  # Print the dimensions of the training and testing datasets
  cat("Training data dimensions:", dim(train_data), "\n")
  cat("Testing data dimensions:", dim(test_data), "\n")
  
  # Extract the selected variables for HMM
  train_df <- data.frame(
    Global_intensity = train_data$Global_intensity,
    Global_active_power = train_data$Global_active_power,
    Sub_metering_3 = train_data$Sub_metering_3
  )
  
  test_df <- data.frame(
    Global_intensity = test_data$Global_intensity,
    Global_active_power = test_data$Global_active_power,
    Sub_metering_3 = test_data$Sub_metering_3
  )
  

  cat("Training features dimensions:", dim(train_df), "\n")
  cat("Test features dimensions:", dim(test_df), "\n")
  cat("Selected variables:", paste(colnames(train_df), collapse=", "), "\n")
  
  # Check for NA values in features
  na_in_train <- any(is.na(train_df))
  na_in_test <- any(is.na(test_df))
  
  if(na_in_train || na_in_test) {
    cat("Warning: Found NA values in feature data. Imputing with means...\n")
    
    # Impute NA values in training data
    if(na_in_train) {
      for(col in names(train_df)) {
        col_na <- is.na(train_df[[col]])
        if(any(col_na)) {
          col_mean <- mean(train_df[[col]], na.rm = TRUE)
          train_df[[col]][col_na] <- col_mean
          cat(paste("Imputed", sum(col_na), "NA values in train", col, "with mean", col_mean, "\n"))
        }
      }
    }
    
    # Impute NA values in test data
    if(na_in_test) {
      for(col in names(test_df)) {
        col_na <- is.na(test_df[[col]])
        if(any(col_na)) {
          col_mean <- mean(test_df[[col]], na.rm = TRUE)
          test_df[[col]][col_na] <- col_mean
          cat(paste("Imputed", sum(col_na), "NA values in test", col, "with mean", col_mean, "\n"))
        }
      }
    }
  }
} else {
  stop("Not enough years in the dataset for proper partitioning.")
}

# Save information about time window and partitioning
cat(paste("Selected time window:", selected_weekday, "from", start_time, "to", end_time), 
    file = "output/time_window_info.txt")
cat(paste("\nTraining years:", paste(training_years, collapse=", ")), 
    file = "output/time_window_info.txt", append = TRUE)
cat(paste("\nTesting year:", testing_year), 
    file = "output/time_window_info.txt", append = TRUE)
cat(paste("\nTraining data dimensions:", paste(dim(train_data), collapse=" x ")), 
    file = "output/time_window_info.txt", append = TRUE)
cat(paste("\nTesting data dimensions:", paste(dim(test_data), collapse=" x ")), 
    file = "output/time_window_info.txt", append = TRUE)


train_multiple_hmms <- function(train_df, test_df, variable, state_range = c(4, 6, 8, 10, 12)) {
  cat("\nStep 4: Training multiple HMM models with different numbers of states...\n")
  
  
  results_df <- data.frame(
    States = integer(),
    LogLikelihood = numeric(),
    NormalizedTrainLL = numeric(),
    BIC = numeric(),
    TestLL = numeric(),
    NormalizedTestLL = numeric()
  )
  
 
  for(num_states in state_range) {
    cat("\nTraining HMM with", num_states, "states for variable:", variable, "\n")
    
    tryCatch({
      
      set.seed(123)
      hmm_mod <- depmix(train_df[[variable]] ~ 1, nstates = num_states, data = train_df)
      
      
      cat("Fitting model...\n")
      hmm_fitted <- fit(hmm_mod, verbose = TRUE)
      
      
      train_ll <- logLik(hmm_fitted)
      train_norm_ll <- train_ll / nrow(train_df)
      bic <- -2 * train_ll + npar(hmm_fitted) * log(nrow(train_df))
      
      #
      test_mod <- depmix(test_df[[variable]] ~ 1, nstates = num_states, data = test_df)
      test_mod <- setpars(test_mod, getpars(hmm_fitted))
      fb <- forwardbackward(test_mod)
      test_ll <- fb$logLike
      test_norm_ll <- test_ll / nrow(test_df)
      

      cat("  Train log-likelihood:", train_ll, "\n")
      cat("  Normalized train log-likelihood:", train_norm_ll, "\n")
      cat("  BIC:", bic, "\n")
      cat("  Test log-likelihood:", test_ll, "\n")
      cat("  Normalized test log-likelihood:", test_norm_ll, "\n")
      

      results_df <- rbind(results_df, data.frame(
        States = num_states,
        LogLikelihood = train_ll,
        NormalizedTrainLL = train_norm_ll,
        BIC = bic,
        TestLL = test_ll,
        NormalizedTestLL = test_norm_ll
      ))
    }, error = function(e) {
      cat("Error in training HMM with", num_states, "states:", e$message, "\n")
    })
  }
  

  return(results_df)
}


results <- train_multiple_hmms(train_df, test_df, "Global_intensity", c(4, 6, 8, 10, 12))


write.csv(results, "output/hmm_model_comparison.csv", row.names = FALSE)


if(nrow(results) > 1) {
  ll_plot <- ggplot(results, aes(x = States, y = LogLikelihood)) +
    geom_point() +
    geom_line() +
    labs(title = "Log-Likelihood vs. Number of States",
         x = "Number of States",
         y = "Log-Likelihood") +
    theme_minimal()
  
  print(ll_plot)
  ggsave("output/log_likelihood_plot.png", ll_plot, width = 8, height = 6)
  

  bic_plot <- ggplot(results, aes(x = States, y = BIC)) +
    geom_point() +
    geom_line() +
    labs(title = "BIC vs. Number of States",
         x = "Number of States",
         y = "BIC") +
    theme_minimal()
  
  print(bic_plot)
  ggsave("output/bic_plot.png", bic_plot, width = 8, height = 6)
  

  norm_ll_data <- data.frame(
    States = rep(results$States, 2),
    Dataset = c(rep("Train", nrow(results)), rep("Test", nrow(results))),
    NormalizedLL = c(results$NormalizedTrainLL, results$NormalizedTestLL)
  )
  
  norm_ll_plot <- ggplot(norm_ll_data, aes(x = States, y = NormalizedLL, color = Dataset)) +
    geom_point() +
    geom_line() +
    labs(title = "Normalized Log-Likelihood vs. Number of States",
         x = "Number of States",
         y = "Normalized Log-Likelihood") +
    theme_minimal()
  
  print(norm_ll_plot)
  ggsave("output/normalized_ll_plot.png", norm_ll_plot, width = 8, height = 6)
}


if(nrow(results) > 0) {

  best_model_idx <- which.min(results$BIC)
  best_states <- results$States[best_model_idx]
  
  cat("\nBest model based on BIC: HMM with", best_states, "states\n")
  cat("BIC:", results$BIC[best_model_idx], "\n")
  cat("Log-likelihood:", results$LogLikelihood[best_model_idx], "\n")
  cat("Normalized train log-likelihood:", results$NormalizedTrainLL[best_model_idx], "\n")
  cat("Normalized test log-likelihood:", results$NormalizedTestLL[best_model_idx], "\n")
  

  cat(paste("Best model based on BIC: HMM with", best_states, "states"), 
      file = "output/best_model_info.txt")
  cat(paste("\nBIC:", results$BIC[best_model_idx]), 
      file = "output/best_model_info.txt", append = TRUE)
  cat(paste("\nLog-likelihood:", results$LogLikelihood[best_model_idx]), 
      file = "output/best_model_info.txt", append = TRUE)
  cat(paste("\nNormalized train log-likelihood:", results$NormalizedTrainLL[best_model_idx]), 
      file = "output/best_model_info.txt", append = TRUE)
  cat(paste("\nNormalized test log-likelihood:", results$NormalizedTestLL[best_model_idx]), 
      file = "output/best_model_info.txt", append = TRUE)
  

  cat("\nRetraining best model for anomaly detection...\n")
  set.seed(123)
  best_mod <- depmix(train_df[["Global_intensity"]] ~ 1, nstates = best_states, data = train_df)
  best_fitted <- fit(best_mod, verbose = TRUE)
  

  best_train_ll <- logLik(best_fitted)
  best_train_norm_ll <- best_train_ll / nrow(train_df)
  

  n_subsets <- 10
  subset_size <- ceiling(nrow(test_df) / n_subsets)
  

  subset_results <- data.frame(
    Subset = integer(),
    Size = integer(),
    LogLikelihood = numeric(),
    NormalizedLL = numeric(),
    Deviation = numeric()
  )
  
  for(i in 1:n_subsets) {
    start_idx <- (i-1) * subset_size + 1
    end_idx <- min(i * subset_size, nrow(test_df))
    
    if (start_idx <= nrow(test_df)) {

      subset_df <- test_df[start_idx:end_idx, , drop = FALSE]
      

      subset_mod <- depmix(subset_df[["Global_intensity"]] ~ 1, nstates = best_states, data = subset_df)
      subset_mod <- setpars(subset_mod, getpars(best_fitted))

      subset_fb <- forwardbackward(subset_mod)
      subset_ll <- subset_fb$logLike
      subset_norm_ll <- subset_ll / nrow(subset_df)
      

      subset_results <- rbind(subset_results, data.frame(
        Subset = i,
        Size = nrow(subset_df),
        LogLikelihood = subset_ll,
        NormalizedLL = subset_norm_ll,
        Deviation = abs(subset_norm_ll - best_train_norm_ll)
      ))
    }
  }
  

  cat("Subset log-likelihood results:\n")
  print(subset_results)
  write.csv(subset_results, "output/subset_results.csv", row.names = FALSE)
  

  subset_ll_plot <- ggplot(subset_results, aes(x = Subset, y = NormalizedLL)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_hline(yintercept = best_train_norm_ll, linetype = "dashed", color = "red") +
    labs(title = paste("Normalized Log-Likelihood for Test Subsets - Global_intensity (", best_states, "states)"),
         subtitle = paste("Red line: Training normalized log-likelihood =", 
                         round(best_train_norm_ll, 4)),
         x = "Subset",
         y = "Normalized Log-Likelihood") +
    theme_minimal()
  
  print(subset_ll_plot)
  ggsave("output/subset_ll_plot.png", subset_ll_plot, width = 10, height = 6)
  

  deviation_plot <- ggplot(subset_results, aes(x = Subset, y = Deviation)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = paste("Deviation from Training Log-Likelihood - Global_intensity (", best_states, "states)"),
         x = "Subset",
         y = "Absolute Deviation") +
    theme_minimal()
  
  print(deviation_plot)
  ggsave("output/deviation_plot.png", deviation_plot, width = 10, height = 6)
  

  threshold <- max(subset_results$Deviation)
  cat("\nThreshold for normal behavior (maximum deviation):", threshold, "\n")
  

  cat(paste("Threshold for normal behavior (maximum deviation):", threshold), 
      file = "output/threshold_info.txt")
  cat(paste("\nTraining normalized log-likelihood:", best_train_norm_ll), 
      file = "output/threshold_info.txt", append = TRUE)
  cat("\nSubset deviations:", file = "output/threshold_info.txt", append = TRUE)
  cat(paste("\n", paste(subset_results$Deviation, collapse = ", ")), 
      file = "output/threshold_info.txt", append = TRUE)
} else {
  cat("\nNo successful models were trained. Please check your data and parameters.\n")
}


cat("\nStep 5: Brief comparison of HMM performance on other variables...\n")


all_vars <- colnames(train_df)
var_results <- data.frame(
  Variable = character(),
  States = integer(),
  LogLikelihood = numeric(),
  NormalizedTrainLL = numeric(),
  BIC = numeric(),
  TestLL = numeric(),
  NormalizedTestLL = numeric(),
  Notes = character()
)

for(var in all_vars) {
  cat("Training model for variable:", var, "\n")
  

  tryCatch({

    set.seed(123)
    var_mod <- depmix(train_df[[var]] ~ 1, nstates = 4, data = train_df)
    var_fitted <- fit(var_mod, verbose = FALSE, emcontrol = em.control(maxit = 200))
    

    train_ll <- logLik(var_fitted)
    train_norm_ll <- train_ll / nrow(train_df)
    bic <- -2 * train_ll + npar(var_fitted) * log(nrow(train_df))
    

    test_mod <- depmix(test_df[[var]] ~ 1, nstates = 4, data = test_df)
    test_mod <- setpars(test_mod, getpars(var_fitted))
    fb <- forwardbackward(test_mod)
    test_ll <- fb$logLike
    test_norm_ll <- test_ll / nrow(test_df)
    

    notes <- ""
    if(train_ll > 0) {
      notes <- "WARNING: Positive log-likelihood indicates model issues"
    }
    

    cat("  Train log-likelihood:", train_ll, "\n")
    cat("  Normalized train log-likelihood:", train_norm_ll, "\n")
    cat("  BIC:", bic, "\n")
    cat("  Test log-likelihood:", test_ll, "\n")
    cat("  Normalized test log-likelihood:", test_norm_ll, "\n")
    if(notes != "") {
      cat("  ", notes, "\n")
    }
    

    var_results <- rbind(var_results, data.frame(
      Variable = var,
      States = 4,
      LogLikelihood = train_ll,
      NormalizedTrainLL = train_norm_ll,
      BIC = bic,
      TestLL = test_ll,
      NormalizedTestLL = test_norm_ll,
      Notes = notes
    ))
  }, error = function(e) {
    cat("Error in training model for", var, ":", e$message, "\n")
    
    var_results <- rbind(var_results, data.frame(
      Variable = var,
      States = 4,
      LogLikelihood = NA,
      NormalizedTrainLL = NA,
      BIC = NA,
      TestLL = NA,
      NormalizedTestLL = NA,
      Notes = paste("ERROR:", e$message)
    ))
  })
}


valid_results <- var_results[var_results$LogLikelihood < 0 & !is.na(var_results$LogLikelihood), ]


cat("\nComparison of HMM performance on different variables (valid models only):\n")
if(nrow(valid_results) > 0) {
  print(valid_results[, c("Variable", "States", "LogLikelihood", "BIC", "NormalizedTestLL")])
  write.csv(valid_results, "output/variable_comparison.csv", row.names = FALSE)
  

  if(nrow(valid_results) < nrow(var_results)) {
    cat("\nWARNING: Some variables showed unusual behavior and were excluded from comparison.\n")
    cat("Full results including problematic variables:\n")
    print(var_results[, c("Variable", "States", "LogLikelihood", "BIC", "Notes")])
    write.csv(var_results, "output/all_variables_comparison.csv", row.names = FALSE)
  }
} else {
  cat("No valid models found. All variables showed unusual behavior.\n")
  print(var_results[, c("Variable", "States", "LogLikelihood", "Notes")])
  write.csv(var_results, "output/all_variables_comparison.csv", row.names = FALSE)
}

cat("\nTerm Project Part II analysis completed.\n")
cat("Check the 'output' directory for all generated results and figures.\n")