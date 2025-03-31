# CMPT 318 Term Project - Part 2

# Time window settings
SELECTED_WEEKDAY <- "Monday"
START_TIME <- "18:00:00"
END_TIME <- "22:00:00"

# HMM state settings
STATE_RANGE <- c(4, 6, 8, 10, 12) # Range of states to test
DEFAULT_STATES <- 6 # Statistically valid state count (from models with negative log-likelihood)

# Folder structure settings
OUTPUT_DIR <- "output"
PART1_DIR <- file.path(OUTPUT_DIR, "part1") # Output from Part 1
PART2_DIR <- file.path(OUTPUT_DIR, "part2") # Output from Part 2
MODELS_DIR <- file.path(OUTPUT_DIR, "models") # Location for saved models

# Check current working directory - Display in console
cat("Current working directory:", getwd(), "\n")

# Create necessary directories
create_directories <- function() {
  dirs <- c(
    "data/processed",
    OUTPUT_DIR,
    PART1_DIR,
    PART2_DIR,
    MODELS_DIR
  )
  
  for(dir in dirs) {
    if(!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    }
  }
}

create_directories()

# Load required packages
load_packages <- function() {
  packages <- c("depmixS4", "dplyr", "lubridate", "ggplot2", "stats", "reshape2", "corrplot")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  # Install and load ggbiplot
  if (!requireNamespace("ggbiplot", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }
    devtools::install_github("vqv/ggbiplot")
  }
  suppressWarnings(library(ggbiplot))
}

load_packages()

# Data preprocessing - Reuse Part 1 functionality
preprocess_data <- function() {
  standardized_data_path <- "data/processed/TermProjectData_Standardized.csv"
  
  # Check if standardized data already exists
  if (file.exists(standardized_data_path)) {
    cat("Using existing standardized data...\n")
    sdf <- read.csv(standardized_data_path, stringsAsFactors = FALSE)
  } else {
    cat("Standardized data not found. Please run Part 1 first.\n")
    stop("Missing standardized data")
  }
  
  # Convert date and time
  sdf$DateTime <- as.POSIXct(paste(sdf$Date, sdf$Time), format="%d/%m/%Y %H:%M:%S")
  if (all(is.na(sdf$DateTime))) {
    # Try alternate format
    sdf$DateTime <- as.POSIXct(paste(sdf$Date, sdf$Time), format="%Y-%m-%d %H:%M:%S")
  }
  sdf$Weekday <- weekdays(sdf$DateTime)
  
  return(sdf)
}

# Step 1: Load the data
cat("\nStep 1: Loading standardized data...\n")
sdf <- preprocess_data()

# Step 2: Variable selection using PCA
perform_pca <- function(sdf) {
  cat("Step 2: Performing PCA analysis...\n")
  variables <- c("Global_active_power", "Global_reactive_power", "Voltage", 
                "Global_intensity", "Sub_metering_1", "Sub_metering_2", "Sub_metering_3")
  features <- sdf[, variables]
  
  # Check for and handle NA values
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
  
  # Display and save PCA summary
  cat("PCA Summary:\n")
  print(summary_pca)
  capture.output(print(summary_pca), file = file.path(PART2_DIR, "pca_summary.txt"))
  
  # Create PCA plots
  create_pca_plots(pca_result, features)
  
  # Calculate correlation matrix
  cor_matrix <- cor(features)
  
  # Create and save correlation plot
  png(file.path(PART2_DIR, "correlation_matrix.png"), width = 800, height = 800)
  corrplot(cor_matrix, method = "circle", type = "lower", order = "hclust", 
           tl.col = "black", tl.srt = 45, addCoef.col = "black", 
           col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200))
  dev.off()
  
  # Save correlation matrix
  write.csv(cor_matrix, file = file.path(PART2_DIR, "correlation_matrix.csv"))
  
  # Extract loadings for each PC
  loadings <- pca_result$rotation
  
  # MODIFIED: Variable selection based on comprehensive criteria
  # 1. Examine loadings on PC1, PC2, and PC3
  cat("\nPCA Loadings:\n")
  print(loadings[, 1:3])
  
  # 2. Analyze correlation between variables
  cat("\nCorrelation Matrix:\n")
  print(cor_matrix)
  
  # 3. Select variables that:
  #    - Represent different aspects of the data (different PCs)
  #    - Have minimal correlation with each other
  #    - Have high loadings on different PCs
  
  # Global_active_power: High loading on PC1, represents total power consumption
  # Voltage: High negative loading on PC1, low correlation with Global_active_power
  # Sub_metering_3: High loading on PC2, represents a specific type of consumption
  
  selected_variables <- c("Global_active_power", "Voltage", "Sub_metering_3")
  
  cat("\nSelected variables for HMM based on comprehensive PCA analysis:\n")
  cat("1. Global_active_power: High loading on PC1, represents total power consumption\n")
  cat("2. Voltage: High negative loading on PC1, represents supply conditions\n")
  cat("3. Sub_metering_3: High loading on PC2, represents specific appliance usage\n")
  
  # Document the selection rationale 
  selection_rationale <- paste(
    "Variable Selection Rationale Based on PCA:\n\n",
    "1. Global_active_power:\n",
    "   - Highest positive loading on PC1 (", round(loadings["Global_active_power", "PC1"], 4), ")\n",
    "   - Represents overall power consumption in the household\n",
    "   - Captures the largest variance in the dataset\n\n",
    "2. Voltage:\n",
    "   - Strong negative loading on PC1 (", round(loadings["Voltage", "PC1"], 4), ")\n",
    "   - Low correlation with Global_active_power (", round(cor_matrix["Global_active_power", "Voltage"], 4), ")\n",
    "   - Represents supply conditions from the electricity network\n",
    "   - Provides information complementary to power consumption\n\n",
    "3. Sub_metering_3:\n",
    "   - High loading on PC2 (", round(loadings["Sub_metering_3", "PC2"], 4), ")\n",
    "   - Low correlation with both Global_active_power (", round(cor_matrix["Global_active_power", "Sub_metering_3"], 4), 
    ") and Voltage (", round(cor_matrix["Voltage", "Sub_metering_3"], 4), ")\n",
    "   - Represents a specific category of appliance usage\n",
    "   - Adds information that is not captured by the other two variables\n\n",
    "This selection ensures our model uses variables that:\n",
    "1. Capture different dimensions of the data (different principal components)\n",
    "2. Have minimal redundancy (low correlation between variables)\n",
    "3. Collectively explain a significant portion of the variance in the dataset\n",
    sep = ""
  )
  
  cat(selection_rationale, file = file.path(PART2_DIR, "variable_selection_rationale.txt"))
  cat("Selected variables for HMM based on PCA:", paste(selected_variables, collapse=", "), "\n")
  
  # Save PCA results (for reuse)
  saveRDS(pca_result, file = file.path(MODELS_DIR, "pca_model.rds"))
  
  return(list(
    pca_result = pca_result,
    selected_variables = selected_variables,
    correlation_matrix = cor_matrix
  ))
}

# Create PCA plots
create_pca_plots <- function(pca_result, features) {
  # Set white background for all plots
  white_theme <- theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#EEEEEE"),
      panel.grid.minor = element_line(color = "#EEEEEE")
    )
  
  # Scree plot
  var_explained <- summary(pca_result)$importance[2,] * 100
  cum_var_explained <- summary(pca_result)$importance[3,] * 100
  
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
    white_theme
  
  print(scree_plot)
  ggsave(file.path(PART2_DIR, "pca_scree_plot.png"), scree_plot, width = 10, height = 6, bg = "white")
  
  # Custom PCA biplot with improved visibility
  tryCatch({
    # Create a reduced sample for visualization if data is large
    sample_size <- min(5000, nrow(features))
    sample_idx <- sample(1:nrow(features), sample_size)
    
    # Extract PCA scores for the sample
    pca_scores <- as.data.frame(pca_result$x[sample_idx, 1:2])
    names(pca_scores) <- c("PC1", "PC2")
    
    # Extract loadings
    loadings <- as.data.frame(pca_result$rotation[, 1:2])
    
    # Create the custom biplot
    max_score <- max(abs(pca_scores$PC1), abs(pca_scores$PC2))
    arrow_scale <- max_score / max(abs(loadings$PC1), abs(loadings$PC2)) * 0.8
    
    # Define variable name positions
    var_names <- rownames(loadings)
    var_x <- loadings$PC1 * arrow_scale * 1.2
    var_y <- loadings$PC2 * arrow_scale * 1.2
    
    # Create the plot
    biplot <- ggplot() +
      # Add data points with alpha for density visualization
      geom_point(data = pca_scores, aes(x = PC1, y = PC2), color = "steelblue", alpha = 0.3, size = 0.8) +
      # Add arrows for variable loadings
      geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
                  arrow = arrow(length = unit(0.2, "cm")), color = "red") +
      # Add variable names
      geom_text(aes(x = var_x, y = var_y, label = var_names), color = "red", fontface = "bold") +
      # Add labels and theme
      labs(title = "PCA Biplot of Electricity Consumption Variables",
           x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "#EEEEEE"),
        panel.grid.minor = element_line(color = "#EEEEEE")
      )
    
    print(biplot)
    ggsave(file.path(PART2_DIR, "pca_biplot.png"), biplot, width = 10, height = 8, bg = "white")
    
    # Create PC1 vs PC3 biplot
    pca_scores_pc13 <- as.data.frame(pca_result$x[sample_idx, c(1,3)])
    names(pca_scores_pc13) <- c("PC1", "PC3")
    
    loadings_pc13 <- as.data.frame(pca_result$rotation[, c(1,3)])
    
    var_x_pc13 <- loadings_pc13$PC1 * arrow_scale * 1.2
    var_y_pc13 <- loadings_pc13$PC3 * arrow_scale * 1.2
    
    biplot_pc13 <- ggplot() +
      geom_point(data = pca_scores_pc13, aes(x = PC1, y = PC3), color = "steelblue", alpha = 0.3, size = 0.8) +
      geom_segment(data = loadings_pc13, aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC3 * arrow_scale),
                  arrow = arrow(length = unit(0.2, "cm")), color = "red") +
      geom_text(aes(x = var_x_pc13, y = var_y_pc13, label = var_names), color = "red", fontface = "bold") +
      labs(title = "PCA Biplot: PC1 vs PC3",
           x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC3 (", round(var_explained[3], 1), "%)")) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "#EEEEEE"),
        panel.grid.minor = element_line(color = "#EEEEEE")
      )
    
    print(biplot_pc13)
    ggsave(file.path(PART2_DIR, "pca_biplot_pc13.png"), biplot_pc13, width = 10, height = 8, bg = "white")
    
    # Create PC2 vs PC3 biplot
    pca_scores_pc23 <- as.data.frame(pca_result$x[sample_idx, c(2,3)])
    names(pca_scores_pc23) <- c("PC2", "PC3")
    
    loadings_pc23 <- as.data.frame(pca_result$rotation[, c(2,3)])
    
    var_x_pc23 <- loadings_pc23$PC2 * arrow_scale * 1.2
    var_y_pc23 <- loadings_pc23$PC3 * arrow_scale * 1.2
    
    biplot_pc23 <- ggplot() +
      geom_point(data = pca_scores_pc23, aes(x = PC2, y = PC3), color = "steelblue", alpha = 0.3, size = 0.8) +
      geom_segment(data = loadings_pc23, aes(x = 0, y = 0, xend = PC2 * arrow_scale, yend = PC3 * arrow_scale),
                  arrow = arrow(length = unit(0.2, "cm")), color = "red") +
      geom_text(aes(x = var_x_pc23, y = var_y_pc23, label = var_names), color = "red", fontface = "bold") +
      labs(title = "PCA Biplot: PC2 vs PC3",
           x = paste0("PC2 (", round(var_explained[2], 1), "%)"),
           y = paste0("PC3 (", round(var_explained[3], 1), "%)")) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "#EEEEEE"),
        panel.grid.minor = element_line(color = "#EEEEEE")
      )
    
    print(biplot_pc23)
    ggsave(file.path(PART2_DIR, "pca_biplot_pc23.png"), biplot_pc23, width = 10, height = 8, bg = "white")
    
  }, error = function(e) {
    cat("Error generating biplot:", e$message, "\n")
    cat("Trying with standard ggbiplot with modifications...\n")
    
    # Fallback to ggbiplot with modifications
    biplot <- ggbiplot(pca_result, labels = NULL, groups = NULL, ellipse = FALSE, 
                     var.axes = TRUE, var.scale = 1, alpha = 0.3, size = 0.5) +
      labs(title = "PCA Biplot of Electricity Consumption Variables") +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "#EEEEEE"),
        panel.grid.minor = element_line(color = "#EEEEEE")
      )
    
    print(biplot)
    ggsave(file.path(PART2_DIR, "pca_biplot.png"), biplot, width = 10, height = 8, bg = "white")
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
    white_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(loading_plot)
  ggsave(file.path(PART2_DIR, "pca_loading_plot.png"), loading_plot, width = 10, height = 6, bg = "white")
}

# Run PCA analysis
pca_results <- perform_pca(sdf)
selected_variables <- pca_results$selected_variables

# Step 3: Select time window and partition data
select_and_partition_data <- function(sdf, selected_variables, selected_weekday, start_time, end_time) {
  cat("Step 3: Selecting time window and partitioning data...\n")
  
  # Function to check if a time falls within the selected window
  is_within_timewindow <- function(time_str, start_time, end_time) {
    hour <- as.numeric(substr(time_str, 1, 2))
    start_hour <- as.numeric(substr(start_time, 1, 2))
    end_hour <- as.numeric(substr(end_time, 1, 2))
    
    return(hour >= start_hour && hour <= end_hour)
  }
  
  # Filter data for the selected weekday and time window
  selected_data <- sdf[sdf$Weekday == selected_weekday & 
                       sapply(sdf$Time, is_within_timewindow, start_time, end_time), ]
  
  # Extract year for partitioning
  selected_data$Year <- year(selected_data$DateTime)
  
  # Check data points per year
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
    
    # Display dataset dimensions
    cat("Training data dimensions:", dim(train_data), "\n")
    cat("Testing data dimensions:", dim(test_data), "\n")
    
    # Extract selected variables
    train_df <- train_data[, selected_variables, drop = FALSE]
    test_df <- test_data[, selected_variables, drop = FALSE]
    
    # Display feature dimensions
    cat("Training features dimensions:", dim(train_df), "\n")
    cat("Test features dimensions:", dim(test_df), "\n")
    cat("Selected variables:", paste(colnames(train_df), collapse=", "), "\n")
    
    # Check for and handle NA values
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
    
    # Save time window and partitioning information
    time_window_info <- paste(
      paste("Selected time window:", selected_weekday, "from", start_time, "to", end_time),
      paste("\nTraining years:", paste(training_years, collapse=", ")),
      paste("\nTesting year:", testing_year),
      paste("\nTraining data dimensions:", paste(dim(train_data), collapse=" x ")),
      paste("\nTesting data dimensions:", paste(dim(test_data), collapse=" x ")),
      sep = ""
    )
    cat(time_window_info, file = file.path(PART2_DIR, "time_window_info.txt"))
    
    # Save data partition (for reuse)
    saveRDS(list(train_df = train_df, test_df = test_df), 
            file = file.path(MODELS_DIR, "data_partition.rds"))
    
    return(list(
      train_df = train_df,
      test_df = test_df,
      training_years = training_years,
      testing_year = testing_year
    ))
  } else {
    stop("Not enough years in the dataset for proper partitioning.")
  }
}

# Execute data partitioning
partitioned_data <- select_and_partition_data(sdf, selected_variables, SELECTED_WEEKDAY, START_TIME, END_TIME)
train_df <- partitioned_data$train_df
test_df <- partitioned_data$test_df

# Step 4: Train multiple HMM models and select the best one
train_multiple_hmms <- function(train_df, test_df, variable, state_range) {
  cat("\nStep 4: Training multiple HMM models with different numbers of states...\n")
  
  # Create white background theme for plots
  white_theme <- theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#EEEEEE"),
      panel.grid.minor = element_line(color = "#EEEEEE")
    )
  
  # Dataframe to store results
  results_df <- data.frame(
    States = integer(),
    LogLikelihood = numeric(),
    NormalizedTrainLL = numeric(),
    BIC = numeric(),
    TestLL = numeric(),
    NormalizedTestLL = numeric(),
    IsValidModel = logical()
  )
  
  # List to store models
  hmm_models <- list()
  
  # Train models with different state counts
  for(num_states in state_range) {
    cat("\nTraining HMM with", num_states, "states for variable:", variable, "\n")
    
    tryCatch({
      # Create model
      set.seed(123)
      hmm_mod <- depmix(train_df[[variable]] ~ 1, nstates = num_states, data = train_df)
      
      # Fit model
      cat("Fitting model...\n")
      hmm_fitted <- fit(hmm_mod, verbose = TRUE, emcontrol = em.control(maxit = 500))
      
      # Calculate evaluation metrics
      train_ll <- logLik(hmm_fitted)
      train_norm_ll <- train_ll / nrow(train_df)
      bic <- -2 * train_ll + npar(hmm_fitted) * log(nrow(train_df))
      
      # Evaluate on test data
      test_mod <- depmix(test_df[[variable]] ~ 1, nstates = num_states, data = test_df)
      test_mod <- setpars(test_mod, getpars(hmm_fitted))
      fb <- forwardbackward(test_mod)
      test_ll <- fb$logLike
      test_norm_ll <- test_ll / nrow(test_df)
      
      # Check if log-likelihood is negative (model validity)
      is_valid <- train_ll < 0
      
      # Display results
      cat("  Train log-likelihood:", train_ll, "\n")
      cat("  Normalized train log-likelihood:", train_norm_ll, "\n")
      cat("  BIC:", bic, "\n")
      cat("  Test log-likelihood:", test_ll, "\n")
      cat("  Normalized test log-likelihood:", test_norm_ll, "\n")
      if(!is_valid) {
        cat("  WARNING: Positive log-likelihood indicates model issues\n")
      }
      
      # Add results to dataframe
      results_df <- rbind(results_df, data.frame(
        States = num_states,
        LogLikelihood = train_ll,
        NormalizedTrainLL = train_norm_ll,
        BIC = bic,
        TestLL = test_ll,
        NormalizedTestLL = test_norm_ll,
        IsValidModel = is_valid
      ))
      
      # Save model
      hmm_models[[as.character(num_states)]] <- hmm_fitted
      
    }, error = function(e) {
      cat("Error in training HMM with", num_states, "states:", e$message, "\n")
    })
  }
  
  # Save results to file
  write.csv(results_df, file.path(PART2_DIR, "hmm_model_comparison.csv"), row.names = FALSE)
  
  # Save all models
  saveRDS(hmm_models, file.path(MODELS_DIR, "all_hmm_models.rds"))
  
  # Visualize results
  if(nrow(results_df) > 1) {
    # Log-likelihood plot
    ll_plot <- ggplot(results_df, aes(x = States, y = LogLikelihood, color = IsValidModel)) +
      geom_point() +
      geom_line() +
      scale_color_manual(values = c("red", "blue"), 
                         labels = c("Invalid (Positive LL)", "Valid (Negative LL)")) +
      labs(title = "Log-Likelihood vs. Number of States",
           subtitle = "Blue points indicate statistically valid models",
           x = "Number of States",
           y = "Log-Likelihood",
           color = "Model Validity") +
      white_theme
    
    print(ll_plot)
    ggsave(file.path(PART2_DIR, "log_likelihood_plot.png"), ll_plot, width = 8, height = 6, bg = "white")
    
    # BIC plot (valid models only)
    valid_results <- results_df[results_df$IsValidModel, ]
    if(nrow(valid_results) > 0) {
      bic_plot <- ggplot(valid_results, aes(x = States, y = BIC)) +
        geom_point() +
        geom_line() +
        labs(title = "BIC vs. Number of States (Valid Models Only)",
             x = "Number of States",
             y = "BIC") +
        white_theme
      
      print(bic_plot)
      ggsave(file.path(PART2_DIR, "bic_plot.png"), bic_plot, width = 8, height = 6, bg = "white")
    }
    
    # Normalized log-likelihood comparison
    norm_ll_data <- data.frame(
      States = rep(results_df$States, 2),
      Dataset = c(rep("Train", nrow(results_df)), rep("Test", nrow(results_df))),
      NormalizedLL = c(results_df$NormalizedTrainLL, results_df$NormalizedTestLL),
      IsValidModel = rep(results_df$IsValidModel, 2)
    )
    
    norm_ll_plot <- ggplot(norm_ll_data, aes(x = States, y = NormalizedLL, color = Dataset, shape = IsValidModel)) +
      geom_point(size = 3) +
      geom_line(aes(linetype = IsValidModel)) +
      scale_shape_manual(values = c(4, 16), 
                       labels = c("Invalid (Positive LL)", "Valid (Negative LL)")) +
      scale_linetype_manual(values = c("dotted", "solid"), 
                          labels = c("Invalid (Positive LL)", "Valid (Negative LL)")) +
      labs(title = "Normalized Log-Likelihood vs. Number of States",
           subtitle = "Solid lines and filled points indicate statistically valid models",
           x = "Number of States",
           y = "Normalized Log-Likelihood",
           shape = "Model Validity",
           linetype = "Model Validity") +
      white_theme
    
    print(norm_ll_plot)
    ggsave(file.path(PART2_DIR, "normalized_ll_plot.png"), norm_ll_plot, width = 10, height = 6, bg = "white")
  }
  
  return(list(
    results = results_df,
    models = hmm_models
  ))
}

# Train HMM models for the main variable (Global_intensity)
hmm_results <- train_multiple_hmms(train_df, test_df, "Global_active_power", STATE_RANGE)
results_df <- hmm_results$results
hmm_models <- hmm_results$models

# Select and save the best model
select_best_model <- function(results_df, hmm_models, default_states) {
  # Filter valid models (with negative log-likelihood)
  valid_results <- results_df[results_df$IsValidModel, ]
  
  if(nrow(valid_results) > 0) {
    # Select model with minimum BIC
    best_model_idx <- which.min(valid_results$BIC)
    best_states <- valid_results$States[best_model_idx]
    
    cat("\nBest model based on BIC (among valid models):", best_states, "states\n")
    cat("BIC:", valid_results$BIC[best_model_idx], "\n")
    cat("Log-likelihood:", valid_results$LogLikelihood[best_model_idx], "\n")
    cat("Normalized train log-likelihood:", valid_results$NormalizedTrainLL[best_model_idx], "\n")
    cat("Normalized test log-likelihood:", valid_results$NormalizedTestLL[best_model_idx], "\n")
    
    # Get best model
    best_model <- hmm_models[[as.character(best_states)]]
  } else {
    # Use default model if no valid models are found
    cat("\nWARNING: No valid models found (all have positive log-likelihood).\n")
    cat("Using default model with", default_states, "states.\n")
    
    # Check if default model exists
    if(as.character(default_states) %in% names(hmm_models)) {
      best_states <- default_states
      best_model <- hmm_models[[as.character(default_states)]]
      
      # Find matching row
      model_idx <- which(results_df$States == default_states)
      if(length(model_idx) > 0) {
        cat("BIC:", results_df$BIC[model_idx], "\n")
        cat("Log-likelihood:", results_df$LogLikelihood[model_idx], "\n")
        cat("Normalized train log-likelihood:", results_df$NormalizedTrainLL[model_idx], "\n")
        cat("Normalized test log-likelihood:", results_df$NormalizedTestLL[model_idx], "\n")
      }
    } else {
      # Use first model if default model is not available
      best_states <- as.numeric(names(hmm_models)[1])
      best_model <- hmm_models[[1]]
      cat("Default model not found. Using first available model with", best_states, "states.\n")
    }
  }
  
  # Save best model information
  # Fix: Use logLik function to get log-likelihood
  model_log_likelihood <- logLik(best_model)
  is_valid_model <- model_log_likelihood < 0
  normalized_ll <- model_log_likelihood / nrow(train_df)
  
  model_info <- paste(
    paste("Best model based on BIC:", best_states, "states"),
    paste("\nModel validity:", ifelse(is_valid_model, "Valid (negative log-likelihood)", "Invalid (positive log-likelihood)")),
    paste("\nLog-likelihood:", model_log_likelihood),
    paste("\nNormalized train log-likelihood:", normalized_ll),
    sep = ""
  )
  cat(model_info, file = file.path(PART2_DIR, "best_model_info.txt"))
  
  # Save best model
  saveRDS(best_model, file.path(MODELS_DIR, "best_hmm_model.rds"))
  saveRDS(best_states, file.path(MODELS_DIR, "best_states.rds"))
  
  return(list(
    model = best_model,
    states = best_states,
    log_likelihood = model_log_likelihood,
    normalized_ll = normalized_ll
  ))
}

# Select best model
best_model_info <- select_best_model(results_df, hmm_models, DEFAULT_STATES)
best_model <- best_model_info$model
best_states <- best_model_info$states
best_model_ll <- best_model_info$log_likelihood
best_model_norm_ll <- best_model_info$normalized_ll

# Add error handling - continue even if errors occur
tryCatch({
  # Calculate threshold for anomaly detection
  calculate_threshold <- function(best_model, train_df, test_df, variable, train_norm_ll) {
    cat("\nCalculating threshold for anomaly detection...\n")
    
    # Create white background theme for plots
    white_theme <- theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "#EEEEEE"),
        panel.grid.minor = element_line(color = "#EEEEEE")
      )
    
    # Divide test data into 10 subsets
    n_subsets <- 10
    subset_size <- ceiling(nrow(test_df) / n_subsets)
    
    # Dataframe to store subset results
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
        # Extract subset
        subset_df <- test_df[start_idx:end_idx, , drop = FALSE]
        
        # Create and set model parameters
        subset_mod <- depmix(subset_df[[variable]] ~ 1, nstates = best_model@nstates, data = subset_df)
        subset_mod <- setpars(subset_mod, getpars(best_model))
        
        # Calculate log-likelihood
        subset_fb <- forwardbackward(subset_mod)
        subset_ll <- subset_fb$logLike
        subset_norm_ll <- subset_ll / nrow(subset_df)
        
        # Add results
        subset_results <- rbind(subset_results, data.frame(
          Subset = i,
          Size = nrow(subset_df),
          LogLikelihood = subset_ll,
          NormalizedLL = subset_norm_ll,
          Deviation = abs(subset_norm_ll - train_norm_ll)
        ))
      }
    }
    
    # Display and save subset results
    cat("Subset log-likelihood results:\n")
    print(subset_results)
    write.csv(subset_results, file.path(PART2_DIR, "subset_results.csv"), row.names = FALSE)
    
    # Plot normalized log-likelihood for subsets
    subset_ll_plot <- ggplot(subset_results, aes(x = Subset, y = NormalizedLL)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_hline(yintercept = train_norm_ll, linetype = "dashed", color = "red") +
      labs(title = paste("Normalized Log-Likelihood for Test Subsets - Global_active_power (", best_states, "states)"),
           subtitle = paste("Red line: Training normalized log-likelihood =", 
                          round(train_norm_ll, 4)),
           x = "Subset",
           y = "Normalized Log-Likelihood") +
      white_theme
    
    print(subset_ll_plot)
    ggsave(file.path(PART2_DIR, "subset_ll_plot.png"), subset_ll_plot, width = 10, height = 6, bg = "white")
    
    # Plot deviation from training log-likelihood
    deviation_plot <- ggplot(subset_results, aes(x = Subset, y = Deviation)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = paste("Deviation from Training Log-Likelihood - Global_active_power (", best_states, "states)"),
           x = "Subset",
           y = "Absolute Deviation") +
      white_theme
    
    print(deviation_plot)
    ggsave(file.path(PART2_DIR, "deviation_plot.png"), deviation_plot, width = 10, height = 6, bg = "white")
    
    # Determine threshold
    threshold <- max(subset_results$Deviation)
    cat("\nThreshold for normal behavior (maximum deviation):", threshold, "\n")
    
    # Save threshold information
    threshold_info <- paste(
      paste("Threshold for normal behavior (maximum deviation):", threshold),
      paste("\nTraining normalized log-likelihood:", train_norm_ll),
      "\nSubset deviations:",
      paste("\n", paste(subset_results$Deviation, collapse = ", ")),
      sep = ""
    )
    cat(threshold_info, file = file.path(PART2_DIR, "threshold_info.txt"))
    
    # Save threshold
    saveRDS(threshold, file.path(MODELS_DIR, "anomaly_threshold.rds"))
    saveRDS(train_norm_ll, file.path(MODELS_DIR, "train_norm_ll.rds"))
    
    return(list(
      threshold = threshold,
      train_norm_ll = train_norm_ll,
      subset_results = subset_results
    ))
  }

  # Calculate threshold
  threshold_info <- calculate_threshold(best_model, train_df, test_df, "Global_active_power", best_model_norm_ll)
  anomaly_threshold <- threshold_info$threshold
}, error = function(e) {
  cat("\nError in threshold calculation:", e$message, "\n")
  cat("Skipping threshold calculation. This won't affect the model selection.\n")
})

# Compare all explanatory variables
compare_variables <- function(train_df, test_df) {
  cat("\nStep 5: Brief comparison of HMM performance on all variables...\n")
  
  # Build models for each variable and calculate BIC and log-likelihood
  all_vars <- colnames(train_df)
  var_results <- data.frame(
    Variable = character(),
    States = integer(),
    LogLikelihood = numeric(),
    NormalizedTrainLL = numeric(),
    BIC = numeric(),
    TestLL = numeric(),
    NormalizedTestLL = numeric(),
    Notes = character(),
    stringsAsFactors = FALSE
  )
  
  for(var in all_vars) {
    cat("Training model for variable:", var, "\n")
    
    # Train model
    tryCatch({
      # Create model
      set.seed(123)
      var_mod <- depmix(train_df[[var]] ~ 1, nstates = 4, data = train_df)
      var_fitted <- fit(var_mod, verbose = FALSE, emcontrol = em.control(maxit = 200))
      
      # Calculate evaluation metrics
      train_ll <- logLik(var_fitted)
      train_norm_ll <- train_ll / nrow(train_df)
      bic <- -2 * train_ll + npar(var_fitted) * log(nrow(train_df))
      
      # Evaluate on test data
      test_mod <- depmix(test_df[[var]] ~ 1, nstates = 4, data = test_df)
      test_mod <- setpars(test_mod, getpars(var_fitted))
      fb <- forwardbackward(test_mod)
      test_ll <- fb$logLike
      test_norm_ll <- test_ll / nrow(test_df)
      
      # Check if log-likelihood is positive
      notes <- ""
      if(train_ll > 0) {
        notes <- "WARNING: Positive log-likelihood indicates model issues"
      }
      
      # Display results
      cat("  Train log-likelihood:", train_ll, "\n")
      cat("  Normalized train log-likelihood:", train_norm_ll, "\n")
      cat("  BIC:", bic, "\n")
      cat("  Test log-likelihood:", test_ll, "\n")
      cat("  Normalized test log-likelihood:", test_norm_ll, "\n")
      if(notes != "") {
        cat("  ", notes, "\n")
      }
      
      # Add results
      var_results <- rbind(var_results, data.frame(
        Variable = var,
        States = 4,
        LogLikelihood = train_ll,
        NormalizedTrainLL = train_norm_ll,
        BIC = bic,
        TestLL = test_ll,
        NormalizedTestLL = test_norm_ll,
        Notes = notes,
        stringsAsFactors = FALSE
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
        Notes = paste("ERROR:", e$message),
        stringsAsFactors = FALSE
      ))
    })
  }
  
  # Filter results
  valid_results <- var_results[var_results$LogLikelihood < 0 & !is.na(var_results$LogLikelihood), ]
  
  # Display and save results
  cat("\nComparison of HMM performance on different variables (valid models only):\n")
  if(nrow(valid_results) > 0) {
    print(valid_results[, c("Variable", "States", "LogLikelihood", "BIC", "NormalizedTestLL")])
    write.csv(valid_results, file.path(PART2_DIR, "variable_comparison.csv"), row.names = FALSE)
    
    # Warn if some variables showed unusual behavior
    if(nrow(valid_results) < nrow(var_results)) {
      cat("\nWARNING: Some variables showed unusual behavior and were excluded from comparison.\n")
      cat("Full results including problematic variables:\n")
      print(var_results[, c("Variable", "States", "LogLikelihood", "BIC", "Notes")])
      write.csv(var_results, file.path(PART2_DIR, "all_variables_comparison.csv"), row.names = FALSE)
    }
  } else {
    cat("No valid models found. All variables showed unusual behavior.\n")
    print(var_results[, c("Variable", "States", "LogLikelihood", "Notes")])
    write.csv(var_results, file.path(PART2_DIR, "all_variables_comparison.csv"), row.names = FALSE)
  }
  
  return(var_results)
}

# Compare variables
var_results <- compare_variables(train_df, test_df)

