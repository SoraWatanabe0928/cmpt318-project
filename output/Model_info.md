# CMPT 318 Term Project - Models and Results

This directory contains the models and results from Part 2 that can be used in Parts 3 and 4.

## Directory Structure
- `output/part1/`: Contains outputs from Part 1 (Feature Scaling)
- `output/part2/`: Contains outputs from Part 2 (PCA, HMM, Threshold)
- `output/models/`: Contains saved models for reuse in Parts 3 and 4

## Saved Models
- `pca_model.rds`: PCA model used for variable selection
- `data_partition.rds`: Training and test data partitions
- `all_hmm_models.rds`: All HMM models with different numbers of states
- `best_hmm_model.rds`: Best HMM model selected based on BIC
- `best_states.rds`: Number of states in the best model
- `anomaly_threshold.rds`: Threshold for anomaly detection
- `train_norm_ll.rds`: Normalized log-likelihood of training data

## How to Load Models in Part 3 and 4
```r
# Load best HMM model
best_model <- readRDS("output/models/best_hmm_model.rds")

# Load anomaly threshold
threshold <- readRDS("output/models/anomaly_threshold.rds")
train_norm_ll <- readRDS("output/models/train_norm_ll.rds")

# Example: Detect anomalies in new data
detect_anomalies <- function(new_data, best_model, threshold, train_norm_ll) {
  # Create a model for the new data
  new_mod <- depmix(new_data$Global_intensity ~ 1, nstates = best_model@nstates, data = new_data)
  new_mod <- setpars(new_mod, getpars(best_model))
  
  # Calculate log-likelihood
  fb <- forwardbackward(new_mod)
  ll <- fb$logLike
  norm_ll <- ll / nrow(new_data)
  
  # Check if it exceeds the threshold
  deviation <- abs(norm_ll - train_norm_ll)
  is_anomaly <- deviation > threshold
  
  return(list(
    norm_ll = norm_ll,
    deviation = deviation,
    is_anomaly = is_anomaly
  ))
}
```
