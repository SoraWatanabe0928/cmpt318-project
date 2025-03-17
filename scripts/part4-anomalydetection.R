
# Check for needed directories
OUTPUT_DIR <- "output"
MODELS_DIR <- file.path(OUTPUT_DIR, "models")  # Where we retrieve trained model
PART4_DIR  <- file.path(OUTPUT_DIR, "part4")   # Outputs for Part 4

if (!dir.exists(MODELS_DIR)) {
  stop("Part 4 requires Trained HMM models from Part 3")
}

dir_list <- c(PART4_DIR)
for (d in dir_list) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
rm(d, dir_list)

# Install packages
packages <- c("depmixS4", "dplyr", "lubridate", "ggplot2")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)  # uses the mirror we set above
  }
  library(pkg, character.only = TRUE)
}
rm(pkg, packages)

# This block is the same as part 3, to load and process the data

# Load scaled data
scaled_data_path <- "data/processed/TermProjectData_Standardized.csv"
if (!file.exists(scaled_data_path)) {
  stop("Scaled data not found at ", scaled_data_path, 
       ". Please run Task 1 first to generate the standardized CSV.")
}

sdf <- read.csv(scaled_data_path, stringsAsFactors = FALSE)
cat("Loaded scaled dataset with", nrow(sdf), "rows and", ncol(sdf), "columns.\n")
rm(scaled_data_path)

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

rm(parsed1, parsed2)

# Check how many remain NA after both attempts
bad_count <- sum(is.na(sdf$DateTime))
if (bad_count > 0) {
  warning("Could not parse DateTime for ", bad_count, " rows. Removing those rows...\n")
  sdf <- sdf[!is.na(sdf$DateTime), ]
  cat("New data size after removing unparsed rows:", nrow(sdf), "\n")
}
rm(bad_count, needs_fallback)

# Read in the best model and its state count
model <- readRDS(file.path(MODELS_DIR, "best_multivariate_hmm.rds"))
states <- readRDS(file.path(MODELS_DIR, "best_states.rds"))

# Same 
selected_vars <- c("Global_active_power", "Voltage", "Sub_metering_3")


cat("Loaded best performing model with", states, "states\n")

# Create a new model using new data
#modNew <- depmix(EventTime~1,data=data2,transition=~Count,nstates=2,
#                 family=multinomial("identity"))
#modNew <- setpars(modNew,getpars(fm))
#modNew <- fit(modNew)
#predStates <- posterior(modNew)
#predStates$state

# Create the train dataset we use to develop our expected logLik value

# TODO FILTER BY YEAR
mutated_df <- sdf

mutated_df$Year <- lubridate::year(mutated_df$DateTime)
# Create week ID values to split per-week
mutated_df$Week_id <- lubridate::year(mutated_df$DateTime) * 100 + lubridate::week(mutated_df$DateTime)

# Learn logLik against the first 3 years, like in part 3
all_years <- sort(unique(mutated_df$Year))
train_years <- all_years[1:3]
test_year   <- all_years[4]
train_data <- mutated_df[mutated_df$Year %in% train_years, ]
test_data  <- mutated_df[mutated_df$Year == test_year, ]

rm(all_years, train_years, test_year)

# Get expected logLik
cat("Learning logLik from training set. This will be slow.\n")
all_weekids <- sort(unique(train_data$Week_id))

logLiks <- list()
for (weekid in all_weekids) {
  week_df <- train_data[train_data$Week_id == weekid, ]
  
  train_df <- week_df[, selected_vars, drop=FALSE]
  
  # Same depmix from Part 3
  newModel <- depmix(
    list(
      train_df[[1]] ~ 1,
      train_df[[2]] ~ 1,
      train_df[[3]] ~ 1
      # Add more lines if you have more selected_vars
    ),
    data    = train_df,
    nstates = states,
    family  = list(gaussian(), gaussian(), gaussian())
  )
  newModel <- setpars(newModel, getpars(model))
  
  lik <- logLik(newModel)
  
  # filter out any positive log-likelihood
  if (lik < 0) {
    logLiks <- c(logLiks, lik)
  }
}
rm(all_weekids, weekid, newModel, week_df, train_df)

rm(train_data)

logLik_mean <- mean(as.numeric(logLiks))

cat("Found the mean logLik to be", logLik_mean, "\n")

# Find the deviation using the last 10 weeks of the last year
deviation_weeks <- tail(sort(unique(test_data$Week_id)), 10)

max_deviation <- 0
for (weekid in deviation_weeks) {
  week_df <- test_data[test_data$Week_id == weekid, ]
  
  test_df <- week_df[, selected_vars, drop=FALSE]
  
  # Same depmix from Part 3
  newModel <- depmix(
    list(
      test_df[[1]] ~ 1,
      test_df[[2]] ~ 1,
      test_df[[3]] ~ 1
      # Add more lines if you have more selected_vars
    ),
    data    = test_df,
    nstates = states,
    family  = list(gaussian(), gaussian(), gaussian())
  )
  newModel <- setpars(newModel, getpars(model))
  
  lik <- as.numeric(logLik(newModel))
  
  deviation <- abs(logLik_mean - lik)
  
  if (deviation > max_deviation) {
    max_deviation <- deviation
  }
}
rm(deviation_weeks, weekid, newModel, deviation, week_df, test_df)

cat("Found the last 10 week's maximum logLik deviation to be", max_deviation, "\n")