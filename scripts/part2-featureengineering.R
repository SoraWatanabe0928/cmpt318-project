# CMPT 318 Term Project - Part II
# Anomaly Detection in Electric Energy Consumption Data
# 設定可能なパラメータ - Part 3、4で変更する場合はここを修正

# 時間枠の設定
SELECTED_WEEKDAY <- "Monday"
START_TIME <- "18:00:00"
END_TIME <- "22:00:00"

# HMMの状態数の設定
STATE_RANGE <- c(4, 6, 8, 10, 12) # 試す状態数の範囲
DEFAULT_STATES <- 6 # 統計的に妥当な状態数（負の対数尤度を持つモデルから）

# フォルダ構造の設定
OUTPUT_DIR <- "output"
PART1_DIR <- file.path(OUTPUT_DIR, "part1") # Part 1の出力
PART2_DIR <- file.path(OUTPUT_DIR, "part2") # Part 2の出力
MODELS_DIR <- file.path(OUTPUT_DIR, "models") # モデルの保存場所

# 作業ディレクトリの確認 - コンソールに表示
cat("Current working directory:", getwd(), "\n")

# 必要なディレクトリの作成
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

# 必要なパッケージのロード
load_packages <- function() {
  packages <- c("depmixS4", "dplyr", "lubridate", "ggplot2", "stats", "reshape2")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  # ggbiplotのインストールとロード
  if (!requireNamespace("ggbiplot", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }
    devtools::install_github("vqv/ggbiplot")
  }
  suppressWarnings(library(ggbiplot))
}

load_packages()

# データの前処理 - Part 1の機能を再利用
preprocess_data <- function() {
  standardized_data_path <- "data/processed/TermProjectData_Standardized.csv"
  
  # 既に標準化されたデータが存在するか確認
  if (file.exists(standardized_data_path)) {
    cat("Using existing standardized data...\n")
    sdf <- read.csv(standardized_data_path, stringsAsFactors = FALSE)
  } else {
    cat("Standardized data not found. Please run Part 1 first.\n")
    stop("Missing standardized data")
  }
  
  # 日付時刻の変換
  sdf$DateTime <- as.POSIXct(paste(sdf$Date, sdf$Time), format="%d/%m/%Y %H:%M:%S")
  if (all(is.na(sdf$DateTime))) {
    # 別の形式を試す
    sdf$DateTime <- as.POSIXct(paste(sdf$Date, sdf$Time), format="%Y-%m-%d %H:%M:%S")
  }
  sdf$Weekday <- weekdays(sdf$DateTime)
  
  return(sdf)
}

# Step 1: データのロード
cat("\nStep 1: Loading standardized data...\n")
sdf <- preprocess_data()

# Step 2: PCA分析による変数選択
perform_pca <- function(sdf) {
  cat("Step 2: Performing PCA analysis...\n")
  variables <- c("Global_active_power", "Global_reactive_power", "Voltage", 
                "Global_intensity", "Sub_metering_1", "Sub_metering_2", "Sub_metering_3")
  features <- sdf[, variables]
  
  # NA値の確認と補完
  na_count <- colSums(is.na(features))
  if(sum(na_count) > 0) {
    cat("Found NA values in the data. Imputing with column means...\n")
    for(col in names(features)) {
      if(na_count[col] > 0) {
        features[[col]][is.na(features[[col]])] <- mean(features[[col]], na.rm = TRUE)
      }
    }
  }
  
  # PCA実行
  pca_result <- prcomp(features, center = TRUE, scale. = TRUE)
  summary_pca <- summary(pca_result)
  
  # PCA要約の表示と保存
  cat("PCA Summary:\n")
  print(summary_pca)
  capture.output(print(summary_pca), file = file.path(PART2_DIR, "pca_summary.txt"))
  
  # PCAプロットの作成
  create_pca_plots(pca_result, features)
  
  # 変数選択
  loadings_pc1 <- abs(pca_result$rotation[, "PC1"])
  selected_variables <- names(loadings_pc1)[order(loadings_pc1, decreasing = TRUE)[1:3]]
  
  cat("Selected variables for HMM based on PCA:", paste(selected_variables, collapse=", "), "\n")
  cat(paste("Selected variables for HMM based on PCA:", paste(selected_variables, collapse=", ")), 
      file = file.path(PART2_DIR, "selected_variables.txt"))
  
  # PCA結果を保存（再利用のため）
  saveRDS(pca_result, file = file.path(MODELS_DIR, "pca_model.rds"))
  
  return(list(
    pca_result = pca_result,
    selected_variables = selected_variables
  ))
}

# PCAプロットの作成
create_pca_plots <- function(pca_result, features) {
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
    theme_minimal()
  
  print(scree_plot)
  ggsave(file.path(PART2_DIR, "pca_scree_plot.png"), scree_plot, width = 10, height = 6)
  
  # Biplot
  tryCatch({
    biplot <- ggbiplot(pca_result, labels = rownames(features), groups = NULL) +
      labs(title = "PCA Biplot of Electricity Consumption Variables")
    
    print(biplot)
    ggsave(file.path(PART2_DIR, "pca_biplot.png"), biplot, width = 10, height = 8)
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
  
  # 最初の3つのコンポーネントに焦点
  loadings_long_subset <- loadings_long[loadings_long$Component %in% c("PC1", "PC2", "PC3"), ]
  
  loading_plot <- ggplot(loadings_long_subset, aes(x = Variable, y = Loading, fill = Component)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "PCA Loading Plot", 
         x = "Variable", 
         y = "Loading") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(loading_plot)
  ggsave(file.path(PART2_DIR, "pca_loading_plot.png"), loading_plot, width = 10, height = 6)
}

# PCA分析の実行
pca_results <- perform_pca(sdf)
selected_variables <- pca_results$selected_variables

# Step 3: データの時間枠選択と分割
select_and_partition_data <- function(sdf, selected_weekday, start_time, end_time) {
  cat("Step 3: Selecting time window and partitioning data...\n")
  
  # 時間枠内かどうかを確認する関数
  is_within_timewindow <- function(time_str, start_time, end_time) {
    hour <- as.numeric(substr(time_str, 1, 2))
    start_hour <- as.numeric(substr(start_time, 1, 2))
    end_hour <- as.numeric(substr(end_time, 1, 2))
    
    return(hour >= start_hour && hour <= end_hour)
  }
  
  # 選択した曜日と時間枠でデータをフィルタリング
  selected_data <- sdf[sdf$Weekday == selected_weekday & 
                       sapply(sdf$Time, is_within_timewindow, start_time, end_time), ]
  
  # 年を抽出して分割に使用
  selected_data$Year <- year(selected_data$DateTime)
  
  # 年別のデータポイント数を確認
  year_count <- table(selected_data$Year)
  cat("Data points by year:\n")
  print(year_count)
  
  # トレーニングデータ（最初の3年）とテストデータ（4年目）に分割
  all_years <- sort(unique(selected_data$Year))
  if(length(all_years) >= 4) {
    training_years <- all_years[1:3]
    testing_year <- all_years[4]
    
    train_data <- selected_data[selected_data$Year %in% training_years, ]
    test_data <- selected_data[selected_data$Year == testing_year, ]
    
    # データセットの次元を表示
    cat("Training data dimensions:", dim(train_data), "\n")
    cat("Testing data dimensions:", dim(test_data), "\n")
    
    # 選択した変数を抽出
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
    
    # 特徴量の次元を表示
    cat("Training features dimensions:", dim(train_df), "\n")
    cat("Test features dimensions:", dim(test_df), "\n")
    cat("Selected variables:", paste(colnames(train_df), collapse=", "), "\n")
    
    # NA値の確認と処理
    na_in_train <- any(is.na(train_df))
    na_in_test <- any(is.na(test_df))
    
    if(na_in_train || na_in_test) {
      cat("Warning: Found NA values in feature data. Imputing with means...\n")
      
      # トレーニングデータのNA値を補完
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
      
      # テストデータのNA値を補完
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
    
    # 時間枠と分割情報を保存
    time_window_info <- paste(
      paste("Selected time window:", selected_weekday, "from", start_time, "to", end_time),
      paste("\nTraining years:", paste(training_years, collapse=", ")),
      paste("\nTesting year:", testing_year),
      paste("\nTraining data dimensions:", paste(dim(train_data), collapse=" x ")),
      paste("\nTesting data dimensions:", paste(dim(test_data), collapse=" x ")),
      sep = ""
    )
    cat(time_window_info, file = file.path(PART2_DIR, "time_window_info.txt"))
    
    # データ分割を保存（再利用のため）
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

# データ分割の実行
partitioned_data <- select_and_partition_data(sdf, SELECTED_WEEKDAY, START_TIME, END_TIME)
train_df <- partitioned_data$train_df
test_df <- partitioned_data$test_df

# Step 4: 複数の状態数でHMMをトレーニングして最適なモデルを選択
train_multiple_hmms <- function(train_df, test_df, variable, state_range) {
  cat("\nStep 4: Training multiple HMM models with different numbers of states...\n")
  
  # 結果を保存するデータフレーム
  results_df <- data.frame(
    States = integer(),
    LogLikelihood = numeric(),
    NormalizedTrainLL = numeric(),
    BIC = numeric(),
    TestLL = numeric(),
    NormalizedTestLL = numeric(),
    IsValidModel = logical()
  )
  
  # モデルを保存するリスト
  hmm_models <- list()
  
  # 異なる状態数でモデルをトレーニング
  for(num_states in state_range) {
    cat("\nTraining HMM with", num_states, "states for variable:", variable, "\n")
    
    tryCatch({
      # モデル作成
      set.seed(123)
      hmm_mod <- depmix(train_df[[variable]] ~ 1, nstates = num_states, data = train_df)
      
      # モデルフィット
      cat("Fitting model...\n")
      hmm_fitted <- fit(hmm_mod, verbose = TRUE)
      
      # 評価指標を計算
      train_ll <- logLik(hmm_fitted)
      train_norm_ll <- train_ll / nrow(train_df)
      bic <- -2 * train_ll + npar(hmm_fitted) * log(nrow(train_df))
      
      # テストデータでの評価
      test_mod <- depmix(test_df[[variable]] ~ 1, nstates = num_states, data = test_df)
      test_mod <- setpars(test_mod, getpars(hmm_fitted))
      fb <- forwardbackward(test_mod)
      test_ll <- fb$logLike
      test_norm_ll <- test_ll / nrow(test_df)
      
      # 対数尤度が負かどうかを確認（モデルの妥当性）
      is_valid <- train_ll < 0
      
      # 結果を表示
      cat("  Train log-likelihood:", train_ll, "\n")
      cat("  Normalized train log-likelihood:", train_norm_ll, "\n")
      cat("  BIC:", bic, "\n")
      cat("  Test log-likelihood:", test_ll, "\n")
      cat("  Normalized test log-likelihood:", test_norm_ll, "\n")
      if(!is_valid) {
        cat("  WARNING: Positive log-likelihood indicates model issues\n")
      }
      
      # 結果をデータフレームに追加
      results_df <- rbind(results_df, data.frame(
        States = num_states,
        LogLikelihood = train_ll,
        NormalizedTrainLL = train_norm_ll,
        BIC = bic,
        TestLL = test_ll,
        NormalizedTestLL = test_norm_ll,
        IsValidModel = is_valid
      ))
      
      # モデルを保存
      hmm_models[[as.character(num_states)]] <- hmm_fitted
      
    }, error = function(e) {
      cat("Error in training HMM with", num_states, "states:", e$message, "\n")
    })
  }
  
  # 結果をファイルに保存
  write.csv(results_df, file.path(PART2_DIR, "hmm_model_comparison.csv"), row.names = FALSE)
  
  # 全モデルを保存
  saveRDS(hmm_models, file.path(MODELS_DIR, "all_hmm_models.rds"))
  
  # 結果の視覚化 - 対数尤度
  if(nrow(results_df) > 1) {
    # 対数尤度プロット
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
      theme_minimal()
    
    print(ll_plot)
    ggsave(file.path(PART2_DIR, "log_likelihood_plot.png"), ll_plot, width = 8, height = 6)
    
    # BICプロット
    valid_results <- results_df[results_df$IsValidModel, ]
    if(nrow(valid_results) > 0) {
      bic_plot <- ggplot(valid_results, aes(x = States, y = BIC)) +
        geom_point() +
        geom_line() +
        labs(title = "BIC vs. Number of States (Valid Models Only)",
             x = "Number of States",
             y = "BIC") +
        theme_minimal()
      
      print(bic_plot)
      ggsave(file.path(PART2_DIR, "bic_plot.png"), bic_plot, width = 8, height = 6)
    }
    
    # 正規化対数尤度の比較
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
      theme_minimal()
    
    print(norm_ll_plot)
    ggsave(file.path(PART2_DIR, "normalized_ll_plot.png"), norm_ll_plot, width = 10, height = 6)
  }
  
  return(list(
    results = results_df,
    models = hmm_models
  ))
}

# 主要変数（Global_intensity）に対してHMMをトレーニング
hmm_results <- train_multiple_hmms(train_df, test_df, "Global_intensity", STATE_RANGE)
results_df <- hmm_results$results
hmm_models <- hmm_results$models

# 最適なモデルの選択と保存
select_best_model <- function(results_df, hmm_models, default_states) {
  # 妥当なモデル（負の対数尤度を持つ）のみをフィルタリング
  valid_results <- results_df[results_df$IsValidModel, ]
  
  if(nrow(valid_results) > 0) {
    # BICが最小のモデルを選択
    best_model_idx <- which.min(valid_results$BIC)
    best_states <- valid_results$States[best_model_idx]
    
    cat("\nBest model based on BIC (among valid models):", best_states, "states\n")
    cat("BIC:", valid_results$BIC[best_model_idx], "\n")
    cat("Log-likelihood:", valid_results$LogLikelihood[best_model_idx], "\n")
    cat("Normalized train log-likelihood:", valid_results$NormalizedTrainLL[best_model_idx], "\n")
    cat("Normalized test log-likelihood:", valid_results$NormalizedTestLL[best_model_idx], "\n")
    
    # 最適モデルを取得
    best_model <- hmm_models[[as.character(best_states)]]
  } else {
    # 妥当なモデルがない場合はデフォルトモデルを使用
    cat("\nWARNING: No valid models found (all have positive log-likelihood).\n")
    cat("Using default model with", default_states, "states.\n")
    
    # デフォルトモデルがあるか確認
    if(as.character(default_states) %in% names(hmm_models)) {
      best_states <- default_states
      best_model <- hmm_models[[as.character(default_states)]]
      
      # 該当する行を検索
      model_idx <- which(results_df$States == default_states)
      if(length(model_idx) > 0) {
        cat("BIC:", results_df$BIC[model_idx], "\n")
        cat("Log-likelihood:", results_df$LogLikelihood[model_idx], "\n")
        cat("Normalized train log-likelihood:", results_df$NormalizedTrainLL[model_idx], "\n")
        cat("Normalized test log-likelihood:", results_df$NormalizedTestLL[model_idx], "\n")
      }
    } else {
      # デフォルトモデルもない場合は、最初のモデルを使用
      best_states <- as.numeric(names(hmm_models)[1])
      best_model <- hmm_models[[1]]
      cat("Default model not found. Using first available model with", best_states, "states.\n")
    }
  }
  
  # 最適モデルの情報を保存
  # エラー修正: logLik関数を使用して対数尤度を取得する
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
  
  # 最適モデルを保存
  saveRDS(best_model, file.path(MODELS_DIR, "best_hmm_model.rds"))
  saveRDS(best_states, file.path(MODELS_DIR, "best_states.rds"))
  
  return(list(
    model = best_model,
    states = best_states,
    log_likelihood = model_log_likelihood,
    normalized_ll = normalized_ll
  ))
}

# 最適モデルの選択
best_model_info <- select_best_model(results_df, hmm_models, DEFAULT_STATES)
best_model <- best_model_info$model
best_states <- best_model_info$states
best_model_ll <- best_model_info$log_likelihood
best_model_norm_ll <- best_model_info$normalized_ll

# エラーハンドリングを追加 - ここでエラーが発生しても続行する
tryCatch({
  # 異常検知のための閾値設定
  calculate_threshold <- function(best_model, train_df, test_df, variable, train_norm_ll) {
    cat("\nCalculating threshold for anomaly detection...\n")
    
    # テストデータを10個のサブセットに分割
    n_subsets <- 10
    subset_size <- ceiling(nrow(test_df) / n_subsets)
    
    # サブセット結果を保存するデータフレーム
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
        # サブセットを抽出
        subset_df <- test_df[start_idx:end_idx, , drop = FALSE]
        
        # モデル作成とパラメータ設定
        subset_mod <- depmix(subset_df[[variable]] ~ 1, nstates = best_model@nstates, data = subset_df)
        subset_mod <- setpars(subset_mod, getpars(best_model))
        
        # 対数尤度計算
        subset_fb <- forwardbackward(subset_mod)
        subset_ll <- subset_fb$logLike
        subset_norm_ll <- subset_ll / nrow(subset_df)
        
        # 結果を追加
        subset_results <- rbind(subset_results, data.frame(
          Subset = i,
          Size = nrow(subset_df),
          LogLikelihood = subset_ll,
          NormalizedLL = subset_norm_ll,
          Deviation = abs(subset_norm_ll - train_norm_ll)
        ))
      }
    }
    
    # サブセット結果を表示・保存
    cat("Subset log-likelihood results:\n")
    print(subset_results)
    write.csv(subset_results, file.path(PART2_DIR, "subset_results.csv"), row.names = FALSE)
    
    # サブセットの正規化対数尤度をプロット
    subset_ll_plot <- ggplot(subset_results, aes(x = Subset, y = NormalizedLL)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_hline(yintercept = train_norm_ll, linetype = "dashed", color = "red") +
      labs(title = paste("Normalized Log-Likelihood for Test Subsets - Global_intensity (", best_states, "states)"),
           subtitle = paste("Red line: Training normalized log-likelihood =", 
                          round(train_norm_ll, 4)),
           x = "Subset",
           y = "Normalized Log-Likelihood") +
      theme_minimal()
    
    print(subset_ll_plot)
    ggsave(file.path(PART2_DIR, "subset_ll_plot.png"), subset_ll_plot, width = 10, height = 6)
    
    # 訓練データからの偏差をプロット
    deviation_plot <- ggplot(subset_results, aes(x = Subset, y = Deviation)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = paste("Deviation from Training Log-Likelihood - Global_intensity (", best_states, "states)"),
           x = "Subset",
           y = "Absolute Deviation") +
      theme_minimal()
    
    print(deviation_plot)
    ggsave(file.path(PART2_DIR, "deviation_plot.png"), deviation_plot, width = 10, height = 6)
    
    # 閾値の決定
    threshold <- max(subset_results$Deviation)
    cat("\nThreshold for normal behavior (maximum deviation):", threshold, "\n")
    
    # 閾値情報を保存
    threshold_info <- paste(
      paste("Threshold for normal behavior (maximum deviation):", threshold),
      paste("\nTraining normalized log-likelihood:", train_norm_ll),
      "\nSubset deviations:",
      paste("\n", paste(subset_results$Deviation, collapse = ", ")),
      sep = ""
    )
    cat(threshold_info, file = file.path(PART2_DIR, "threshold_info.txt"))
    
    # 閾値を保存
    saveRDS(threshold, file.path(MODELS_DIR, "anomaly_threshold.rds"))
    saveRDS(train_norm_ll, file.path(MODELS_DIR, "train_norm_ll.rds"))
    
    return(list(
      threshold = threshold,
      train_norm_ll = train_norm_ll,
      subset_results = subset_results
    ))
  }

  # 閾値計算
  threshold_info <- calculate_threshold(best_model, train_df, test_df, "Global_intensity", best_model_norm_ll)
  anomaly_threshold <- threshold_info$threshold
}, error = function(e) {
  cat("\nError in threshold calculation:", e$message, "\n")
  cat("Skipping threshold calculation. This won't affect the model selection.\n")
})

# すべての説明変数に対する簡単な比較
compare_variables <- function(train_df, test_df) {
  cat("\nStep 5: Brief comparison of HMM performance on all variables...\n")
  
  # 各変数のモデルを構築し、BICとログ尤度を計算
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
    
    # トレーニングモデル
    tryCatch({
      # モデル作成
      set.seed(123)
      var_mod <- depmix(train_df[[var]] ~ 1, nstates = 4, data = train_df)
      var_fitted <- fit(var_mod, verbose = FALSE, emcontrol = em.control(maxit = 200))
      
      # 評価指標を計算
      train_ll <- logLik(var_fitted)
      train_norm_ll <- train_ll / nrow(train_df)
      bic <- -2 * train_ll + npar(var_fitted) * log(nrow(train_df))
      
      # テストデータでの評価
      test_mod <- depmix(test_df[[var]] ~ 1, nstates = 4, data = test_df)
      test_mod <- setpars(test_mod, getpars(var_fitted))
      fb <- forwardbackward(test_mod)
      test_ll <- fb$logLike
      test_norm_ll <- test_ll / nrow(test_df)
      
      # 対数尤度が正か負かの確認
      notes <- ""
      if(train_ll > 0) {
        notes <- "WARNING: Positive log-likelihood indicates model issues"
      }
      
      # 結果を表示
      cat("  Train log-likelihood:", train_ll, "\n")
      cat("  Normalized train log-likelihood:", train_norm_ll, "\n")
      cat("  BIC:", bic, "\n")
      cat("  Test log-likelihood:", test_ll, "\n")
      cat("  Normalized test log-likelihood:", test_norm_ll, "\n")
      if(notes != "") {
        cat("  ", notes, "\n")
      }
      
      # 結果を追加
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
  
  # 結果をフィルタリング
  valid_results <- var_results[var_results$LogLikelihood < 0 & !is.na(var_results$LogLikelihood), ]
  
  # 結果を表示・保存
  cat("\nComparison of HMM performance on different variables (valid models only):\n")
  if(nrow(valid_results) > 0) {
    print(valid_results[, c("Variable", "States", "LogLikelihood", "BIC", "NormalizedTestLL")])
    write.csv(valid_results, file.path(PART2_DIR, "variable_comparison.csv"), row.names = FALSE)
    
    # 異常を示した変数があれば警告
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

# 変数比較
var_results <- compare_variables(train_df, test_df)

# Part 3, 4での利用のためのREADMEファイルを作成
create_readme <- function() {
  readme_content <- paste(
    "# CMPT 318 Term Project - Models and Results\n\n",
    "This directory contains the models and results from Part 2 that can be used in Parts 3 and 4.\n\n",
    "## Directory Structure\n",
    "- `output/part1/`: Contains outputs from Part 1 (Feature Scaling)\n",
    "- `output/part2/`: Contains outputs from Part 2 (PCA, HMM, Threshold)\n",
    "- `output/models/`: Contains saved models for reuse in Parts 3 and 4\n\n",
    "## Saved Models\n",
    "- `pca_model.rds`: PCA model used for variable selection\n",
    "- `data_partition.rds`: Training and test data partitions\n",
    "- `all_hmm_models.rds`: All HMM models with different numbers of states\n",
    "- `best_hmm_model.rds`: Best HMM model selected based on BIC\n",
    "- `best_states.rds`: Number of states in the best model\n",
    "- `anomaly_threshold.rds`: Threshold for anomaly detection\n",
    "- `train_norm_ll.rds`: Normalized log-likelihood of training data\n\n",
    "## How to Load Models in Part 3 and 4\n",
    "```r\n",
    "# Load best HMM model\n",
    "best_model <- readRDS(\"output/models/best_hmm_model.rds\")\n",
    "\n",
    "# Load anomaly threshold\n",
    "threshold <- readRDS(\"output/models/anomaly_threshold.rds\")\n",
    "train_norm_ll <- readRDS(\"output/models/train_norm_ll.rds\")\n",
    "\n",
    "# Example: Detect anomalies in new data\n",
    "detect_anomalies <- function(new_data, best_model, threshold, train_norm_ll) {\n",
    "  # Create a model for the new data\n",
    "  new_mod <- depmix(new_data$Global_intensity ~ 1, nstates = best_model@nstates, data = new_data)\n",
    "  new_mod <- setpars(new_mod, getpars(best_model))\n",
    "  \n",
    "  # Calculate log-likelihood\n",
    "  fb <- forwardbackward(new_mod)\n",
    "  ll <- fb$logLike\n",
    "  norm_ll <- ll / nrow(new_data)\n",
    "  \n",
    "  # Check if it exceeds the threshold\n",
    "  deviation <- abs(norm_ll - train_norm_ll)\n",
    "  is_anomaly <- deviation > threshold\n",
    "  \n",
    "  return(list(\n",
    "    norm_ll = norm_ll,\n",
    "    deviation = deviation,\n",
    "    is_anomaly = is_anomaly\n",
    "  ))\n",
    "}\n",
    "```\n",
    sep = ""
  )
  
  cat(readme_content, file = file.path(OUTPUT_DIR, "README.md"))
}

create_readme()

cat("\nTerm Project Part II analysis completed.\n")
cat("Check the following directories for all generated results and figures:\n")
cat("  - ", PART2_DIR, " (Analysis results and plots)\n")
cat("  - ", MODELS_DIR, " (Saved models for Part 3 and 4)\n")