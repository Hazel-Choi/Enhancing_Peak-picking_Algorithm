library(DBI)
library(RSQLite)
library(dplyr)
library(ggplot2)
library(gridExtra)


# load both peak data
peaks_conventional <- read.csv("/nobackup/xhnh76/dissertation/dataset/conventional/detected_peaks_with_groups.csv")
peaks_enhanced <- read.csv("/nobackup/xhnh76/dissertation/dataset/new_pipeline/final_peak_table.csv")

con <- dbConnect(SQLite(), dbname = "peak_data.db")

dbWriteTable(con, "peaks_conventional", peaks_conventional, overwrite = TRUE)
dbWriteTable(con, "peaks_enhanced", peaks_enhanced, overwrite = TRUE)

# Check if the tables are in the database
dbListTables(con)

# Query some rows from 'peaks_conventional' to verify
dbGetQuery(con, "SELECT * FROM peaks_conventional LIMIT 5")
# Query some rows from 'peaks_enhanced' to verify
dbGetQuery(con, "SELECT * FROM peaks_enhanced LIMIT 5")


# The SQL query to merge the datasets based on mz and rt tolerances
tolerance_mz <- 10  # ppm tolerance for mz
tolerance_rt <- 10  # seconds tolerance for rt

query <- paste0("
  SELECT c.sample AS sample_conventional, c.sample_group AS sample_group_conventional, 
         c.mz AS mz_conventional, c.rt AS rt_conventional, 
         c.maxo AS maxo_conventional, c.\"into\" AS into_conventional, 
         e.sample AS sample_enhanced, e.sample_group AS sample_group_enhanced, 
         e.mz AS mz_enhanced, e.rt AS rt_enhanced, 
         e.maxo AS maxo_enhanced, e.\"into\" AS into_enhanced
  FROM peaks_conventional c
  JOIN peaks_enhanced e
  ON c.sample = e.sample AND c.sample_group = e.sample_group
  WHERE ABS(c.mz - e.mz) <= ", tolerance_mz, " 
  AND ABS(c.rt - e.rt) <= ", tolerance_rt, "
")

# Execute the query
merged_peaks <- dbGetQuery(con, query) #took around 6 hours to run


library(dplyr)
library(ggplot2)

# Manual quantile normalization function
quantile_normalize <- function(x, target) {
  ranks <- rank(x, ties.method = "min")
  target_quantiles <- quantile(target, probs = (ranks - 1) / (length(x) - 1))
  return(target_quantiles)
}

# Define the function to normalize intensity, perform statistical comparison, and calculate FDR, sensitivity, specificity
perform_statistical_comparison <- function(peaks_conventional, peaks_enhanced, merged_peaks, replicate_threshold = 3) {

  # Manual quantile normalization for maxo and into values
  normalized_merged_peaks <- merged_peaks %>%
    mutate(maxo_norm_conventional = quantile_normalize(maxo_conventional, maxo_enhanced),
           into_norm_conventional = quantile_normalize(into_conventional, into_enhanced),
           maxo_norm_enhanced = quantile_normalize(maxo_enhanced, maxo_conventional),
           into_norm_enhanced = quantile_normalize(into_enhanced, into_conventional))

  # Peak Intensity Comparison (maxo)
  plot_maxo <- ggplot(normalized_merged_peaks, aes(x = maxo_norm_conventional, y = maxo_norm_enhanced)) +
    geom_point() +
    labs(title = "Peak Intensity Comparison (maxo)", 
         x = "Conventional Method (maxo, normalized)", 
         y = "Enhanced Method (maxo, normalized)") +
    theme_minimal()

  # Pearson correlation for maxo
  cor_maxo <- cor(normalized_merged_peaks$maxo_norm_conventional, normalized_merged_peaks$maxo_norm_enhanced, use = "complete.obs")

  # Save plot_maxo as PDF
  pdf("plot_maxo_comparison.pdf")
  print(plot_maxo)
  dev.off()

  # Peak Intensity Comparison (into)
  plot_into <- ggplot(normalized_merged_peaks, aes(x = into_norm_conventional, y = into_norm_enhanced)) +
    geom_point() +
    labs(title = "Peak Intensity Comparison (into)", 
         x = "Conventional Method (into, normalized)", 
         y = "Enhanced Method (into, normalized)") +
    theme_minimal()

  # Pearson correlation for into
  cor_into <- cor(normalized_merged_peaks$into_norm_conventional, normalized_merged_peaks$into_norm_enhanced, use = "complete.obs")

  # Save plot_into as PDF
  pdf("plot_into_comparison.pdf")
  print(plot_into)
  dev.off()

  # SNR Comparison
  plot_snr <- ggplot(normalized_merged_peaks, aes(x = maxo_conventional, y = maxo_enhanced)) +
    geom_point() +
    labs(title = "SNR Comparison", 
         x = "Conventional Method (SNR)", 
         y = "Enhanced Method (SNR)") +
    theme_minimal()

  # Save plot_snr as PDF
  pdf("plot_snr_comparison.pdf")
  print(plot_snr)
  dev.off()

  #### Consistency and False Discovery Rate (FDR), Sensitivity, Specificity ####

  # Consistency across replicates
  consistency_conventional <- peaks_conventional %>%
    group_by(mz_conventional, rt_conventional) %>%
    summarize(rep_count_conventional = n_distinct(sample_conventional))

  consistency_enhanced <- peaks_enhanced %>%
    group_by(mz_enhanced, rt_enhanced) %>%
    summarize(rep_count_enhanced = n_distinct(sample_enhanced))

  # Merge consistency data
  replicate_consistency <- full_join(consistency_conventional, consistency_enhanced, 
                                     by = c("mz_conventional" = "mz_enhanced", 
                                            "rt_conventional" = "rt_enhanced"))

  # Classify peaks as TP, FP, FN based on replicate counts
  replicate_consistency <- replicate_consistency %>%
    mutate(TP = ifelse(rep_count_conventional >= replicate_threshold & rep_count_enhanced >= replicate_threshold, 1, 0),
           FP_conventional = ifelse(rep_count_conventional >= replicate_threshold & rep_count_enhanced < replicate_threshold, 1, 0),
           FP_enhanced = ifelse(rep_count_enhanced >= replicate_threshold & rep_count_conventional < replicate_threshold, 1, 0),
           FN_conventional = ifelse(rep_count_conventional < replicate_threshold & rep_count_enhanced >= replicate_threshold, 1, 0),
           FN_enhanced = ifelse(rep_count_enhanced < replicate_threshold & rep_count_conventional >= replicate_threshold, 1, 0))

  # Sum up true positives, false positives, and false negatives
  TP <- sum(replicate_consistency$TP)
  FP_conventional <- sum(replicate_consistency$FP_conventional)
  FP_enhanced <- sum(replicate_consistency$FP_enhanced)
  FN_conventional <- sum(replicate_consistency$FN_conventional)
  FN_enhanced <- sum(replicate_consistency$FN_enhanced)

  # Assume that the total possible number of negatives (for TN) is the sum of all possible peaks minus detected peaks
  total_peaks <- nrow(replicate_consistency)
  TN_conventional <- total_peaks - (TP + FP_conventional + FN_conventional)
  TN_enhanced <- total_peaks - (TP + FP_enhanced + FN_enhanced)

  # Conventional method metrics
  FDR_conventional <- FP_conventional / (FP_conventional + TP)
  sensitivity_conventional <- TP / (TP + FN_conventional)
  specificity_conventional <- TN_conventional / (TN_conventional + FP_conventional)

  # Enhanced method metrics
  FDR_enhanced <- FP_enhanced / (FP_enhanced + TP)
  sensitivity_enhanced <- TP / (TP + FN_enhanced)
  specificity_enhanced <- TN_enhanced / (TN_enhanced + FP_enhanced)

  # Output FDR, Sensitivity, and Specificity for both methods
  print(paste("Conventional Method - FDR:", FDR_conventional, 
              "Sensitivity:", sensitivity_conventional, 
              "Specificity:", specificity_conventional))

  print(paste("Enhanced Method - FDR:", FDR_enhanced, 
              "Sensitivity:", sensitivity_enhanced, 
              "Specificity:", specificity_enhanced))

  # Return correlations and other results
  return(list(cor_maxo = cor_maxo, cor_into = cor_into, 
              FDR_conventional = FDR_conventional, FDR_enhanced = FDR_enhanced,
              sensitivity_conventional = sensitivity_conventional, sensitivity_enhanced = sensitivity_enhanced,
              specificity_conventional = specificity_conventional, specificity_enhanced = specificity_enhanced))
}

# Run the comparison function
comparison_results <- perform_statistical_comparison(peaks_conventional, peaks_enhanced, merged_peaks)


