# Load necessary libraries
library(xcms)
library(MSnbase)
library(dplyr)
library(BiocParallel)
library(progress)  # Progress bar package


bpparam <- MulticoreParam(workers = 4)  # Adjust 'workers' based on available cores

# Step 1: Import mzXML Files from Directory and Define Sample Groups
import_mzXML_files_with_groups <- function(directory) {
  mzxml_files <- list.files(path = directory, pattern = "*.mzXML", full.names = TRUE)
  
  sample_names <- basename(mzxml_files)
  sample_groups <- rep(NA, length(sample_names))
  
  sample_groups[grep("UWA479|UWA420|UWA506|UWA531|UWA431|UWA623", sample_names)] <- "Control"
  sample_groups[grep("UWA626|UWA596|UWA419|UWA422|UWA556|UWA558", sample_names)] <- "AsymAD"
  sample_groups[grep("UWA576|UWA579|UWA580|UWA612|UWA614|UWA530", sample_names)] <- "AD"
  
  return(list(files = mzxml_files, groups = sample_groups))
}

# Step 2: Calculate Signal Density and Set Dynamic Windows with Progress Bar and Error Handling
calculate_density_and_noise <- function(msdata, window_size = 50) {
  rt <- rtime(msdata)
  intensity_list <- intensity(msdata)
  
  # Convert intensity list to a numeric vector
  intensity <- unlist(intensity_list)
  
  # Define sliding windows
  windows <- seq(min(rt), max(rt), by = window_size)
  
  # Create a progress bar
  pb <- progress_bar$new(
    format = "  Calculating signal density [:bar] :percent in :elapsed",
    total = length(windows),
    clear = FALSE
  )
  
  # Calculate signal density in each window
  signal_density <- sapply(windows, function(w) {
    pb$tick()  # Update progress bar
    
    idx <- which(rt >= w & rt < (w + window_size))
    if (length(idx) > 0) {
      return(sum(intensity[idx]))  # Sum the intensities within the window
    } else {
      return(0)  # Return 0 for empty windows
    }
  })
  
  # Estimate noise level (median signal intensity in low-density regions)
  noise_level <- median(intensity[signal_density < quantile(signal_density, 0.1)])
  
  return(list(signal_density = signal_density, noise_level = noise_level, windows = windows))
}

# Step 3: Adjust Window Sizes Dynamically
adjust_window_sizes <- function(signal_density, threshold = 0.75, small_window = 10, large_window = 50) {
  # Use small windows for dense areas (density above threshold) and large windows for sparse areas
  window_sizes <- ifelse(signal_density > quantile(signal_density, threshold), small_window, large_window)
  return(window_sizes)
}

# Step 4: Peak Detection and Dynamic Parameter Selection with Error Handling and Progress Bar
test_peak_detection_with_best_params <- function(msdata, windows, window_sizes, param_grid, bpparam) {
  results <- list()

  # Initialize progress bar
  pb <- progress_bar$new(
    format = "  Processing windows [:bar] :percent in :elapsed",
    total = length(windows),
    clear = FALSE
  )

  # Loop through windows
  for (i in seq_along(windows)) {
    window_rt <- windows[i]
    window_size <- window_sizes[i]
    pb$tick()

    result <- tryCatch({
      msdata_window <- filterRt(msdata, rt = c(window_rt, window_rt + window_size))

      window_results <- bplapply(seq_along(param_grid), function(j) {
        params <- param_grid[[j]]
        start_time <- Sys.time()

        # Create CentWaveParam with the current parameter set
        cwp <- CentWaveParam(ppm = params$ppm, peakwidth = params$peakwidth, snthresh = params$snthresh, prefilter = params$prefilter)

        chrom_peaks <- tryCatch({
          findChromPeaks(msdata_window, param = cwp)
        }, error = function(e) {
          message(sprintf("Error detecting peaks with parameter set %d: %s", j, e$message))
          return(NULL)
        })

        end_time <- Sys.time()
        time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))

        # If peaks are detected
        if (!is.null(chrom_peaks) && nrow(chromPeaks(chrom_peaks)) > 0) {
          num_peaks <- nrow(chromPeaks(chrom_peaks))
          snr <- calculate_snr(chrom_peaks)
          peak_width_stats <- calculate_peak_width_stats(chrom_peaks)
          intensity_stats <- calculate_intensity_stats(chrom_peaks)

          # Scoring based on peak detection metrics
          score <- (0.3 * num_peaks) + 
                   (0.2 * snr) + 
                   (0.2 * -peak_width_stats$median_width) +  
                   (0.2 * intensity_stats$mean_intensity) +  
                   (0.1 * -time_taken)

          return(list(
            params = params,
            num_peaks = num_peaks,
            time_taken = time_taken,
            snr = snr,
            peak_width_stats = peak_width_stats,
            intensity_stats = intensity_stats,
            score = score,
            detected_peaks = chrom_peaks
          ))
        } else {
          return(NULL)
        }
      }, BPPARAM = bpparam)

      return(window_results)

    }, error = function(e) {
      message(sprintf("Error detected in window %d: %s", i, e$message))
      return(NULL)
    })

    # Store results if they exist
    if (!is.null(result)) {
      results[[paste0("Window_", i)]] <- result
    }
  }

  # Flatten results to a single list
  all_results <- do.call(c, results)
  
  # Filter out NULL values
  all_results <- Filter(Negate(is.null), all_results)
  
  if (length(all_results) == 0) {
    message("No valid results found. Returning default parameters.")
    return(param_grid[[1]])  # Return default parameter set
  }

  # Find the best result by highest score
  best_result <- all_results[[which.max(sapply(all_results, function(x) x$score))]]

  return(best_result$params)  # Return only the best parameter set
}


# Function to calculate SNR
calculate_snr <- function(chrom_peaks) {
  peak_data <- chromPeaks(chrom_peaks)  # Access peak data matrix
  if ("sn" %in% colnames(peak_data)) {  # Check if 'sn' column exists
    snr <- peak_data[, "sn"]
    return(median(snr, na.rm = TRUE))  # Use mean SNR
  } else {
    return(NA)  # Return NA if 'sn' column is not present
  }
}

# Function to calculate peak width stats
calculate_peak_width_stats <- function(chrom_peaks) {
  peak_data <- chromPeaks(chrom_peaks)  # Access peak data matrix
  peak_width <- peak_data[, "rtmax"] - peak_data[, "rtmin"]  # Calculate peak width
  return(list(
    median_width = median(peak_width, na.rm = TRUE), 
    width_range = range(peak_width, na.rm = TRUE)
  ))
}

# Function to calculate intensity stats
calculate_intensity_stats <- function(chrom_peaks) {
  peak_data <- chromPeaks(chrom_peaks)  # Access peak data matrix
  intensity <- peak_data[, "into"]  # Access 'into' (integrated intensity) column
  return(list(
    mean_intensity = mean(intensity, na.rm = TRUE), 
    median_intensity = median(intensity, na.rm = TRUE),
    intensity_range = range(intensity, na.rm = TRUE)
  ))
}


run_peak_detection_with_best_params <- function(directory, param_grid) {
  mzxml_info <- import_mzXML_files_with_groups(directory)
  mzxml_files <- mzxml_info$files
  sample_groups <- mzxml_info$groups

  # Create progress bar for file processing
  pb_files <- progress_bar$new(
    format = "  Processing files [:bar] :percent in :elapsed",
    total = length(mzxml_files),
    clear = FALSE
  )

  peak_detection_results <- list()

  # Process each file (replicate)
  for (file_idx in seq_along(mzxml_files)) {
    pb_files$tick()  # Update file processing progress bar
    
    # Read the mzXML file and handle any errors that occur during reading
    raw_data <- tryCatch({
      readMSData(files = mzxml_files[file_idx], mode = "onDisk")
    }, error = function(e) {
      message(sprintf("Error reading file %s: %s", mzxml_files[file_idx], e$message))
      return(NULL)
    })

    if (is.null(raw_data)) next  # Skip if reading the file failed

    # Add sample group labels to the raw_data object
    pData(raw_data)$sample_group <- sample_groups[file_idx]

    # Step 1: Calculate signal density and set dynamic windows
    density_info <- tryCatch({
      calculate_density_and_noise(raw_data, window_size = 50)
    }, error = function(e) {
      message(sprintf("Error calculating density for file %s: %s", mzxml_files[file_idx], e$message))
      return(NULL)
    })

    if (is.null(density_info)) next  # Skip if signal density calculation failed

    # Step 2: Adjust window sizes based on signal density
    window_sizes <- tryCatch({
      adjust_window_sizes(density_info$signal_density)
    }, error = function(e) {
      message(sprintf("Error adjusting windows for file %s: %s", mzxml_files[file_idx], e$message))
      return(NULL)
    })

    if (is.null(window_sizes)) next  # Skip if window adjustment failed

    # Step 3: Find the best parameter set for peak detection
    best_params <- tryCatch({
      test_peak_detection_with_best_params(raw_data, density_info$windows, window_sizes, param_grid, bpparam = bpparam)
    }, error = function(e) {
      message(sprintf("Error detecting peaks for file %s: %s", mzxml_files[file_idx], e$message))
      return(NULL)
    })

    default_params_set <- list(
      ppm = 10,
      peakwidth = c(5, 20),
      snthresh = 5,
      prefilter = c(3, 200)
    )

    # Initialize variables to store the best score and the corresponding params
    best_score <- -Inf
    best_params_set <- NULL

    # Loop through the first 3 sublists in best_params to compare scores
    for (i in 1:3) {
      current_score <- best_params[[i]]$score  # Access the score of the current sublist
  
      # Check if the current score is higher than the best score found so far
      if (current_score > best_score) {
        best_score <- current_score            # Update the best score
        best_params_set <- best_params[[i]]$params  # Store the params of the best sublist
      }
    }

    if (is.null(best_params_set)) {
      message("No valid best parameters found, using default parameter set.")
      best_params_set <- default_params_set
    }  # Skip further processing if peak detection failed
    

    # Step 4: Apply the best parameter set to detect peaks
    cwp <- CentWaveParam(ppm = best_params_set$ppm, 
                         peakwidth = c(best_params_set$peakwidth[1], best_params_set$peakwidth[2]), 
                         snthresh = best_params_set$snthresh, 
                         prefilter = c(best_params_set$prefilter[1], best_params_set$prefilter[2])
    )

    detected_peaks <- tryCatch({
      findChromPeaks(raw_data, param = cwp)
    }, error = function(e) {
      message(sprintf("Error detecting final peaks for file %s: %s", mzxml_files[file_idx], e$message))
      return(NULL)
    })

    if (!is.null(detected_peaks)) {
      peak_detection_results[[file_idx]] <- detected_peaks
    }
  }

  return(peak_detection_results)  # Return the list of detected peaks for further processing
}
  


run_peak_alignment_grouping <- function(peak_detection_results, sample_groups) {
  if (length(peak_detection_results) == 0) {
    message("No peaks detected. Cannot proceed with alignment and grouping.")
    return(NULL)
  }

  # Step 1: Combine the peak detection results into a single XCMSnExp object
  combined_peaks <- tryCatch({
    do.call(c, peak_detection_results)  # Combine all XCMSnExp objects
  }, error = function(e) {
    message("Error combining peak detection results: ", e$message)
    return(NULL)
  })

  if (is.null(combined_peaks)) return(NULL)  # Exit if combining failed

  # Step 2: Align chromatographic peaks across samples
  aligned_peaks <- tryCatch({
    adjustRtime(combined_peaks, param = ObiwarpParam(binSize = 0.6, response = 1, distFun = "cor"))
  }, error = function(e) {
    message("Error aligning peaks: ", e$message)
    return(NULL)
  })

  if (is.null(aligned_peaks)) return(NULL)  # Exit if peak alignment failed

  # Step 3: Group peaks across samples using sample group information
  grouped_aligned_peaks <- tryCatch({
    groupChromPeaks(aligned_peaks, param = PeakDensityParam(sampleGroups = sample_groups))
  }, error = function(e) {
    message("Error grouping peaks: ", e$message)
    return(NULL)
  })

  if (is.null(grouped_aligned_peaks)) return(NULL)  # Exit if peak grouping failed

  # Step 4: Fill missing peaks (optional)
  filled_aligned_peaks <- tryCatch({
    fillChromPeaks(grouped_aligned_peaks, BPPARAM = SerialParam())
  }, error = function(e) {
    message("Error filling peaks: ", e$message)
    return(NULL)
  })

  if (is.null(filled_aligned_peaks)) return(NULL)  # Exit if filling peaks failed

  # Step 5: Extract the final peak table with sample group information
  peak_table <- chromPeaks(filled_aligned_peaks)
  peak_data <- as.data.frame(peak_table)

  # Step 6: Map sample index to sample group and add to peak table
  peak_sample_indices <- peak_data$sample
  peak_sample_groups <- sample_groups[peak_sample_indices]
  peak_data$sample_group <- peak_sample_groups

  return(peak_data)  # Return the final peak table with sample group information
}




## USAGE ##
directory <- "/nobackup/xhnh76/dissertation/dataset"
param_grid <- list(
  list(ppm = 5, peakwidth = c(2, 10), snthresh = 5, prefilter = c(3, 100)),
  list(ppm = 10, peakwidth = c(5, 20), snthresh = 3, prefilter = c(3, 200)),
  list(ppm = 15, peakwidth = c(10, 30), snthresh = 8, prefilter = c(5, 500))
)

peak_detection_results <- run_peak_detection_with_best_params(directory, param_grid)

sample_groups <- import_mzXML_files_with_groups(directory)$groups
final_peak_table <- run_peak_alignment_grouping(peak_detection_results, sample_groups)

# Save final peak table
if (!is.null(final_peak_table)) {
  write.csv(final_peak_table, file = "final_peak_table.csv", row.names = FALSE)
}





## With Error log ##

# Open the error log file at the start of the script
error_log <- file("error_log.txt", open = "wt")

run_peak_detection_with_best_params <- function(directory, param_grid) {
  mzxml_info <- import_mzXML_files_with_groups(directory)
  mzxml_files <- mzxml_info$files
  sample_groups <- mzxml_info$groups

  # Create progress bar for file processing
  pb_files <- progress_bar$new(
    format = "  Processing files [:bar] :percent in :elapsed",
    total = length(mzxml_files),
    clear = FALSE
  )

  peak_detection_results <- list()

  # Process each file (replicate)
  for (file_idx in seq_along(mzxml_files)) {
    pb_files$tick()  # Update file processing progress bar
    
    # Read the mzXML file and handle any errors that occur during reading
    raw_data <- tryCatch({
      readMSData(files = mzxml_files[file_idx], mode = "onDisk")
    }, error = function(e) {
      message <- sprintf("Error reading file %s: %s\n", mzxml_files[file_idx], e$message)
      write(message, error_log)
      return(NULL)
    })

    if (is.null(raw_data)) next  # Skip if reading the file failed

    # Add sample group labels to the raw_data object
    pData(raw_data)$sample_group <- sample_groups[file_idx]

    # Step 1: Calculate signal density and set dynamic windows
    density_info <- tryCatch({
      calculate_density_and_noise(raw_data, window_size = 50)
    }, error = function(e) {
      message <- sprintf("Error calculating density for file %s: %s\n", mzxml_files[file_idx], e$message)
      write(message, error_log)
      return(NULL)
    })

    if (is.null(density_info)) next  # Skip if signal density calculation failed

    # Step 2: Adjust window sizes based on signal density
    window_sizes <- tryCatch({
      adjust_window_sizes(density_info$signal_density)
    }, error = function(e) {
      message <- sprintf("Error adjusting windows for file %s: %s\n", mzxml_files[file_idx], e$message)
      write(message, error_log)
      return(NULL)
    })

    if (is.null(window_sizes)) next  # Skip if window adjustment failed

    # Step 3: Find the best parameter set for peak detection
    best_params <- tryCatch({
      test_peak_detection_with_best_params(raw_data, density_info$windows, window_sizes, param_grid, bpparam = bpparam)
    }, error = function(e) {
      message(sprintf("Error detecting peaks for file %s: %s", mzxml_files[file_idx], e$message))
      return(NULL)
    })

    default_params_set <- list(
      ppm = 10,
      peakwidth = c(5, 20),
      snthresh = 5,
      prefilter = c(3, 200)
    )

    # Initialize variables to store the best score and the corresponding params
    best_score <- -Inf
    best_params_set <- NULL

    # Loop through the first 3 sublists in best_params to compare scores
    for (i in 1:3) {
      if (!is.null(best_params[[i]]) && !is.null(best_params[[i]]$score)) {
        current_score <- best_params[[i]]$score  # Access the score of the current sublist

        if (current_score > best_score) {
          best_score <- current_score
          best_params_set <- best_params[[i]]$params
        }
      }
    }

    if (is.null(best_params_set)) {
      message("No valid best parameters found, using default parameter set.")
      best_params_set <- default_params_set
    }  
    
    # Step 4: Apply the best parameter set to detect peaks
    cwp <- CentWaveParam(ppm = best_params$ppm, 
                         peakwidth = c(best_params$peakwidth[1], best_params$peakwidth[2]), 
                         snthresh = best_params$snthresh, 
                         prefilter = c(best_params$prefilter[1], best_params$prefilter[2])
    )

    detected_peaks <- tryCatch({
      findChromPeaks(raw_data, param = cwp)
    }, error = function(e) {
      message <- sprintf("Error detecting peaks for file %s: %s\n", mzxml_files[file_idx], e$message)
      write(message, error_log)
      return(NULL)
    })

    if (!is.null(detected_peaks)) {
      peak_detection_results[[file_idx]] <- detected_peaks
    }
  }

  return(peak_detection_results)  # Return the list of detected peaks for further processing
}

# Close the error log file at the end of the script
close(error_log)


