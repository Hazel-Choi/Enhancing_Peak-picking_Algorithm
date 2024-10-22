# Doesn't work since hashes are not identical

library(digest)
library(dplyr)

# load both peak data
peaks_conventional <- read.csv("/nobackup/xhnh76/dissertation/dataset/conventional/detected_peaks_with_groups.csv")
peaks_enhanced <- read.csv("/nobackup/xhnh76/dissertation/dataset/new_pipeline/final_peak_table.csv")


# Generate hashes
# Function to create a hash for each peak based on mz, rt, and intensity ('into')
generate_peak_hash <- function(mz, rt, into) {
  return(digest(paste(mz, rt, into, sep = "_")))
}

# Generate hashes for each peak in both datasets
peaks_conventional <- peaks_conventional %>%
  mutate(hash = mapply(generate_peak_hash, mz, rt, into))

peaks_enhanced <- peaks_enhanced %>%
  mutate(hash = mapply(generate_peak_hash, mz, rt, into))

# Check the hashes
head(peaks_conventional)
head(peaks_enhanced)


# Compare the hashes
# Find common peaks (shared hashes)
common_hashes <- intersect(peaks_conventional$hash, peaks_enhanced$hash)

# Find unique peaks in each dataset
unique_to_conventional <- setdiff(peaks_conventional$hash, peaks_enhanced$hash)
unique_to_enhanced <- setdiff(peaks_enhanced$hash, peaks_conventional$hash)

# Summary of the comparison
cat("Number of common peaks:", length(common_hashes), "\n")
cat("Number of unique peaks in conventional:", length(unique_to_conventional), "\n")
cat("Number of unique peaks in enhanced:", length(unique_to_enhanced), "\n")


# Save hash results
# Save the common and unique peaks
write.csv(common_hashes, "common_peaks.csv", row.names = FALSE)
write.csv(unique_to_conventional, "unique_peaks_conventional.csv", row.names = FALSE)
write.csv(unique_to_enhanced, "unique_peaks_enhanced.csv", row.names = FALSE)


# Visualise the comparison
# Load the necessary library for Venn diagram
library(VennDiagram)

# Draw a Venn diagram to visualize the comparison
venn.plot <- venn.diagram(
  x = list(
    "Conventional Peaks" = peaks_conventional$hash,
    "Enhanced Peaks" = peaks_enhanced$hash
  ),
  category.names = c("Conventional", "Enhanced"),
  filename = NULL,
  fill = c("blue", "red"),
  alpha = 0.5
)

# Display the Venn diagram
grid.draw(venn.plot)
