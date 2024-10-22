library(ggplot2)
library(dplyr)
library(gridExtra)

# Load your CSV files if not loaded yet
# peaks_conventional <- read.csv("peaks_conventional.csv")
# peaks_enhanced <- read.csv("peaks_enhanced.csv")

### 1. Peak Intensity Comparison
pdf("peak_intensity_comparison.pdf", width = 10, height = 8)
# Plot maximum intensity comparison (maxo) with different line types, colors, and transparency
plot_maxo <- ggplot() +
  geom_density(aes(x = peaks_conventional$maxo, color = "Conventional"), size = 1.2, linetype = "solid", alpha = 0.7) +
  geom_density(aes(x = peaks_enhanced$maxo, color = "Enhanced"), size = 1.2, linetype = "dashed", alpha = 0.7) +
  labs(title = "Maximum Intensity (maxo) Comparison",
       x = "Max Intensity (maxo)", y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("Conventional" = "blue", "Enhanced" = "red"))

# Plot integrated intensity comparison (into) with different line types, colors, and transparency
plot_into <- ggplot() +
  geom_density(aes(x = peaks_conventional$into, color = "Conventional"), size = 1.2, linetype = "solid", alpha = 0.7) +
  geom_density(aes(x = peaks_enhanced$into, color = "Enhanced"), size = 1.2, linetype = "dashed", alpha = 0.7) +
  labs(title = "Integrated Intensity (into) Comparison",
       x = "Integrated Intensity (into)", y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("Conventional" = "blue", "Enhanced" = "red"))

# Arrange both plots side by side
grid.arrange(plot_maxo, plot_into, ncol = 2)
dev.off()

### 2. Peak Width and Shape Comparison
pdf("peak_width_comparison.pdf", width = 10, height = 8)
# Peak width comparison with distinguishable color and transparency
plot_width <- ggplot() +
  geom_density(aes(x = peaks_conventional$rt_width, fill = "Conventional"), alpha = 0.4, color = "blue", size = 1.2) +
  geom_density(aes(x = peaks_enhanced$rt_width, fill = "Enhanced"), alpha = 0.4, color = "red", size = 1.2) +
  labs(title = "Peak Width Comparison (RT)",
       x = "Retention Time Width", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Conventional" = "blue", "Enhanced" = "red"))

print(plot_width)
dev.off()

### 3. Signal-to-Noise Ratio (SNR) Comparison
pdf("snr_comparison.pdf", width = 10, height = 8)
# SNR comparison with better visibility (different colors, line types)
plot_snr <- ggplot() +
  geom_density(aes(x = peaks_conventional$sn, color = "Conventional"), size = 1.2, linetype = "solid", alpha = 0.7) +
  geom_density(aes(x = peaks_enhanced$sn, color = "Enhanced"), size = 1.2, linetype = "dashed", alpha = 0.7) +
  labs(title = "Signal-to-Noise Ratio (SNR) Comparison",
       x = "SNR", y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("Conventional" = "blue", "Enhanced" = "red"))

print(plot_snr)
dev.off()

### 4. Peak Overlap Based on Retention Time and m/z
pdf("peak_overlap_comparison.pdf", width = 10, height = 8)
# Plot overlap in m/z with distinct colors and transparency
plot_overlap <- ggplot(merged_peaks, aes(x = mz_conventional, y = mz_enhanced, color = sample_group_conventional)) +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(title = "Peak Overlap Based on m/z",
       x = "Conventional m/z", y = "Enhanced m/z") +
  theme_minimal() +
  scale_color_manual(values = c("Conventional" = "blue", "Enhanced" = "red"))

print(plot_overlap)
dev.off()

### 5. Area Under the Curve (AUC) Comparison
pdf("auc_comparison.pdf", width = 10, height = 8)
# AUC (into) comparison with different line types and transparency
plot_auc <- ggplot() +
  geom_density(aes(x = peaks_conventional$into, color = "Conventional"), size = 1.2, linetype = "solid", alpha = 0.7) +
  geom_density(aes(x = peaks_enhanced$into, color = "Enhanced"), size = 1.2, linetype = "dashed", alpha = 0.7) +
  labs(title = "Area Under the Curve (AUC) Comparison (into)",
       x = "Area Under the Curve (into)", y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("Conventional" = "blue", "Enhanced" = "red"))

print(plot_auc)
dev.off()
