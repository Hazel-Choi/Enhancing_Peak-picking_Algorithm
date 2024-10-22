# Load required libraries
library(xcms)
library(CAMERA)
library(MSnbase)

# Define the directory where the mzXML files are located
mzxml_dir <- getwd()

# List all mzXML files
mzxml_files <- list.files(mzxml_dir, pattern = "*.mzXML", full.names = TRUE)

# Define sample group labels corresponding to your samples
# Control group: ["UWA479", "UWA420", "UWA506", "UWA531", "UWA431", "UWA623"]
# AsymAD group: ["UWA626", "UWA596", "UWA419", "UWA422", "UWA556", "UWA558"]
# AD group: ["UWA576", "UWA579", "UWA580", "UWA612", "UWA614", "UWA530"]
sample_names <- basename(mzxml_files)

# Create a vector of sample group labels (Control, AsymAD, AD)
sample_groups <- rep(NA, length(sample_names))
sample_groups[grep("UWA479|UWA420|UWA506|UWA531|UWA431|UWA623", sample_names)] <- "Control"
sample_groups[grep("UWA626|UWA596|UWA419|UWA422|UWA556|UWA558", sample_names)] <- "AsymAD"
sample_groups[grep("UWA576|UWA579|UWA580|UWA612|UWA614|UWA530", sample_names)] <- "AD"


# Set parameters for the CentWave peak detection method
peak_params <- CentWaveParam(ppm = 15, 
                             peakwidth = c(10, 60), 
                             snthresh = 10, 
                             noise = 1000, 
                             prefilter = c(3, 100), 
                             mzCenterFun = "wMean", 
                             integrate = 1, 
                             mzdiff = -0.001)

# Read the mzXML files and create an xcmsSet object
raw_data <- readMSData(files = mzxml_files, mode = "onDisk")

# Add the sample group labels to the raw_data object
pData(raw_data)$sample_group <- sample_groups

# Perform the peak detection using the CentWave method
detected_peaks <- findChromPeaks(raw_data, param = peak_params)

#Peak Alignment
# Align the chromatographic peaks across samples
aligned_peaks <- adjustRtime(detected_peaks, param = ObiwarpParam(binSize = 0.6, response = 1, distFun = "cor"))

# Group peaks across samples with the sample group information
grouped_aligned_peaks <- groupChromPeaks(aligned_peaks, 
                                         param = PeakDensityParam(sampleGroups = pData(raw_data)$sample_group))

# Fill in missing peaks (optional)
filled_aligned_peaks <- fillChromPeaks(grouped_aligned_peaks, BPPARAM = SerialParam())


# Extract peak table with sample group information
peak_table <- chromPeaks(filled_aligned_peaks)
peak_data <- as.data.frame(peak_table)

# Extract the sample index from the peak table
peak_sample_indices <- peak_data$sample

# Map the sample index to sample group
peak_sample_groups <- pData(filled_aligned_peaks)$sample_group[peak_sample_indices]

# Add sample group information to the peak table
peak_data$sample_group <- peak_sample_groups


# Save peak table to CSV with group information
write.csv(peak_data, file = "detected_peaks_with_groups.csv")


#Peak Quantification
# Extract peak data including intensity values
peak_quant_table <- featureValues(filled_aligned_peaks, method = "maxint")
# Convert the peak data to a dataframe
peak_quant_df <- as.data.frame(peak_quant_table)

#Handle NA values
#sum(is.na(peak_quant_df)) = 454926
peak_quant_df[is.na(peak_quant_df)] <- 0


#Normalization
# Perform normalization using Total Ion Current (TIC)
tic_values <- apply(peak_quant_df, 2, sum)
normalized_peak_quant_df <- sweep(peak_quant_df, 2, tic_values, FUN = "/")


#Differential Expression - Statistical Analysis
library(multcomp)
library(reshape2)
# Convert the peak intensity data to a long format for statistical analysis
normalized_peak_quant_df$peak_id <- rownames(normalized_peak_quant_df)  # Adding peak IDs for reference
long_peak_data <- melt(normalized_peak_quant_df, id.vars = "peak_id", variable.name = "sample", value.name = "intensity")
long_peak_data$group <- rep(sample_groups, each = nrow(normalized_peak_quant_df))

# Perform ANOVA for each peak
# Get rid of 'peak_id' column from normalized_peak_quant_df
normalized_peak_quant_df_no_id <- normalized_peak_quant_df[, -which(colnames(normalized_peak_quant_df) == "peak_id")]

anova_results <- apply(normalized_peak_quant_df_no_id, 1, function(x) {
  temp_data <- data.frame(intensity = x, group = sample_groups)
  aov(intensity ~ group, data = temp_data)
})

#filter ANOVA results
significant_anova_results <- anova_results[sapply(anova_results, function(a) summary(a)[[1]][["Pr(>F)"]][1]) < 0.05]

# Perform Tukey's HSD test for pairwise comparisons
tukey_results <- lapply(anova_results, TukeyHSD)


# WGCNA (Weighted Gene Co-expression Network Analysis)
library(WGCNA)

# Prepare data for WGCNA (transpose to have peaks as columns and samples as rows)
datExpr <- t(normalized_peak_quant_df_no_id)  # Make sure peak IDs are removed
datExpr <- as.data.frame(datExpr)  # Convert to data frame

# Check for missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Choose the soft thresholding power for network construction
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the scale-free topology fit index
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")

# Select the appropriate soft-thresholding power
softPower <- sft$powerEstimate

# Construct the network
net <- blockwiseModules(datExpr, power = softPower, TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,
                        pamRespectsDendro = FALSE, verbose = 3)

# Convert labels to colors for visualization
moduleColors <- labels2colors(net$colors)

# Plot the dendrogram with module colors
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module Colors",
                    dendroLabels = FALSE, hang = 0.03)

# Relate modules to clinical traits (e.g., control, AD, AsymAD)
# Assume you have clinical data in a vector called 'sampleGroups'
MEs <- net$MEs  # Module eigengenes

# Create sample groups matrix for heatmap
sample_groups_factor <- factor(sample_groups, levels = c("Control", "AsymAD", "AD"))
trait_matrix <- model.matrix(~ sample_groups_factor - 1)
colnames(trait_matrix) <- c("Control", "AsymAD", "AD")

# Correlate module eigengenes (MEs) with the trait matrix
moduleTraitCor <- cor(MEs, trait_matrix, use = "p")

# Plot the heatmap with proper xLabels (traits) and yLabels (modules)
labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = colnames(trait_matrix),  # Correctly labeled traits (Control, AsymAD, AD)
               yLabels = names(MEs),              # Module eigengenes
               ySymbols = names(MEs), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50))







#Annotation
# Use MetaboAnalystR for annotation
library(MetaboAnalystR)

# extract information from peak_data for annotation
peak_data_for_annotation <- peak_data[, c("mz", "rt")]
write.csv(peak_data_for_annotation, file = "peak_data_for_annotation.csv")

# Assuming you have 'mz' and 'rt' columns for annotation
mSet <- InitDataObjects("pktable", "stat", FALSE)

# Load the peak data for annotation (m/z and rt)
mSet <- Read.TextData(mSet, "peak_data_for_annotation.csv", "rowu", "disc")

# Perform compound identification (annotation) based on m/z values
mSet <- PerformPeakMapping(mSet, "mummichog", 5)

# View identified compounds
View(mSet$peak_annotation)




#Functional Enrichment Analysis
library(clusterProfiler)
library(org.Hs.eg.db)

# Perform GO enrichment analysis
enrich_res <- enrichGO(gene = rownames(datExpr), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)

# View the results
head(enrich_res)




#Clinical Correlation
# Extract cognitive scores information from xlsx file
library(readxl)
file_path <- "/nobackup/xhnh76/dissertation/dataset/Higginbotham_MembraneProteomeManuscript_TableS1_FINAL.xlsx"
data <- read_excel(file_path, sheet = "Table S1")
sample_data <- data[6:nrow(data), c(1,9)]
sample_data <- sample_data[1:(nrow(sample_data)-5),]
sample_data <- sample_data[-c(7,8,15,16),]
# Rename columns
colnames(sample_data) <- c("Sample", "CASI")

# Convert data types
sample_data$CASI <- as.numeric(sample_data$CASI)

# Match file names for MEs
rownames_cleaned <- gsub("_.*", "", rownames(MEs))
# Calculate means of the replicates for each sample
MEs_aggregated <- aggregate(MEs, by = list(Sample = rownames_cleaned), FUN = mean)

# Match sample orders
sample_data_ordered <- sample_data[match(MEs_aggregated$Sample, sample_data$Sample), ]
MEs_aggregated_values <- MEs_aggregated[, -1]

cognitive_scores <- sample_data_ordered$

# Correlate the modules with cognitive scores
if (nrow(MEs_aggregated_values) == length(cognitive_scores)) {
  
  # Calculate the correlation between module eigengenes and cognitive scores
  module_cognitive_cor <- cor(MEs_aggregated_values, cognitive_scores, use = "pairwise.complete.obs")
  
  # View the correlation results
  print(module_cognitive_cor)
  
} else {
  stop("The number of rows in MEs_aggregated_values and the length of cognitive_scores do not match.")
}

# Plot the correlation
cor_values <- as.numeric(module_cognitive_cor)  # Convert to a numeric vector

barplot(cor_values, 
        names.arg = colnames(MEs_aggregated_values), 
        las = 2, 
        col = ifelse(cor_values > 0, "blue", "red"),
        ylab = "Correlation with Cognitive Score", 
        main = "Module-Cognitive Score Correlations",
        ylim = c(-1, 1))


#Data visualization
library(ggplot2)

# Volcano plot for differential expression
create_volcano_data <- function(group1, group2, normalized_data, sample_groups) {
  
  # Select the relevant samples for the two groups
  group_selection <- sample_groups %in% c(group1, group2)
  data_subset <- normalized_data[, group_selection]
  group_labels <- sample_groups[group_selection]
  
  # Calculate log fold change (logFC) and p-values for each peak
  volcano_data <- data.frame(
    logFC = apply(data_subset, 1, function(x) {
      mean(x[group_labels == group1]) - mean(x[group_labels == group2])
    }),
    pvalue = sapply(anova_results, function(x) summary(x)[[1]][["Pr(>F)"]][1])  # Extract p-value from ANOVA
  )
  
  # Add adjusted p-values
  volcano_data$adj_pvalue <- p.adjust(volcano_data$pvalue, method = "BH")
  
  return(volcano_data)
}

volcano_data_AD_Control <- create_volcano_data("AD", "Control", normalized_peak_quant_df_no_id, sample_groups)
volcano_data_AD_AsymAD <- create_volcano_data("AD", "AsymAD", normalized_peak_quant_df_no_id, sample_groups)
volcano_data_AsymAD_Control <- create_volcano_data("AsymAD", "Control", normalized_peak_quant_df_no_id, sample_groups)


plot_volcano <- function(volcano_data, comparison_label) {
  ggplot(volcano_data, aes(x = logFC, y = -log10(pvalue))) +
    geom_point(aes(color = adj_pvalue < 0.05)) +  # Color by significance (adjusted p-value)
    theme_minimal() +
    xlab(paste("Log Fold Change (", comparison_label, ")", sep = "")) + 
    ylab("-log10(p-value)") +
    ggtitle(paste("Volcano Plot of Differential Expression: ", comparison_label)) +
    scale_color_manual(values = c("black", "red"), guide = "none")  # Red for significant, black for non-significant
}

# Plot AD vs Control
plot_volcano(volcano_data_AD_Control, "AD vs Control")

# Plot AD vs AsymAD
plot_volcano(volcano_data_AD_AsymAD, "AD vs AsymAD")

# Plot AsymAD vs Control
plot_volcano(volcano_data_AsymAD_Control, "AsymAD vs Control")









### Annotate Manually ###

##BASH - download NCBI RefSeq file
#wget https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/complete.wp_protein.50.protein.faa.gz
#gunzip complete.wp_protein.50.protein.faa.gz
library(Biostrings)
library(dplyr)
refseq_fasta <- readAAStringSet("./complete.wp_protein.50.protein.faa")
head(refseq_fasta)

# Function to calculate the molecular weight of a protein sequence
calculate_mw <- function(sequence) {
  # Amino acid molecular weights (excluding water)
  aa_weights <- c(A = 71.08, R = 156.19, N = 114.10, D = 115.09, C = 103.15,
                  E = 129.12, Q = 128.13, G = 57.05, H = 137.14, I = 113.16,
                  L = 113.16, K = 128.17, M = 131.19, F = 147.18, P = 97.12,
                  S = 87.08, T = 101.10, W = 186.21, Y = 163.18, V = 99.13)
  
  # Split the protein sequence into individual amino acids
  aa_comp <- alphabetFrequency(AAString(sequence), baseOnly = TRUE)
  
  # Calculate molecular weight by summing the weights of individual amino acids
  mw <- sum(aa_weights[names(aa_comp)] * aa_comp)
  
  # Add the weight of water (18.015 Da) for the complete protein
  mw <- mw + 18.015
  
  return(mw)
}

# Function to match detected m/z values to molecular weights of proteins in the database
match_peaks_to_refseq <- function(peak_data, refseq_data, tolerance = 0.5) {
  # Calculate the molecular weights for each protein in RefSeq
  refseq_mw <- sapply(as.character(refseq_data), calculate_mw)
  
  # Initialize an empty list to store matched proteins
  matched_peaks <- list()
  
  for (i in 1:nrow(peak_data)) {
    mz_value <- peak_data$mz[i]
    
    # Find protein entries whose molecular weight matches within the specified tolerance
    matched_indices <- which(abs(refseq_mw - mz_value) <= tolerance)
    
    if (length(matched_indices) > 0) {
      # Store matches in a list, preserving all columns from peak_data
      matched_peaks[[i]] <- cbind(
        peak_data[i, ],
        protein_id = names(refseq_data)[matched_indices],
        protein_sequence = as.character(refseq_data[matched_indices]),
        molecular_weight = refseq_mw[matched_indices]
      )
    } else {
      # Preserve all columns from peak_data even when no match is found
      matched_peaks[[i]] <- cbind(
        peak_data[i, ],
        protein_id = NA,
        protein_sequence = NA,
        molecular_weight = NA
      )
    }
  }
  
  # Combine all matched data frames into one
  matched_peaks_df <- do.call(rbind, matched_peaks)
  
  return(matched_peaks_df)
}

# Perform annotation on the peak data
annotated_peak_data <- match_peaks_to_refseq(peak_data, refseq_fasta, tolerance = 0.5)
print(annotated_peak_data)

# Save the annotated peak data
write.csv(annotated_peak_data, file = "annotated_peaks_refseq.csv")

#########################################################################

# Run CAMERA for peak annotation (optional)
annotated_peaks <- xsAnnotate(filled_peaks)
annotated_peaks <- groupFWHM(annotated_peaks)
annotated_peaks <- findIsotopes(annotated_peaks)
annotated_peaks <- findAdducts(annotated_peaks)

# Save annotated peaks
saveRDS(annotated_peaks, file = "annotated_peaks_with_groups.rds")

