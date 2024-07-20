#Save unique sample names in R

mzxml_dir <- "/Users/hyeseon/Documents/Durham Academics/2024_Dissertation/mzXML" #local path where .mzXML files are stored
mzxml_files <- list.files(mzxml_dir, pattern = "\\.mzXML$", full.names = TRUE)
samples <- sub("_[A-E]\\.mzXML$", "", basename(mzxml_files))
names <- unique(samples)
lapply(names, write, "unique_sample_names.txt", append=TRUE, ncolumns=length(names)) #save list of unique names into a text file