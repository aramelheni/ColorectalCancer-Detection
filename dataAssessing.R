# ==================== DADA2 Pipeline for Microbiome Analysis ====================
# This script processes paired-end 16S rRNA sequencing data using the DADA2 pipeline
#
# KAGGLE USAGE:
# 1. Upload this script as part of your Kaggle dataset or notebook
# 2. Upload raw FASTQ files as a Kaggle dataset
# 3. Upload SILVA reference database as a Kaggle dataset:
#    - silva_nr99_v138.1_train_set.fa
#    - silva_species_assignment_v138.1.fa
# 4. Run test.py which will call this R script with appropriate paths
#
# LOCAL USAGE:
# - Called by test.py with command line arguments
# - Requires DADA2 package (Bioconductor): installed automatically if missing
# ================================================================================

args <- commandArgs(trailingOnly = TRUE)

# Input arguments from Python
forward_file <- args[1]
reverse_file <- args[2]
reads_folder <- args[3]
plots_folder <- args[4]
results_folder <- args[5]
silva_path <- args[6]

# Install DADA2 if not already installed (Bioconductor package)
cat("📦 Checking for DADA2 package...\n")
if (!requireNamespace("dada2", quietly = TRUE)) {
  cat("⚙️  DADA2 not found. Installing from Bioconductor...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
  }
  BiocManager::install("dada2", update = FALSE, ask = FALSE)
  cat("✅ DADA2 installation complete!\n")
}

cat("📦 Loading DADA2 package...\n")
library(dada2)
cat("✅ DADA2 loaded\n\n")

# Construct filtered file paths
filtFs <- file.path(reads_folder, "SRR29923448_1_F_filt.fastq")
filtRs <- file.path(reads_folder, "SRR29923448_2_R_filt.fastq")

# 1. Filtering and trimming
cat("🔬 Step 1: Filtering and trimming reads...\n")
out <- filterAndTrim(
  c(forward_file), c(filtFs),
  c(reverse_file), c(filtRs),
  truncLen = c(240, 230),
  maxN = 0, maxEE = c(2, 2), truncQ = 2,
  rm.phix = TRUE, compress = FALSE, multithread = FALSE
)

write.table(out, file.path(results_folder, "filtering_summary.txt"), sep = "\t", quote = FALSE, col.names = NA)
cat("✅ Filtering complete!\n")
cat("   Reads in:", out[1], "| Reads out:", out[2], "\n\n")

# 2. Plot quality profiles
cat("📊 Step 2: Generating quality profiles...\n")
pdf(file.path(plots_folder, "quality_profiles_forward.pdf"))
plotQualityProfile(filtFs)
dev.off()

pdf(file.path(plots_folder, "quality_profiles_reverse.pdf"))
plotQualityProfile(filtRs)
dev.off()
cat("✅ Quality profiles saved\n\n")

# 3. Learn error rates
cat("📈 Step 3: Learning error rates...\n")
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
cat("✅ Error rates learned\n\n")

# 4. Denoising and merging
cat("🧬 Step 4: Denoising sequences...\n")
dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)

cat("🔗 Step 5: Merging paired reads...\n")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
cat("✅ Merging complete\n\n")

# Export merged sequences as FASTQ for comparison
cat("💾 Exporting merged sequences as FASTQ...\n")
dada2_output_fastq <- file.path(results_folder, "dada2_output.fastq")

# Write FASTQ with sequences and quality scores from merged data
# mergers is a data.frame for single samples
fastq_conn <- file(dada2_output_fastq, "w")
for (i in seq_len(nrow(mergers))) {
  seq_id <- paste0("seq", i)
  sequence <- mergers$sequence[i]
  # Use the quality scores from merged data (average quality)
  qual_scores <- mergers$qual[i,]
  # Convert numeric quality to ASCII (Phred+33)
  qual_string <- intToUtf8(qual_scores + 33, multiple = FALSE)
  
  # Write FASTQ format: @id, sequence, +, quality
  writeLines(paste0("@", seq_id), fastq_conn)
  writeLines(sequence, fastq_conn)
  writeLines("+", fastq_conn)
  writeLines(qual_string, fastq_conn)
}
close(fastq_conn)
cat("✅ FASTQ export complete:", dada2_output_fastq, "\n\n")

# 5. Chimera removal and taxonomy assignment
cat("🧪 Step 6: Creating sequence table...\n")
seqtab <- makeSequenceTable(mergers)
cat("   Total ASVs:", ncol(seqtab), "\n")

cat("🧹 Step 7: Removing chimeras...\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE)
cat("✅ Chimera removal complete\n")
cat("   ASVs remaining:", ncol(seqtab_nochim), "\n\n")

# Assign taxonomy
cat("🏷️  Step 8: Assigning taxonomy...\n")
silva_train_set <- file.path(silva_path, "silva_nr99_v138.1_train_set.fa")
silva_species <- file.path(silva_path, "silva_species_assignment_v138.1.fa")

cat("   Using SILVA database from:", silva_path, "\n")
taxa <- assignTaxonomy(seqtab_nochim, silva_train_set, multithread = FALSE)
cat("   Genus-level classification complete\n")
taxa <- addSpecies(taxa, silva_species)
cat("✅ Species-level classification complete\n\n")

# Set ASV sequences as row names in the taxonomy table
rownames(taxa) <- colnames(seqtab_nochim)

# Save as CSV (compatible with Python/Kaggle)
cat("💾 Step 9: Saving results...\n")
seqtab_output <- file.path(results_folder, "seqtab_nochim.csv")
taxa_output <- file.path(results_folder, "taxonomy_table_species.csv")

write.csv(seqtab_nochim, seqtab_output, row.names = TRUE)
write.csv(taxa, taxa_output, row.names = TRUE)

cat("✅ DADA2 pipeline complete!\n")
cat("📄 Output files:\n")
cat("   -", seqtab_output, "\n")
cat("   -", taxa_output, "\n")
cat("\n🎉 Ready for machine learning analysis!\n")
