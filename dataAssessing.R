args <- commandArgs(trailingOnly = TRUE)

# Input arguments from Python
forward_file <- args[1]
reverse_file <- args[2]
reads_folder <- args[3]
plots_folder <- args[4]
results_folder <- args[5]
silva_path <- args[6]

# Load required package
library(dada2)

# Construct filtered file paths
filtFs <- file.path(reads_folder, "SRR29923448_1_F_filt.fastq")
filtRs <- file.path(reads_folder, "SRR29923448_2_R_filt.fastq")

# 1. Filtering and trimming
out <- filterAndTrim(
  c(forward_file), c(filtFs),
  c(reverse_file), c(filtRs),
  truncLen = c(240, 230),
  maxN = 0, maxEE = c(2, 2), truncQ = 2,
  rm.phix = TRUE, compress = FALSE, multithread = FALSE
)

write.table(out, file.path(results_folder, "filtering_summary.txt"), sep = "\t", quote = FALSE, col.names = NA)
print("✅ Filtering complete!")

# 2. Plot quality profiles
pdf(file.path(plots_folder, "quality_profiles_forward.pdf"))
plotQualityProfile(filtFs)
dev.off()

pdf(file.path(plots_folder, "quality_profiles_reverse.pdf"))
plotQualityProfile(filtRs)
dev.off()
print("📊 Quality profiles saved.")

# 3. Learn error rates
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
print("📈 Learned error rates.")

# 4. Denoising and merging
dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
print("🔗 Merging complete.")

# 5. Chimera removal and taxonomy assignment
seqtab <- makeSequenceTable(mergers)
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE)
print("🧬 Chimera removal complete.")

# Assign taxonomy
silva_train_set <- file.path(silva_path, "silva_nr99_v138.1_train_set.fa")
silva_species <- file.path(silva_path, "silva_species_assignment_v138.1.fa")

taxa <- assignTaxonomy(seqtab_nochim, silva_train_set, multithread = FALSE)
taxa <- addSpecies(taxa, silva_species)

# ✅ Set ASV sequences as row names in the taxonomy table
rownames(taxa) <- colnames(seqtab_nochim)

# ✅ Save as CSV only (no RDS)
write.csv(seqtab_nochim, file.path(results_folder, "seqtab_nochim.csv"), row.names = TRUE)
write.csv(taxa, file.path(results_folder, "taxa.csv"), row.names = TRUE)

print("✅ DADA2 pipeline complete — CSVs saved, no RDS files generated.")
