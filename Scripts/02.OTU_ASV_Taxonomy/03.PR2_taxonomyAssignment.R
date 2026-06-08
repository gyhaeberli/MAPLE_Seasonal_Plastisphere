# ============================================
# 18S rRNA Taxonomy Assignment Script
# MAPLE Seasonla Plastisphere data
# ============================================

# This script performs taxonomic assignment on LULU-curated 18S ASVs
# using PR2 database with both DADA2 and DECIPHER IDTAXA algorithms.
# 
# Processes BOTH nopool and pool datasets in a single run.
#


library(dada2)
library(Biostrings)
library(DECIPHER)


# ============================================
# PATH CONFIGURATION
# ============================================

BASE_DIR     <- "/home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001"
DATA_DIR     <- "/data/glennsdata/MAPLE/18S/EUK_data"

RAW_DATA_DIR   <- file.path(BASE_DIR, "01.RawData")
RWORKSPACE_DIR <- file.path(DATA_DIR, "R_workspaces")
TAXA_DB_DIR    <- file.path(RAW_DATA_DIR, "taxa_db")
TAXONOMY_PR2   <- file.path(TAXA_DB_DIR, "taxonomy_PR2")
TAXONOMY_18S   <- file.path(TAXONOMY_PR2, "taxonomy_18S")
PROCESSED_DIR  <- file.path(DATA_DIR, "Processed_data", "taxonomy_PR2")

## Create directories if they don't exist
invisible(lapply(
  c(RWORKSPACE_DIR, TAXA_DB_DIR, TAXONOMY_PR2, TAXONOMY_18S, PROCESSED_DIR),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

# ============================================
# CONFIGURATION 
# ============================================

setwd(RAW_DATA_DIR)

# LULU OTU tables

nP.lulu  = readRDS(file.path(RWORKSPACE_DIR, "nP_lulu.rds"))
psP.lulu = readRDS(file.path(RWORKSPACE_DIR, "psP_lulu.rds"))
fP.lulu  = readRDS(file.path(RWORKSPACE_DIR, "fP_lulu.rds"))

# SWARM OTU tables

load(file.path(RWORKSPACE_DIR, "swarm_results.RData"))

nP.sw = nP_result$otu_table
colnames(nP.sw) = nP_result$representative_sequences[colnames(nP.sw)]

psP.sw = psP_result$otu_table
colnames(psP.sw) = psP_result$representative_sequences[colnames(psP.sw)]

fP.sw = fP_result$otu_table
colnames(fP.sw) = fP_result$representative_sequences[colnames(fP.sw)]

# Database paths - download from https://github.com/pr2database/pr2database/releases
db_dir <- TAXA_DB_DIR

# PR2 database files (download these)
pr2_dada2_db  <- file.path(db_dir, "pr2_version_5.1.1_SSU_dada2.fasta.gz")
pr2_idtaxa_db <- file.path(db_dir, "pr2_version_5.1.0_SSU.decipher.trained.rds")


# --- Output directory ---
output_dir <- TAXONOMY_PR2
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================
# PR2 v5.x TAXONOMIC LEVELS
# ============================================
# PR2 uses 9 taxonomic levels (different from default DADA2):
PR2_TAX_LEVELS <- c("Domain", "Supergroup", "Division", "Subdivision", 
                    "Class", "Order", "Family", "Genus", "Species")


# ============================================
# STORE ALL 6 DATASETS
# ============================================

datasets <- list(
  nP_lulu  = nP.lulu,
  psP_lulu = psP.lulu,
  fP_lulu  = fP.lulu,
  nP_sw    = nP.sw,
  psP_sw   = psP.sw,
  fP_sw    = fP.sw
)


# ============================================
# GET ALL UNIQUE SEQUENCES ACROSS ALL 6 DATASETS
# ============================================
# Why do this? Running assignTaxonomy once on all unique sequences
# is much faster than running it 6 times separately. Many sequences
# are shared across pooling methods.

cat("\n========================================\n")
cat("Preparing sequences from all 6 datasets\n")
cat("========================================\n")

all_seqs_list <- lapply(datasets, colnames)

# Report per-dataset counts
for (name in names(datasets)) {
  cat(sprintf("  %-10s: %d sequences\n", name, length(all_seqs_list[[name]])))
}

# Combine and deduplicate
all_unique_seqs <- unique(unlist(all_seqs_list))
cat("\nTotal unique sequences across all datasets:", length(all_unique_seqs), "\n")

# Report sharing between LULU and SWARM
lulu_seqs  <- unique(unlist(all_seqs_list[c("nP_lulu", "psP_lulu", "fP_lulu")]))
swarm_seqs <- unique(unlist(all_seqs_list[c("nP_sw", "psP_sw", "fP_sw")]))
cat("Unique LULU sequences:", length(lulu_seqs), "\n")
cat("Unique SWARM sequences:", length(swarm_seqs), "\n")
cat("Shared between LULU and SWARM:", length(intersect(lulu_seqs, swarm_seqs)), "\n")

# Create DNAStringSet for DECIPHER
seqs_dna <- DNAStringSet(all_unique_seqs)
names(seqs_dna) <- paste0("ASV", seq_along(seqs_dna))

# Create sequence -> ASV ID mapping (useful for tracking)
seq_to_asv <- setNames(names(seqs_dna), all_unique_seqs)

cat("\nSequence length range:", min(width(seqs_dna)), "-", max(width(seqs_dna)), "bp\n")
cat("Median length:", median(width(seqs_dna)), "bp\n")

# ============================================
# STEP 1: DADA2 assignTaxonomy with PR2
# Run ONCE on all unique sequences
# ============================================

cat("\n========================================\n")
cat("STEP 1: DADA2 assignTaxonomy (PR2 v5.1.1)\n")
cat("========================================\n")

tax_dada2_all  <- NULL
boot_dada2_all <- NULL

if (!file.exists(pr2_dada2_db)) {
  stop("PR2 DADA2 database not found at: ", pr2_dada2_db,
       "\nDownload from: https://github.com/pr2database/pr2database/releases")
}

cat("Database:", basename(pr2_dada2_db), "\n")
cat("Sequences:", length(all_unique_seqs), "\n")
cat("This may take 30-90 minutes...\n\n")

start_time <- Sys.time()

taxa_dada2 <- assignTaxonomy(
  all_unique_seqs,
  refFasta      = pr2_dada2_db,
  minBoot       = 50,      # minimum bootstrap confidence to keep assignment
  tryRC         = TRUE,    # try reverse complement if no match found
  outputBootstraps = TRUE, # also return confidence scores
  taxLevels     = PR2_TAX_LEVELS,
  multithread   = TRUE,
  verbose       = TRUE
)

end_time <- Sys.time()
cat("\nDone! Time:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")

# Pull out taxonomy and bootstrap tables
tax_dada2_all  <- taxa_dada2$tax
boot_dada2_all <- taxa_dada2$boot

# Use sequences as rownames (makes subsetting per-dataset simple later)
rownames(tax_dada2_all)  <- all_unique_seqs
rownames(boot_dada2_all) <- all_unique_seqs

# Report classification success per rank
cat("\nClassification success by rank:\n")
for (rank in colnames(tax_dada2_all)) {
  n  <- sum(!is.na(tax_dada2_all[, rank]))
  pct <- round(n / nrow(tax_dada2_all) * 100, 1)
  cat(sprintf("  %-12s: %5d / %5d (%5.1f%%)\n", rank, n, nrow(tax_dada2_all), pct))
}

# Save master DADA2 result immediately (in case later steps fail)
saveRDS(list(tax = tax_dada2_all, boot = boot_dada2_all),
        file.path(output_dir, "dada2_pr2_taxonomy_ALL_sequences.rds"))
cat("\nSaved: dada2_pr2_taxonomy_ALL_sequences.rds\n")

save.image("dada2_PR2_taxonomyAssign.RData")

# ============================================
# STEP 2: DECIPHER IDTAXA with PR2
# Run ONCE on all unique sequences
# ============================================

cat("\n========================================\n")
cat("STEP 2: DECIPHER IDTAXA (PR2 v5.1.0)\n")
cat("========================================\n")

tax_idtaxa_all  <- NULL
conf_idtaxa_all <- NULL
ids_all         <- NULL

if (!file.exists(pr2_idtaxa_db)) {
  cat("WARNING: IDTAXA database not found at:", pr2_idtaxa_db, "\n")
  cat("Skipping IDTAXA. Download from: https://github.com/pr2database/pr2database/releases\n")
} else {
  
  cat("Loading trained classifier...\n")
  trainingSet <- readRDS(pr2_idtaxa_db)
  cat("Loaded PR2 DECIPHER training set\n\n")
  
  cat("Sequences:", length(seqs_dna), "\n")
  cat("This may take 10-40 minutes...\n\n")
  
  start_time <- Sys.time()
  
  ids_all <- IdTaxa(
    seqs_dna,
    trainingSet,
    type      = "extended",
    strand    = "top",
    threshold = 50,       # minimum confidence threshold (same as minBoot above)
    processors = NULL,    # NULL = use all available cores
    verbose   = TRUE
  )
  
  end_time <- Sys.time()
  cat("\nDone! Time:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")
  
  # --- Extract results ---
  # PR2 v5.1.0 IDTAXA output uses positional taxonomy (no rank labels).
  # Position 1 = Root (skip), positions 2-10 = the 9 PR2 levels.
  
  cat("Extracting IDTAXA results...\n")
  
  tax_idtaxa_all <- t(sapply(ids_all, function(x) {
    taxa <- rep(NA_character_, 9)
    n    <- length(x$taxon)
    if (n > 1) {
      available       <- min(n - 1, 9)
      taxa[1:available] <- x$taxon[2:(available + 1)]
    }
    # Remove "unclassified_" entries (treat as NA)
    taxa[!is.na(taxa) & startsWith(taxa, "unclassified_")] <- NA
    taxa
  }))
  
  conf_idtaxa_all <- t(sapply(ids_all, function(x) {
    conf <- rep(NA_real_, 9)
    n    <- length(x$confidence)
    if (n > 1) {
      available       <- min(n - 1, 9)
      conf[1:available] <- x$confidence[2:(available + 1)]
    }
    conf
  }))
  
  colnames(tax_idtaxa_all)  <- PR2_TAX_LEVELS
  colnames(conf_idtaxa_all) <- PR2_TAX_LEVELS
  rownames(tax_idtaxa_all)  <- all_unique_seqs
  rownames(conf_idtaxa_all) <- all_unique_seqs
  
  # Report classification success per rank
  cat("\nClassification success by rank:\n")
  for (rank in colnames(tax_idtaxa_all)) {
    n        <- sum(!is.na(tax_idtaxa_all[, rank]))
    pct      <- round(n / nrow(tax_idtaxa_all) * 100, 1)
    mean_conf <- round(mean(conf_idtaxa_all[!is.na(tax_idtaxa_all[, rank]), rank], na.rm = TRUE), 1)
    cat(sprintf("  %-12s: %5d / %5d (%5.1f%%) [mean confidence: %5.1f%%]\n",
                rank, n, nrow(tax_idtaxa_all), pct, mean_conf))
  }
  
  # Save master IDTAXA result
  saveRDS(list(tax = tax_idtaxa_all, conf = conf_idtaxa_all),
          file.path(output_dir, "idtaxa_pr2_taxonomy_ALL_sequences.rds"))
  cat("\nSaved: idtaxa_pr2_taxonomy_ALL_sequences.rds\n")
  
  rm(trainingSet)
  gc()
}

# ============================================
# STEP 3: Map taxonomy back to each dataset
# ============================================

cat("\n========================================\n")
cat("STEP 3: Mapping taxonomy to each dataset\n")
cat("========================================\n")

results <- list()

for (dataset_name in names(datasets)) {
  
  cat("\n--- Dataset:", dataset_name, "---\n")
  
  dataset_seqs <- colnames(datasets[[dataset_name]])
  dataset_outdir <- file.path(output_dir, dataset_name)
  dir.create(dataset_outdir, showWarnings = FALSE, recursive = TRUE)
  
  results[[dataset_name]] <- list()
  
  # Subset DADA2 results to this dataset's sequences
  if (!is.null(tax_dada2_all)) {
    tax_dada2  <- tax_dada2_all[dataset_seqs, , drop = FALSE]
    boot_dada2 <- boot_dada2_all[dataset_seqs, , drop = FALSE]
    
    results[[dataset_name]]$dada2 <- list(tax = tax_dada2, boot = boot_dada2)
    
    saveRDS(list(tax = tax_dada2, boot = boot_dada2),
            file.path(dataset_outdir, "dada2_taxonomy.rds"))
    
    cat("  Saved DADA2 results:", nrow(tax_dada2), "sequences\n")
  }
  
  # Subset IDTAXA results to this dataset's sequences
  if (!is.null(tax_idtaxa_all)) {
    tax_idtaxa  <- tax_idtaxa_all[dataset_seqs, , drop = FALSE]
    conf_idtaxa <- conf_idtaxa_all[dataset_seqs, , drop = FALSE]
    
    results[[dataset_name]]$idtaxa <- list(tax = tax_idtaxa, conf = conf_idtaxa)
    
    saveRDS(list(tax = tax_idtaxa, conf = conf_idtaxa),
            file.path(dataset_outdir, "idtaxa_taxonomy.rds"))
    
    cat("  Saved IDTAXA results:", nrow(tax_idtaxa), "sequences\n")
  }
  
  # Save the DADA2 result as the primary tax_table for phyloseq
  # (you can switch to IDTAXA or a consensus later)
  tax_for_phyloseq <- if (!is.null(tax_dada2_all)) tax_dada2 else tax_idtaxa
  results[[dataset_name]]$tax_for_phyloseq <- tax_for_phyloseq
  saveRDS(tax_for_phyloseq, file.path(dataset_outdir, "tax_table_for_phyloseq.rds"))
  cat("  Saved phyloseq tax_table\n")
}

# ============================================
# STEP 4: Algorithm comparison summary
# ============================================

if (!is.null(tax_dada2_all) && !is.null(tax_idtaxa_all)) {
  
  cat("\n========================================\n")
  cat("STEP 4: DADA2 vs IDTAXA Comparison\n")
  cat("========================================\n")
  
  for (dataset_name in names(datasets)) {
    cat("\n---", dataset_name, "---\n")
    
    td <- results[[dataset_name]]$dada2$tax
    ti <- results[[dataset_name]]$idtaxa$tax
    
    for (rank in c("Division", "Class", "Order", "Family", "Genus")) {
      d_assigned <- sum(!is.na(td[, rank]))
      i_assigned <- sum(!is.na(ti[, rank]))
      both       <- sum(!is.na(td[, rank]) & !is.na(ti[, rank]))
      n          <- nrow(td)
      cat(sprintf("  %-10s  DADA2: %4d (%4.1f%%)  IDTAXA: %4d (%4.1f%%)  Both: %4d (%4.1f%%)\n",
                  rank,
                  d_assigned, d_assigned/n*100,
                  i_assigned, i_assigned/n*100,
                  both, both/n*100))
    }
  }
}

# ============================================
# STEP 5: Save complete workspace
# ============================================

cat("\n========================================\n")
cat("STEP 5: Saving complete workspace\n")
cat("========================================\n")

saveRDS(list(
  dada2  = list(tax = tax_dada2_all, boot = boot_dada2_all),
  idtaxa = list(tax = tax_idtaxa_all, conf = conf_idtaxa_all),
  seq_to_asv = seq_to_asv,
  all_unique_seqs = all_unique_seqs
), file.path(output_dir, "all_taxonomy_master.rds"))

saveRDS(results, file.path(output_dir, "all_taxonomy_by_dataset.rds"))

save(datasets, results,
     tax_dada2_all, boot_dada2_all,
     tax_idtaxa_all, conf_idtaxa_all,
     all_unique_seqs, seq_to_asv, ids_all,
     file = file.path(output_dir, "postTaxonomy_PR2_MAPLE.RData"))

cat("\nSaved:\n")
cat("  all_taxonomy_master.rds\n")
cat("  all_taxonomy_by_dataset.rds\n")
cat("  postTaxonomy_PR2_MAPLE.RData\n")

cat("\n========================================\n")
cat("TAXONOMY ASSIGNMENT COMPLETE\n")
cat("========================================\n")

cat("\nOutput structure:\n")
cat("  ", output_dir, "/\n")
for (name in names(datasets)) {
  cat("    ├──", name, "/ → dada2_taxonomy.rds, idtaxa_taxonomy.rds, tax_table_for_phyloseq.rds\n")
}
cat("    ├── dada2_pr2_taxonomy_ALL_sequences.rds\n")
cat("    ├── idtaxa_pr2_taxonomy_ALL_sequences.rds\n")
cat("    ├── all_taxonomy_master.rds\n")
cat("    ├── all_taxonomy_by_dataset.rds\n")
cat("    └── postTaxonomy_PR2_MAPLE.RData\n")




## prepare fasta files 

load(file.path(output_dir, "postTaxonomy_PR2_MAPLE.RData"))

# Write all unique sequences to FASTA
library(Biostrings)
dna <- DNAStringSet(all_unique_seqs)
names(dna) <- paste0("ASV", seq_along(all_unique_seqs))
writeXStringSet(dna, file.path(output_dir, "all_ASVs.fasta"))




## BLAST + LCA:
# 
# 
# python3 03_BLAST_LCA_18S_v2.py \
# -i /data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001/01.RawData/taxa_db/taxonomy_PR2/all_ASVs.fasta \
# -d /data/bigdata/ncbi_nt_db/nt \
# -o /data/glennsdata/MAPLE/18S/EUK_data/X204SC25074012-Z01-F001/01.RawData/taxa_db/taxonomy_PR2/blast_lca_results.csv \
# -t 12




# ============================================
# ENSEMBLE CODE - optimizing taxa assignement
# ============================================
#
# Strategy:
#   1. BLAST ALL ASVs against NCBI (with rRNA filter)
#   2. High-confidence BLAST (≥99% identity, single species, informative name) overrides PR2
#   3. Otherwise use ensemble of DADA2 + IDTAXA + BLAST/LCA
#
# ============================================



library(dplyr)
library(tidyr)
library(seqinr)

setwd(output_dir)

# ============================================
# CONFIGURATION
# ============================================

taxonomy_dir <- TAXONOMY_18S
blast_all_file <- file.path(taxonomy_dir, "blast_lca_results.csv")
output_dir <- taxonomy_dir

PR2_RANKS <- c("Domain", "Supergroup", "Division", "Subdivision",
               "Class", "Order", "Family", "Genus", "Species")

DADA2_MIN_BOOT <- 50
IDTAXA_MIN_CONF <- 50

# BLAST override thresholds
BLAST_HIGH_IDENTITY <- 99.0    # Minimum identity for override

# Uninformative species patterns to exclude from override
UNINFORMATIVE_PATTERNS <- c("uncultured", "environmental", "unidentified", 
                            "unknown", "metagenome", "clone", "\\bsp\\.$",
                            "\\bsp\\b", "\\bcf\\.", "\\baff\\.")

# ============================================
# VALIDATION FUNCTIONS
# ============================================

validate_taxonomy_table <- function(df, method_name) {
  cat(sprintf("\n[VALIDATION] %s taxonomy table\n", method_name))
  
  cat(sprintf("  Rows: %d\n", nrow(df)))
  
  # Check for empty strings
  empty_counts <- sapply(PR2_RANKS, function(rank) {
    if (rank %in% colnames(df)) {
      sum(df[[rank]] == "", na.rm = TRUE)
    } else {
      NA
    }
  })
  
  if (any(empty_counts > 0, na.rm = TRUE)) {
    cat("  ✗ WARNING: Empty strings found:\n")
    for (rank in names(empty_counts)) {
      if (!is.na(empty_counts[rank]) && empty_counts[rank] > 0) {
        cat(sprintf("      %s: %d\n", rank, empty_counts[rank]))
      }
    }
  }
  
  # Assignment rates
  cat("  Assignment rates:\n")
  prev_rate <- 100
  for (rank in PR2_RANKS) {
    if (rank %in% colnames(df)) {
      rate <- mean(!is.na(df[[rank]]) & df[[rank]] != "") * 100
      indicator <- ifelse(rate <= prev_rate + 5, "✓", "⚠")
      cat(sprintf("    %s %-12s: %5.1f%%\n", indicator, rank, rate))
      prev_rate <- rate
    }
  }
  
  invisible(TRUE)
}

# ============================================
# LOAD DATA
# ============================================

cat("\n========================================\n")
cat("Loading Taxonomy Data\n")
cat("========================================\n")

load(file.path(TAXONOMY_PR2, "postTaxonomy_PR2_MAPLE.RData"))

cat("Loaded postTaxonomy_18S.RData\n")
cat(sprintf("  Unique sequences: %d\n", length(all_unique_seqs)))

# Load BLAST ALL results
if (!file.exists(blast_all_file)) {
  stop(paste("\nBLAST results not found:", blast_all_file,
             "\n\nRun BLAST on ALL ASVs first."))
}

blast_all <- read.csv(blast_all_file, stringsAsFactors = FALSE)
cat(sprintf("  BLAST ALL results: %d ASVs\n", nrow(blast_all)))

# Convert empty strings to NA
for (rank in PR2_RANKS) {
  if (rank %in% colnames(blast_all)) {
    blast_all[[rank]][blast_all[[rank]] == ""] <- NA
    blast_all[[rank]][grepl("^\\s*$", blast_all[[rank]])] <- NA
  }
}

# ============================================
# IDENTIFY HIGH-CONFIDENCE BLAST HITS
# ============================================

cat("\n========================================\n")
cat("Identifying High-Confidence BLAST Hits\n")
cat("========================================\n")

# Check for informative species names (exclude uncultured, environmental, etc.)
blast_all$informative_species <- !is.na(blast_all$Species) & 
  !grepl(paste(UNINFORMATIVE_PATTERNS, collapse = "|"), 
         blast_all$Species, ignore.case = TRUE)

cat(sprintf("  Informative species names: %d / %d (%.1f%%)\n",
            sum(blast_all$informative_species),
            sum(!is.na(blast_all$Species)),
            100 * sum(blast_all$informative_species) / sum(!is.na(blast_all$Species))))

# High confidence = high identity AND single taxid AND informative species
blast_all$high_confidence <- FALSE
blast_all$high_confidence[
  !is.na(blast_all$best_identity) &
    blast_all$best_identity >= BLAST_HIGH_IDENTITY &
    blast_all$n_taxids == 1 &
    blast_all$informative_species
] <- TRUE

n_high_conf <- sum(blast_all$high_confidence)
cat(sprintf("  Threshold: ≥%.1f%% identity, single species, informative name\n", BLAST_HIGH_IDENTITY))
cat(sprintf("  High-confidence hits: %d / %d (%.1f%%)\n", 
            n_high_conf, nrow(blast_all), 100 * n_high_conf / nrow(blast_all)))

# Show some examples
if (n_high_conf > 0) {
  cat("\n  Example high-confidence assignments:\n")
  examples <- blast_all %>%
    filter(high_confidence) %>%
    select(ASV_ID, best_identity, Species) %>%
    head(10)
  print(examples)
}

# Show what was filtered out
n_filtered <- sum(!is.na(blast_all$best_identity) & 
                    blast_all$best_identity >= BLAST_HIGH_IDENTITY &
                    blast_all$n_taxids == 1 &
                    !is.na(blast_all$Species) &
                    !blast_all$informative_species)
cat(sprintf("\n  Filtered out (uninformative): %d\n", n_filtered))

if (n_filtered > 0) {
  cat("  Examples of filtered uninformative species:\n")
  filtered_examples <- blast_all %>%
    filter(!is.na(best_identity) & 
             best_identity >= BLAST_HIGH_IDENTITY &
             n_taxids == 1 &
             !is.na(Species) &
             !informative_species) %>%
    select(ASV_ID, best_identity, Species) %>%
    head(5)
  print(filtered_examples)
}

# Also identify "good" BLAST hits (lower threshold for LCA contribution)
blast_all$good_hit <- !is.na(blast_all$best_identity) & blast_all$best_identity >= 80

cat(sprintf("\n  Good BLAST hits (≥80%%): %d / %d (%.1f%%)\n",
            sum(blast_all$good_hit), nrow(blast_all), 
            100 * sum(blast_all$good_hit) / nrow(blast_all)))

# ============================================
# PREPARE TAXONOMY TABLES
# ============================================

cat("\n========================================\n")
cat("Preparing Taxonomy Tables\n")
cat("========================================\n")

prepare_tax_table <- function(tax_matrix, conf_matrix = NULL, threshold = 50, 
                              method_name = "method") {
  if (!is.null(conf_matrix)) {
    tax_filtered <- tax_matrix
    tax_filtered[conf_matrix < threshold] <- NA
  } else {
    tax_filtered <- tax_matrix
  }
  
  tax_df <- as.data.frame(tax_filtered, stringsAsFactors = FALSE)
  tax_df$sequence <- rownames(tax_filtered)
  
  for (rank in PR2_RANKS) {
    if (!rank %in% colnames(tax_df)) {
      tax_df[[rank]] <- NA_character_
    }
    tax_df[[rank]][tax_df[[rank]] == ""] <- NA
  }
  
  tax_df <- tax_df[, c("sequence", PR2_RANKS)]
  cat(sprintf("  %s: %d sequences\n", method_name, nrow(tax_df)))
  return(tax_df)
}

cat("\nPreparing DADA2...\n")
tax_dada2_df <- prepare_tax_table(tax_dada2_all, boot_dada2_all, 
                                  threshold = DADA2_MIN_BOOT, method_name = "DADA2")
validate_taxonomy_table(tax_dada2_df, "DADA2")

cat("\nPreparing IDTAXA...\n")
tax_idtaxa_df <- prepare_tax_table(tax_idtaxa_all, conf_idtaxa_all, 
                                   threshold = IDTAXA_MIN_CONF, method_name = "IDTAXA")
validate_taxonomy_table(tax_idtaxa_df, "IDTAXA")

cat("\nPreparing BLAST (all ASVs)...\n")
tax_blast_df <- blast_all %>%
  select(sequence, all_of(PR2_RANKS))
validate_taxonomy_table(tax_blast_df, "BLAST")

# 
# Blast has some inconsistencies in rank that root from the way NCBI classifies. the assignment rates are not decreasing with rank!
#   -> The gaps at Supergroup and Subdivision are a known limitation of mapping NCBI→PR2
# Not necesserally an issue, since most analysis use the OTU tbale, but if this becomes an issue: 
#   
#   # When you need a complete lineage string, fill gaps explicitly
#   ensemble_tax_lineage <- ensemble_tax %>%
#   mutate(across(all_of(PR2_RANKS), ~replace_na(.x, "unclassified")))
# 

# ============================================
# COMBINE METHODS
# ============================================

cat("\n========================================\n")
cat("Combining Methods\n")
cat("========================================\n")

master_seqs <- data.frame(sequence = all_unique_seqs, stringsAsFactors = FALSE)

combined <- master_seqs %>%
  left_join(tax_dada2_df %>% 
              rename_with(~paste0(., ".dada2"), all_of(PR2_RANKS)),
            by = "sequence") %>%
  left_join(tax_idtaxa_df %>% 
              rename_with(~paste0(., ".idtaxa"), all_of(PR2_RANKS)),
            by = "sequence") %>%
  left_join(tax_blast_df %>%
              rename_with(~paste0(., ".blast"), all_of(PR2_RANKS)),
            by = "sequence")

# Add BLAST metadata
combined <- combined %>%
  left_join(blast_all %>% select(sequence, best_identity, n_taxids, high_confidence, informative_species),
            by = "sequence")

cat(sprintf("Combined: %d sequences x %d columns\n", nrow(combined), ncol(combined)))

# ============================================
# ENSEMBLE WITH BLAST OVERRIDE
# ============================================

cat("\n========================================\n")
cat("Computing Ensemble with BLAST Override\n")
cat("========================================\n")

compute_ensemble_with_override <- function(combined_df) {
  
  ensemble_tax <- data.frame(sequence = combined_df$sequence, stringsAsFactors = FALSE)
  ensemble_tax$assignment_source <- NA_character_
  
  for (rank in PR2_RANKS) {
    ensemble_tax[[rank]] <- NA_character_
  }
  
  for (i in 1:nrow(combined_df)) {
    row <- combined_df[i, ]
    
    # RULE 1: High-confidence BLAST (with informative species) overrides everything
    if (!is.na(row$high_confidence) && row$high_confidence == TRUE) {
      for (rank in PR2_RANKS) {
        blast_val <- row[[paste0(rank, ".blast")]]
        if (!is.na(blast_val) && blast_val != "") {
          ensemble_tax[[rank]][i] <- blast_val
        }
      }
      ensemble_tax$assignment_source[i] <- "BLAST_override"
      next
    }
    
    # RULE 2: For other sequences, use consensus
    for (rank in PR2_RANKS) {
      vals <- c(
        row[[paste0(rank, ".dada2")]],
        row[[paste0(rank, ".idtaxa")]],
        row[[paste0(rank, ".blast")]]
      )
      vals <- vals[!is.na(vals) & vals != ""]
      
      if (length(vals) == 0) {
        next
      }
      
      # Majority vote
      tab <- table(vals)
      if (max(tab) >= 2) {
        ensemble_tax[[rank]][i] <- names(tab)[which.max(tab)]
        if (is.na(ensemble_tax$assignment_source[i])) {
          ensemble_tax$assignment_source[i] <- "consensus"
        }
      } else {
        # Tiebreaker: IDTAXA > BLAST > DADA2 (IDTAXA most conservative for protists)
        # But if BLAST has good identity, prefer BLAST for non-protists
        blast_val <- row[[paste0(rank, ".blast")]]
        idtaxa_val <- row[[paste0(rank, ".idtaxa")]]
        dada2_val <- row[[paste0(rank, ".dada2")]]
        
        # Check if likely non-protist (Metazoa, Fungi, Viridiplantae)
        supergroup <- coalesce(row$Supergroup.blast, row$Supergroup.dada2, row$Supergroup.idtaxa)
        is_metazoan <- !is.na(supergroup) && supergroup %in% c("Obazoa", "Metazoa", "Fungi", "Archaeplastida", "Viridiplantae")
        
        # Also check if BLAST species is informative before preferring it
        blast_informative <- !is.na(row$informative_species) && row$informative_species == TRUE
        
        if (is_metazoan && !is.na(row$best_identity) && row$best_identity >= 95 && 
            !is.na(blast_val) && blast_val != "" && blast_informative) {
          # Prefer BLAST for metazoans with good hits AND informative species
          ensemble_tax[[rank]][i] <- blast_val
          if (is.na(ensemble_tax$assignment_source[i])) {
            ensemble_tax$assignment_source[i] <- "BLAST_preferred"
          }
        } else if (!is.na(idtaxa_val) && idtaxa_val != "") {
          ensemble_tax[[rank]][i] <- idtaxa_val
          if (is.na(ensemble_tax$assignment_source[i])) {
            ensemble_tax$assignment_source[i] <- "IDTAXA_tiebreak"
          }
        } else if (!is.na(blast_val) && blast_val != "") {
          ensemble_tax[[rank]][i] <- blast_val
          if (is.na(ensemble_tax$assignment_source[i])) {
            ensemble_tax$assignment_source[i] <- "BLAST_tiebreak"
          }
        } else if (!is.na(dada2_val) && dada2_val != "") {
          ensemble_tax[[rank]][i] <- dada2_val
          if (is.na(ensemble_tax$assignment_source[i])) {
            ensemble_tax$assignment_source[i] <- "DADA2_tiebreak"
          }
        }
      }
    }
    
    if (is.na(ensemble_tax$assignment_source[i])) {
      ensemble_tax$assignment_source[i] <- "unassigned"
    }
  }
  
  return(ensemble_tax)
}

cat("Processing...\n")
ensemble_tax <- compute_ensemble_with_override(combined)

# Summary of assignment sources
cat("\nAssignment sources:\n")
print(table(ensemble_tax$assignment_source, useNA = "ifany"))

validate_taxonomy_table(ensemble_tax, "ENSEMBLE")

# ============================================
# AGREEMENT STATISTICS
# ============================================

cat("\n========================================\n")
cat("Method Agreement\n")
cat("========================================\n")

calc_agreement <- function(combined_df, rank, method1, method2) {
  col1 <- paste0(rank, ".", method1)
  col2 <- paste0(rank, ".", method2)
  
  if (!all(c(col1, col2) %in% colnames(combined_df))) return(NA)
  
  v1 <- combined_df[[col1]]
  v2 <- combined_df[[col2]]
  
  both_assigned <- !is.na(v1) & !is.na(v2) & v1 != "" & v2 != ""
  if (sum(both_assigned) == 0) return(NA)
  
  agreement <- sum(v1[both_assigned] == v2[both_assigned]) / sum(both_assigned)
  return(round(agreement * 100, 1))
}

cat("\nDADA2 vs IDTAXA:\n")
for (rank in PR2_RANKS) {
  agreement <- calc_agreement(combined, rank, "dada2", "idtaxa")
  if (!is.na(agreement)) cat(sprintf("  %-12s: %5.1f%%\n", rank, agreement))
}

cat("\nDADA2 vs BLAST:\n")
for (rank in PR2_RANKS) {
  agreement <- calc_agreement(combined, rank, "dada2", "blast")
  if (!is.na(agreement)) cat(sprintf("  %-12s: %5.1f%%\n", rank, agreement))
}

cat("\nIDTAXA vs BLAST:\n")
for (rank in PR2_RANKS) {
  agreement <- calc_agreement(combined, rank, "idtaxa", "blast")
  if (!is.na(agreement)) cat(sprintf("  %-12s: %5.1f%%\n", rank, agreement))
}

# ============================================
# SUMMARY TABLE
# ============================================

cat("\n========================================\n")
cat("Classification Summary\n")
cat("========================================\n")

cat(sprintf("\n%-12s | %7s | %7s | %7s | %10s\n", 
            "Rank", "DADA2", "IDTAXA", "BLAST", "Ensemble"))
cat(paste(rep("-", 55), collapse = ""), "\n")

for (rank in PR2_RANKS) {
  dada2_pct <- round(mean(!is.na(combined[[paste0(rank, ".dada2")]]) & 
                            combined[[paste0(rank, ".dada2")]] != "") * 100, 1)
  idtaxa_pct <- round(mean(!is.na(combined[[paste0(rank, ".idtaxa")]]) & 
                             combined[[paste0(rank, ".idtaxa")]] != "") * 100, 1)
  blast_pct <- round(mean(!is.na(combined[[paste0(rank, ".blast")]]) & 
                            combined[[paste0(rank, ".blast")]] != "") * 100, 1)
  ensemble_pct <- round(mean(!is.na(ensemble_tax[[rank]]) & 
                               ensemble_tax[[rank]] != "") * 100, 1)
  
  cat(sprintf("%-12s | %6.1f%% | %6.1f%% | %6.1f%% | %9.1f%%\n", 
              rank, dada2_pct, idtaxa_pct, blast_pct, ensemble_pct))
}

# ============================================
# MAP TO DATASETS
# ============================================

cat("\n========================================\n")
cat("Mapping to Datasets\n")
cat("========================================\n")

ensemble_matrix <- as.matrix(ensemble_tax[, PR2_RANKS])
rownames(ensemble_matrix) <- ensemble_tax$sequence

# Also keep assignment source
assignment_source <- ensemble_tax$assignment_source
names(assignment_source) <- ensemble_tax$sequence

ensemble_results <- list()

for (dataset_name in names(datasets)) {
  cat(sprintf("\n--- %s ---\n", toupper(dataset_name)))
  
  seqtab <- datasets[[dataset_name]]
  dataset_seqs <- colnames(seqtab)
  
  tax_ensemble <- ensemble_matrix[dataset_seqs, , drop = FALSE]
  source_vector <- assignment_source[dataset_seqs]
  
  tax_ensemble_asv <- tax_ensemble
  rownames(tax_ensemble_asv) <- paste0("ASV", seq_len(nrow(tax_ensemble)))
  
  ensemble_results[[dataset_name]] <- list(
    tax = tax_ensemble,
    tax_asv = tax_ensemble_asv,
    source = source_vector
  )
  
  dataset_outdir <- file.path(output_dir, dataset_name)
  if (!dir.exists(dataset_outdir)) dir.create(dataset_outdir, recursive = TRUE)
  
  saveRDS(tax_ensemble, file.path(dataset_outdir, "ensemble_taxonomy_v2.rds"))
  write.csv(tax_ensemble_asv, file.path(dataset_outdir, "ensemble_taxonomy_v2_assignments.csv"))
  saveRDS(tax_ensemble, file.path(dataset_outdir, "tax_table_for_phyloseq.rds"))
  
  cat("  ✓ Saved\n")
  
  cat("  Classification:\n")
  for (rank in PR2_RANKS) {
    n_assigned <- sum(!is.na(tax_ensemble[, rank]) & tax_ensemble[, rank] != "")
    pct <- round(n_assigned / nrow(tax_ensemble) * 100, 1)
    cat(sprintf("    %-12s: %4d / %4d (%5.1f%%)\n", 
                rank, n_assigned, nrow(tax_ensemble), pct))
  }
  
  cat("\n  Assignment sources:\n")
  print(table(source_vector, useNA = "ifany"))
}

# ============================================
# SAVE RESULTS
# ============================================

cat("\n========================================\n")
cat("Saving Results\n")
cat("========================================\n")

saveRDS(ensemble_matrix, file.path(output_dir, "ensemble_taxonomy_v2_master.rds"))
write.csv(ensemble_tax, file.path(output_dir, "ensemble_taxonomy_v2_all.csv"), row.names = FALSE)
cat("✓ ensemble_taxonomy_v2_master.rds\n")
cat("✓ ensemble_taxonomy_v2_all.csv\n")

write.csv(combined, file.path(output_dir, "all_methods_combined_v2.csv"), row.names = FALSE)
cat("✓ all_methods_combined_v2.csv\n")

results$ensemble_v2 <- ensemble_results
saveRDS(results, file.path(output_dir, "all_taxonomy_by_dataset_v2.rds"))
cat("✓ all_taxonomy_by_dataset_v2.rds\n")

save(datasets, results, ensemble_results, ensemble_matrix, combined, ensemble_tax,
     tax_dada2_all, boot_dada2_all, tax_idtaxa_all, conf_idtaxa_all,
     blast_all, all_unique_seqs, seq_to_asv,
     file = file.path(output_dir, "postTaxonomy_18S_ensemble_v2.RData"))
cat("✓ postTaxonomy_18S_ensemble_v2.RData\n")

# ============================================
# CREATE PHYLOSEQ OBJECTS
# ============================================

# PHYLOSEQ CREATED AFTER TAX CLEANING

# 
# if (requireNamespace("phyloseq", quietly = TRUE)) {
#   
#   cat("\n========================================\n")
#   cat("Creating Phyloseq Objects\n")
#   cat("========================================\n")
#   
#   library(phyloseq)
#   
#   for (dataset_name in names(datasets)) {
#     cat(sprintf("\n%s:\n", toupper(dataset_name)))
#     
#     seqtab <- datasets[[dataset_name]]
#     tax <- ensemble_results[[dataset_name]]$tax
#     
#     seqs <- colnames(seqtab)
#     asv_names <- paste0("ASV", seq_along(seqs))
#     
#     seqtab_renamed <- seqtab
#     colnames(seqtab_renamed) <- asv_names
#     
#     tax_renamed <- tax
#     rownames(tax_renamed) <- asv_names
#     
#     otu <- otu_table(seqtab_renamed, taxa_are_rows = FALSE)
#     tax_ps <- tax_table(tax_renamed)
#     
#     ps <- phyloseq(otu, tax_ps)
#     
#     dna_seqs <- Biostrings::DNAStringSet(seqs)
#     names(dna_seqs) <- asv_names
#     ps <- merge_phyloseq(ps, refseq(dna_seqs))
#     
#     saveRDS(ps, file.path(output_dir, dataset_name, 
#                           paste0("phyloseq_", dataset_name, "_v2.rds")))
#     cat(sprintf("  ✓ Saved: %d samples, %d ASVs\n", nsamples(ps), ntaxa(ps)))
#   }
# }


cat("\n========================================\n")
cat("COMPLETE\n")
cat("========================================\n")
cat("  - BLAST run on ALL ASVs with rRNA filter\n")
cat("  - Identity-aware LCA (species only assigned with 100% identity)\n")
cat("  - Uninformative species filtered from override (uncultured, environmental, etc.)\n")
cat("  - High-confidence BLAST (≥99%, single species, informative) overrides PR2\n")
cat("  - Metazoans with good BLAST (≥95%) prefer BLAST over IDTAXA\n")
cat("  - Assignment source tracked for each ASV\n")




## Post summary
# in Astbury

library(dplyr)

load(file.path(PROCESSED_DIR, "postTaxonomy_18S_ensemble_v2.RData"))

PR2_RANKS <- c("Domain", "Supergroup", "Division", "Subdivision", 
               "Class", "Order", "Family", "Genus", "Species")

blast_all <- read.csv(file.path(TAXONOMY_18S, "blast_lca_results.csv"),
                      stringsAsFactors = FALSE)

summary_table <- data.frame(Rank = PR2_RANKS)

for (rank in PR2_RANKS) {
  
  # DADA2: use bootstrap threshold of 50
  dada2_assigned <- !is.na(tax_dada2_all[, rank]) & 
    boot_dada2_all[, rank] >= 50
  summary_table$DADA2[summary_table$Rank == rank] <- 
    round(mean(dada2_assigned) * 100, 1)
  
  # IDTAXA: use confidence threshold of 50
  idtaxa_assigned <- !is.na(tax_idtaxa_all[, rank]) & 
    conf_idtaxa_all[, rank] >= 50
  summary_table$IDTAXA[summary_table$Rank == rank] <- 
    round(mean(idtaxa_assigned) * 100, 1)
  
  # BLAST
  if (rank %in% colnames(blast_all)) {
    blast_assigned <- !is.na(blast_all[[rank]]) & blast_all[[rank]] != ""
    summary_table$BLAST[summary_table$Rank == rank] <- 
      round(mean(blast_assigned) * 100, 1)
  } else {
    summary_table$BLAST[summary_table$Rank == rank] <- NA
  }
}

print(summary_table)


# Assignement %

ensemble_fP_lulu <- readRDS(file.path(TAXONOMY_18S, "fP_lulu", "ensemble_taxonomy_v2.rds"))

# Assignment rate per rank
for (rank in PR2_RANKS) {
  n_assigned <- sum(!is.na(ensemble_fP_lulu[, rank]) & ensemble_fP_lulu[, rank] != "")
  pct <- round(n_assigned / nrow(ensemble_fP_lulu) * 100, 1)
  cat(sprintf("  %-12s: %5.1f%%\n", rank, pct))
}

# Overall: ASVs assigned at least to Domain level
n_any <- sum(!is.na(ensemble_fP_lulu[, "Domain"]) & ensemble_fP_lulu[, "Domain"] != "")
cat(sprintf("\nASVs assigned at any level: %d / %d (%.1f%%)\n", 
            n_any, nrow(ensemble_fP_lulu), n_any / nrow(ensemble_fP_lulu) * 100))

