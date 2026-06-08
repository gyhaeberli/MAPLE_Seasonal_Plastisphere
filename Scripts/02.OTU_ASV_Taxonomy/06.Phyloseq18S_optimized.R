################################################################################
############## PHYLOSEQ PIPELINE - v4
############# reorganized from v3.1 (optimized by ClaudeAI Sonnet 4.6)
################################################################################
#
# PART 1 — BUILD RAW OBJECTS
#   1.1  Load metadata, OTU tables, and taxonomy
#   1.2  Build master ASV ID map
#   1.3  Create 18 raw global phyloseq objects and save
#   1.4  Split each object by experiment and save
#   1.5  Assemble and save raw master list
#
# PART 2 — CLEAN AND FINALIZE
#   2.1  Check experimental balance (before cleaning)
#   2.2  Identify and remove invalid samples
#   2.3  Remove duplicate and zero-read samples
#   2.x  Quality filtering (read counts, OTU counts per experimental design)
#   2.4  Remove contaminant taxa
#   2.5  Check experimental balance (after cleaning)
#   2.6  QC overview across all 18 cleaned objects
#   2.7  Save all 18 cleaned global objects
#   2.8  Split cleaned objects by experiment and save
#   2.9  Create PR2 no-metazoan versions (ST + LT only) and save
#   2.10 Assemble and save final master list
#
################################################################################

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(readr)
library(dada2)
library(dplyr)
library(patchwork)
library(DECIPHER)
library(phangorn)
library(Biostrings)
library(purrr)

################################################################################
########################## PATH CONFIGURATION ##################################
################################################################################

GITHUB_DIR   <- "~/GitHub/MAPLE_Seasonal_Plastisphere"
SCRIPTS_DIR  <- file.path(GITHUB_DIR, "Scripts", "02.OTU_ASV_Taxonomy")
PROC_DIR     <- file.path(GITHUB_DIR, "Processed_data")
RAW_DIR      <- file.path(GITHUB_DIR, "Raw_data")

OTU_FILT_DIR <- file.path(PROC_DIR, "OTU_tables", "filtered")
TAXONOMY_DIR <- file.path(PROC_DIR, "Taxonomy")
METADATA_DIR <- file.path(RAW_DIR, "Metadata")
ps_out       <- file.path(PROC_DIR, "Phyloseq_objects")

## Create directories if they don't exist
invisible(lapply(
  c(ps_out),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

setwd(SCRIPTS_DIR)

################################################################################
################################################################################
##                                                                            ##
##                       PART 1 — BUILD RAW OBJECTS                          ##
##                                                                            ##
################################################################################
################################################################################

################################################################################
# 1.1  LOAD INPUT DATA
################################################################################

# ---- Metadata ----

metadata_list <- list(
  nP.lulu   = read_rds(file.path(METADATA_DIR, "nP.lulu_metadata.rds")),
  psP.lulu  = read_rds(file.path(METADATA_DIR, "psP.lulu_metadata.rds")),
  fP.lulu   = read_rds(file.path(METADATA_DIR, "fP.lulu_metadata.rds")),
  nP.swarm  = read_rds(file.path(METADATA_DIR, "nP.swarm_metadata.rds")),
  psP.swarm = read_rds(file.path(METADATA_DIR, "psP.swarm_metadata.rds")),
  fP.swarm  = read_rds(file.path(METADATA_DIR, "fP.swarm_metadata.rds"))
)

compare_metadata <- function(metadata_list) {
  sample_comparison <- lapply(metadata_list, function(x) sort(rownames(x)))
  all_same <- length(unique(sample_comparison)) == 1
  cat("All methods have same samples:", all_same, "\n\n")
  
  if (!all_same) {
    all_samples     <- unique(unlist(sample_comparison))
    presence_matrix <- sapply(metadata_list, function(meta) all_samples %in% rownames(meta))
    rownames(presence_matrix) <- all_samples
    cat("Samples NOT in all methods:\n")
    print(all_samples[rowSums(presence_matrix) < ncol(presence_matrix)])
  }
  
  key_cols <- c("site", "substrate", "season", "experiment", "bio_replicate")
  for (col in key_cols) {
    if (all(sapply(metadata_list, function(x) col %in% colnames(x)))) {
      first_meta       <- metadata_list[[1]][, col]
      identical_values <- all(sapply(metadata_list[-1], function(meta) {
        identical(meta[rownames(metadata_list[[1]]), col], first_meta)
      }))
      cat(col, "identical across methods:", identical_values, "\n")
    }
  }
  
  cat("\nRead counts/OTU richness (should differ by method):\n")
  depth_summary <- do.call(rbind, lapply(names(metadata_list), function(method) {
    data.frame(
      method       = method,
      mean_depth   = mean(metadata_list[[method]]$totseq, na.rm = TRUE),
      median_depth = median(metadata_list[[method]]$totseq, na.rm = TRUE),
      mean_OTUs    = mean(metadata_list[[method]]$OTUs, na.rm = TRUE)
    )
  }))
  print(depth_summary)
  return(invisible(list(all_same_samples = all_same, depth_summary = depth_summary)))
}

metadata_comparison <- compare_metadata(metadata_list)

collapse_metadata <- function(meta_df) {
  meta_df$sample_base <- sub("\\_rep[0-9]+$", "", rownames(meta_df))
  meta_collapsed      <- meta_df[!duplicated(meta_df$sample_base), ]
  rownames(meta_collapsed) <- meta_collapsed$sample_base
  meta_collapsed$sample_base <- NULL
  return(meta_collapsed)
}

metadata_collapsed_list <- lapply(metadata_list, collapse_metadata)

# ---- OTU tables ----

# Uncollapsed
otu_rep50_list <- list(
  nP.lulu   = read_rds(file.path(OTU_FILT_DIR, "LULU_nP_repfilt50.rds")),
  psP.lulu  = read_rds(file.path(OTU_FILT_DIR, "LULU_psP_repfilt50.rds")),
  fP.lulu   = read_rds(file.path(OTU_FILT_DIR, "LULU_fP_repfilt50.rds")),
  nP.swarm  = read_rds(file.path(OTU_FILT_DIR, "SWARM_nP_repfilt50.rds")),
  psP.swarm = read_rds(file.path(OTU_FILT_DIR, "SWARM_psP_repfilt50.rds")),
  fP.swarm  = read_rds(file.path(OTU_FILT_DIR, "SWARM_fP_repfilt50.rds")))

otu_rep75_list <- list(
  nP.lulu   = read_rds(file.path(OTU_FILT_DIR, "LULU_nP_repfilt75.rds")),
  psP.lulu  = read_rds(file.path(OTU_FILT_DIR, "LULU_psP_repfilt75.rds")),
  fP.lulu   = read_rds(file.path(OTU_FILT_DIR, "LULU_fP_repfilt75.rds")),
  nP.swarm  = read_rds(file.path(OTU_FILT_DIR, "SWARM_nP_repfilt75.rds")),
  psP.swarm = read_rds(file.path(OTU_FILT_DIR, "SWARM_psP_repfilt75.rds")),
  fP.swarm  = read_rds(file.path(OTU_FILT_DIR, "SWARM_fP_repfilt75.rds")))

otu_norep_list <- list(
  nP.lulu   = read_rds(file.path(OTU_FILT_DIR, "LULU_nP_norepfilt.rds")),
  psP.lulu  = read_rds(file.path(OTU_FILT_DIR, "LULU_psP_norepfilt.rds")),
  fP.lulu   = read_rds(file.path(OTU_FILT_DIR, "LULU_fP_norepfilt.rds")),
  nP.swarm  = read_rds(file.path(OTU_FILT_DIR, "SWARM_nP_norepfilt.rds")),
  psP.swarm = read_rds(file.path(OTU_FILT_DIR, "SWARM_psP_norepfilt.rds")),
  fP.swarm  = read_rds(file.path(OTU_FILT_DIR, "SWARM_fP_norepfilt.rds")))

# Collapsed
otu_rep50_collap_list <- list(
  nP.lulu   = read_rds(file.path(OTU_FILT_DIR, "LULU_nP_CollapsedRepFilt50.rds")),
  psP.lulu  = read_rds(file.path(OTU_FILT_DIR, "LULU_psP_CollapsedRepFilt50.rds")),
  fP.lulu   = read_rds(file.path(OTU_FILT_DIR, "LULU_fP_CollapsedRepFilt50.rds")),
  nP.swarm  = read_rds(file.path(OTU_FILT_DIR, "SWARM_nP_CollapsedRepFilt50.rds")),
  psP.swarm = read_rds(file.path(OTU_FILT_DIR, "SWARM_psP_CollapsedRepFilt50.rds")),
  fP.swarm  = read_rds(file.path(OTU_FILT_DIR, "SWARM_fP_CollapsedRepFilt50.rds")))

otu_rep75_collap_list <- list(
  nP.lulu   = read_rds(file.path(OTU_FILT_DIR, "LULU_nP_CollapsedRepFilt75.rds")),
  psP.lulu  = read_rds(file.path(OTU_FILT_DIR, "LULU_psP_CollapsedRepFilt75.rds")),
  fP.lulu   = read_rds(file.path(OTU_FILT_DIR, "LULU_fP_CollapsedRepFilt75.rds")),
  nP.swarm  = read_rds(file.path(OTU_FILT_DIR, "SWARM_nP_CollapsedRepFilt75.rds")),
  psP.swarm = read_rds(file.path(OTU_FILT_DIR, "SWARM_psP_CollapsedRepFilt75.rds")),
  fP.swarm  = read_rds(file.path(OTU_FILT_DIR, "SWARM_fP_CollapsedRepFilt75.rds")))

otu_norep_collap_list <- list(
  nP.lulu   = read_rds(file.path(OTU_FILT_DIR, "LULU_nP_collap_norepfilt.rds")),
  psP.lulu  = read_rds(file.path(OTU_FILT_DIR, "LULU_psP_collap_norepfilt.rds")),
  fP.lulu   = read_rds(file.path(OTU_FILT_DIR, "LULU_fP_collap_norepfilt.rds")),
  nP.swarm  = read_rds(file.path(OTU_FILT_DIR, "SWARM_nP_collap_norepfilt.rds")),
  psP.swarm = read_rds(file.path(OTU_FILT_DIR, "SWARM_psP_collap_norepfilt.rds")),
  fP.swarm  = read_rds(file.path(OTU_FILT_DIR, "SWARM_fP_collap_norepfilt.rds")))

# ---- Taxonomy ----

rename_tax_list <- function(rds_path, pattern, replacement) {
  obj        <- readRDS(rds_path)
  names(obj) <- sub(pattern, replacement, names(obj))
  obj
}

silva_taxonomy_norep_list <- c(
  rename_tax_list(file.path(TAXONOMY_DIR, "LULU_silva_taxonomy_filtered_norepfilt.rds"),  "silva132_lulu_(.*)", "\\1.lulu"),
  rename_tax_list(file.path(TAXONOMY_DIR, "SWARM_silva_taxonomy_filtered_norepfilt.rds"), "silva132_sw_(.*)",   "\\1.swarm"))

silva_taxonomy_rep50_list <- c(
  rename_tax_list(file.path(TAXONOMY_DIR, "LULU_silva_taxonomy_filtered_rep50.rds"),  "silva132_lulu_(.*)", "\\1.lulu"),
  rename_tax_list(file.path(TAXONOMY_DIR, "SWARM_silva_taxonomy_filtered_rep50.rds"), "silva132_sw_(.*)",   "\\1.swarm"))

silva_taxonomy_rep75_list <- c(
  rename_tax_list(file.path(TAXONOMY_DIR, "LULU_silva_taxonomy_filtered_rep75.rds"),  "silva132_lulu_(.*)", "\\1.lulu"),
  rename_tax_list(file.path(TAXONOMY_DIR, "SWARM_silva_taxonomy_filtered_rep75.rds"), "silva132_sw_(.*)",   "\\1.swarm"))

silva_taxonomy_collap_norep_list <- c(
  rename_tax_list(file.path(TAXONOMY_DIR, "LULU_silva_taxonomy_filtered_collapsed_norepfilt.rds"),  "silva132_lulu_(.*)", "\\1.lulu"),
  rename_tax_list(file.path(TAXONOMY_DIR, "SWARM_silva_taxonomy_filtered_collapsed_norepfilt.rds"), "silva132_sw_(.*)",   "\\1.swarm"))

silva_taxonomy_collap_rep50_list <- c(
  rename_tax_list(file.path(TAXONOMY_DIR, "LULU_silva_taxonomy_filtered_collapsed_rep50.rds"),  "silva132_lulu_(.*)", "\\1.lulu"),
  rename_tax_list(file.path(TAXONOMY_DIR, "SWARM_silva_taxonomy_filtered_collapsed_rep50.rds"), "silva132_sw_(.*)",   "\\1.swarm"))

silva_taxonomy_collap_rep75_list <- c(
  rename_tax_list(file.path(TAXONOMY_DIR, "LULU_silva_taxonomy_filtered_collapsed_rep75.rds"),  "silva132_lulu_(.*)", "\\1.lulu"),
  rename_tax_list(file.path(TAXONOMY_DIR, "SWARM_silva_taxonomy_filtered_collapsed_rep75.rds"), "silva132_sw_(.*)",   "\\1.swarm"))

pr2_taxonomy_norep_list <- c(
  rename_tax_list(file.path(TAXONOMY_DIR, "LULU_pr2_taxonomy_filtered_norepfilt.rds"),  "pr2_lulu_(.*)", "\\1.lulu"),
  rename_tax_list(file.path(TAXONOMY_DIR, "SWARM_pr2_taxonomy_filtered_norepfilt.rds"), "pr2_sw_(.*)",   "\\1.swarm"))

pr2_taxonomy_rep50_list <- c(
  rename_tax_list(file.path(TAXONOMY_DIR, "LULU_pr2_taxonomy_filtered_rep50.rds"),  "pr2_lulu_(.*)", "\\1.lulu"),
  rename_tax_list(file.path(TAXONOMY_DIR, "SWARM_pr2_taxonomy_filtered_rep50.rds"), "pr2_sw_(.*)",   "\\1.swarm"))

pr2_taxonomy_rep75_list <- c(
  rename_tax_list(file.path(TAXONOMY_DIR, "LULU_pr2_taxonomy_filtered_rep75.rds"),  "pr2_lulu_(.*)", "\\1.lulu"),
  rename_tax_list(file.path(TAXONOMY_DIR, "SWARM_pr2_taxonomy_filtered_rep75.rds"), "pr2_sw_(.*)",   "\\1.swarm"))

pr2_taxonomy_collap_norep_list <- c(
  rename_tax_list(file.path(TAXONOMY_DIR, "LULU_pr2_taxonomy_filtered_collapsed_norepfilt.rds"),  "pr2_lulu_(.*)", "\\1.lulu"),
  rename_tax_list(file.path(TAXONOMY_DIR, "SWARM_pr2_taxonomy_filtered_collapsed_norepfilt.rds"), "pr2_sw_(.*)",   "\\1.swarm"))

pr2_taxonomy_collap_rep50_list <- c(
  rename_tax_list(file.path(TAXONOMY_DIR, "LULU_pr2_taxonomy_filtered_collapsed_rep50.rds"),  "pr2_lulu_(.*)", "\\1.lulu"),
  rename_tax_list(file.path(TAXONOMY_DIR, "SWARM_pr2_taxonomy_filtered_collapsed_rep50.rds"), "pr2_sw_(.*)",   "\\1.swarm"))

pr2_taxonomy_collap_rep75_list <- c(
  rename_tax_list(file.path(TAXONOMY_DIR, "LULU_pr2_taxonomy_filtered_collapsed_rep75.rds"),  "pr2_lulu_(.*)", "\\1.lulu"),
  rename_tax_list(file.path(TAXONOMY_DIR, "SWARM_pr2_taxonomy_filtered_collapsed_rep75.rds"), "pr2_sw_(.*)",   "\\1.swarm"))


################################################################################
# 1.2  MASTER ASV ID MAP
# Built from the norep OTU tables (the largest possible sequence set).
# All rep-filtered sets are subsets, so every sequence gets a stable unique ID.
################################################################################

all_seqs          <- unique(unlist(lapply(otu_norep_list, colnames)))
all_seqs_sorted   <- sort(all_seqs)
master_seq_to_asv <- setNames(paste0("ASV_", seq_along(all_seqs_sorted)), all_seqs_sorted)
cat(sprintf("Master ASV map: %d unique sequences\n", length(master_seq_to_asv)))

################################################################################
# 1.3  OBJECT SPECIFICATION TABLE + BUILD FUNCTIONS
################################################################################

phylo_specs <- list(
  list(var  = "phyloGlobal_silva_norep",
       otu  = "otu_norep_list",        tax  = "silva_taxonomy_norep_list",
       meta = "metadata_list",         prefix = "SILVA - All OTUs",
       file = "phyloGlobal_silva_norep.rds"),
  list(var  = "phyloGlobal_silva_rep50",
       otu  = "otu_rep50_list",        tax  = "silva_taxonomy_rep50_list",
       meta = "metadata_list",         prefix = "SILVA - in 2+ rep",
       file = "phyloGlobal_silva_rep50.rds"),
  list(var  = "phyloGlobal_silva_rep75",
       otu  = "otu_rep75_list",        tax  = "silva_taxonomy_rep75_list",
       meta = "metadata_list",         prefix = "SILVA - in 3+ rep",
       file = "phyloGlobal_silva_rep75.rds"),
  list(var  = "phyloGlobal_silva_collap_norep",
       otu  = "otu_norep_collap_list", tax  = "silva_taxonomy_collap_norep_list",
       meta = "metadata_collapsed_list", prefix = "SILVA - All OTUs (collapsed)",
       file = "phyloGlobal_silva_collap_norep.rds"),
  list(var  = "phyloGlobal_silva_collap_rep50",
       otu  = "otu_rep50_collap_list", tax  = "silva_taxonomy_collap_rep50_list",
       meta = "metadata_collapsed_list", prefix = "SILVA - in 2+ rep (collapsed)",
       file = "phyloGlobal_silva_collap_rep50.rds"),
  list(var  = "phyloGlobal_silva_collap_rep75",
       otu  = "otu_rep75_collap_list", tax  = "silva_taxonomy_collap_rep75_list",
       meta = "metadata_collapsed_list", prefix = "SILVA - in 3+ rep (collapsed)",
       file = "phyloGlobal_silva_collap_rep75.rds"),
  list(var  = "phyloGlobal_pr2_norep",
       otu  = "otu_norep_list",        tax  = "pr2_taxonomy_norep_list",
       meta = "metadata_list",         prefix = "PR2 - All OTUs",
       file = "phyloGlobal_pr2_norep.rds"),
  list(var  = "phyloGlobal_pr2_rep50",
       otu  = "otu_rep50_list",        tax  = "pr2_taxonomy_rep50_list",
       meta = "metadata_list",         prefix = "PR2 - in 2+ rep",
       file = "phyloGlobal_pr2_rep50.rds"),
  list(var  = "phyloGlobal_pr2_rep75",
       otu  = "otu_rep75_list",        tax  = "pr2_taxonomy_rep75_list",
       meta = "metadata_list",         prefix = "PR2 - in 3+ rep",
       file = "phyloGlobal_pr2_rep75.rds"),
  list(var  = "phyloGlobal_pr2_collap_norep",
       otu  = "otu_norep_collap_list", tax  = "pr2_taxonomy_collap_norep_list",
       meta = "metadata_collapsed_list", prefix = "PR2 - All OTUs (collapsed)",
       file = "phyloGlobal_pr2_collap_norep.rds"),
  list(var  = "phyloGlobal_pr2_collap_rep50",
       otu  = "otu_rep50_collap_list", tax  = "pr2_taxonomy_collap_rep50_list",
       meta = "metadata_collapsed_list", prefix = "PR2 - in 2+ rep (collapsed)",
       file = "phyloGlobal_pr2_collap_rep50.rds"),
  list(var  = "phyloGlobal_pr2_collap_rep75",
       otu  = "otu_rep75_collap_list", tax  = "pr2_taxonomy_collap_rep75_list",
       meta = "metadata_collapsed_list", prefix = "PR2 - in 3+ rep (collapsed)",
       file = "phyloGlobal_pr2_collap_rep75.rds")
)

global_var_names <- sapply(phylo_specs, `[[`, "var")

create_phyloseq_list <- function(otu_list, tax_list, meta_list,
                                 seq_to_asv_map = NULL, name_prefix) {
  cat("\n=== Creating", name_prefix, "===\n")
  
  ps_list <- lapply(names(otu_list), function(method) {
    cat("  Processing:", method, "\n")
    
    otu_mat <- otu_list[[method]]
    tax_mat <- tax_list[[method]]
    meta_df <- meta_list[[method]]
    
    common_taxa <- intersect(colnames(otu_mat), rownames(tax_mat))
    if (length(common_taxa) == 0) { warning("No common taxa for ", method, "!"); return(NULL) }
    cat("    Matched taxa:", length(common_taxa), "/", ncol(otu_mat), "\n")
    
    otu_mat <- otu_mat[, common_taxa, drop = FALSE]
    tax_mat <- tax_mat[common_taxa, , drop = FALSE]
    
    common_samples <- intersect(rownames(otu_mat), rownames(meta_df))
    if (length(common_samples) == 0) { warning("No common samples for ", method, "!"); return(NULL) }
    
    otu_mat <- otu_mat[common_samples, , drop = FALSE]
    meta_df <- meta_df[common_samples, , drop = FALSE]
    cat("    Final:", nrow(otu_mat), "samples x", ncol(otu_mat), "taxa\n")
    
    original_sequences <- colnames(otu_mat)
    
    if (!is.null(seq_to_asv_map)) {
      asv_ids <- seq_to_asv_map[original_sequences]
      if (any(is.na(asv_ids))) warning(sum(is.na(asv_ids)), " sequences not in master map for ", method)
    } else {
      unique_seqs <- sort(unique(original_sequences))
      local_map   <- setNames(paste0("ASV_", seq_along(unique_seqs)), unique_seqs)
      asv_ids     <- local_map[original_sequences]
    }
    
    colnames(otu_mat)    <- asv_ids
    rownames(tax_mat)    <- asv_ids
    sequences_dna        <- DNAStringSet(original_sequences)
    names(sequences_dna) <- asv_ids
    
    phyloseq(
      otu_table(t(otu_mat), taxa_are_rows = TRUE),
      tax_table(as.matrix(tax_mat)),
      sample_data(meta_df),
      sequences_dna
    )
  })
  
  names(ps_list) <- names(otu_list)
  ps_list[!sapply(ps_list, is.null)]
}

split_by_experiment <- function(phylo_list, exp_col = "experiment") {
  result <- list()
  for (iteration_name in names(phylo_list)) {
    ps <- phylo_list[[iteration_name]]
    if (!inherits(ps, "phyloseq")) { warning(iteration_name, " is not phyloseq, skipping"); next }
    meta        <- as.data.frame(sample_data(ps))
    experiments <- unique(meta[[exp_col]])
    result[[iteration_name]] <- list()
    for (exp in experiments) {
      samples_to_keep <- rownames(meta)[meta[[exp_col]] == exp]
      ps_exp          <- prune_samples(samples_to_keep, ps)
      ps_exp          <- prune_taxa(taxa_sums(ps_exp) > 0, ps_exp)
      result[[iteration_name]][[exp]] <- ps_exp
      cat("Split:", iteration_name, "-", exp,
          "- Samples:", nsamples(ps_exp), "- Taxa:", ntaxa(ps_exp), "\n")
    }
  }
  return(result)
}

################################################################################
# 1.4  CREATE AND SAVE RAW GLOBAL PHYLOSEQ OBJECTS
################################################################################

cat("\n### CREATING RAW PHYLOSEQ OBJECTS ###\n\n")

for (spec in phylo_specs) {
  ps_obj <- create_phyloseq_list(
    otu_list       = get(spec$otu),
    tax_list       = get(spec$tax),
    meta_list      = get(spec$meta),
    seq_to_asv_map = master_seq_to_asv,
    name_prefix    = spec$prefix
  )
  saveRDS(ps_obj, file.path(ps_out, spec$file))
  assign(spec$var, ps_obj)
  rm(ps_obj); gc()
}

cat("\nAll 18 raw global phyloseq objects saved.\n")

################################################################################
# 1.5  SAVE RAW GLOBAL OBJECTS + MASTER LIST
################################################################################

cat("\n### SAVING RAW GLOBAL OBJECTS ###\n\n")

for (var_name in global_var_names) {
  saveRDS(get(var_name), file.path(ps_out, paste0(var_name, ".rds")))
  cat("Saved ->", paste0(var_name, ".rds"), "\n")
}

phyloGlobal_all_raw <- setNames(
  lapply(global_var_names, function(nm) readRDS(file.path(ps_out, paste0(nm, ".rds")))),
  sub("phyloGlobal_", "", global_var_names)
)

saveRDS(phyloGlobal_all_raw, file.path(ps_out, "FINAL_PHYLOSEQ_OBJECTS_RAW.rds"))
cat("\nRaw master list saved -> FINAL_PHYLOSEQ_OBJECTS_RAW.rds\n")


################################################################################
################################################################################
##                                                                            ##
##                       PART 2 — CLEAN AND FINALIZE                          ##
##                                                                            ##
################################################################################
################################################################################

################################################################################
# CLEANING FUNCTIONS
################################################################################

check_balance_simple <- function(phyloseq_obj) {
  meta <- as.data.frame(sample_data(phyloseq_obj))
  cat("\n=================================================================\n")
  cat("EXPERIMENTAL BALANCE CHECK - Problems Only\n")
  cat("=================================================================\n\n")
  
  sites <- unique(meta$site)
  sites <- sites[!is.na(sites) & sites != "LAB"]
  
  for (site in sites) {
    cat("\n#################################################################\n")
    cat("###  SITE:", site, "\n")
    cat("#################################################################\n\n")
    
    meta_site <- meta[meta$site == site, ]
    
    cat("ST EXPERIMENT: Season x Substrate\n")
    cat("-----------------------------------------------------------------\n")
    meta_ST <- meta_site[meta_site$experiment == "ST", ]
    if (nrow(meta_ST) > 0) {
      print(table(meta_ST$season, meta_ST$substrate))
      cat("\nExpected: 12 samples per cell (3 bio_reps x 4 PCR_reps)\n\n")
      
      cat("CELLS WITH WRONG NUMBER OF BIOLOGICAL REPLICATES:\n")
      st_bio           <- table(meta_ST$season, meta_ST$substrate, meta_ST$bio_replicate)
      found_bio_problem <- FALSE
      for (season in dimnames(st_bio)[[1]]) {
        for (substrate in dimnames(st_bio)[[2]]) {
          n_bio <- sum(st_bio[season, substrate, ] > 0)
          if (n_bio != 3) {
            cat(sprintf("  %s x %s: %d biological replicates (expected 3)\n", season, substrate, n_bio))
            found_bio_problem <- TRUE
          }
        }
      }
      if (!found_bio_problem) cat("  All cells have 3 biological replicates\n")
      
      cat("\nBIOLOGICAL REPLICATES WITH WRONG NUMBER OF PCR REPLICATES:\n")
      pcr_counts <- meta_ST %>%
        group_by(season, substrate, bio_replicate) %>%
        summarise(n_PCR = n(), .groups = "drop") %>%
        filter(n_PCR != 4)
      if (nrow(pcr_counts) > 0) print(as.data.frame(pcr_counts)) else cat("  All biological replicates have 4 PCR replicates\n")
    } else { cat("No ST samples found\n") }
    
    cat("\nLT EXPERIMENT: Sampling Event x Substrate\n")
    cat("-----------------------------------------------------------------\n")
    meta_LT <- meta_site[meta_site$experiment == "LT", ]
    if (nrow(meta_LT) > 0) {
      print(table(meta_LT$sampling_event, meta_LT$substrate))
      cat("\nExpected: 12 samples per cell\n\n")
      
      cat("CELLS WITH WRONG NUMBER OF BIOLOGICAL REPLICATES:\n")
      lt_bio            <- table(meta_LT$sampling_event, meta_LT$substrate, meta_LT$bio_replicate)
      found_bio_problem <- FALSE
      for (event in dimnames(lt_bio)[[1]]) {
        for (substrate in dimnames(lt_bio)[[2]]) {
          n_bio <- sum(lt_bio[event, substrate, ] > 0)
          if (n_bio != 3) {
            cat(sprintf("  %s x %s: %d biological replicates (expected 3)\n", event, substrate, n_bio))
            found_bio_problem <- TRUE
          }
        }
      }
      if (!found_bio_problem) cat("  All cells have 3 biological replicates\n")
      
      cat("\nBIOLOGICAL REPLICATES WITH WRONG NUMBER OF PCR REPLICATES:\n")
      pcr_counts <- meta_LT %>%
        group_by(sampling_event, substrate, bio_replicate) %>%
        summarise(n_PCR = n(), .groups = "drop") %>%
        filter(n_PCR != 4)
      if (nrow(pcr_counts) > 0) print(as.data.frame(pcr_counts)) else cat("  All biological replicates have 4 PCR replicates\n")
    } else { cat("No LT samples found\n") }
    
    cat("\nWATER EXPERIMENT: Sampling Event\n")
    cat("-----------------------------------------------------------------\n")
    meta_Water <- meta_site[meta_site$experiment == "Water", ]
    if (nrow(meta_Water) > 0) {
      print(table(meta_Water$sampling_event))
      cat("\nExpected: 12 samples per event\n\n")
      
      cat("SAMPLING EVENTS WITH WRONG NUMBER OF BIOLOGICAL REPLICATES:\n")
      water_bio         <- table(meta_Water$sampling_event, meta_Water$bio_replicate)
      found_bio_problem <- FALSE
      for (event in dimnames(water_bio)[[1]]) {
        n_bio <- sum(water_bio[event, ] > 0)
        if (n_bio != 3) {
          cat(sprintf("  %s: %d biological replicates (expected 3)\n", event, n_bio))
          found_bio_problem <- TRUE
        }
      }
      if (!found_bio_problem) cat("  All sampling events have 3 biological replicates\n")
      
      cat("\nBIOLOGICAL REPLICATES WITH WRONG NUMBER OF PCR REPLICATES:\n")
      pcr_counts <- meta_Water %>%
        group_by(sampling_event, bio_replicate) %>%
        summarise(n_PCR = n(), .groups = "drop") %>%
        filter(n_PCR != 4)
      if (nrow(pcr_counts) > 0) print(as.data.frame(pcr_counts)) else cat("  All biological replicates have 4 PCR replicates\n")
    } else { cat("No Water samples found\n") }
  }
  cat("=================================================================\n\n")
}

valid_design <- data.frame(
  library = c("GAB_EUK1","GAB_EUK1","GAB_EUK2","GAB_EUK2","GAB_EUK3","GAB_EUK3","GAB_EUK4"),
  plate   = c("P1","P2","P3","P4","P5","P6","P7"),
  valid_adapters = c(
    "A1,B2,C3,D4", "E5,F6,G7,H8",
    "A1,B2,C3,D4", "A1,F6,G7,D4",
    "A1,B2,C3,D4", "A1,B2,C3,D4",
    "A1,B2,C3,D4"),
  stringsAsFactors = FALSE)

find_invalid_samples <- function(phyloseq_obj, design_table) {
  meta              <- as.data.frame(sample_data(phyloseq_obj))
  rownames_meta     <- rownames(meta)
  meta              <- data.frame(meta, stringsAsFactors = FALSE)
  rownames(meta)    <- rownames_meta
  meta$read_count   <- sample_sums(phyloseq_obj)
  
  cat("Total samples before filtering:", nrow(meta), "\n")
  meta <- meta[meta$experiment != "Degradation", ]
  cat("Total samples after excluding Degradation:", nrow(meta), "\n\n")
  if (nrow(meta) == 0) { cat("No samples remaining!\n"); return(NULL) }
  
  invalid_samples  <- list()
  lib_plate_combos <- unique(data.frame(
    library = as.character(meta$library),
    plate   = as.character(meta$plate),
    stringsAsFactors = FALSE))
  
  lib_plate_combos <- lib_plate_combos[
    !is.na(lib_plate_combos$library) & lib_plate_combos$library != "NA" &
      !is.na(lib_plate_combos$plate)   & lib_plate_combos$plate   != "NA", ]
  
  cat("Library/plate combinations to check:", nrow(lib_plate_combos), "\n\n")
  
  for (i in 1:nrow(lib_plate_combos)) {
    lib        <- lib_plate_combos$library[i]
    plt        <- lib_plate_combos$plate[i]
    design_row <- design_table[design_table$library == lib & design_table$plate == plt, ]
    
    if (nrow(design_row) == 0) {
      cat("WARNING: No design entry for", lib, plt, "\n")
      next
    }
    
    valid_adapters <- trimws(strsplit(design_row$valid_adapters, ",")[[1]])
    combo_samples  <- meta[
      !is.na(meta$library) & !is.na(meta$plate) &
        meta$library == lib & meta$plate == plt, ]
    if (nrow(combo_samples) == 0) next
    
    invalid <- combo_samples[!combo_samples$adapt_primer %in% valid_adapters, ]
    if (nrow(invalid) > 0) {
      cols_to_show <- c("clean_sample_names","bio_replicate","PCR_rep",
                        "adapt_primer","library","plate","read_count",
                        "site","substrate","season","experiment")
      cols_to_show    <- cols_to_show[cols_to_show %in% colnames(invalid)]
      invalid_display <- data.frame(invalid[, cols_to_show, drop = FALSE],
                                    stringsAsFactors = FALSE)
      cat("\nINVALID SAMPLES - Library:", lib, "| Plate:", plt, "\n")
      print(invalid_display)
      invalid_samples[[paste(lib, plt, sep = "_")]] <- invalid_display
    }
  }
  
  # Keep only dataframe elements, handle empty list
  invalid_samples <- Filter(is.data.frame, invalid_samples)
  
  if (length(invalid_samples) == 0) {
    cat("\nTotal samples checked:", nrow(meta), "\n")
    cat("Total INVALID samples found: 0\n")
    return(NULL)
  }
  
  total_invalid <- sum(sapply(invalid_samples, nrow))
  cat("\nTotal samples checked:", nrow(meta), "\n")
  cat("Total INVALID samples found:", total_invalid, "\n")
  
  if (total_invalid > 0) {
    all_invalid <- do.call(rbind, invalid_samples)
    all_invalid$lib_plate_combo <- rep(names(invalid_samples),
                                       sapply(invalid_samples, nrow))
    write.csv(all_invalid, "invalid_samples_to_remove.csv", row.names = TRUE)
    cat("Full list saved to: invalid_samples_to_remove.csv\n")
  }
  
  return(if (total_invalid == 0) NULL else invalid_samples)
}

filter_invalid_samples <- function(phyloseq_input, invalid_samples_list, input_name = NULL) {
  is_list <- is.list(phyloseq_input) && !inherits(phyloseq_input, "phyloseq")
  
  if (is_list) {
    cat("\n### FILTERING PHYLOSEQ LIST:", ifelse(is.null(input_name), "unnamed", input_name), "###\n\n")
    cleaned_list <- list()
    removed_list <- list()
    for (i in seq_along(phyloseq_input)) {
      ps_name                <- names(phyloseq_input)[i]
      result                 <- filter_invalid_samples(phyloseq_input[[i]], invalid_samples_list, ps_name)
      cleaned_list[[ps_name]] <- result$cleaned
      removed_list[[ps_name]] <- result$removed
    }
    return(list(cleaned = cleaned_list, removed = removed_list))
    
  } else {
    if (is.null(invalid_samples_list) || length(invalid_samples_list) == 0) {
      cat("No invalid samples to remove.\n")
      return(list(cleaned = phyloseq_input, removed = NULL))
    }
    
    invalid_rownames <- unlist(lapply(invalid_samples_list, rownames))
    names(invalid_rownames) <- NULL
    all_samples   <- sample_names(phyloseq_input)
    invalid_in_ps <- invalid_rownames[invalid_rownames %in% all_samples]
    
    if (length(invalid_in_ps) == 0) {
      cat("  No matching invalid samples in", input_name, "\n")
      return(list(cleaned = phyloseq_input, removed = NULL))
    }
    
    ps_removed <- prune_samples(invalid_in_ps, phyloseq_input)
    ps_cleaned <- prune_samples(setdiff(all_samples, invalid_in_ps), phyloseq_input)
    ps_cleaned <- prune_taxa(taxa_sums(ps_cleaned) > 0, ps_cleaned)
    
    cat("  ", input_name, "- Removed:", length(invalid_in_ps),
        "| Remaining samples:", nsamples(ps_cleaned),
        "| Remaining taxa:", ntaxa(ps_cleaned), "\n")
    
    return(list(cleaned = ps_cleaned, removed = ps_removed,
                removed_sample_names = invalid_in_ps))
  }
}

remove_dups_and_zeros_list <- function(ps_list) {
  lapply(ps_list, function(ps) {
    all_samples <- sample_names(ps)
    dup_samples <- all_samples[grepl("dup", all_samples, ignore.case = TRUE)]
    if (length(dup_samples) > 0) {
      cat("  - Removing", length(dup_samples), "duplicate samples\n")
      ps <- prune_samples(setdiff(all_samples, dup_samples), ps)
    }
    sample_depths <- sample_sums(ps)
    zero_samples  <- names(sample_depths[sample_depths == 0])
    if (length(zero_samples) > 0) {
      cat("  - Removing", length(zero_samples), "samples with 0 reads\n")
      ps <- prune_samples(sample_depths > 0, ps)
    }
    prune_taxa(taxa_sums(ps) > 0, ps)
  })
}

remove_contaminant_taxa <- function(ps_list, list_name = "unnamed") {
  contaminant_taxa <- c(
    "Diptera","Hymenoptera","Coleoptera","Lepidoptera",
    "Insecta","Mammalia","Aves","Reptilia","Amphibia",
    "Embryophyta","Pinophyta","Magnoliopsida","Liliopsida", "Craniata")
  
  method_logs  <- list()
  cleaned_list <- list()
  
  for (method in names(ps_list)) {
    ps     <- ps_list[[method]]
    tax_df <- as.data.frame(as(tax_table(ps), "matrix"))
    
    flagged <- rownames(tax_df)[
      tax_df$Genus %in% contaminant_taxa |
        tax_df$Order %in% contaminant_taxa |
        tax_df$Class %in% contaminant_taxa
    ]
    
    cat(method, "- Contaminant ASVs removed:", length(flagged),
        "| ASVs remaining:", ntaxa(ps) - length(flagged), "\n")
    
    if (length(flagged) > 0) {
      method_logs[[method]] <- tax_df[flagged, ]
      ps <- prune_taxa(!taxa_names(ps) %in% flagged, ps)
    }
    cleaned_list[[method]] <- ps
  }
  
  assign(paste0("contaminant_log_", list_name), method_logs, envir = globalenv())
  cat("  -> Log saved as: contaminant_log_", list_name, "\n\n", sep = "")
  return(cleaned_list)
}

remove_contaminant_taxa_pr2 <- function(ps_list, list_name = "unnamed") {
  contaminant_taxa <- c(
    "Diptera","Hymenoptera","Coleoptera","Lepidoptera",
    "Insecta","Mammalia","Aves","Reptilia","Amphibia",
    "Embryophyta","Pinophyta","Magnoliopsida","Liliopsida","Hexapoda", "Craniata")
  
  method_logs  <- list()
  cleaned_list <- list()
  
  for (method in names(ps_list)) {
    ps     <- ps_list[[method]]
    tax_df <- as.data.frame(as(tax_table(ps), "matrix"))
    
    keep_euk      <- rownames(tax_df)[grepl("^Eukaryota", tax_df$Domain)]
    n_removed_euk <- ntaxa(ps) - length(keep_euk)
    ps     <- prune_taxa(keep_euk, ps)
    ps     <- prune_samples(sample_sums(ps) > 0, ps)
    tax_df <- as.data.frame(as(tax_table(ps), "matrix"))
    cat(method, "- Non-eukaryote ASVs removed:", n_removed_euk, "\n")
    
    ranks_to_check <- intersect(c("Division","Subdivision","Class","Order","Family"), colnames(tax_df))
    flagged <- rownames(tax_df)[
      apply(tax_df[, ranks_to_check, drop = FALSE], 1, function(row) {
        any(row %in% contaminant_taxa, na.rm = TRUE)
      })
    ]
    
    cat(method, "- Contaminant ASVs removed:", length(flagged),
        "| ASVs remaining:", ntaxa(ps) - length(flagged), "\n")
    
    if (length(flagged) > 0) {
      method_logs[[method]] <- tax_df[flagged, ]
      ps <- prune_taxa(!taxa_names(ps) %in% flagged, ps)
    }
    cleaned_list[[method]] <- ps
  }
  
  assign(paste0("contaminant_log_", list_name), method_logs, envir = globalenv())
  cat("  -> Log saved as: contaminant_log_", list_name, "\n\n", sep = "")
  return(cleaned_list)
}

remove_metazoa_from_ps <- function(ps) {
  tax_df  <- as.data.frame(as(tax_table(ps), "matrix"))
  flagged <- rownames(tax_df)[!is.na(tax_df$Subdivision) & tax_df$Subdivision == "Metazoa"]
  if (length(flagged) > 0) {
    ps <- prune_taxa(!taxa_names(ps) %in% flagged, ps)
    ps <- prune_samples(sample_sums(ps) > 0, ps)
  }
  cat("  Metazoan ASVs removed:", length(flagged), "| Remaining:", ntaxa(ps), "\n")
  return(ps)
}

apply_noMetazoa_to_biofilm_exps <- function(split_list, list_name) {
  # Only ST and LT get metazoa removed; Water is untouched (planktonic metazoa are legitimate there)
  biofilm_exps <- c("ST", "LT")
  lapply(names(split_list), function(method) {
    exp_list <- split_list[[method]]
    lapply(names(exp_list), function(exp) {
      if (exp %in% biofilm_exps) {
        cat(list_name, "-", method, "-", exp, ": removing metazoa\n")
        remove_metazoa_from_ps(exp_list[[exp]])
      } else {
        cat(list_name, "-", method, "-", exp, ": kept intact\n")
        exp_list[[exp]]
      }
    }) %>% setNames(names(exp_list))
  }) %>% setNames(names(split_list))
}

assess_phyloseq_list <- function(ps_list, list_name) {
  cat("\n===", list_name, "===\n")
  summary_df <- map_df(names(ps_list), function(method) {
    ps        <- ps_list[[method]]
    samp_data <- as(sample_data(ps), "data.frame")
    data.frame(
      method        = method,
      total_samples = nsamples(ps),
      total_taxa    = ntaxa(ps),
      n_ST          = sum(samp_data$experiment == "ST",          na.rm = TRUE),
      n_LT          = sum(samp_data$experiment == "LT",          na.rm = TRUE),
      n_DEG         = sum(samp_data$experiment == "Degradation", na.rm = TRUE),
      median_depth  = median(sample_sums(ps)),
      min_depth     = min(sample_sums(ps)),
      max_depth     = max(sample_sums(ps)),
      total_reads   = sum(sample_sums(ps)))
  })
  print(summary_df)
  return(summary_df)
}

################################################################################
# 2.1  BALANCE CHECK BEFORE CLEANING
################################################################################

cat("\n### BALANCE CHECK - BEFORE CLEANING ###\n")
check_balance_simple(phyloGlobal_silva_norep$nP.lulu)

################################################################################
# 2.2  IDENTIFY AND REMOVE INVALID SAMPLES
# Invalid samples are defined by metadata (adapter/plate combinations),
# which is identical across all taxonomy versions and filtering levels,
# so deriving them from one object is sufficient.
################################################################################

invalid_samples <- find_invalid_samples(phyloGlobal_silva_norep$nP.lulu, valid_design)

cat("\n### REMOVING INVALID SAMPLES FROM ALL OBJECTS ###\n\n")

for (var_name in global_var_names) {
  filtered     <- filter_invalid_samples(get(var_name), invalid_samples, var_name)
  cleaned_name <- paste0(var_name, "_cleaned")
  assign(cleaned_name, filtered$cleaned)
  if (var_name == "phyloGlobal_silva_norep") removed_samples_phyloseq <- filtered$removed
  rm(filtered); gc()
}

################################################################################
# 2.3  REMOVE DUPLICATE 
################################################################################

cat("\n### REMOVING DUPLICATE SAMPLES ###\n\n")

cleaned_var_names <- paste0(global_var_names, "_cleaned")

for (obj_name in cleaned_var_names) {
  cat("Processing:", obj_name, "\n")
  assign(obj_name, lapply(get(obj_name), function(ps) {
    all_samples <- sample_names(ps)
    dup_samples <- all_samples[grepl("dup", all_samples, ignore.case = TRUE)]
    if (length(dup_samples) > 0) {
      cat("  - Removing", length(dup_samples), "duplicate samples\n")
      ps <- prune_samples(setdiff(all_samples, dup_samples), ps)
      ps <- prune_taxa(taxa_sums(ps) > 0, ps)
    }
    ps
  }))
  gc()
}


################################################################################
# 2.4  REMOVE CONTAMINANT TAXA
################################################################################

cat("\n### REMOVING CONTAMINANT TAXA ###\n\n")

contam_specs <- lapply(phylo_specs, function(spec) {
  list(
    name  = paste0(spec$var, "_cleaned"),
    label = sub("phyloGlobal_", "", spec$var),
    pr2   = grepl("^phyloGlobal_pr2", spec$var)
  )
})

for (cs in contam_specs) {
  if (cs$pr2) {
    assign(cs$name, remove_contaminant_taxa_pr2(get(cs$name), cs$label))
  } else {
    assign(cs$name, remove_contaminant_taxa(get(cs$name), cs$label))
  }
  gc()
}


################################################################################
# ZERO-DEPTH SAMPLE REMOVAL — SANITY CHECK
# Quality filtering was handled comprehensively in Script 1.
# This is a safety net only — flags anything unexpected that slipped through
# after invalid sample and duplicate removal.
################################################################################

cat("\n### ZERO-DEPTH SANITY CHECK ###\n\n")

for (obj_name in cleaned_var_names) {
  any_removed <- FALSE
  assign(obj_name, lapply(names(get(obj_name)), function(method) {
    ps   <- get(obj_name)[[method]]
    zero <- sample_names(ps)[sample_sums(ps) == 0]
    if (length(zero) > 0) {
      cat(sprintf("  REMOVED: %s | %s — %d zero-depth sample(s): %s\n",
                  obj_name, method, length(zero), paste(zero, collapse = ", ")))
      any_removed <<- TRUE
      ps <- prune_samples(sample_sums(ps) > 0, ps)
      ps <- prune_taxa(taxa_sums(ps) > 0, ps)
    }
    ps
  }) %>% setNames(names(get(obj_name))))
  if (!any_removed) cat(sprintf("  %s — no zero-depth samples found\n", obj_name))
  gc()
}



################################################################################
# 2.5  BALANCE CHECK AFTER CLEANING
################################################################################

cat("\n### BALANCE CHECK - AFTER CLEANING ###\n")
check_balance_simple(phyloGlobal_silva_norep_cleaned$nP.lulu)
check_balance_simple(phyloGlobal_silva_norep_cleaned$nP.swarm)

################################################################################
# 2.6  QC OVERVIEW
################################################################################

cat("\n### QC OVERVIEW ###\n")
assess_phyloseq_list(phyloGlobal_silva_norep_cleaned,        "SILVA norep")
assess_phyloseq_list(phyloGlobal_silva_rep50_cleaned,        "SILVA rep50")
assess_phyloseq_list(phyloGlobal_silva_rep75_cleaned,        "SILVA rep75")
assess_phyloseq_list(phyloGlobal_silva_collap_norep_cleaned, "SILVA collap norep")
assess_phyloseq_list(phyloGlobal_silva_collap_rep50_cleaned, "SILVA collap rep50")
assess_phyloseq_list(phyloGlobal_silva_collap_rep75_cleaned, "SILVA collap rep75")
assess_phyloseq_list(phyloGlobal_pr2_norep_cleaned,          "PR2 norep")
assess_phyloseq_list(phyloGlobal_pr2_rep50_cleaned,          "PR2 rep50")
assess_phyloseq_list(phyloGlobal_pr2_rep75_cleaned,          "PR2 rep75")
assess_phyloseq_list(phyloGlobal_pr2_collap_norep_cleaned,   "PR2 collap norep")
assess_phyloseq_list(phyloGlobal_pr2_collap_rep50_cleaned,   "PR2 collap rep50")
assess_phyloseq_list(phyloGlobal_pr2_collap_rep75_cleaned,   "PR2 collap rep75")

################################################################################
# 2.7  SAVE ALL 18 CLEANED GLOBAL OBJECTS
################################################################################

cat("\n### SAVING CLEANED GLOBAL OBJECTS ###\n\n")

for (var_name in global_var_names) {
  cleaned_name <- paste0(var_name, "_cleaned")
  rds_name     <- paste0(cleaned_name, ".rds")
  saveRDS(get(cleaned_name), file.path(ps_out, rds_name))
  cat("Saved ->", rds_name, "\n")
}

saveRDS(removed_samples_phyloseq, file.path(ps_out, "removed_samples.rds"))
saveRDS(invalid_samples,          file.path(ps_out, "invalid_samples_list.rds"))

################################################################################
# 2.8  SPLIT CLEANED OBJECTS BY EXPERIMENT AND SAVE
################################################################################

cat("\n### SPLITTING CLEANED OBJECTS BY EXPERIMENT ###\n\n")

for (var_name in global_var_names) {
  cleaned_name <- paste0(var_name, "_cleaned")
  split_name   <- sub("phyloGlobal_", "phyloExpSplit_", cleaned_name)
  
  split_obj <- split_by_experiment(get(cleaned_name))
  assign(split_name, split_obj)
  saveRDS(split_obj, file.path(ps_out, paste0(split_name, ".rds")))
  cat("Saved ->", paste0(split_name, ".rds"), "\n")
  rm(split_obj); gc()
}

################################################################################
# 2.9  PR2 NO-METAZOAN VERSIONS (ST + LT only)
################################################################################

cat("\n### CREATING NO-METAZOAN PR2 VERSIONS (ST + LT only) ###\n\n")

pr2_split_names <- c(
  "phyloExpSplit_pr2_norep_cleaned",
  "phyloExpSplit_pr2_rep50_cleaned",
  "phyloExpSplit_pr2_rep75_cleaned",
  "phyloExpSplit_pr2_collap_norep_cleaned",
  "phyloExpSplit_pr2_collap_rep50_cleaned",
  "phyloExpSplit_pr2_collap_rep75_cleaned"
)

for (split_name in pr2_split_names) {
  nometa_name <- paste0(split_name, "_noMeta")
  nometa_obj  <- apply_noMetazoa_to_biofilm_exps(get(split_name), split_name)
  assign(nometa_name, nometa_obj)
  saveRDS(nometa_obj, file.path(ps_out, paste0(nometa_name, ".rds")))
  cat("Saved ->", paste0(nometa_name, ".rds"), "\n")
  rm(nometa_obj); gc()
}

################################################################################
# 2.10 ASSEMBLE AND SAVE FINAL MASTER LIST
################################################################################

cat("\nAssembling FINAL_PHYLOSEQ_OBJECTS master list...\n")

all_cleaned_split_names <- paste0(
  sub("phyloGlobal_", "phyloExpSplit_", global_var_names), "_cleaned"
)
all_nometa_names <- paste0(pr2_split_names, "_noMeta")

all_final_names <- c(all_cleaned_split_names, all_nometa_names)

phyloExpSplit_all <- setNames(
  lapply(all_final_names, function(nm) readRDS(file.path(ps_out, paste0(nm, ".rds")))),
  sub("phyloExpSplit_", "", sub("_cleaned", "", all_final_names))
)

saveRDS(phyloExpSplit_all, file.path(ps_out, "FINAL_PHYLOSEQ_OBJECTS.rds"))
save.image("all_phylo_objects_v4.RData")

cat("\n### PIPELINE COMPLETE ###\n")
cat("Master ASV map:           master_seq_to_asv\n")
cat("Raw phyloseq objects:     FINAL_PHYLOSEQ_OBJECTS_RAW.rds\n")
cat("Cleaned phyloseq objects: FINAL_PHYLOSEQ_OBJECTS.rds\n")
cat("Output directory:        ", ps_out, "\n")

