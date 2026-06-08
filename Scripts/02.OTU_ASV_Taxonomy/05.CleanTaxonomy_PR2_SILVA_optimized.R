################################################################################
############## TAXONOMY WITH CLEANED OTU TABLES ###############################

# 3 pooling methods: nP, psP, fP
# 2 replicate filtering levels: 2/4 and 3/4
# 2 clustering methods: LULU vs SWARM
# 2 taxonomy sources: SILVA (DADA2) vs PR2
# (collapsed and uncollapsed)

# OTU tables: rownames are sample names
# OTU tables: colnames are sequences (the actual ASV sequences)
# SILVA taxonomy: rownames are sequences (matching the colnames of OTU tables)

library("dplyr")

################################################################################
########################## PATH CONFIGURATION ##################################
################################################################################

GITHUB_DIR    <- "~/GitHub/MAPLE_Seasonal_Plastisphere"
SCRIPTS_DIR   <- file.path(GITHUB_DIR, "Scripts", "02.OTU_ASV_Taxonomy")
PROC_DIR      <- file.path(GITHUB_DIR, "Processed_data")
RAW_DIR       <- file.path(GITHUB_DIR, "Raw_data")
RESULTS_DIR   <- file.path(GITHUB_DIR, "Results")

OTU_FILT_DIR  <- file.path(PROC_DIR, "OTU_tables", "filtered")
TAXONOMY_DIR  <- file.path(PROC_DIR, "Taxonomy")
SILVA_TAX_DIR <- file.path(RAW_DIR, "TaxonomyAssignements", "silvaDB_dada2")
PR2_TAX_DIR   <- file.path(RAW_DIR, "TaxonomyAssignements", "taxonomy_PR2")
TRACKING_DIR  <- file.path(RESULTS_DIR, "SanityChecks", "Read_tracking")

## Create directories if they don't exist
invisible(lapply(
  c(TAXONOMY_DIR, TRACKING_DIR),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

setwd(SCRIPTS_DIR)

################################################################################
########################## Load data ###########################################
################################################################################

### Cleaned OTU tables (after control filtering)

clean_lulu_OTUtables <- list(
  nP.lulu_rep50_otu  = readRDS(file.path(OTU_FILT_DIR, "LULU_nP_repfilt50.rds")),
  psP.lulu_rep50_otu = readRDS(file.path(OTU_FILT_DIR, "LULU_psP_repfilt50.rds")),
  fP.lulu_rep50_otu  = readRDS(file.path(OTU_FILT_DIR, "LULU_fP_repfilt50.rds")),
  nP.lulu_rep75_otu  = readRDS(file.path(OTU_FILT_DIR, "LULU_nP_repfilt75.rds")),
  psP.lulu_rep75_otu = readRDS(file.path(OTU_FILT_DIR, "LULU_psP_repfilt75.rds")),
  fP.lulu_rep75_otu  = readRDS(file.path(OTU_FILT_DIR, "LULU_fP_repfilt75.rds")),
  nP.lulu_norep_otu  = readRDS(file.path(OTU_FILT_DIR, "LULU_nP_norepfilt.rds")),
  psP.lulu_norep_otu = readRDS(file.path(OTU_FILT_DIR, "LULU_psP_norepfilt.rds")),
  fP.lulu_norep_otu  = readRDS(file.path(OTU_FILT_DIR, "LULU_fP_norepfilt.rds"))
)

clean_swarm_OTUtables <- list(
  nP.swarm_rep50_otu  = readRDS(file.path(OTU_FILT_DIR, "SWARM_nP_repfilt50.rds")),
  psP.swarm_rep50_otu = readRDS(file.path(OTU_FILT_DIR, "SWARM_psP_repfilt50.rds")),
  fP.swarm_rep50_otu  = readRDS(file.path(OTU_FILT_DIR, "SWARM_fP_repfilt50.rds")),
  nP.swarm_rep75_otu  = readRDS(file.path(OTU_FILT_DIR, "SWARM_nP_repfilt75.rds")),
  psP.swarm_rep75_otu = readRDS(file.path(OTU_FILT_DIR, "SWARM_psP_repfilt75.rds")),
  fP.swarm_rep75_otu  = readRDS(file.path(OTU_FILT_DIR, "SWARM_fP_repfilt75.rds")),
  nP.swarm_norep_otu  = readRDS(file.path(OTU_FILT_DIR, "SWARM_nP_norepfilt.rds")),
  psP.swarm_norep_otu = readRDS(file.path(OTU_FILT_DIR, "SWARM_psP_norepfilt.rds")),
  fP.swarm_norep_otu  = readRDS(file.path(OTU_FILT_DIR, "SWARM_fP_norepfilt.rds"))
)

### Collapsed OTU tables

clean_collap_lulu_OTUtables <- list(
  nP.lulu_collapsed_rep50_otu  = readRDS(file.path(OTU_FILT_DIR, "LULU_nP_CollapsedRepFilt50.rds")),
  psP.lulu_collapsed_rep50_otu = readRDS(file.path(OTU_FILT_DIR, "LULU_psP_CollapsedRepFilt50.rds")),
  fP.lulu_collapsed_rep50_otu  = readRDS(file.path(OTU_FILT_DIR, "LULU_fP_CollapsedRepFilt50.rds")),
  nP.lulu_collapsed_rep75_otu  = readRDS(file.path(OTU_FILT_DIR, "LULU_nP_CollapsedRepFilt75.rds")),
  psP.lulu_collapsed_rep75_otu = readRDS(file.path(OTU_FILT_DIR, "LULU_psP_CollapsedRepFilt75.rds")),
  fP.lulu_collapsed_rep75_otu  = readRDS(file.path(OTU_FILT_DIR, "LULU_fP_CollapsedRepFilt75.rds")),
  nP.lulu_collapsed_norep_otu  = readRDS(file.path(OTU_FILT_DIR, "LULU_nP_collap_norepfilt.rds")),
  psP.lulu_collapsed_norep_otu = readRDS(file.path(OTU_FILT_DIR, "LULU_psP_collap_norepfilt.rds")),
  fP.lulu_collapsed_norep_otu  = readRDS(file.path(OTU_FILT_DIR, "LULU_fP_collap_norepfilt.rds"))
)

clean_collap_swarm_OTUtables <- list(
  nP.sw_collapsed_rep50_otu  = readRDS(file.path(OTU_FILT_DIR, "SWARM_nP_CollapsedRepFilt50.rds")),
  psP.sw_collapsed_rep50_otu = readRDS(file.path(OTU_FILT_DIR, "SWARM_psP_CollapsedRepFilt50.rds")),
  fP.sw_collapsed_rep50_otu  = readRDS(file.path(OTU_FILT_DIR, "SWARM_fP_CollapsedRepFilt50.rds")),
  nP.sw_collapsed_rep75_otu  = readRDS(file.path(OTU_FILT_DIR, "SWARM_nP_CollapsedRepFilt75.rds")),
  psP.sw_collapsed_rep75_otu = readRDS(file.path(OTU_FILT_DIR, "SWARM_psP_CollapsedRepFilt75.rds")),
  fP.sw_collapsed_rep75_otu  = readRDS(file.path(OTU_FILT_DIR, "SWARM_fP_CollapsedRepFilt75.rds")),
  nP.sw_collapsed_norep_otu  = readRDS(file.path(OTU_FILT_DIR, "SWARM_nP_collap_norepfilt.rds")),
  psP.sw_collapsed_norep_otu = readRDS(file.path(OTU_FILT_DIR, "SWARM_psP_collap_norepfilt.rds")),
  fP.sw_collapsed_norep_otu  = readRDS(file.path(OTU_FILT_DIR, "SWARM_fP_collap_norepfilt.rds"))
)

################################################################################
### SILVA taxonomy assignments
################################################################################

SILVA_TAX_LEVELS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

clean_silva_taxonomy <- function(tax_mat) {
  tax_mat[, 1:7, drop = FALSE]  # Keep only columns 1-7 (Kingdom through Species)
}

silva_tax_lulu <- lapply(
  list(
    silva132_lulu_nP  = readRDS(file.path(SILVA_TAX_DIR, "taxa132_nP.rds")),
    silva132_lulu_psP = readRDS(file.path(SILVA_TAX_DIR, "taxa132_psP.rds")),
    silva132_lulu_fP  = readRDS(file.path(SILVA_TAX_DIR, "taxa132_fP.rds"))
  ),
  clean_silva_taxonomy
)

silva_tax_swarm <- lapply(
  list(
    silva132_sw_nP  = readRDS(file.path(SILVA_TAX_DIR, "taxa132_swarmed_nP.rds")),
    silva132_sw_psP = readRDS(file.path(SILVA_TAX_DIR, "taxa132_swarmed_psP.rds")),
    silva132_sw_fP  = readRDS(file.path(SILVA_TAX_DIR, "taxa132_swarmed_fP.rds"))
  ),
  clean_silva_taxonomy
)

cat("SILVA rank coverage (lulu nP):\n")
for (rank in SILVA_TAX_LEVELS) {
  n   <- sum(!is.na(silva_tax_lulu$silva132_lulu_nP[, rank]))
  pct <- round(n / nrow(silva_tax_lulu$silva132_lulu_nP) * 100, 1)
  cat(sprintf("  %-12s: %5d / %5d (%5.1f%%)\n", rank, n, nrow(silva_tax_lulu$silva132_lulu_nP), pct))
}

################################################################################
### PR2 taxonomy assignments
################################################################################

PR2_TAX_LEVELS <- c("Domain", "Supergroup", "Division", "Subdivision",
                    "Class", "Order", "Family", "Genus", "Species")

pr2_master <- readRDS(file.path(PR2_TAX_DIR, "all_taxonomy_by_dataset_v2.rds"))

pr2_tax_lulu <- list(
  pr2_lulu_nP  = pr2_master$nP_lulu$dada2$tax,
  pr2_lulu_psP = pr2_master$psP_lulu$dada2$tax,
  pr2_lulu_fP  = pr2_master$fP_lulu$dada2$tax
)

pr2_tax_swarm <- list(
  pr2_sw_nP  = pr2_master$nP_sw$dada2$tax,
  pr2_sw_psP = pr2_master$psP_sw$dada2$tax,
  pr2_sw_fP  = pr2_master$fP_sw$dada2$tax
)

cat("PR2 LULU nP - first rowname (should be a DNA sequence):\n")
cat(substr(rownames(pr2_tax_lulu$pr2_lulu_nP)[1], 1, 50), "...\n\n")

cat("PR2 rank coverage:\n")
for (rank in PR2_TAX_LEVELS) {
  n   <- sum(!is.na(pr2_tax_lulu$pr2_lulu_nP[, rank]))
  pct <- round(n / nrow(pr2_tax_lulu$pr2_lulu_nP) * 100, 1)
  cat(sprintf("  %-12s: %5d / %5d (%5.1f%%)\n", rank, n, nrow(pr2_tax_lulu$pr2_lulu_nP), pct))
}

################################################################################
################## FILTER FUNCTION ############################################
################################################################################

# SILVA and PR2 both use sequences as rownames, so the same function handles both

filter_silva_taxonomy <- function(silva_tax, otu_table) {
  kept_seq     <- colnames(otu_table)
  filtered_tax <- silva_tax[rownames(silva_tax) %in% kept_seq, , drop = FALSE]
  
  cat(sprintf("  Original OTUs: %d\n", nrow(silva_tax)))
  cat(sprintf("  Kept OTUs:     %d\n", nrow(filtered_tax)))
  cat(sprintf("  Removed:       %d (%.1f%%)\n",
              nrow(silva_tax) - nrow(filtered_tax),
              (1 - nrow(filtered_tax) / nrow(silva_tax)) * 100))
  
  return(filtered_tax)
}

################################################################################
################## FILTERING CONFIG TABLE #####################################
################################################################################

# This table drives all filtering. Each row specifies:
#   tax_list        : the taxonomy list (e.g. silva_tax_lulu)
#   tax_name_prefix : prefix used in tax_list names (e.g. "silva132_lulu_")
#   otu_list        : the OTU table list to draw from
#   otu_suffix      : suffix to build the OTU table name (e.g. ".lulu_rep50_otu")
#   output_name     : the R object to assign results to
#   save_path       : where to save the RDS file

filter_configs <- list(
  
  ## ---- SILVA LULU UNCOLLAPSED ----
  list(tax_list = silva_tax_lulu, tax_prefix = "silva132_lulu_",
       otu_list = clean_lulu_OTUtables,        otu_suffix = ".lulu_rep50_otu",
       output_name = "lulu_silva_tax_filtered_rep50",
       save_path   = file.path(TAXONOMY_DIR, "LULU_silva_taxonomy_filtered_rep50.rds")),
  
  list(tax_list = silva_tax_lulu, tax_prefix = "silva132_lulu_",
       otu_list = clean_lulu_OTUtables,        otu_suffix = ".lulu_rep75_otu",
       output_name = "lulu_silva_tax_filtered_rep75",
       save_path   = file.path(TAXONOMY_DIR, "LULU_silva_taxonomy_filtered_rep75.rds")),
  
  list(tax_list = silva_tax_lulu, tax_prefix = "silva132_lulu_",
       otu_list = clean_lulu_OTUtables,        otu_suffix = ".lulu_norep_otu",
       output_name = "lulu_silva_tax_filtered_norepfilt",
       save_path   = file.path(TAXONOMY_DIR, "LULU_silva_taxonomy_filtered_norepfilt.rds")),
  
  ## ---- SILVA LULU COLLAPSED ----
  list(tax_list = silva_tax_lulu, tax_prefix = "silva132_lulu_",
       otu_list = clean_collap_lulu_OTUtables, otu_suffix = ".lulu_collapsed_rep50_otu",
       output_name = "lulu_silva_tax_filtcollapsed_rep50",
       save_path   = file.path(TAXONOMY_DIR, "LULU_silva_taxonomy_filtered_collapsed_rep50.rds")),
  
  list(tax_list = silva_tax_lulu, tax_prefix = "silva132_lulu_",
       otu_list = clean_collap_lulu_OTUtables, otu_suffix = ".lulu_collapsed_rep75_otu",
       output_name = "lulu_silva_tax_filtcollapsed_rep75",
       save_path   = file.path(TAXONOMY_DIR, "LULU_silva_taxonomy_filtered_collapsed_rep75.rds")),
  
  list(tax_list = silva_tax_lulu, tax_prefix = "silva132_lulu_",
       otu_list = clean_collap_lulu_OTUtables, otu_suffix = ".lulu_collapsed_norep_otu",
       output_name = "lulu_silva_tax_filtcollapsed_norepfilt",
       save_path   = file.path(TAXONOMY_DIR, "LULU_silva_taxonomy_filtered_collapsed_norepfilt.rds")),
  
  ## ---- SILVA SWARM UNCOLLAPSED ----
  list(tax_list = silva_tax_swarm, tax_prefix = "silva132_sw_",
       otu_list = clean_swarm_OTUtables,        otu_suffix = ".swarm_rep50_otu",
       output_name = "swarm_silva_tax_filtered_rep50",
       save_path   = file.path(TAXONOMY_DIR, "SWARM_silva_taxonomy_filtered_rep50.rds")),
  
  list(tax_list = silva_tax_swarm, tax_prefix = "silva132_sw_",
       otu_list = clean_swarm_OTUtables,        otu_suffix = ".swarm_rep75_otu",
       output_name = "swarm_silva_tax_filtered_rep75",
       save_path   = file.path(TAXONOMY_DIR, "SWARM_silva_taxonomy_filtered_rep75.rds")),
  
  list(tax_list = silva_tax_swarm, tax_prefix = "silva132_sw_",
       otu_list = clean_swarm_OTUtables,        otu_suffix = ".swarm_norep_otu",
       output_name = "swarm_silva_tax_filtered_norepfilt",
       save_path   = file.path(TAXONOMY_DIR, "SWARM_silva_taxonomy_filtered_norepfilt.rds")),
  
  ## ---- SILVA SWARM COLLAPSED ----
  list(tax_list = silva_tax_swarm, tax_prefix = "silva132_sw_",
       otu_list = clean_collap_swarm_OTUtables, otu_suffix = ".sw_collapsed_rep50_otu",
       output_name = "swarm_silva_tax_filtcollapsed_rep50",
       save_path   = file.path(TAXONOMY_DIR, "SWARM_silva_taxonomy_filtered_collapsed_rep50.rds")),
  
  list(tax_list = silva_tax_swarm, tax_prefix = "silva132_sw_",
       otu_list = clean_collap_swarm_OTUtables, otu_suffix = ".sw_collapsed_rep75_otu",
       output_name = "swarm_silva_tax_filtcollapsed_rep75",
       save_path   = file.path(TAXONOMY_DIR, "SWARM_silva_taxonomy_filtered_collapsed_rep75.rds")),
  
  list(tax_list = silva_tax_swarm, tax_prefix = "silva132_sw_",
       otu_list = clean_collap_swarm_OTUtables, otu_suffix = ".sw_collapsed_norep_otu",
       output_name = "swarm_silva_tax_filtcollapsed_norepfilt",
       save_path   = file.path(TAXONOMY_DIR, "SWARM_silva_taxonomy_filtered_collapsed_norepfilt.rds")),
  
  ## ---- PR2 LULU UNCOLLAPSED ----
  list(tax_list = pr2_tax_lulu, tax_prefix = "pr2_lulu_",
       otu_list = clean_lulu_OTUtables,        otu_suffix = ".lulu_rep50_otu",
       output_name = "lulu_pr2_tax_filtered_rep50",
       save_path   = file.path(TAXONOMY_DIR, "LULU_pr2_taxonomy_filtered_rep50.rds")),
  
  list(tax_list = pr2_tax_lulu, tax_prefix = "pr2_lulu_",
       otu_list = clean_lulu_OTUtables,        otu_suffix = ".lulu_rep75_otu",
       output_name = "lulu_pr2_tax_filtered_rep75",
       save_path   = file.path(TAXONOMY_DIR, "LULU_pr2_taxonomy_filtered_rep75.rds")),
  
  list(tax_list = pr2_tax_lulu, tax_prefix = "pr2_lulu_",
       otu_list = clean_lulu_OTUtables,        otu_suffix = ".lulu_norep_otu",
       output_name = "lulu_pr2_tax_filtered_norepfilt",
       save_path   = file.path(TAXONOMY_DIR, "LULU_pr2_taxonomy_filtered_norepfilt.rds")),
  
  ## ---- PR2 LULU COLLAPSED ----
  list(tax_list = pr2_tax_lulu, tax_prefix = "pr2_lulu_",
       otu_list = clean_collap_lulu_OTUtables, otu_suffix = ".lulu_collapsed_rep50_otu",
       output_name = "lulu_pr2_tax_filtcollapsed_rep50",
       save_path   = file.path(TAXONOMY_DIR, "LULU_pr2_taxonomy_filtered_collapsed_rep50.rds")),
  
  list(tax_list = pr2_tax_lulu, tax_prefix = "pr2_lulu_",
       otu_list = clean_collap_lulu_OTUtables, otu_suffix = ".lulu_collapsed_rep75_otu",
       output_name = "lulu_pr2_tax_filtcollapsed_rep75",
       save_path   = file.path(TAXONOMY_DIR, "LULU_pr2_taxonomy_filtered_collapsed_rep75.rds")),
  
  list(tax_list = pr2_tax_lulu, tax_prefix = "pr2_lulu_",
       otu_list = clean_collap_lulu_OTUtables, otu_suffix = ".lulu_collapsed_norep_otu",
       output_name = "lulu_pr2_tax_filtcollapsed_norepfilt",
       save_path   = file.path(TAXONOMY_DIR, "LULU_pr2_taxonomy_filtered_collapsed_norepfilt.rds")),
  
  ## ---- PR2 SWARM UNCOLLAPSED ----
  list(tax_list = pr2_tax_swarm, tax_prefix = "pr2_sw_",
       otu_list = clean_swarm_OTUtables,        otu_suffix = ".swarm_rep50_otu",
       output_name = "swarm_pr2_tax_filtered_rep50",
       save_path   = file.path(TAXONOMY_DIR, "SWARM_pr2_taxonomy_filtered_rep50.rds")),
  
  list(tax_list = pr2_tax_swarm, tax_prefix = "pr2_sw_",
       otu_list = clean_swarm_OTUtables,        otu_suffix = ".swarm_rep75_otu",
       output_name = "swarm_pr2_tax_filtered_rep75",
       save_path   = file.path(TAXONOMY_DIR, "SWARM_pr2_taxonomy_filtered_rep75.rds")),
  
  list(tax_list = pr2_tax_swarm, tax_prefix = "pr2_sw_",
       otu_list = clean_swarm_OTUtables,        otu_suffix = ".swarm_norep_otu",
       output_name = "swarm_pr2_tax_filtered_norepfilt",
       save_path   = file.path(TAXONOMY_DIR, "SWARM_pr2_taxonomy_filtered_norepfilt.rds")),
  
  ## ---- PR2 SWARM COLLAPSED ----
  list(tax_list = pr2_tax_swarm, tax_prefix = "pr2_sw_",
       otu_list = clean_collap_swarm_OTUtables, otu_suffix = ".sw_collapsed_rep50_otu",
       output_name = "swarm_pr2_tax_filtcollapsed_rep50",
       save_path   = file.path(TAXONOMY_DIR, "SWARM_pr2_taxonomy_filtered_collapsed_rep50.rds")),
  
  list(tax_list = pr2_tax_swarm, tax_prefix = "pr2_sw_",
       otu_list = clean_collap_swarm_OTUtables, otu_suffix = ".sw_collapsed_rep75_otu",
       output_name = "swarm_pr2_tax_filtcollapsed_rep75",
       save_path   = file.path(TAXONOMY_DIR, "SWARM_pr2_taxonomy_filtered_collapsed_rep75.rds")),
  
  list(tax_list = pr2_tax_swarm, tax_prefix = "pr2_sw_",
       otu_list = clean_collap_swarm_OTUtables, otu_suffix = ".sw_collapsed_norep_otu",
       output_name = "swarm_pr2_tax_filtcollapsed_norepfilt",
       save_path   = file.path(TAXONOMY_DIR, "SWARM_pr2_taxonomy_filtered_collapsed_norepfilt.rds"))
)

################################################################################
################## RUN FILTERING ##############################################
################################################################################

# This single loop replaces ~300 lines of near-identical lapply blocks.
# For each config entry it:
#   1. Filters each pooling method's taxonomy against the correct OTU table
#   2. Assigns the result to the named object in the global environment
#   3. Saves the result to disk

for (cfg in filter_configs) {
  
  pooling_methods <- c("nP", "psP", "fP")
  
  result <- lapply(pooling_methods, function(pooling) {
    
    tax_name  <- paste0(cfg$tax_prefix, pooling)
    otu_name  <- paste0(pooling, cfg$otu_suffix)
    otu_table <- cfg$otu_list[[otu_name]]
    
    if (is.null(otu_table)) stop("OTU table not found: ", otu_name)
    
    cat(tax_name, ":\n")
    filter_silva_taxonomy(cfg$tax_list[[tax_name]], otu_table)
  })
  
  # Restore original names (e.g. "silva132_lulu_nP", "pr2_sw_fP", ...)
  names(result) <- paste0(cfg$tax_prefix, pooling_methods)
  
  assign(cfg$output_name, result, envir = .GlobalEnv)
  saveRDS(result, cfg$save_path)
  
  cat("  -> Saved:", cfg$output_name, "\n\n")
}

save.image("TaxAssignment_ready4Phyloseq.RData")

################################################################################
################### TAXONOMY FILTERING DIAGNOSTICS #############################
################################################################################

library(dplyr)
library(tidyr)

cat("\n========================================\n")
cat("TAXONOMY FILTERING PIPELINE TRACKING\n")
cat("========================================\n\n")

create_taxonomy_summary <- function(method_name, original_tax, filtered_tax,
                                    otu_table, tax_source) {
  n_original <- nrow(original_tax)
  n_filtered <- nrow(filtered_tax)
  n_lost     <- n_original - n_filtered
  
  data.frame(
    Method            = method_name,
    TaxSource         = tax_source,
    Samples           = nrow(otu_table),
    OTUs_in_Table     = ncol(otu_table),
    Total_Reads       = sum(otu_table),
    Original_Taxa     = n_original,
    Filtered_Taxa     = n_filtered,
    Taxa_Lost         = n_lost,
    Taxa_Lost_Pct     = round((n_lost / n_original) * 100, 2),
    Taxa_Retained_Pct = round((n_filtered / n_original) * 100, 2),
    stringsAsFactors  = FALSE
  )
}

# Diagnostic config reuses the same filter_configs, adding a tax_source label

diag_source_map <- list(
  "silva132_lulu_" = "SILVA",
  "silva132_sw_"   = "SILVA",
  "pr2_lulu_"      = "PR2",
  "pr2_sw_"        = "PR2"
)

taxonomy_tracking <- list()

for (cfg in filter_configs) {
  
  tax_source    <- diag_source_map[[cfg$tax_prefix]]
  result_list   <- get(cfg$output_name, envir = .GlobalEnv)
  pooling_methods <- c("nP", "psP", "fP")
  
  for (pooling in pooling_methods) {
    
    tax_name  <- paste0(cfg$tax_prefix, pooling)
    otu_name  <- paste0(pooling, cfg$otu_suffix)
    otu_table <- cfg$otu_list[[otu_name]]
    
    method_label <- paste0(cfg$output_name, "_", pooling)
    
    taxonomy_tracking[[method_label]] <- create_taxonomy_summary(
      method_name  = method_label,
      original_tax = cfg$tax_list[[tax_name]],
      filtered_tax = result_list[[tax_name]],
      otu_table    = otu_table,
      tax_source   = tax_source
    )
    
    cat(sprintf("  %s: %d -> %d taxa (%.1f%% retained)\n",
                method_label,
                nrow(cfg$tax_list[[tax_name]]),
                nrow(result_list[[tax_name]]),
                taxonomy_tracking[[method_label]]$Taxa_Retained_Pct))
  }
}

cat("\n========================================\n")
cat("EXPORTING TAXONOMY SUMMARY TABLES\n")
cat("========================================\n\n")

taxonomy_summary_all <- bind_rows(taxonomy_tracking)

write.csv(taxonomy_summary_all,
          file.path(TRACKING_DIR, "taxonomy_filtering_summary_ALL.csv"),
          row.names = FALSE)

saveRDS(taxonomy_tracking,
        file.path(TRACKING_DIR, "taxonomy_tracking_object.rds"))

cat("Files saved:\n")
cat("  - taxonomy_filtering_summary_ALL.csv\n")
cat("  - taxonomy_tracking_object.rds\n")
cat("\n========================================\n")
cat("TAXONOMY FILTERING DIAGNOSTICS COMPLETE!\n")
cat("========================================\n")


### Positive control removed from script 04.Control filtering
# Join taxonomy onto the log — run after Script 2 has been run
tax_raw <- readRDS(file.path(SILVA_TAX_DIR, "taxa132_nP.rds"))

pos_filter_log <- read.csv(
  file.path(RESULTS_DIR, "SanityChecks", "positive_control_filter_log.csv"),
  stringsAsFactors = FALSE
)

pos_filter_log_tax <- pos_filter_log[pos_filter_log$method == "nP.lulu", ]
pos_filter_log_tax <- cbind(
  pos_filter_log_tax,
  tax_raw[pos_filter_log_tax$otu, , drop = FALSE]
)

# Quick taxonomy summary of what was removed
cat("Phylum breakdown of removed OTUs:\n")
print(table(pos_filter_log_tax$Genus, useNA = "ifany"))

