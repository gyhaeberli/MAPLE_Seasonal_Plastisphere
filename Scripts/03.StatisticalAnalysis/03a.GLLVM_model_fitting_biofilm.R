################################################################################
# GLLVM FINAL FITTING — FOULING TRACK, PA MODEL (M2 ALL-2-WAY, 2 LV)
#
# Model was selected in GLLVM_models_fit_v2.R
#
# PURPOSE:
#   Self-contained script that prepares the 18S community data and fits
#   the selected final model:
#     - Track:   full dataset
#     - Family:  binomial (presence/absence incidence matrix)
#     - Formula: M2 — all 2-way interactions: (site + season + substrate)^2
#     - LV:      num_lv = 2
#     - Method:  VA, probit link
#
#   All filtering parameters are hardcoded to the values used in the v2
#   screening pipeline. No offset is used (PA model on binary response).
#
# FILTERING STEPS (identical to prepare_gllvm_data_v2):
#   Step 1   — Library size filter (>= 2000 reads)
#   Step 2   — OTU prevalence (>= 5%) / top-500 abundance filter
#   Step 2b  — OTU degeneracy filter (>= 3 OTUs per PCR replicate)
#   Step 2c  — Per-factor-level presence filter (>= 2 levels of season & site)
#   Step 3   — Factor level assignment with reference levels
#   Step 4   — Readable OTU labels
#   Step 5   — Random effect structure (PCR replicate grouping)
#
################################################################################

library(gllvm)
library(phyloseq)
library(dplyr)

cat("=== GLLVM FINAL SCRIPT — FOULING / PA / M2 / LV2 ===\n\n")

setwd("~/GitHub/MAPLE_Seasonal_Plastisphere/Scripts/03.StatisticalAnalysis")


SAVE_DIR_BASE  <- "~/Github/MAPLE_Seasonal_Plastisphere/Processed_data/gllvm_models"
OUT_DIR <- "~/Github/MAPLE_Seasonal_Plastisphere/Results/gllvm_results"

# ==============================================================================
# HARDCODED PARAMETERS
# ==============================================================================

# --- Filtering ---
MIN_LIB_SIZE        <- 2000   # Step 1:  minimum reads per sample
PREVALENCE_THRESH   <- 0.05   # Step 2:  minimum prevalence across all samples
TOP_N_OTUS          <- 500    # Step 2:  top OTUs by total abundance to retain
MIN_BIO_SAMPLES     <- 2      # Step 2:  min biological samples with detection
MIN_OTUS_PER_SAMPLE <- 3      # Step 2b: min OTUs per PCR replicate
FILTER_FACTORS      <- c("season", "site", "substrate")  # Step 2c: factors to check
MIN_LEVELS_PRESENT  <- 2      # Step 2c: min distinct factor levels with detection

# --- Response matrix ---
N_OTUS <- 200   # top OTUs by total abundance used in model

# --- Model ---
SEED   <- 123
NUM_LV <- 2


# ==============================================================================
# SECTION 1: LOAD DATA
# ==============================================================================

cat("--- Loading phyloseq objects ---\n")

ps_all <- readRDS(
  "~/GitHub/MAPLE_Seasonal_Plastisphere/Processed_data/Phyloseq_objects/FINAL_PHYLOSEQ_OBJECTS.rds"
)

# fouling track: full eukaryotic community, PCR replicates uncollapsed
ps <- ps_all$pr2_rep50$fP.lulu$ST
rm(ps_all); gc()

cat("Loaded phyloseq object for fouling track.\n\n")


# ==============================================================================
# SECTION 2: DATA PREPARATION
# Replicates prepare_gllvm_data_v2() exactly, steps 1-5
# ==============================================================================

# ── Step 1: Library size filter ───────────────────────────────────────────────

cat("--- Step 1: Library size filter (>=", MIN_LIB_SIZE, "reads) ---\n")

otu_mat  <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)

metadata <- as.data.frame(sample_data(ps))
class(metadata) <- "data.frame"
tax_mat  <- as.data.frame(as(tax_table(ps), "matrix"))

cat("  Original:", nrow(otu_mat), "samples x", ncol(otu_mat), "OTUs\n")

library_sizes <- rowSums(otu_mat)
cat("  Library size range:", min(library_sizes), "-", max(library_sizes), "\n")
cat("  Mean:", round(mean(library_sizes)),
    "| CV:", round(sd(library_sizes) / mean(library_sizes), 3), "\n")

keep <- library_sizes >= MIN_LIB_SIZE
cat("  Removed (<", MIN_LIB_SIZE, "reads):", sum(!keep),
    "| Retained:", sum(keep), "\n\n")

otu_mat       <- otu_mat[keep, ]
metadata      <- metadata[keep, ]
library_sizes <- library_sizes[keep]

otu_mat  <- otu_mat[, colSums(otu_mat) > 0]
tax_mat  <- tax_mat[colnames(otu_mat), , drop = FALSE]
keep_s   <- rowSums(otu_mat) > 0
otu_mat  <- otu_mat[keep_s, ]
metadata <- metadata[keep_s, ]
library_sizes <- library_sizes[keep_s]

# ── Step 2: OTU prevalence / abundance filter ─────────────────────────────────

cat("--- Step 2: OTU prevalence / abundance filter ---\n")

prevalence <- colSums(otu_mat > 0) / nrow(otu_mat)
abund_otus <- prevalence >= PREVALENCE_THRESH

bio_present <- sapply(colnames(otu_mat), function(otu) {
  length(unique(metadata$clean_sample_names[otu_mat[, otu] > 0]))
})

top_abund <- names(sort(colSums(otu_mat), decreasing = TRUE)[1:TOP_N_OTUS])
keep_otus <- abund_otus | (colnames(otu_mat) %in% top_abund & bio_present >= MIN_BIO_SAMPLES)

otu_filt <- otu_mat[, keep_otus]
cat("  OTUs retained:", sum(keep_otus), "| Removed:", sum(!keep_otus), "\n")

zero_s <- rowSums(otu_filt) == 0
if (sum(zero_s) > 0) {
  cat("  WARNING:", sum(zero_s), "zero-sum samples removed\n")
  otu_filt      <- otu_filt[!zero_s, ]
  metadata      <- metadata[!zero_s, ]
  library_sizes <- library_sizes[!zero_s]
}

# ── Step 2b: OTU degeneracy filter ───────────────────────────────────────────

cat("--- Step 2b: OTU degeneracy filter (threshold =", MIN_OTUS_PER_SAMPLE, ") ---\n")

n_otus_per_pcr <- rowSums(otu_filt > 0)



# Show the distribution BEFORE removing anything
cat("  Distribution of OTU counts per PCR replicate:\n")
print(table(n_otus_per_pcr))

cat("\n  Replicates with < 3 OTUs:", sum(n_otus_per_pcr < MIN_OTUS_PER_SAMPLE),
    "out of", length(n_otus_per_pcr),
    paste0("(", round(100 * mean(n_otus_per_pcr < MIN_OTUS_PER_SAMPLE), 1), "%)\n"))

degenerate <- n_otus_per_pcr < MIN_OTUS_PER_SAMPLE


if (sum(degenerate) > 0) {
  cat("  Removed", sum(degenerate), "PCR replicates with <",
      MIN_OTUS_PER_SAMPLE, "OTUs:\n")
  print(data.frame(
    sample = rownames(otu_filt)[degenerate],
    n_otus = n_otus_per_pcr[degenerate],
    reads  = rowSums(otu_filt)[degenerate]
  ))
  otu_filt      <- otu_filt[!degenerate, ]
  metadata      <- metadata[!degenerate, ]
  library_sizes <- library_sizes[!degenerate]
} else {
  cat("  No degenerate samples found.\n")
}

# ── Step 2c: Per-factor-level presence filter (v2 addition) ──────────────────

cat("--- Step 2c: Per-factor-level presence filter ---\n")
cat("  Factors checked:", paste(FILTER_FACTORS, collapse = ", "),
    "| Min levels required:", MIN_LEVELS_PRESENT, "\n")

keep_level <- rep(TRUE, ncol(otu_filt))
names(keep_level) <- colnames(otu_filt)

for (fac in FILTER_FACTORS) {
  if (!fac %in% colnames(metadata)) {
    cat("  WARNING: factor '", fac, "' not found in metadata — skipping\n")
    next
  }
  fac_vec <- metadata[[fac]]
  n_levels_present <- apply(otu_filt, 2, function(x) {
    length(unique(fac_vec[x > 0]))
  })
  fails_this_factor <- n_levels_present < MIN_LEVELS_PRESENT
  cat("  Factor '", fac, "': removing", sum(fails_this_factor),
      "OTUs present in <", MIN_LEVELS_PRESENT, "levels\n")
  keep_level <- keep_level & !fails_this_factor
}

cat("  OTUs retained after Step 2c:", sum(keep_level),
    "| Total removed:", sum(!keep_level), "\n")

removed_step2c <- colnames(otu_filt)[!keep_level]
if (length(removed_step2c) > 0) {
  cat("  Removed taxa (candidates for indicator species analysis):\n")
  cat("   ", paste(head(removed_step2c, 20), collapse = ", "))
  if (length(removed_step2c) > 20)
    cat(" ... and", length(removed_step2c) - 20, "more")
  cat("\n")
}

otu_filt <- otu_filt[, keep_level]

zero_s2c <- rowSums(otu_filt) == 0
if (sum(zero_s2c) > 0) {
  cat("  WARNING:", sum(zero_s2c),
      "zero-sum samples created by Step 2c — removed\n")
  otu_filt      <- otu_filt[!zero_s2c, ]
  metadata      <- metadata[!zero_s2c, ]
  library_sizes <- library_sizes[!zero_s2c]
}

lib_sizes_filt <- rowSums(otu_filt)
cat("  Final matrix:", nrow(otu_filt), "samples x", ncol(otu_filt), "OTUs\n")
cat("  Library size range (post-filter):",
    min(lib_sizes_filt), "-", max(lib_sizes_filt), "\n\n")

# ── Step 3: Factor levels ─────────────────────────────────────────────────────

cat("--- Step 3: Factor levels ---\n")

metadata$season    <- factor(metadata$season,
                             levels = c("Winter", "Spring", "Summer", "Fall", "Winter2"))
metadata$substrate <- factor(metadata$substrate,
                             levels = c("Glass", "PE", "Weathered_PE", "PET", "Weathered_PET"))
metadata$site      <- factor(metadata$site, levels = c("SELVA", "TBS"))

cat("  Season ref: Winter | Substrate ref: Glass | Site ref: SELVA\n\n")

# ── Step 4: Readable OTU labels ──────────────────────────────────────────────

cat("--- Step 4: OTU labels ---\n")

make_label <- function(id, tm) {
  for (rank in c("Genus", "Family", "Class")) {
    v <- tm[id, rank]
    if (!is.na(v) && nchar(v) > 0 && !v %in% c("uncultured", "Unknown", ""))
      return(paste0(v, "_", substr(id, 1, 12)))
  }
  paste0("OTU_", substr(id, 1, 12))
}

otu_id_to_label <- setNames(
  sapply(colnames(otu_filt), make_label, tm = tax_mat),
  colnames(otu_filt)
)
cat("  Labels created:", length(otu_id_to_label), "\n\n")

# ── Step 5: Random effect structure ──────────────────────────────────────────

cat("--- Step 5: Random effect structure ---\n")

metadata$clean_sample_names <- as.factor(metadata$clean_sample_names)
study_design <- data.frame(clean_sample_names = metadata$clean_sample_names)

cat("  Groups:", nlevels(metadata$clean_sample_names), "\n")
cat("  Reps per group:",
    paste(range(table(metadata$clean_sample_names)), collapse = "-"), "\n\n")

# ==============================================================================
# SECTION 3: RESPONSE MATRIX
# Replicates prepare_response_matrices_v2() PA matrix exactly
# ==============================================================================

cat("--- Building PA response matrix (top", N_OTUS, "OTUs by abundance) ---\n")

top_otus  <- names(sort(colSums(otu_filt), decreasing = TRUE)[1:N_OTUS])
count_mat <- otu_filt[, top_otus]
colnames(count_mat) <- otu_id_to_label[top_otus]
pa_mat <- (count_mat > 0) * 1L

cat("  PA matrix:", nrow(pa_mat), "x", ncol(pa_mat),
    "| zeros:", round(mean(pa_mat == 0) * 100, 1), "%\n\n")



# ==============================================================================
# SECTION 4: MODEL FITTING
#   Family:  binomial (PA)
#   Formula: M2 — (site + season + substrate)^2  (all 2-way interactions)
#   LV:      num_lv = 2
#   Method:  VA, probit link
#   No offset (binary response)
# ==============================================================================

cat("--- Fitting final model: PA / M2 (all 2-way) / LV2 ---\n\n")

X <- metadata[, c("site", "season", "substrate")]

set.seed(SEED)

fit_final <- gllvm(
  y           = pa_mat,
  X           = X,
  formula     = ~ (site + season + substrate)^2,
  family      = "binomial",
  num.lv      = NUM_LV,
  studyDesign = study_design,
  row.eff     = ~(1 | clean_sample_names),
  method      = "VA",
  link        = "probit",
  control     = list(reltol = 1e-6),
  seed        = SEED
)

saveRDS(fit_final, file.path(SAVE_DIR_BASE, "GLLVM_final_incM2.rds"))

fit_final <- readRDS(file.path(SAVE_DIR_BASE, "GLLVM_final_incM2.rds"))

cat("\n--- Best model fitted successfully ---\n")
cat("  AIC:    ", round(AIC(fit_final), 1), "\n")
cat("  BIC:    ", round(BIC(fit_final), 1), "\n")
cat("  AICc:   ", round(summary(fit_final)$AICc, 1), "\n")
cat("  logLik: ", round(as.numeric(logLik(fit_final)), 1), "\n")
cat("  Xcoef range:", paste(round(range(coef(fit_final)$Xcoef), 2), collapse = " to "), "\n\n")


# ==============================================================================
# SECTION 5: NULL / UNCONSTRAINED MODEL (M7)
#   No fixed effects — latent variables only
#   Used as baseline for variance partitioning and to assess how much
#   community structure is explained by the fixed effects in M2
# ==============================================================================

cat("--- Fitting null model: PA / M7 (no fixed effects) / LV2 ---\n\n")

set.seed(SEED)

fit_null <- gllvm(
  y           = pa_mat,
  family      = "binomial",
  num.lv      = NUM_LV,
  studyDesign = study_design,
  row.eff     = ~(1 | clean_sample_names),
  method      = "VA",
  link        = "probit",
  control     = list(reltol = 1e-6),
  seed        = SEED
)

saveRDS(fit_null, file.path(SAVE_DIR_BASE, "GLLVM_null_incM7.rds"))

fit_null  <- readRDS(file.path(SAVE_DIR_BASE, "GLLVM_null_incM7.rds"))

cat("\n--- Null model fitted successfully ---\n")
cat("  AIC:    ", round(AIC(fit_null), 1), "\n")
cat("  BIC:    ", round(BIC(fit_null), 1), "\n")
cat("  logLik: ", round(as.numeric(logLik(fit_null)), 1), "\n\n")


# Sanity check — should match best_inc_f_v2 / uncons_inc_f_v2 results
lv_null <- sum(apply(getLV(fit_null), 2, var))
lv_best <- sum(apply(getLV(fit_final), 2, var))
cat("--- Sanity check ---\n")
cat("  LV variance — null:", round(lv_null, 3),
    "| best:", round(lv_best, 3),
    "| reduction:", round(100 * (lv_null - lv_best) / lv_null, 1), "%\n\n")



#save.image("gllvm_final_fit.RData")

#load("gllvm_final_fit.RData")




# ==============================================================================
#  DOES SITE, SEASON AND SUBSTRATE STRUCTURE THE COMMUNITY?
# ==============================================================================

cat("=== Q1: Community structuring — variance partitioning ===\n\n")


### FROM GLLVM PACKAGE

vp_gllvm = VP(fit_final)
vp_gllvm



#TABLE
# ── Build VP table ────────────────────────────────────────────────────────────
vp_values <- colMeans(vp_gllvm$PropExplainedVarSp) * 100

vp_display <- data.frame(
  Component = c(
    "Site",
    "Season",
    "Substrate",
    "Site × Season",
    "Site × Substrate",
    "Season × Substrate",
    "Latent variable 1",
    "Latent variable 2",
    "Sample random effect"
  ),
  Type = c(
    "Main effect",
    "Main effect",
    "Main effect",
    "Interaction",
    "Interaction",
    "Interaction",
    "Residual",
    "Residual",
    "Residual"
  ),
  Variance_explained = round(as.numeric(vp_values), 1),
  stringsAsFactors = FALSE
)

colnames(vp_display) <- c("Component", "Type", "Variance explained (%)")

# ── Render ────────────────────────────────────────────────────────────────────
kbl(vp_display,
    booktabs = TRUE,
    align    = c("l", "l", "r"),
    na       = "") %>%
  kable_classic(full_width = FALSE, html_font = "Arial") %>%
  pack_rows("Main effects", 1, 3) %>%
  pack_rows("Interactions", 4, 6) %>%
  pack_rows("Residual", 7, 9) %>%
  row_spec(7:9, italic = TRUE, color = "gray") 
#row_spec(4:6, bold = TRUE) %>%



ft <- flextable(vp_display) %>%
  autofit()

read_docx() %>%
  body_add_flextable(ft) %>%
  print(target = file.path(OUT_DIR, "Tables/vp_gllvm.docx"))

print(doc, target = "vp_gllvm.docx")




# ==============================================================================
# SECTION 6: SAVE DATA OBJECT FOR ANALYSIS SCRIPT
# ==============================================================================

model_data_final <- list(
  otu_filt        = otu_filt,
  tax_mat         = tax_mat,
  metadata        = metadata,
  otu_id_to_label = otu_id_to_label,
  study_design    = study_design,
  design_factors  = c("site", "season", "substrate"),
  pa_mat          = pa_mat,
  top_otus        = top_otus,
  count_mat       = count_mat,
  removed_step2c  = removed_step2c
)

saveRDS(model_data_final,
        file.path(SAVE_DIR_BASE, "data_for_gllvm_FINAL.rds"))
cat("Model is fitted and saved.\n")

