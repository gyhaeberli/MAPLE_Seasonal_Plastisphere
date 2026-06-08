################################################################################
# GLLVM FINAL FITTING — WATER TRACK, PA MODEL (M2 SITE × SEASON, 2 LV)
#
# PURPOSE:
#   Self-contained script that prepares the 18S water community data and fits
#   a GLLVM on water samples alone.
#
#   KEY DIFFERENCE FROM BIOFILM MODEL:
#   Water samples have no substrate variable, so the formula uses only
#   site and season: ~ (site + season)
#
# DATA SOURCE:
#   ps_all$pr2_rep50$fP.lulu$Water  (PCR replicates, uncollapsed)
#   Filtered to seasonal sampling events (3, 5, 7, 9, 11) only.
#
# FILTERING STEPS (identical logic to biofilm script):
#   Step 1   — Library size filter (>= 2000 reads)
#   Step 2   — OTU prevalence (>= 5%) / top-500 abundance filter
#   Step 2b  — OTU degeneracy filter (>= 3 OTUs per PCR replicate)
#   Step 2c  — Per-factor-level presence filter (>= 2 levels of season & site)
#   Step 3   — Factor level assignment with reference levels
#   Step 4   — Readable OTU labels
#   Step 5   — Random effect structure (PCR replicate grouping)
#
# MODEL:
#   Family:  binomial (presence/absence)
#   Formula: ~ (site + season)^2
#   LV:      num_lv = 2
#   Method:  VA, probit link
#   Random:  row.eff = ~(1 | clean_sample_names)
#
################################################################################

library(gllvm)
library(phyloseq)
library(dplyr)

cat("=== GLLVM FINAL — WATER TRACK / PA / SITE × SEASON / LV2 ===\n\n")

setwd("~/GitHub/MAPLE_Seasonal_Plastisphere/Scripts/03.StatisticalAnalysis")

SAVE_DIR_BASE <- "~/Github/MAPLE_Seasonal_Plastisphere/Processed_data/gllvm_models"

# ==============================================================================
# HARDCODED PARAMETERS
# ==============================================================================

# --- Filtering ---
MIN_LIB_SIZE        <- 2000
PREVALENCE_THRESH   <- 0.05
TOP_N_OTUS          <- 500
MIN_BIO_SAMPLES     <- 2
MIN_OTUS_PER_SAMPLE <- 3
FILTER_FACTORS      <- c("season", "site")   # no substrate for water model
MIN_LEVELS_PRESENT  <- 2

# --- Response matrix ---
N_OTUS <- 200

# --- Seasonal sampling events ---
# All ST water events (deployment + retrieval)
st_events <- c(1, 3, 5, 6, 7, 8, 9, 10, 11, 12)

# Season lookup covering all ST events
season_lookup <- c(
  "1"  = "Winter",  "3"  = "Winter",
  "5"  = "Spring",  "6"  = "Spring",
  "7"  = "Summer",  "8"  = "Summer",
  "9"  = "Fall",    "10" = "Fall",
  "11" = "Winter2", "12" = "Winter2"
)

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

ps_wf_raw <- ps_all$pr2_rep50$fP.lulu$Water

rm(ps_all); gc()

# Filter to seasonal sampling events
wf_meta        <- as.data.frame(sample_data(ps_wf_raw))
class(wf_meta) <- "data.frame"
samples_seasonal <- rownames(wf_meta)[wf_meta$sampling_event %in% st_events]
ps_wf          <- prune_samples(samples_seasonal, ps_wf_raw)
ps_wf          <- prune_taxa(taxa_sums(ps_wf) > 0, ps_wf)

cat("  WF (seasonal):", nsamples(ps_wf), "samples,", ntaxa(ps_wf), "taxa\n")

# Add season column (needed downstream)
sd_wf         <- as.data.frame(sample_data(ps_wf))
class(sd_wf)  <- "data.frame"

# Season already in metadata from index.info() — just verify
cat("  Season distribution:\n")
print(table(sd_wf$season, useNA = "always"))
cat("  Role distribution:\n")
print(table(sd_wf$st_role, useNA = "always"))

cat("  Season distribution:\n")
print(table(sd_wf$season))
cat("\n")


# ==============================================================================
# SECTION 2: DATA PREPARATION
# ==============================================================================

otu_mat  <- as(otu_table(ps_wf), "matrix")
if (taxa_are_rows(ps_wf)) otu_mat <- t(otu_mat)

metadata <- as.data.frame(sample_data(ps_wf))
class(metadata) <- "data.frame"
tax_mat  <- as.data.frame(as(tax_table(ps_wf), "matrix"))

cat("  Original:", nrow(otu_mat), "samples x", ncol(otu_mat), "OTUs\n\n")


# ── Step 1: Library size filter ───────────────────────────────────────────────

cat("--- Step 1: Library size filter (>=", MIN_LIB_SIZE, "reads) ---\n")

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

top_abund <- names(sort(colSums(otu_mat), decreasing = TRUE)[1:min(TOP_N_OTUS, ncol(otu_mat))])
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
degenerate     <- n_otus_per_pcr < MIN_OTUS_PER_SAMPLE

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


# ── Step 2c: Per-factor-level presence filter ─────────────────────────────────

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
  cat("  Removed taxa:\n")
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

metadata$season <- factor(metadata$season,
                          levels = c("Winter", "Spring", "Summer",
                                     "Fall", "Winter2"))
metadata$site   <- factor(metadata$site, levels = c("SELVA", "TBS"))

cat("  Season ref: Winter | Site ref: SELVA\n\n")


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
# ==============================================================================

cat("--- Building PA response matrix (top", N_OTUS, "OTUs by abundance) ---\n")

n_select  <- min(N_OTUS, ncol(otu_filt))
top_otus  <- names(sort(colSums(otu_filt), decreasing = TRUE)[1:n_select])
count_mat <- otu_filt[, top_otus]
colnames(count_mat) <- otu_id_to_label[top_otus]
pa_mat <- (count_mat > 0) * 1L

cat("  PA matrix:", nrow(pa_mat), "x", ncol(pa_mat),
    "| zeros:", round(mean(pa_mat == 0) * 100, 1), "%\n\n")


# ==============================================================================
# SECTION 4: SAVE ENVIRONMENT FOR HPC FITTING
# ==============================================================================

cat("--- Saving environment for HPC fitting ---\n")

env_water <- list(
  # Matrices
  pa_mat          = pa_mat,
  count_mat       = count_mat,
  otu_filt        = otu_filt,
  
  # Metadata and design
  metadata        = metadata,
  X               = metadata[, c("site", "season")],
  study_design    = study_design,
  
  # Taxonomy
  tax_mat         = tax_mat,
  otu_id_to_label = otu_id_to_label,
  top_otus        = top_otus,
  
  # Filtering artefacts
  removed_step2c  = removed_step2c,
  
  # Parameters (for reference)
  params = list(
    NUM_LV              = NUM_LV,
    SEED                = SEED,
    N_OTUS              = N_OTUS,
    MIN_LIB_SIZE        = MIN_LIB_SIZE,
    PREVALENCE_THRESH   = PREVALENCE_THRESH,
    TOP_N_OTUS          = TOP_N_OTUS,
    MIN_OTUS_PER_SAMPLE = MIN_OTUS_PER_SAMPLE,
    FILTER_FACTORS      = FILTER_FACTORS,
    MIN_LEVELS_PRESENT  = MIN_LEVELS_PRESENT
  )
)

saveRDS(env_water,
        file.path(SAVE_DIR_BASE, "env_water_WF_for_HPC.rds"))

cat("  Saved: env_water_WF_for_HPC.rds\n\n")


# ==============================================================================
# SECTION 5: GLLVM FITTING
# Run this section on HPC after loading the saved environment above.
# To resume on HPC:
# env <- readRDS(".../env_water_WF_for_HPC.rds")
 env = readRDS(file.path(SAVE_DIR_BASE, "env_water_WF_for_HPC.rds"))
list2env(env, envir = .GlobalEnv)
 list2env(params, envir = .GlobalEnv)
# Then run the gllvm() calls below.
# ==============================================================================

cat("--- Fitting final model: PA / site + season (additive) / LV2 ---\n\n")

set.seed(SEED)

fit_water <- gllvm(
  y           = pa_mat,
  X           = X,
  formula     = ~ site + season,
  family      = "binomial",
  num.lv      = NUM_LV,
  studyDesign = study_design,
  row.eff     = ~(1 | clean_sample_names),
  method      = "VA",
  link        = "probit",
  control     = list(reltol = 1e-6),
  seed        = SEED
)

saveRDS(fit_water, "GLLVM_final_water_WF_add.rds")

cat("\n--- Water model fitted successfully ---\n")
cat("  AIC:    ", round(AIC(fit_water), 1), "\n")
cat("  BIC:    ", round(BIC(fit_water), 1), "\n")
cat("  AICc:   ", round(summary(fit_water)$AICc, 1), "\n")
cat("  logLik: ", round(as.numeric(logLik(fit_water)), 1), "\n")
cat("  Xcoef range:",
    paste(round(range(coef(fit_water)$Xcoef), 2), collapse = " to "), "\n\n")


# ── Null model ────────────────────────────────────────────────────────────────

cat("--- Fitting null model: PA / no fixed effects / LV2 ---\n\n")

set.seed(SEED)

fit_water_null <- gllvm(
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

saveRDS(fit_water_null,"GLLVM_null_water_WF.rds")

cat("\n--- Null model fitted successfully ---\n")
cat("  AIC:    ", round(AIC(fit_water_null), 1), "\n")
cat("  BIC:    ", round(BIC(fit_water_null), 1), "\n")
cat("  logLik: ", round(as.numeric(logLik(fit_water_null)), 1), "\n\n")


# ── Sanity check ──────────────────────────────────────────────────────────────

lv_null <- sum(apply(getLV(fit_water_null), 2, var))
lv_best <- sum(apply(getLV(fit_water),      2, var))
cat("--- Sanity check ---\n")
cat("  LV variance — null:", round(lv_null, 3),
    "| best:", round(lv_best, 3),
    "| reduction:", round(100 * (lv_null - lv_best) / lv_null, 1), "%\n\n")

cat("Script complete.\n")