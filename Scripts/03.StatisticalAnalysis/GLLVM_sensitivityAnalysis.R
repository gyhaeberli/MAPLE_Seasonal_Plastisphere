################################################################################
# GLLVM SENSITIVITY ANALYSIS — TOP-N OTU THRESHOLD
#
# PURPOSE:
#   Refits the FINAL selected model (M2, all 2-way interactions, 2 LV,
#   binomial/probit/VA) using different numbers of top OTUs (by total
#   abundance) as the response matrix, to check whether results (variance
#   partitioning, model fit) are sensitive to the choice of N_OTUS = 200
#   used in the main analysis.
#
#   Values tested: 50, 100, 300, 400
#   (200 is your original final model — rerun separately / already saved)
#
#   Data preparation (Steps 1-2c) is IDENTICAL to your final script and is
#   only run ONCE, since it does not depend on N_OTUS. N_OTUS only affects
#   Section 3 (building the response matrix) onward.
#
################################################################################

library(gllvm)
library(phyloseq)
library(dplyr)

cat("=== GLLVM SENSITIVITY ANALYSIS — TOP-N OTU THRESHOLD ===\n\n")

setwd("~/GitHub/MAPLE_Seasonal_Plastisphere/Scripts/03.StatisticalAnalysis")

SAVE_DIR_BASE <- "~/Github/MAPLE_Seasonal_Plastisphere/Processed_data/gllvm_models"
OUT_DIR       <- "~/Github/MAPLE_Seasonal_Plastisphere/Results/gllvm_results/sensitivity_analysis"

# Create a subfolder just for sensitivity outputs, so they don't mix with
# your final model files
SENS_DIR <- file.path(SAVE_DIR_BASE, "sensitivity_topN")
dir.create(SENS_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# HARDCODED PARAMETERS (same as final script)
# ==============================================================================

# --- Filtering ---
MIN_LIB_SIZE        <- 2000
PREVALENCE_THRESH   <- 0.05
TOP_N_OTUS          <- 500     # candidate pool size (Step 2) -- NOT the same as N_OTUS below
MIN_BIO_SAMPLES     <- 2
MIN_OTUS_PER_SAMPLE <- 3
FILTER_FACTORS      <- c("season", "site", "substrate")
MIN_LEVELS_PRESENT  <- 2

# --- Sensitivity: list of top-N OTU thresholds to test ---
N_OTUS_LIST <- c(50, 100, 300, 400)

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

ps <- ps_all$pr2_rep50$fP.lulu$ST
rm(ps_all); gc()

cat("Loaded phyloseq object for fouling track.\n\n")


# ==============================================================================
# SECTION 2: DATA PREPARATION (identical to final script, run once)
# ==============================================================================

# ── Step 1: Library size filter ──────────────────────────────────────────────

cat("--- Step 1: Library size filter (>=", MIN_LIB_SIZE, "reads) ---\n")

otu_mat  <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)

metadata <- as.data.frame(sample_data(ps))
class(metadata) <- "data.frame"
tax_mat  <- as.data.frame(as(tax_table(ps), "matrix"))

library_sizes <- rowSums(otu_mat)
keep <- library_sizes >= MIN_LIB_SIZE
cat("  Removed (<", MIN_LIB_SIZE, "reads):", sum(!keep), "| Retained:", sum(keep), "\n\n")

otu_mat       <- otu_mat[keep, ]
metadata      <- metadata[keep, ]
library_sizes <- library_sizes[keep]

otu_mat  <- otu_mat[, colSums(otu_mat) > 0]
tax_mat  <- tax_mat[colnames(otu_mat), , drop = FALSE]
keep_s   <- rowSums(otu_mat) > 0
otu_mat  <- otu_mat[keep_s, ]
metadata <- metadata[keep_s, ]
library_sizes <- library_sizes[keep_s]

# ── Step 2: OTU prevalence / abundance filter ────────────────────────────────

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
  otu_filt      <- otu_filt[!zero_s, ]
  metadata      <- metadata[!zero_s, ]
  library_sizes <- library_sizes[!zero_s]
}

# ── Step 2b: OTU degeneracy filter ───────────────────────────────────────────

cat("--- Step 2b: OTU degeneracy filter (threshold =", MIN_OTUS_PER_SAMPLE, ") ---\n")

n_otus_per_pcr <- rowSums(otu_filt > 0)
degenerate <- n_otus_per_pcr < MIN_OTUS_PER_SAMPLE

if (sum(degenerate) > 0) {
  cat("  Removed", sum(degenerate), "PCR replicates with <", MIN_OTUS_PER_SAMPLE, "OTUs\n")
  otu_filt      <- otu_filt[!degenerate, ]
  metadata      <- metadata[!degenerate, ]
  library_sizes <- library_sizes[!degenerate]
} else {
  cat("  No degenerate samples found.\n")
}

# ── Step 2c: Per-factor-level presence filter ────────────────────────────────

cat("--- Step 2c: Per-factor-level presence filter ---\n")

keep_level <- rep(TRUE, ncol(otu_filt))
names(keep_level) <- colnames(otu_filt)

for (fac in FILTER_FACTORS) {
  if (!fac %in% colnames(metadata)) next
  fac_vec <- metadata[[fac]]
  n_levels_present <- apply(otu_filt, 2, function(x) length(unique(fac_vec[x > 0])))
  fails_this_factor <- n_levels_present < MIN_LEVELS_PRESENT
  cat("  Factor '", fac, "': removing", sum(fails_this_factor), "OTUs\n")
  keep_level <- keep_level & !fails_this_factor
}

otu_filt <- otu_filt[, keep_level]

zero_s2c <- rowSums(otu_filt) == 0
if (sum(zero_s2c) > 0) {
  otu_filt      <- otu_filt[!zero_s2c, ]
  metadata      <- metadata[!zero_s2c, ]
  library_sizes <- library_sizes[!zero_s2c]
}

cat("  Final filtered matrix:", nrow(otu_filt), "samples x", ncol(otu_filt), "OTUs\n\n")

# ── Step 3: Factor levels ────────────────────────────────────────────────────

metadata$season    <- factor(metadata$season,
                             levels = c("Winter", "Spring", "Summer", "Fall", "Winter2"))
metadata$substrate <- factor(metadata$substrate,
                             levels = c("Glass", "PE", "Weathered_PE", "PET", "Weathered_PET"))
metadata$site      <- factor(metadata$site, levels = c("SELVA", "TBS"))

# ── Step 4: Readable OTU labels ──────────────────────────────────────────────

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

# ── Step 5: Random effect structure ──────────────────────────────────────────

metadata$clean_sample_names <- as.factor(metadata$clean_sample_names)
study_design <- data.frame(clean_sample_names = metadata$clean_sample_names)

X <- metadata[, c("site", "season", "substrate")]

cat("Data preparation complete. This filtered dataset is reused for every\n")
cat("N_OTUS value below — only the response matrix (top-N columns) changes.\n\n")


save.image("gllvm_sensivity_check_dataprep.RData")


# ==============================================================================
# SECTION 3: SENSITIVITY FUNCTION — REFIT MODEL FOR A GIVEN N_OTUS
# ==============================================================================

library(gllvm)
library(phyloseq)
library(dplyr)

cat("=== GLLVM SENSITIVITY ANALYSIS — TOP-N OTU THRESHOLD - on ASTBURY HPC ===\n\n")

load("gllvm_sensivity_check_dataprep.RData")

SAVE_DIR_BASE <- "/home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/STATISTICS"
OUT_DIR       <- "/home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/STATISTICS/sensitivity_analysis"

# Create a subfolder just for sensitivity outputs, so they don't mix with
# your final model files
SENS_DIR <- file.path(SAVE_DIR_BASE, "sensitivity_topN")
dir.create(SENS_DIR, showWarnings = FALSE, recursive = TRUE)


# fit_sensitivity_model():
#   Takes a single N_OTUS value, builds the top-N presence/absence matrix
#   from the already-filtered otu_filt, fits the same M2 model, saves the
#   fitted object to disk, and returns a list of summary stats (AIC, BIC,
#   logLik, variance partitioning). Returns NULL if the fit fails, so a
#   failed call doesn't stop the others.
fit_sensitivity_model <- function(n) {
  
  cat("\n================================================================\n")
  cat(" Fitting model with top", n, "OTUs\n")
  cat("================================================================\n\n")
  
  # --- Safety check: make sure we have enough OTUs to select from ---
  if (n > ncol(otu_filt)) {
    cat("  !! Requested N_OTUS =", n, "exceeds available OTUs (",
        ncol(otu_filt), ") after filtering. Skipping.\n")
    return(NULL)
  }
  
  # --- Build PA response matrix for this N ---
  top_otus_n  <- names(sort(colSums(otu_filt), decreasing = TRUE)[1:n])
  count_mat_n <- otu_filt[, top_otus_n]
  colnames(count_mat_n) <- otu_id_to_label[top_otus_n]
  pa_mat_n <- (count_mat_n > 0) * 1L
  
  cat("  PA matrix:", nrow(pa_mat_n), "x", ncol(pa_mat_n),
      "| zeros:", round(mean(pa_mat_n == 0) * 100, 1), "%\n\n")
  
  # --- Fit the model (same formula/settings as the final model) ---
  set.seed(SEED)
  
  fit_n <- tryCatch({
    gllvm(
      y           = pa_mat_n,
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
  }, error = function(e) {
    cat("  !! Model fitting FAILED for N_OTUS =", n, ":", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(fit_n)) return(NULL)
  
  # --- Save the fitted model object ---
  fit_path <- file.path(SENS_DIR, paste0("GLLVM_sensitivity_top", n, ".rds"))
  saveRDS(fit_n, fit_path)
  cat("  Saved:", fit_path, "\n")
  
  # --- Variance partitioning for this fit ---
  vp_n <- tryCatch(VP(fit_n), error = function(e) {
    cat("  !! VP() failed for N_OTUS =", n, ":", conditionMessage(e), "\n")
    NULL
  })
  
  vp_values_n <- if (!is.null(vp_n)) colMeans(vp_n$PropExplainedVarSp) * 100 else NA
  
  cat("  AIC:", round(AIC(fit_n), 1),
      "| BIC:", round(BIC(fit_n), 1),
      "| logLik:", round(as.numeric(logLik(fit_n)), 1), "\n")
  
  # --- Return summary stats ---
  list(
    N_OTUS = n,
    n_taxa = ncol(pa_mat_n),
    AIC    = AIC(fit_n),
    BIC    = BIC(fit_n),
    logLik = as.numeric(logLik(fit_n)),
    VP     = vp_values_n
  )
}

# --- Call the function once per N_OTUS value (explicit calls, no loop) ---
result_50  <- fit_sensitivity_model(50)
result_100 <- fit_sensitivity_model(100)
result_300 <- fit_sensitivity_model(300)
result_400 <- fit_sensitivity_model(400)

# Collect into a named list for the summary section below.
# NULL entries (failed fits) are dropped automatically.
sensitivity_results <- Filter(
  Negate(is.null),
  list("50" = result_50, "100" = result_100, "300" = result_300, "400" = result_400)
)


# ==============================================================================
# SECTION 4: BUILD COMPARISON TABLE ACROSS N_OTUS VALUES
# ==============================================================================

cat("\n\n=== SENSITIVITY SUMMARY ===\n\n")

# Basic fit-statistics table (AIC/BIC/logLik per N_OTUS)
fit_summary <- do.call(rbind, lapply(sensitivity_results, function(r) {
  data.frame(
    N_OTUS = r$N_OTUS,
    n_taxa = r$n_taxa,
    AIC    = round(r$AIC, 1),
    BIC    = round(r$BIC, 1),
    logLik = round(r$logLik, 1)
  )
}))
rownames(fit_summary) <- NULL
print(fit_summary)

# Variance-partitioning comparison table (% variance per component, per N_OTUS)
# Each row = a VP component (site, season, substrate, interactions, LVs, random effect)
# Each column = one N_OTUS value tested
vp_table <- do.call(cbind, lapply(sensitivity_results, function(r) r$VP))
colnames(vp_table) <- paste0("N", sapply(sensitivity_results, function(r) r$N_OTUS))

cat("\n--- Variance explained (%) by component, across N_OTUS ---\n")
print(round(vp_table, 1))

# Save both tables
write.csv(fit_summary, file.path(OUT_DIR, "Tables/sensitivity_fit_summary.csv"), row.names = FALSE)
write.csv(vp_table,    file.path(OUT_DIR, "Tables/sensitivity_VP_comparison.csv"))

cat("\nSaved summary tables to:\n")
cat(" ", file.path(OUT_DIR, "Tables/sensitivity_fit_summary.csv"), "\n")
cat(" ", file.path(OUT_DIR, "Tables/sensitivity_VP_comparison.csv"), "\n")

cat("\n=== SENSITIVITY ANALYSIS COMPLETE ===\n")
