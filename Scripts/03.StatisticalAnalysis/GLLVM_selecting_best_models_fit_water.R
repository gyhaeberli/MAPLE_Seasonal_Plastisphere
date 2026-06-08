################################################################################
# GLLVM WATER TRACK — LV SELECTION + MODEL SCREENING + FINAL FIT
#
# PURPOSE:
#   Prepares the 18S water community data, then runs:
#
#   SECTION 5 — LV SELECTION
#     Tests num_lv = 0:3 on the M2 formula (~ (site + season)^2) using the
#     PA (binomial, probit) family. Prints an AIC/BIC diagnostic table and
#     saves the results as a CSV.
#
#   SECTION 6 — MODEL SCREENING (M1–M7)
#     Fits all seven covariate structures at the LV chosen in Section 5.
#     Because the water model has NO substrate variable, several biofilm
#     formulas collapse:
#
#       M1  site * season * substrate  →  site * season        (= M2, M4)
#       M2  (site+season+substrate)²   →  (site + season)^2    (= site*season)
#       M3  site + season * substrate  →  site + season        (= M5, M6)
#       M4  site * season + substrate  →  site * season        (= M1, M2)
#       M5  site * substrate + season  →  site + season        (= M3, M6)
#       M6  site + season + substrate  →  site + season        (= M3, M5)
#       M7  null (no fixed effects)    →  unchanged
#
#     All seven labels are kept for traceability. Models with identical
#     formulas will produce identical fits; the diagnostic table makes this
#     visible through identical AIC/BIC values.
#
#   SECTION 7 — FINAL FIT
#     Refits the best model (lowest AICc) at the selected LV, adds a null
#     model for LV-variance sanity check, and saves both.
#
# DATA SOURCE:
#   ps_all$pr2_rep50$fP.lulu$Water  (PCR replicates, uncollapsed)
#   Filtered to seasonal sampling events (3, 5, 7, 9, 11) only.
#
# FILTERING STEPS (identical to biofilm script):
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

cat("=== GLLVM WATER TRACK — LV SELECTION + MODEL SCREENING + FINAL FIT ===\n\n")

setwd("~/GitHub/MAPLE_Seasonal_Plastisphere/Scripts/03.StatisticalAnalysis")

SAVE_DIR_BASE      <- "~/Github/MAPLE_Seasonal_Plastisphere/Processed_data/gllvm_models"
SAVE_DIR_SELECTION <- file.path(SAVE_DIR_BASE, "water_selection")

dir.create(SAVE_DIR_SELECTION, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# HARDCODED PARAMETERS
# ==============================================================================

# --- Filtering ---
MIN_LIB_SIZE        <- 2000
PREVALENCE_THRESH   <- 0.05
TOP_N_OTUS          <- 500
MIN_BIO_SAMPLES     <- 2
MIN_OTUS_PER_SAMPLE <- 3
FILTER_FACTORS      <- c("season", "site")   # no substrate for water
MIN_LEVELS_PRESENT  <- 2

# --- Response matrix ---
N_OTUS <- 200

# --- Seasonal sampling events ---
st_events <- c(3, 5, 7, 9, 11)

season_lookup <- c(
  "3"  = "Winter",
  "5"  = "Spring",
  "7"  = "Summer",
  "9"  = "Fall",
  "11" = "Winter2"
)

# --- LV selection ---
LV_MAX  <- 3    # tests num_lv = 0, 1, 2, 3
SEED    <- 123

# --- Final model (set after inspecting Section 5 output) ---
# Update NUM_LV_FINAL before running Section 7.
# Default is 2 (matching the original water script); change if selection
# suggests otherwise.
NUM_LV_FINAL <- 2


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

# Add season column
sd_wf         <- as.data.frame(sample_data(ps_wf))
class(sd_wf)  <- "data.frame"
sd_wf$season  <- season_lookup[as.character(sd_wf$sampling_event)]
sample_data(ps_wf) <- sample_data(sd_wf)

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

# Convenience shorthand used by Sections 5–7
X <- metadata[, c("site", "season")]


# ==============================================================================
# SECTION 4: SAVE ENVIRONMENT FOR HPC FITTING
# ==============================================================================

cat("--- Saving environment for HPC fitting ---\n")

save.image("gllvm_water_model_selection.RData")





load("gllvm_water_model_selection.RData")
library(gllvm)
library(phyloseq)
library(dplyr)




# ==============================================================================
# SECTION 5: LV SELECTION
#
# Tests num_lv = 0, 1, 2, 3 on the M2 formula (~ (site + season)^2).
# Why M2? It is the intended final formula, so LV selection is done in the
# context of that fixed-effect structure rather than in a null model.
# Using the null model (no fixed effects) for LV selection would inflate
# the optimal LV count because the LVs would absorb variance that the
# covariates will explain in the final model.
#
# AICc is the primary criterion (corrects AIC for small sample sizes).
# BIC is reported as a secondary, more conservative criterion.
#
# num_lv = 0 is included as a baseline: a model with no latent variables
# at all. If AICc favours lv=0, the data do not support residual community
# structure beyond what the fixed effects explain.
# ==============================================================================

cat("=================================================================\n")
cat(" SECTION 5: LV SELECTION (num_lv = 0 to", LV_MAX, ")\n")
cat(" Formula: ~ (site + season)^2   Family: binomial (probit)\n")
cat("=================================================================\n\n")

load("gllvm_water_model_selection.RData")
library(gllvm)
library(phyloseq)
library(dplyr)

lv_results <- data.frame(
  num_lv     = 0:LV_MAX,
  AICc       = NA_real_,
  AIC        = NA_real_,
  BIC        = NA_real_,
  logLik     = NA_real_,
  converged  = FALSE
)

for (i in seq_along(lv_results$num_lv)) {
  
  lv_i <- lv_results$num_lv[i]
  cat("  Fitting num_lv =", lv_i, "...")
  
  set.seed(SEED)
  
  fit_lv <- tryCatch(
    gllvm(
      y           = pa_mat,
      X           = X,
      formula     = ~ (site + season)^2,
      family      = "binomial",
      num.lv      = lv_i,
      studyDesign = study_design,
      row.eff     = ~(1 | clean_sample_names),
      method      = "VA",
      link        = "probit",
      control     = list(reltol = 1e-6),
      sd.errors   = FALSE,      # skip SEs to save time during selection
      seed        = SEED
    ),
    error = function(e) {
      cat(" ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (!is.null(fit_lv)) {
    lv_results$AICc[i]      <- tryCatch(summary(fit_lv)$AICc, error = function(e) NA)
    lv_results$AIC[i]       <- tryCatch(AIC(fit_lv),           error = function(e) NA)
    lv_results$BIC[i]       <- tryCatch(BIC(fit_lv),           error = function(e) NA)
    lv_results$logLik[i]    <- tryCatch(as.numeric(logLik(fit_lv)), error = function(e) NA)
    lv_results$converged[i] <- TRUE
    cat(" AICc =", round(lv_results$AICc[i], 1), "\n")
    
    # Save each fitted LV model (useful for diagnostics later)
    saveRDS(fit_lv,
            file.path(SAVE_DIR_SELECTION,
                      paste0("lv_select_water_lv", lv_i, ".rds")))
    rm(fit_lv); gc()
  }
}



# --------


# Delta columns (distance from best)
lv_results$delta_AICc <- lv_results$AICc - min(lv_results$AICc, na.rm = TRUE)
lv_results$delta_BIC  <- lv_results$BIC  - min(lv_results$BIC,  na.rm = TRUE)

# Best choices
best_lv_AICc <- lv_results$num_lv[which.min(lv_results$AICc)]
best_lv_BIC  <- lv_results$num_lv[which.min(lv_results$BIC)]

cat("\n--- LV SELECTION RESULTS ---\n")
print(lv_results[, c("num_lv", "AICc", "delta_AICc", "BIC", "delta_BIC",
                     "logLik", "converged")])
cat("\n  Best by AICc: num_lv =", best_lv_AICc, "\n")
cat("  Best by BIC:  num_lv =", best_lv_BIC,  "\n")
cat("\n  NOTE: NUM_LV_FINAL is currently set to", NUM_LV_FINAL, ".\n")
cat("  If the selection above suggests a different value, update\n")
cat("  NUM_LV_FINAL at the top of this script before running Section 7.\n\n")

write.csv(lv_results, "lv_selection_water_PA.csv")
cat("  Saved: lv_selection_water_PA.csv\n\n")

# ==============================================================================
# SECTION 6: MODEL SCREENING (M1–M7) AT NUM_LV_FINAL
#
# Because the water model has no substrate variable, several biofilm formulas
# collapse to the same expression:
#
#   M1, M2, M4  →  ~ (site + season)^2   (labelled "interaction")
#   M3, M5, M6  →  ~ site + season        (labelled "additive")
#   M7          →  null (no fixed effects)
#
# Only the three unique formulas are fitted (one model run each). Results are
# then expanded back into the full M1–M7 table by copying the fit statistics
# to every label that shares the same formula. This keeps the output
# comparable to biofilm results while avoiding redundant HPC computation.
#
# sd.errors = FALSE is used here for speed. SEs are computed only for the
# final model in Section 7.
# ==============================================================================

cat("=================================================================\n")
cat(" SECTION 6: MODEL SCREENING (M1–M7)\n")
cat(" num_lv =", NUM_LV_FINAL, "  Family: binomial (probit)\n")
cat("=================================================================\n\n")

# ── Full label → formula mapping (for the output table) ──────────────────────
#   NULL signals "no fixed effects" (M7).
water_formulas <- list(
  M1 = ~ (site + season)^2,   # site * season * substrate  →  site * season
  M2 = ~ (site + season)^2,   # (site+season+substrate)^2  →  (site+season)^2
  M3 = ~ site + season,        # site + season * substrate  →  site + season
  M4 = ~ (site + season)^2,   # site * season + substrate  →  site * season
  M5 = ~ site + season,        # site * substrate + season  →  site + season
  M6 = ~ site + season,        # site + season + substrate  →  site + season
  M7 = NULL                    # null (no fixed effects)    →  unchanged
)

# ── Unique formulas to actually fit ──────────────────────────────────────────
#   Key  = short label used for the RDS filename and progress messages.
#   Value = formula (NULL for the null model).
unique_formulas <- list(
  interaction = ~ (site + season)^2,
  additive    = ~ site + season,
  null        = NULL
)

# Maps each M-label to the unique_formulas key whose formula it matches.
formula_group <- c(
  M1 = "interaction",
  M2 = "interaction",
  M3 = "additive",
  M4 = "interaction",
  M5 = "additive",
  M6 = "additive",
  M7 = "null"
)

# ── Fit the three unique models ───────────────────────────────────────────────

unique_fits <- list()   # stores fit statistics keyed by unique_formulas name

for (grp in names(unique_formulas)) {
  
  formula_i <- unique_formulas[[grp]]
  include_X <- !is.null(formula_i)
  
  cat("  Fitting [", grp, "]",
      if (include_X) paste(":", deparse(formula_i)) else ": null", "...\n")
  
  set.seed(SEED)
  
  args <- list(
    y           = pa_mat,
    family      = "binomial",
    num.lv      = NUM_LV_FINAL,
    studyDesign = study_design,
    row.eff     = ~(1 | clean_sample_names),
    method      = "VA",
    link        = "probit",
    control     = list(reltol = 1e-6),
    sd.errors   = FALSE,
    seed        = SEED
  )
  if (include_X) {
    args$X       <- X
    args$formula <- formula_i
  }
  
  fit_m <- tryCatch(
    do.call(gllvm, args),
    error = function(e) {
      cat("    ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (!is.null(fit_m)) {
    unique_fits[[grp]] <- list(
      AICc      = tryCatch(summary(fit_m)$AICc,          error = function(e) NA),
      AIC       = tryCatch(AIC(fit_m),                    error = function(e) NA),
      BIC       = tryCatch(BIC(fit_m),                    error = function(e) NA),
      logLik    = tryCatch(as.numeric(logLik(fit_m)),     error = function(e) NA),
      converged = TRUE
    )
    saveRDS(fit_m,
            file.path(SAVE_DIR_SELECTION,
                      paste0("screen_", grp, "_water_PA_lv", NUM_LV_FINAL, ".rds")))
    cat("    AICc =", round(unique_fits[[grp]]$AICc, 1),
        "| BIC =",  round(unique_fits[[grp]]$BIC, 1), "\n")
    rm(fit_m); gc()
  } else {
    unique_fits[[grp]] <- list(
      AICc = NA, AIC = NA, BIC = NA, logLik = NA, converged = FALSE
    )
  }
}

# ── Expand results into the full M1–M7 table ─────────────────────────────────

screening_results <- data.frame(
  model     = names(water_formulas),
  group     = formula_group[names(water_formulas)],
  formula   = sapply(water_formulas, function(f)
    if (is.null(f)) "null" else deparse(f)),
  AICc      = sapply(formula_group[names(water_formulas)],
                     function(g) unique_fits[[g]]$AICc),
  AIC       = sapply(formula_group[names(water_formulas)],
                     function(g) unique_fits[[g]]$AIC),
  BIC       = sapply(formula_group[names(water_formulas)],
                     function(g) unique_fits[[g]]$BIC),
  logLik    = sapply(formula_group[names(water_formulas)],
                     function(g) unique_fits[[g]]$logLik),
  converged = sapply(formula_group[names(water_formulas)],
                     function(g) unique_fits[[g]]$converged),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Delta columns and rank (computed on unique values to avoid ties distorting rank)
screening_results$delta_AICc <- screening_results$AICc -
  min(screening_results$AICc, na.rm = TRUE)
screening_results$delta_BIC  <- screening_results$BIC  -
  min(screening_results$BIC,  na.rm = TRUE)

# Rank applied to the group-level AICc so duplicate rows share the same rank
group_rank <- rank(
  sapply(unique_formulas, function(f)
    unique_fits[[names(unique_formulas)[
      sapply(unique_formulas, function(uf)
        identical(deparse(uf), deparse(f)))
    ]]]$AICc),
  na.last = "keep"
)
screening_results$rank_AICc <- group_rank[formula_group[names(water_formulas)]]

cat("\n--- MODEL SCREENING RESULTS (sorted by AICc) ---\n")
cat("  (M1=M2=M4 and M3=M5=M6 share a formula; only 3 models were fitted)\n\n")
print(screening_results[order(screening_results$AICc, na.last = TRUE),
                        c("model", "group", "formula", "AICc", "delta_AICc",
                          "BIC", "delta_BIC", "logLik", "rank_AICc",
                          "converged")])

# Best: take the winning group and report the lowest M-label in that group
best_group      <- names(which.min(sapply(unique_fits, `[[`, "AICc")))
best_model_name <- names(formula_group)[formula_group == best_group][1]
best_formula    <- unique_formulas[[best_group]]

cat("\n  Best group by AICc:", best_group,
    if (!is.null(best_formula)) paste("—", deparse(best_formula)) else "— null")
cat("\n  Representative label:", best_model_name, "\n\n")

write.csv(screening_results,
          file.path(SAVE_DIR_SELECTION,
                    paste0("model_screening_water_PA_lv", NUM_LV_FINAL, ".csv")),
          row.names = FALSE)
cat("  Saved: model_screening_water_PA_lv", NUM_LV_FINAL, ".csv\n\n", sep = "")