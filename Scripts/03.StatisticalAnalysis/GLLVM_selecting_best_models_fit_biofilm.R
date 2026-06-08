################################################################################
# GLLVM ANALYSIS — SCRIPT 1: DATA PREPARATION & MODEL FITTING (VERSION 2)
#
# PURPOSE:
#   Prepares data and fits all GLLVM models. Designed to run on a server.
#   Script 2 (gllvm_02_analysis.R) handles all post-fitting analysis and
#   should be run locally after transferring the .rds files back.
#
# TRACKS:
#   fouling : full eukaryotic community (includes metazoa)
#   biofilm : microeukaryotes only (metazoa removed)
#
# MODEL FAMILIES (four parallel arms):
#   NB      : Negative binomial on raw counts          (VA,  log link)
#   PA      : Binomial presence/absence on 0/1 matrix  (VA,  probit link)
#   ZINB    : Zero-inflated NB on raw counts            (VA,  probit link)
#   betaH   : Hurdle beta on proportions [0,1)          (EVA, logit link)
#
#
# OUTPUTS (all in SAVE_DIR_BASE):
#   fouling_for_gllvm.rds / biofilm_for_gllvm.rds  — prepared data objects
#   gllvm_NB/     — NB models
#   gllvm_PA/     — PA models
#   gllvm_ZINB/   — ZINB models
#   gllvm_betaH/  — betaH models
#   gllvm_lv/     — LV selection CSVs
#
# REFERENCES:
#   Niku et al. (2019). Methods Ecol Evol. doi:10.1111/2041-210X.13303
#   Korhonen et al. (2024). Methods Ecol Evol. doi:10.1111/2041-210X.14437
#   Tang et al. (2025). JDSSV. doi:10.52933/jdssv.v5i6.133
#   Smithson & Verkuilen (2006). Psychological Methods 11:54-71
#   Burnham & Anderson (2002). Model Selection and Multimodel Inference.
################################################################################
library(gllvm)
library(phyloseq)
library(dplyr)

cat("=== GLLVM FITTING SCRIPT ===\n\n")

setwd("~/GitHub/MAPLE_Seasonal_Plastisphere/Scripts/03.StatisticalAnalysis")


# V2 APPENDIX — UPDATED FILTERING + REFITTING
#
# PURPOSE:
#   This section re-runs the full data preparation and model fitting pipeline
#   with two key changes relative to v1:
#
#   CHANGE 1 — Per-factor-level presence filter (Step 2c, new)
#     v1 filtered OTUs by overall prevalence (>= 5% of samples) and total
#     abundance (top 500). This does not prevent retention of taxa whose
#     detections are concentrated in a single factor level — e.g. a taxon
#     present in 12% of samples but 99% of reads in Fall only. Such taxa
#     produce inestimable species-specific coefficients for season (the taxon's
#     baseline and the season effect are fully confounded), causing the
#     optimizer to push estimates toward ±infinity on the log scale.
#     Diagnostic output from v1 confirmed this problem: range(Xcoef) = -57/+68
#     for the NB fouling model, with hundreds of taxa exceeding |10| even in
#     the additive M6 model.
#     Step 2c removes any OTU present in fewer than min_levels_present distinct
#     levels of each factor in filter_factors. Default: >= 2 levels of season
#     and >= 2 levels of site. Substrate is excluded from the default because
#     substrate-level selectivity may reflect genuine biology.
#    
#
#
# DESIGN DECISIONS:
#   - Same phyloseq objects as v1 (uncollapsed PCR replicates) for direct
#     comparison. This isolates the effect of the filtering changes.
#   - Same LV strategy: num_lv = 2 for NB, PA, ZINB; num_lv = 1 for betaH.
#     betaH lv=1 was confirmed as best by AICc and AIC in v1 LV selection.
#   - Full M1-M7 series refit for all four families.
#   - Outputs saved to *_v2 directories; filenames unchanged within those
#     directories so downstream loading code only needs the directory path
#     updated.
#   - Nothing from v1 is overwritten.
#
# TAXA REMOVED BY STEP 2C — NOTE FOR INTERPRETATION:
#   Taxa removed by the per-factor-level filter are genuine season or site
#   specialists. They are not lost — they are candidates for a separate
#   indicator species analysis (De Caceres & Legendre, 2009, Ecology,
#   90(12), 3566-3574) which is better suited to identifying taxa that
#   characterise specific factor level combinations. The GLLVM here answers
#   the community-level question (how does overall composition shift across
#   season, substrate, and site); indicator species analysis answers the
#   taxon-level question (which organisms define specific seasonal or
#   substrate communities).
#

# ==============================================================================
# ==============================================================================

N_OTUS              <- 200      # top OTUs by total abundance
LV_SCREEN           <- 2        # num_lv for the M1-M7 screening pass
LV_BEST             <- 5        # num_lv for refitting the best model
LV_MAX              <- 5        # max num_lv tested in select_num_lv()
MIN_OTUS_PER_SAMPLE <- 3        # remove PCR reps below this after OTU filter
SEED                <- 123



# ── V2 save directories ───────────────────────────────────────────────────────
SAVE_DIR_BASE  <- "~/GitHub/MAPLE_Seasonal_Plastisphere/Processed_data/gllvm_models"

SAVE_DIR_NB_V2    <- file.path(SAVE_DIR_BASE, "gllvm_NB_v2")
SAVE_DIR_PA_V2    <- file.path(SAVE_DIR_BASE, "gllvm_PA_v2")
SAVE_DIR_ZINB_V2  <- file.path(SAVE_DIR_BASE, "gllvm_ZINB_v2")
SAVE_DIR_BETAH_V2 <- file.path(SAVE_DIR_BASE, "gllvm_betaH_v2")

for (d in c(SAVE_DIR_NB_V2, SAVE_DIR_PA_V2, SAVE_DIR_ZINB_V2, SAVE_DIR_BETAH_V2)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── V2 LV settings ────────────────────────────────────────────────────────────
# num_lv = 2 for NB, PA, ZINB  — consistent with v1 screening pass
# num_lv = 1 for betaH         — confirmed best by AICc and AIC in v1 LV
#                                 selection for both fouling and biofilm

LV_SCREEN_NB_PA_ZINB <- 2
LV_SCREEN_BETAH      <- 1


# ==============================================================================
# SECTION V2-1: UPDATED DATA PREPARATION FUNCTION
#
# prepare_gllvm_data_v2() is identical to prepare_gllvm_data() with the
# addition of Step 2c (per-factor-level presence filter). All original
# steps are preserved and numbered identically so diffs are easy to read.
#
# New arguments vs prepare_gllvm_data():
#   filter_factors     — character vector of metadata columns to check.
#                        Default: c("season", "site").
#   min_levels_present — integer. OTUs present in fewer than this many
#                        distinct levels of ANY factor in filter_factors
#                        are removed. Default: 2.
# ==============================================================================

prepare_gllvm_data_v2 <- function(ps,
                                  track,
                                  min_lib_size        = 2000,
                                  prevalence_thresh   = 0.05,
                                  top_n_otus          = 500,
                                  min_bio_samples     = 2,
                                  min_otus_per_sample = MIN_OTUS_PER_SAMPLE,
                                  filter_factors      = c("season", "site"),
                                  min_levels_present  = 2,
                                  fig_dir             = "Figures") {
  
  cat("\n======================================================\n")
  cat("  PREPARING DATA V2 FOR TRACK:", toupper(track), "\n")
  cat("======================================================\n\n")
  
  # ── Step 1: extract and filter by library size ────────────────────────────
  cat("Step 1: Library size filter...\n")
  
  otu_mat  <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  
  metadata <- as.data.frame(sample_data(ps))
  class(metadata) <- "data.frame"
  tax_mat  <- as.data.frame(as(tax_table(ps), "matrix"))
  
  cat("  Original:", nrow(otu_mat), "samples ×", ncol(otu_mat), "OTUs\n")
  
  library_sizes <- rowSums(otu_mat)
  cat("  Library size range:", min(library_sizes), "–", max(library_sizes), "\n")
  cat("  Mean:", round(mean(library_sizes)), "| CV:",
      round(sd(library_sizes) / mean(library_sizes), 3), "\n")
  
  keep <- library_sizes >= min_lib_size
  cat("  Removed (<", min_lib_size, "reads):", sum(!keep),
      "| Retained:", sum(keep), "\n\n")
  
  otu_mat       <- otu_mat[keep, ]
  metadata      <- metadata[keep, ]
  library_sizes <- library_sizes[keep]
  
  # remove empty OTUs/samples
  otu_mat  <- otu_mat[, colSums(otu_mat) > 0]
  tax_mat  <- tax_mat[colnames(otu_mat), , drop = FALSE]
  keep_s   <- rowSums(otu_mat) > 0
  otu_mat  <- otu_mat[keep_s, ]
  metadata <- metadata[keep_s, ]
  library_sizes <- library_sizes[keep_s]
  
  # save library size distribution plot
  png(file.path(fig_dir, paste0("library_size_v2_", track, ".png")),
      width = 10, height = 5, units = "in", res = 200)
  par(mfrow = c(1, 2))
  hist(library_sizes, breaks = 50,
       main = paste("Library sizes v2 —", track),
       xlab = "Reads", col = "steelblue", border = "white")
  abline(v = mean(library_sizes), col = "red",  lwd = 2, lty = 2)
  abline(v = min_lib_size,        col = "blue", lwd = 2, lty = 3)
  boxplot(library_sizes ~ metadata$season,
          main = "by season", ylab = "Reads", las = 2,
          col = rainbow(length(unique(metadata$season))))
  abline(h = min_lib_size, col = "blue", lwd = 2, lty = 3)
  par(mfrow = c(1, 1))
  dev.off()
  
  # ── Step 2: OTU prevalence/abundance filter ───────────────────────────────
  # Unchanged from v1. Removes OTUs below overall prevalence threshold or
  # outside the top_n_otus by abundance with < min_bio_samples detections.
  # Limitation retained intentionally: overall prevalence does not capture
  # whether detections are distributed across factor levels. Step 2c below
  # addresses this gap.
  cat("Step 2: OTU prevalence/abundance filter...\n")
  
  prevalence <- colSums(otu_mat > 0) / nrow(otu_mat)
  abund_otus <- prevalence >= prevalence_thresh
  
  bio_present <- sapply(colnames(otu_mat), function(otu) {
    length(unique(metadata$clean_sample_names[otu_mat[, otu] > 0]))
  })
  
  top_abund <- names(sort(colSums(otu_mat), decreasing = TRUE)[1:top_n_otus])
  keep_otus <- abund_otus | (colnames(otu_mat) %in% top_abund & bio_present >= min_bio_samples)
  
  otu_filt <- otu_mat[, keep_otus]
  cat("  OTUs retained:", sum(keep_otus), "| Removed:", sum(!keep_otus), "\n")
  
  zero_s <- rowSums(otu_filt) == 0
  if (sum(zero_s) > 0) {
    cat("  WARNING:", sum(zero_s), "zero-sum samples removed\n")
    otu_filt      <- otu_filt[!zero_s, ]
    metadata      <- metadata[!zero_s, ]
    library_sizes <- library_sizes[!zero_s]
  }
  
  # ── Step 2b: OTU degeneracy filter ───────────────────────────────────────
  # Unchanged from v1. Removes PCR replicates reduced to < min_otus_per_sample
  # OTUs after the prevalence filter. These carry negligible co-occurrence
  # information for latent variable estimation.
  cat("Step 2b: OTU degeneracy filter (threshold =", min_otus_per_sample, ")...\n")
  
  n_otus_per_pcr <- rowSums(otu_filt > 0)
  degenerate     <- n_otus_per_pcr < min_otus_per_sample
  
  if (sum(degenerate) > 0) {
    cat("  Removed", sum(degenerate), "PCR replicates with <",
        min_otus_per_sample, "OTUs:\n")
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
  
  # ── Step 2c: per-factor-level presence filter (NEW in v2) ─────────────────
  # Removes OTUs that are detectable in fewer than min_levels_present distinct
  # levels of any factor in filter_factors.
  #
  # Why Steps 2 and 2b are insufficient:
  #   Overall prevalence filtering retains taxa whose detections are
  #   concentrated in a single factor level. When a taxon is only ever present
  #   in (e.g.) Fall, its species-specific season coefficient is unidentifiable:
  #   the taxon's baseline abundance and the Fall effect are fully confounded.
  #   The GLLVM optimizer resolves this by pushing the coefficient toward
  #   ±infinity. v1 diagnostic confirmed: range(Xcoef) = -57/+68 for NB
  #   fouling M1, with 226 taxa > |10| even in the additive M6 model.
  #
  # Factor selection:
  #   filter_factors defaults to c("season", "site"). Season is the most
  #   problematic dimension (confirmed by Ephelota_ASV_2990 diagnostic: 99%
  #   of reads in Fall). Site (2 levels) is included because a taxon detected
  #   at only one site cannot support a stable site coefficient.
  #   Substrate is excluded from the default: substrate-level selectivity
  #   (e.g. a taxon only on PE) may be genuine colonisation biology rather
  #   than a data sparsity artefact, and applying the filter to substrate
  #   would remove ecologically informative patterns.
  #   Users can add "substrate" to filter_factors if needed.
  #
  # Implementation:
  #   For each factor, count distinct levels with >= 1 detection per OTU.
  #   Remove OTUs failing the threshold for ANY factor.
  #
  # Estimability rationale: Warton et al. (2015), MEE, 6(12), 1395-1404.
  #   https://doi.org/10.1111/2041-210X.12414
  #   Separation geometry: Heinze & Schemper (2002), Stat Med, 21, 2409-2419.
  #   https://doi.org/10.1002/sim.1047
  #
  # NOTE: taxa removed here are candidates for indicator species analysis
  # (De Caceres & Legendre 2009) — they are season/site specialists that
  # are better characterised by fidelity/specificity metrics than by GLLVM
  # fixed-effect coefficients.
  
  cat("Step 2c: Per-factor-level presence filter (v2 addition)...\n")
  cat("  Factors checked:", paste(filter_factors, collapse = ", "),
      "| Min levels required:", min_levels_present, "\n")
  
  keep_level <- rep(TRUE, ncol(otu_filt))
  names(keep_level) <- colnames(otu_filt)
  
  for (fac in filter_factors) {
    
    if (!fac %in% colnames(metadata)) {
      cat("  WARNING: factor '", fac, "' not found in metadata — skipping\n")
      next
    }
    
    fac_vec <- metadata[[fac]]
    
    # count distinct levels with at least one detection per OTU
    n_levels_present <- apply(otu_filt, 2, function(x) {
      length(unique(fac_vec[x > 0]))
    })
    
    fails_this_factor <- n_levels_present < min_levels_present
    cat("  Factor '", fac, "': removing", sum(fails_this_factor),
        "OTUs present in <", min_levels_present, "levels\n")
    
    keep_level <- keep_level & !fails_this_factor
  }
  
  cat("  OTUs retained after Step 2c:", sum(keep_level),
      "| Total removed:", sum(!keep_level), "\n")
  
  # report removed taxa for downstream indicator species analysis
  removed_2c <- colnames(otu_filt)[!keep_level]
  if (length(removed_2c) > 0) {
    cat("  Removed taxa (candidates for indicator species analysis):\n")
    cat("   ", paste(head(removed_2c, 20), collapse = ", "))
    if (length(removed_2c) > 20)
      cat(" ... and", length(removed_2c) - 20, "more")
    cat("\n")
  }
  
  otu_filt <- otu_filt[, keep_level]
  
  # remove any zero-sum samples created by Step 2c
  zero_s2c <- rowSums(otu_filt) == 0
  if (sum(zero_s2c) > 0) {
    cat("  WARNING:", sum(zero_s2c),
        "zero-sum samples created by Step 2c — removed\n")
    otu_filt      <- otu_filt[!zero_s2c, ]
    metadata      <- metadata[!zero_s2c, ]
    library_sizes <- library_sizes[!zero_s2c]
  }
  
  lib_sizes_filt <- rowSums(otu_filt)
  lib_offset     <- log(lib_sizes_filt)
  
  cat("  Final matrix:", nrow(otu_filt), "samples ×", ncol(otu_filt), "OTUs\n")
  cat("  Library size range (post-filter):",
      min(lib_sizes_filt), "–", max(lib_sizes_filt), "\n\n")
  
  # ── Step 3: factor levels ─────────────────────────────────────────────────
  cat("Step 3: Factor levels...\n")
  metadata$season    <- factor(metadata$season,
                               levels = c("Winter","Spring","Summer","Fall","Winter2"))
  metadata$substrate <- factor(metadata$substrate,
                               levels = c("Glass","PE","Weathered_PE","PET","Weathered_PET"))
  metadata$site      <- factor(metadata$site, levels = c("SELVA","TBS"))
  cat("  Season ref: Winter | Substrate ref: Glass | Site ref: SELVA\n\n")
  
  # ── Step 4: readable OTU labels ──────────────────────────────────────────
  cat("Step 4: OTU labels...\n")
  make_label <- function(id, tm) {
    for (rank in c("Genus","Family","Class")) {
      v <- tm[id, rank]
      if (!is.na(v) && nchar(v) > 0 && !v %in% c("uncultured","Unknown",""))
        return(paste0(v, "_", substr(id, 1, 12)))
    }
    paste0("OTU_", substr(id, 1, 12))
  }
  otu_id_to_label <- setNames(
    sapply(colnames(otu_filt), make_label, tm = tax_mat),
    colnames(otu_filt)
  )
  cat("  Labels created:", length(otu_id_to_label), "\n\n")
  
  # ── Step 5: random effect structure ──────────────────────────────────────
  cat("Step 5: Random effect structure...\n")
  metadata$clean_sample_names <- as.factor(metadata$clean_sample_names)
  study_design <- data.frame(clean_sample_names = metadata$clean_sample_names)
  cat("  Groups:", nlevels(metadata$clean_sample_names), "\n")
  cat("  Reps per group:",
      paste(range(table(metadata$clean_sample_names)), collapse = "–"), "\n\n")
  
  cat("  v2 preparation complete for", track, "\n\n")
  
  list(
    track              = track,
    otu_filt           = otu_filt,
    metadata           = metadata,
    tax_mat            = tax_mat,
    lib_offset         = lib_offset,
    lib_sizes_filt     = lib_sizes_filt,
    lib_sizes_original = library_sizes,
    study_design       = study_design,
    otu_id_to_label    = otu_id_to_label,
    design_factors     = c("site", "season", "substrate"),
    removed_step2c     = removed_2c        # passed to indicator species script
  )
}


# ==============================================================================
# SECTION V2-2: UPDATED RESPONSE MATRIX FUNCTION
# ==============================================================================

prepare_response_matrices_v2 <- function(d, n_otus = N_OTUS) {
  
  cat("\n--- Response matrices v2 for track:", toupper(d$track), "---\n")
  
  # ── OTU selection: prevalence-first ranking ───────────────────────────────
  top_otus <- names(sort(colSums(d$otu_filt), decreasing = TRUE)[1:n_otus])
  
  cat("  OTU selection: total abundance ranking\n")
  cat("  Abundance range of selected OTUs:",
      min(colSums(d$otu_filt)[top_otus]), "–",
      max(colSums(d$otu_filt)[top_otus]), "reads\n")
  
  # Count matrix
  count_mat <- d$otu_filt[, top_otus]
  colnames(count_mat) <- d$otu_id_to_label[top_otus]
  
  # Presence/absence matrix
  pa_mat <- (count_mat > 0) * 1L
  
  # Proportion matrix — denominator = original library size
  prop_mat <- d$otu_filt[, top_otus] / d$lib_sizes_original
  colnames(prop_mat) <- d$otu_id_to_label[top_otus]
  
  # Exact-1 check and nudge for betaH
  n_ones <- sum(prop_mat == 1)
  if (n_ones == 0) {
    cat("  Exact 1s: 0 — betaH applicable without adjustment\n")
  } else {
    cat("  Exact 1s:", n_ones, "— applying nudge (1 - .Machine$double.eps)\n")
    prop_mat[prop_mat == 1] <- 1 - .Machine$double.eps
    cat("  Verification — exact 1s remaining:", sum(prop_mat == 1), "\n")
  }
  
  cat("  Count:", nrow(count_mat), "×", ncol(count_mat),
      "| zeros:", round(mean(count_mat == 0) * 100, 1), "%\n")
  cat("  Prop: max =", format(max(prop_mat), digits = 6),
      "| zeros:", round(mean(prop_mat == 0) * 100, 1), "%\n\n")
  
  list(count = count_mat, pa = pa_mat, prop = prop_mat)
}





# ==============================================================================
# SECTION V2-3: DATA PREPARATION — run from raw phyloseq objects
#
# Uses the same phyloseq objects as v1 (uncollapsed PCR replicates) so that
# any differences in model output are attributable to the filtering changes
# alone, not to a change in the input data structure.
# ==============================================================================

cat("\n=== V2: DATA PREPARATION ===\n")

ps_all_v2 <- readRDS(
  "~/GitHub/MAPLE_Seasonal_Plastisphere/Processed_data/Phyloseq_objects/FINAL_PHYLOSEQ_OBJECTS.rds"
)

data_fouling_v2 <- prepare_gllvm_data_v2(
  ps    = ps_all_v2$pr2_rep50$fP.lulu$ST,
  track = "fouling"
)
saveRDS(data_fouling_v2,
        file.path(SAVE_DIR_BASE, "fouling_for_gllvm_v2.rds"))

data_biofilm_v2 <- prepare_gllvm_data_v2(
  ps    = ps_all_v2$pr2_rep50_noMeta$fP.lulu$ST,
  track = "biofilm"
)
saveRDS(data_biofilm_v2,
        file.path(SAVE_DIR_BASE, "biofilm_for_gllvm_v2.rds"))

rm(ps_all_v2); gc()

# ── V2 response matrices ──────────────────────────────────────────────────────

mats_fouling_v2 <- prepare_response_matrices_v2(data_fouling_v2)
mats_biofilm_v2 <- prepare_response_matrices_v2(data_biofilm_v2)






################################################################################
# RUN ON ASTBURY
################################################################################



#save.image("gllvm_data_to_fit_v2.RData")

library(gllvm)
library(phyloseq)
library(dplyr)
load("gllvm_data_to_fit_v2.RData")

SAVE_DIR_BASE <- "/data/glennsdata/MAPLE/18S/STATISTICS/GLLVM/models_v2"

SAVE_DIR_NB_V2    <- file.path(SAVE_DIR_BASE, "gllvm_NB_v2")
SAVE_DIR_PA_V2    <- file.path(SAVE_DIR_BASE, "gllvm_PA_v2")
SAVE_DIR_ZINB_V2  <- file.path(SAVE_DIR_BASE, "gllvm_ZINB_v2")
SAVE_DIR_BETAH_V2 <- file.path(SAVE_DIR_BASE, "gllvm_betaH_v2")
SAVE_DIR_LV    <- file.path("GLLVM/gllvm_lv")

for (d in c(SAVE_DIR_NB_v2, SAVE_DIR_PA_v2, SAVE_DIR_ZINB_v2, SAVE_DIR_BETAH_v2, SAVE_DIR_LV)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

dir.create("Figures", showWarnings = FALSE)

# ==============================================================================
# SECTION 3: LV SELECTION (optional — run if num_lv is uncertain)
#
# Fits M6 additive with num_lv = 1:LV_MAX for each family.
# Read AICc and BIC output before setting LV_SCREEN and LV_BEST above.
#
# NOTE: LV selection for NB showed lv=5 as best by AICc, which is unusual.
# The screening pass (lv=2) is still run for all M1-M7 to identify the best
# covariate structure first; then the best model is refit with lv=5 to assess
# whether the extra LVs affect inference.
# ==============================================================================

select_num_lv <- function(d, Y, family, method, link = "probit",
                          use_offset = TRUE, formula = ~ site + season + substrate,
                          max_lv = LV_MAX, seed = SEED, arm_label, out_dir) {
  
  cat("\n--- LV selection:", arm_label, "|", toupper(d$track), "---\n")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ctrl <- if (family == "betaH") list(reltol = 1e-12) else list(reltol = 1e-6)
  
  results <- data.frame(num_lv = 1:max_lv, AICc = NA, AIC = NA, BIC = NA)
  
  for (i in 1:max_lv) {
    cat("  lv =", i, "...")
    set.seed(seed)
    
    args <- list(
      y           = Y,
      X           = d$metadata[, d$design_factors],
      family      = family,
      num.lv      = i,
      formula     = formula,
      studyDesign = d$study_design,
      row.eff     = ~(1 | clean_sample_names),
      method      = method,
      link        = link,
      control     = ctrl,
      sd.errors   = FALSE,
      seed        = seed
    )
    if (use_offset) args$offset <- d$lib_offset
    
    fit <- tryCatch(
      do.call(gllvm, args),
      error = function(e) { cat(" ERROR:", conditionMessage(e), "\n"); NULL }
    )
    
    if (!is.null(fit)) {
      results$AICc[i] <- tryCatch(summary(fit)$AICc, error = function(e) NA)
      results$AIC[i]  <- tryCatch(AIC(fit),           error = function(e) NA)
      results$BIC[i]  <- tryCatch(BIC(fit),           error = function(e) NA)
      cat(" AICc =", round(results$AICc[i], 1), "\n")
      rm(fit); gc()
    }
  }
  
  results$delta_AICc <- results$AICc - min(results$AICc, na.rm = TRUE)
  results$delta_BIC  <- results$BIC  - min(results$BIC,  na.rm = TRUE)
  
  cat("\n"); print(results)
  cat("Best by AICc: lv =", results$num_lv[which.min(results$AICc)], "\n")
  cat("Best by BIC:  lv =", results$num_lv[which.min(results$BIC)],  "\n\n")
  
  fname <- file.path(out_dir, paste0("lv_selection_", arm_label, "_", d$track, ".csv"))
  write.csv(results, fname, row.names = FALSE)
  cat("Saved →", fname, "\n")
  
  invisible(results)
}

# Uncomment to run LV selection (server only — slow):
lv_nb_f    <- select_num_lv(data_fouling_v2, mats_fouling_v2$count, "negative.binomial", "VA",  arm_label="NB",    out_dir=SAVE_DIR_LV)
lv_nb_b    <- select_num_lv(data_biofilm_v2, mats_biofilm_v2$count, "negative.binomial", "VA",  arm_label="NB",    out_dir=SAVE_DIR_LV)
lv_pa_f    <- select_num_lv(data_fouling_v2, mats_fouling_v2$pa,    "binomial",          "VA",  use_offset=FALSE,  arm_label="PA",    out_dir=SAVE_DIR_LV)
lv_pa_b    <- select_num_lv(data_biofilm_v2, mats_biofilm_v2$pa,    "binomial",          "VA",  use_offset=FALSE,  arm_label="PA",    out_dir=SAVE_DIR_LV)
lv_zinb_f  <- select_num_lv(data_fouling_v2, mats_fouling_v2$count, "ZINB",              "VA",  arm_label="ZINB",  out_dir=SAVE_DIR_LV)
lv_zinb_b  <- select_num_lv(data_biofilm_v2, mats_biofilm_v2$count, "ZINB",              "VA",  arm_label="ZINB",  out_dir=SAVE_DIR_LV)
lv_bh_f    <- select_num_lv(data_fouling_v2, mats_fouling_v2$prop,  "betaH",             "EVA", use_offset=FALSE, link="logit", arm_label="betaH", out_dir=SAVE_DIR_LV)
lv_bh_b    <- select_num_lv(data_biofilm_v2, mats_biofilm_v2$prop,  "betaH",             "EVA", use_offset=FALSE, link="logit", arm_label="betaH", out_dir=SAVE_DIR_LV)




# ==============================================================================
# SECTION V2-4: MODEL FITTING — PHASE 1 SCREENING PASS
#
# fit_gllvm_arm() fits M1–M7 for any family/method combination.
# Each model is saved as an .rds immediately after fitting to protect
# against server crashes mid-run.
#
# MODEL SERIES (M1–M7):
#   M1: site × season × substrate   (full 3-way)
#   M2: (site + season + substrate)² (all 2-way interactions)
#   M3: site + season × substrate
#   M4: site × season + substrate
#   M5: site × substrate + season
#   M6: site + season + substrate    (additive)
#   M7: no fixed effects             (unconstrained — null for LV variance)
#
# FAMILY-SPECIFIC NOTES:
#   NB    — family="negative.binomial", method="VA", with offset
#   PA    — family="binomial",          method="VA", no offset (binary)
#   ZINB  — family="ZINB",              method="VA", with offset
#           Tang et al. (2025): ZINB outperforms NB at lower sparsity
#   betaH — family="betaH",             method="EVA", no offset (proportions)
#           reltol=1e-12 required for EVA convergence (Korhonen vignette 2025)
# ==============================================================================

fit_gllvm_arm <- function(d, Y, family, method, link = "probit",
                          use_offset = TRUE, num_lv = LV_SCREEN,
                          seed = SEED, save_dir, arm_label) {
  
  cat("\n======================================================\n")
  cat(" ", arm_label, "— FITTING:", toupper(d$track),
      "| num_lv =", num_lv, "\n")
  cat("======================================================\n\n")
  
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  X            <- d$metadata[, d$design_factors]
  study_design <- d$study_design
  ctrl         <- if (family == "betaH") list(reltol = 1e-12) else list(reltol = 1e-6)
  
  run_one <- function(formula, label, include_X = TRUE) {
    cat("  Fitting", label, "...\n")
    set.seed(seed)
    
    args <- list(
      y           = Y,
      family      = family,
      num.lv      = num_lv,
      studyDesign = study_design,
      row.eff     = ~(1 | clean_sample_names),
      method      = method,
      link        = link,
      control     = ctrl,
      seed        = seed
    )
    if (use_offset)  args$offset  <- d$lib_offset
    if (include_X) { args$X <- X; args$formula <- formula }
    
    fit <- tryCatch(
      do.call(gllvm, args),
      error = function(e) { cat("    ERROR:", conditionMessage(e), "\n"); NULL }
    )
    
    if (!is.null(fit)) {
      fname <- file.path(save_dir,
                         paste0(label, "_", arm_label, "_", d$track,
                                "_lv", num_lv, ".rds"))
      saveRDS(fit, fname)
      cat("    Saved →", fname, "\n")
    }
    invisible(fit)
  }
  
  models <- list(
    M1 = run_one(~ site * season * substrate,     "M1_3way"),
    M2 = run_one(~ (site + season + substrate)^2, "M2_all2way"),
    M3 = run_one(~ site + season * substrate,      "M3_seasonXsubstr"),
    M4 = run_one(~ site * season + substrate,      "M4_siteXseason"),
    M5 = run_one(~ site * substrate + season,      "M5_siteXsubstr"),
    M6 = run_one(~ site + season + substrate,      "M6_additive"),
    M7 = run_one(NULL, "M7_unconstrained", include_X = FALSE)
  )
  
  # Free memory before returning
  rm(X, study_design, ctrl); gc()
  
  invisible(models)
}


cat("\n=== V2: PHASE 1 SCREENING PASS ===\n")
cat("  NB, PA, ZINB: num_lv =", LV_SCREEN_NB_PA_ZINB, "\n")
cat("  betaH:        num_lv =", LV_SCREEN_BETAH, "\n\n")

# ── NB ────────────────────────────────────────────────────────────────────────
fit_gllvm_arm(data_fouling_v2, mats_fouling_v2$count, "negative.binomial", "VA",
              use_offset = TRUE,  num_lv = LV_SCREEN_NB_PA_ZINB,
              save_dir = SAVE_DIR_NB_V2, arm_label = "NB")

fit_gllvm_arm(data_biofilm_v2, mats_biofilm_v2$count, "negative.binomial", "VA",
              use_offset = TRUE,  num_lv = LV_SCREEN_NB_PA_ZINB,
              save_dir = SAVE_DIR_NB_V2, arm_label = "NB")

# ── PA ────────────────────────────────────────────────────────────────────────
fit_gllvm_arm(data_fouling_v2, mats_fouling_v2$pa, "binomial", "VA",
              use_offset = FALSE, num_lv = LV_SCREEN_NB_PA_ZINB,
              save_dir = SAVE_DIR_PA_V2, arm_label = "PA")

fit_gllvm_arm(data_biofilm_v2, mats_biofilm_v2$pa, "binomial", "VA",
              use_offset = FALSE, num_lv = LV_SCREEN_NB_PA_ZINB,
              save_dir = SAVE_DIR_PA_V2, arm_label = "PA")

# ── ZINB ──────────────────────────────────────────────────────────────────────
fit_gllvm_arm(data_fouling_v2, mats_fouling_v2$count, "ZINB", "VA",
              use_offset = TRUE,  num_lv = LV_SCREEN_NB_PA_ZINB,
              save_dir = SAVE_DIR_ZINB_V2, arm_label = "ZINB")

fit_gllvm_arm(data_biofilm_v2, mats_biofilm_v2$count, "ZINB", "VA",
              use_offset = TRUE,  num_lv = LV_SCREEN_NB_PA_ZINB,
              save_dir = SAVE_DIR_ZINB_V2, arm_label = "ZINB")

# ── betaH ─────────────────────────────────────────────────────────────────────
fit_gllvm_arm(data_fouling_v2, mats_fouling_v2$prop, "betaH", "EVA",
              use_offset = FALSE, num_lv = LV_SCREEN_BETAH, link = "logit",
              save_dir = SAVE_DIR_BETAH_V2, arm_label = "betaH")

fit_gllvm_arm(data_biofilm_v2, mats_biofilm_v2$prop, "betaH", "EVA",
              use_offset = FALSE, num_lv = LV_SCREEN_BETAH, link = "logit",
              save_dir = SAVE_DIR_BETAH_V2, arm_label = "betaH")













