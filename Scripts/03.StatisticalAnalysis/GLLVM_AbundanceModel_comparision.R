################################################################################
# GLLVM SENSITIVITY CHECK — FOULING TRACK, ABUNDANCE MODEL (NEGATIVE BINOMIAL)
#
# PURPOSE:
#   Sensitivity check against the original binomial presence/absence (PA)
#   model. Same data filtering pipeline, but:
#     - Response: raw OTU counts (abundance), NOT presence/absence
#     - Family:   negative.binomial (instead of binomial)
#     - Link:     default (log) for negative binomial — "probit" removed,
#                 since it only applies to binary/binomial models
#     - Offset:   log(library size) added, since abundance models normally
#                 need to account for differing sequencing depth per sample.
#                 ASSUMPTION FLAGGED: your original PA model had no offset
#                 because presence/absence doesn't need one. If you want a
#                 stricter apples-to-apples comparison, remove the offset
#                 argument in Section 4 (search "OFFSET ASSUMPTION").
#
#   This script is a direct merge of:
#     1. GLLVM_final_fit.R      (data prep + model fitting)      -> Sections 1-6
#     2. GLLVM_analysis.R       (residual co-occurrence analysis) -> Sections 7+
#
#   All filtering parameters (Steps 1-2c) are identical to the original v2
#   pipeline, since the question here is only "does family/response type
#   change the conclusions", not "does filtering change them".
#
################################################################################

library(gllvm)
library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(knitr)
library(scales)
library(vegan)
library(RColorBrewer)
library(ggsci)
library(kableExtra)
library(ggrepel)
library(igraph)
library(ggraph)
library(tidygraph)
library(cluster)
library(circlize)
library(flextable)
library(officer)

cat("=== GLLVM SENSITIVITY CHECK — FOULING / NEGATIVE BINOMIAL / ABUNDANCE / M2 / LV2 ===\n\n")

setwd("~/GitHub/MAPLE_Seasonal_Plastisphere/Scripts/03.StatisticalAnalysis")

SAVE_DIR_BASE <- "~/Github/MAPLE_Seasonal_Plastisphere/Processed_data/gllvm_models"
OUT_DIR       <- "~/Github/MAPLE_Seasonal_Plastisphere/Results/gllvm_results/NB_sensitivity"
TABLE_DIR     <- file.path(OUT_DIR, "Tables")
FIG_DIR       <- file.path(OUT_DIR, "Figures")
for (d in c(OUT_DIR, FIG_DIR, TABLE_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# HARDCODED PARAMETERS (identical to original v2 pipeline)
# ==============================================================================

# --- Filtering ---
MIN_LIB_SIZE        <- 2000
PREVALENCE_THRESH   <- 0.05
TOP_N_OTUS          <- 500
MIN_BIO_SAMPLES     <- 2
MIN_OTUS_PER_SAMPLE <- 3
FILTER_FACTORS      <- c("season", "site", "substrate")
MIN_LEVELS_PRESENT  <- 2

# --- Response matrix ---
N_OTUS <- 200   # top OTUs by total abundance used in model (same as PA model)

# --- Model ---
SEED   <- 123
NUM_LV <- 2

# --- Downstream analysis thresholds (same as original analysis script) ---
R_STRONG   <- 0.5
N_CLUSTERS <- 5


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
# SECTION 2: DATA PREPARATION (identical to original v2 pipeline, Steps 1-5)
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
# SECTION 3: RESPONSE MATRIX — ABUNDANCE (COUNTS), NOT PRESENCE/ABSENCE
#
# This is the key change from the original PA pipeline: we keep the raw
# counts for the top N_OTUS OTUs instead of converting to 0/1.
# ==============================================================================

cat("--- Building abundance (count) response matrix (top", N_OTUS, "OTUs) ---\n")

top_otus  <- names(sort(colSums(otu_filt), decreasing = TRUE)[1:N_OTUS])
count_mat <- otu_filt[, top_otus]
colnames(count_mat) <- otu_id_to_label[top_otus]

# Library sizes recomputed on the final filtered matrix, used as the model offset
library_sizes_final <- rowSums(otu_filt)

cat("  Count matrix:", nrow(count_mat), "x", ncol(count_mat),
    "| zeros:", round(mean(count_mat == 0) * 100, 1), "%\n\n")


save.image("gllvm_abundance_sensitivity_dataprep.RData")


# ==============================================================================
# SECTION 4: MODEL FITTING
#   Family:  negative.binomial (abundance)   [was: binomial, PA]
#   Formula: M2 — (site + season + substrate)^2 (all 2-way interactions)
#   LV:      num_lv = 2
#   Method:  VA
#   Link:    default (log) — "probit" is not valid for negative binomial
#            and has been removed.
#   Offset:  log(library size)  <-- OFFSET ASSUMPTION, see header note.
#            Remove the "offset =" line below if you want a model with no
#            offset, matching the original PA script exactly.
# ==============================================================================

load("gllvm_abundance_sensitivity_dataprep.RData")


library(gllvm)
library(phyloseq)
library(dplyr)

cat("=== GLLVM SENSITIVITY ANALYSIS — TOP-N OTU THRESHOLD - on ASTBURY HPC ===\n\n")

SEED <- 123
SAVE_DIR_BASE <- "/home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/STATISTICS"
OUT_DIR       <- "/home/ahomew/g/gabrihae/data/glennsdata/MAPLE/18S/STATISTICS/sensitivity_analysis"

# Create a subfolder just for sensitivity outputs, so they don't mix with
# your final model files
SENS_DIR <- file.path(SAVE_DIR_BASE, "sensitivity_abundance")
dir.create(SENS_DIR, showWarnings = FALSE, recursive = TRUE)



cat("--- Fitting final model: Abundance / NB / M2 (all 2-way) / LV2 ---\n\n")

X <- metadata[, c("site", "season", "substrate")]

set.seed(SEED)

fit_final <- gllvm(
  y           = count_mat,
  X           = X,
  formula     = ~ (site + season + substrate)^2,
  family      = "negative.binomial",
  num.lv      = NUM_LV,
  studyDesign = study_design,
  row.eff     = ~(1 | clean_sample_names),
  offset      = log(library_sizes_final),   # OFFSET ASSUMPTION — see header
  method      = "VA",
  control     = list(reltol = 1e-6),
  seed        = SEED
)

saveRDS(fit_final, file.path(SAVE_DIR_BASE, "GLLVM_final_incM2_NB.rds"))

fit_final <- readRDS(file.path(SAVE_DIR_BASE, "GLLVM_final_incM2_NB.rds"))

cat("\n--- Best model fitted successfully ---\n")
cat("  AIC:    ", round(AIC(fit_final), 1), "\n")
cat("  BIC:    ", round(BIC(fit_final), 1), "\n")
cat("  AICc:   ", round(summary(fit_final)$AICc, 1), "\n")
cat("  logLik: ", round(as.numeric(logLik(fit_final)), 1), "\n")
cat("  Xcoef range:", paste(round(range(coef(fit_final)$Xcoef), 2), collapse = " to "), "\n\n")


# ==============================================================================
# SECTION 5: NULL / UNCONSTRAINED MODEL
#   No fixed effects — latent variables only. Same offset logic as above.
# ==============================================================================

cat("--- Fitting null model: Abundance / NB / M7 (no fixed effects) / LV2 ---\n\n")

set.seed(SEED)

fit_null <- gllvm(
  y           = count_mat,
  family      = "negative.binomial",
  num.lv      = NUM_LV,
  studyDesign = study_design,
  row.eff     = ~(1 | clean_sample_names),
  offset      = log(library_sizes_final),   # OFFSET ASSUMPTION — see header
  method      = "VA",
  control     = list(reltol = 1e-6),
  seed        = SEED
)

saveRDS(fit_null, file.path(SAVE_DIR_BASE, "GLLVM_null_incM7_NB.rds"))

fit_null <- readRDS(file.path(SAVE_DIR_BASE, "GLLVM_null_incM7_NB.rds"))

cat("\n--- Null model fitted successfully ---\n")
cat("  AIC:    ", round(AIC(fit_null), 1), "\n")
cat("  BIC:    ", round(BIC(fit_null), 1), "\n")
cat("  logLik: ", round(as.numeric(logLik(fit_null)), 1), "\n\n")

lv_null <- sum(apply(getLV(fit_null), 2, var))
lv_best <- sum(apply(getLV(fit_final), 2, var))
cat("--- Sanity check ---\n")
cat("  LV variance — null:", round(lv_null, 3),
    "| best:", round(lv_best, 3),
    "| reduction:", round(100 * (lv_null - lv_best) / lv_null, 1), "%\n\n")


# ==============================================================================
# SECTION 6: VARIANCE PARTITIONING (Q1: does site/season/substrate structure
# the community, under the abundance/NB model?)
# ==============================================================================

cat("=== Q1: Community structuring — variance partitioning (NB model) ===\n\n")

vp_gllvm <- VP(fit_final)
vp_gllvm

vp_values <- colMeans(vp_gllvm$PropExplainedVarSp) * 100

vp_display <- data.frame(
  Component = c(
    "Site", "Season", "Substrate",
    "Site x Season", "Site x Substrate", "Season x Substrate",
    "Latent variable 1", "Latent variable 2", "Sample random effect"
  ),
  Type = c(
    "Main effect", "Main effect", "Main effect",
    "Interaction", "Interaction", "Interaction",
    "Residual", "Residual", "Residual"
  ),
  Variance_explained = round(as.numeric(vp_values), 1),
  stringsAsFactors = FALSE
)
colnames(vp_display) <- c("Component", "Type", "Variance explained (%)")

ft_vp <- flextable(vp_display) %>% autofit()
save_as_docx(ft_vp, path = file.path(TABLE_DIR, "vp_gllvm_NB.docx"))

# ── Save data object for downstream analysis (mirrors data_for_gllvm_FINAL.rds)

model_data_final <- list(
  otu_filt          = otu_filt,
  tax_mat           = tax_mat,
  metadata          = metadata,
  otu_id_to_label   = otu_id_to_label,
  study_design      = study_design,
  design_factors    = c("site", "season", "substrate"),
  count_mat         = count_mat,
  library_sizes     = library_sizes_final,
  top_otus          = top_otus,
  removed_step2c    = removed_step2c
)

saveRDS(model_data_final, file.path(SAVE_DIR_BASE, "data_for_gllvm_NB_FINAL.rds"))
cat("NB model fitted and data object saved.\n\n")


################################################################################
################################################################################
#
# FROM HERE ON: RESIDUAL CO-OCCURRENCE ANALYSIS
# (merged and adapted from GLLVM_analysis.R, repointed at the NB abundance
#  model and count-based data instead of the binomial PA model)
#
################################################################################
################################################################################

cat("=== RESIDUAL CO-OCCURRENCE ANALYSIS (NB abundance model) ===\n\n")

model_data <- model_data_final
metadata   <- model_data$metadata

best_model <- fit_final
null_model <- fit_null

cat("Models loaded successfully.\n\n")

# ==============================================================================
# CHECK FOR RESIDUAL STRUCTURE (LV ordination)
# ==============================================================================

lv_scores <- getLV(best_model)
colnames(lv_scores) <- paste0("LV", seq_len(ncol(lv_scores)))

stopifnot(nrow(lv_scores) == nrow(metadata))
lv_df <- data.frame(lv_scores, metadata)

p_lv_season <- ggplot(lv_df, aes(x = LV1, y = LV2, colour = season, shape = site)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(group = season), linewidth = 0.4, alpha = 0.5) +
  coord_equal() +
  labs(x = "Latent variable 1", y = "Latent variable 2",
       title = "Residual (LV) ordination by season and site — NB abundance model") +
  theme_bw(base_size = 11)

print(p_lv_season)
ggsave(file.path(FIG_DIR, "S3A_LV_ordination_season_site_NB.png"),
       p_lv_season, width = 7, height = 6, dpi = 200)


# ==============================================================================
# 3-A: POLLOCK DECOMPOSITION
# ==============================================================================

cat("--- 3-A: Pollock decomposition ---\n\n")

B     <- coef(best_model)$Xcoef
X_raw <- best_model$X

X_numeric <- model.matrix(
  ~ site + season + substrate +
    site:season + site:substrate + season:substrate,
  data = X_raw
)[, -1]

if (!identical(colnames(X_numeric), colnames(B))) {
  cat("Mismatch — only in X_numeric:", setdiff(colnames(X_numeric), colnames(B)), "\n")
  cat("Only in B:",                    setdiff(colnames(B), colnames(X_numeric)), "\n")
  stop("Design matrix columns do not match coefficient matrix. Check factor levels.")
}

fitted_env <- X_numeric %*% t(B)

env_cor <- cor(fitted_env)
diag(env_cor) <- NA

res_cor <- getResidualCor(best_model)
diag(res_cor) <- NA

rcov_best <- getResidualCov(best_model, adjust = 0)
rcov_null <- getResidualCov(null_model, adjust = 0)
trace_reduction_pct <- round((1 - rcov_best$trace / rcov_null$trace) * 100, 1)

cat("Residual covariance trace — null model:", round(rcov_null$trace, 3), "\n")
cat("Residual covariance trace — best model:", round(rcov_best$trace, 3), "\n")
cat("Reduction by site + season + substrate: ", trace_reduction_pct, "%\n\n")

otu_labels  <- colnames(res_cor)
label_to_id <- setNames(names(model_data$otu_id_to_label),
                        model_data$otu_id_to_label)
otu_ids <- label_to_id[otu_labels]

tax_cols_use <- intersect(
  c("Division", "Class", "Order", "Family", "Genus", "Species"),
  colnames(model_data$tax_mat)
)

tax_lookup <- model_data$tax_mat[otu_ids, tax_cols_use, drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID") %>%
  mutate(OTU_Label = model_data$otu_id_to_label[OTU_ID])

pairs <- which(upper.tri(res_cor), arr.ind = TRUE)

pair_df <- data.frame(
  OTU_j   = otu_labels[pairs[, 1]],
  OTU_k   = otu_labels[pairs[, 2]],
  res_cor  = res_cor[pairs],
  env_cor  = env_cor[pairs],
  stringsAsFactors = FALSE
) %>%
  left_join(tax_lookup %>% select(OTU_Label, Class_j = Class, Division_j = Division),
            by = c("OTU_j" = "OTU_Label")) %>%
  left_join(tax_lookup %>% select(OTU_Label, Class_k = Class, Division_k = Division),
            by = c("OTU_k" = "OTU_Label")) %>%
  mutate(
    same_class    = !is.na(Class_j) & !is.na(Class_k) & Class_j == Class_k,
    same_division = !is.na(Division_j) & !is.na(Division_k) & Division_j == Division_k,
    quadrant = case_when(
      env_cor >= 0 & res_cor >= 0 ~ "Q1",
      env_cor >= 0 & res_cor <  0 ~ "Q2",
      env_cor <  0 & res_cor >= 0 ~ "Q3",
      env_cor <  0 & res_cor <  0 ~ "Q4"
    ),
    surprise   = abs(res_cor) - abs(env_cor),
    strong_res = abs(res_cor) > R_STRONG,
    concordant = sign(res_cor) == sign(env_cor)
  ) %>%
  filter(!is.na(res_cor), !is.na(env_cor))

cat("Total OTU pairs:", nrow(pair_df), "\n")
cat("Strong residual pairs (|r| >", R_STRONG, "):", sum(pair_df$strong_res),
    sprintf("(%.1f%%)\n", sum(pair_df$strong_res) / nrow(pair_df) * 100))

cat("Quadrant breakdown:\n")
print(table(pair_df$quadrant))
cat("\n")

write.csv(pair_df, file.path(TABLE_DIR, "S3A_pollock_pairs_NB.csv"), row.names = FALSE)

all_pollock <- ggplot(pair_df, aes(x = env_cor, y = res_cor)) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40") +
  geom_point(data = filter(pair_df, !strong_res),
             colour = "grey85", alpha = 0.3, size = 1.2) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_point(data = filter(pair_df, strong_res),
             aes(colour = quadrant),
             alpha = 0.75, size = 1.2) +
  scale_colour_manual(
    name   = "Quadrant",
    values = c("Q1" = "#2166AC", "Q2" = "#F4A582", "Q3" = "#2AB5A0", "Q4" = "#B2182B")
  ) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(x = "Design-based correlation\n", y = "Residual correlation\n") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

ggsave(file.path(FIG_DIR, "S3A_pollock_all_highlighted_NB.png"),
       all_pollock, width = 8, height = 8, dpi = 200)
print(all_pollock)

# ── Quadrant summary table ─────────────────────────────────────────────────────

quadrant_summary <- pair_df %>%
  filter(!is.na(quadrant)) %>%
  group_by(quadrant) %>%
  summarise(
    n_pairs        = n(),
    pct_pairs      = round(n() / nrow(pair_df) * 100, 1),
    n_strong       = sum(strong_res),
    pct_strong     = round(sum(strong_res) / sum(pair_df$strong_res) * 100, 1),
    mean_res_cor   = round(mean(res_cor),  3),
    mean_env_cor   = round(mean(env_cor),  3),
    mean_surprise  = round(mean(surprise), 3),
    n_same_class   = sum(same_class),
    pct_same_class = round(sum(same_class) / sum(pair_df$same_class) * 100, 1),
    .groups        = "drop"
  )

write.csv(quadrant_summary, file.path(TABLE_DIR, "S3A_quadrant_summary_NB.csv"), row.names = FALSE)

quad_sum <- flextable(quadrant_summary) %>%
  set_header_labels(
    quadrant = "Quadrant", n_pairs = "N pairs", pct_pairs = "% pairs",
    n_strong = "N strong", pct_strong = "% strong",
    mean_res_cor = "Mean res. cor.", mean_env_cor = "Mean env. cor.",
    mean_surprise = "Mean CoOc Score", n_same_class = "N same class",
    pct_same_class = "% same class"
  ) %>%
  theme_vanilla() %>%
  autofit()
save_as_docx(quad_sum, path = file.path(TABLE_DIR, "quadrant_summary_biof_NB.docx"))

# ── Top surprising pairs ───────────────────────────────────────────────────────

top_surprise <- pair_df %>%
  filter(strong_res) %>%
  arrange(desc(surprise)) %>%
  slice_head(n = 100) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write.csv(top_surprise, file.path(TABLE_DIR, "S3A_top_surprising_pairs_NB.csv"), row.names = FALSE)


# ==============================================================================
# 3-B: CLUSTERING OF RESIDUAL CORRELATION MATRIX
# ==============================================================================

cat("--- 3-B: Residual correlation clustering ---\n\n")

dist_mat   <- as.dist(1 - getResidualCor(best_model))
hclust_res <- hclust(dist_mat, method = "ward.D2")

sil_widths <- sapply(2:10, function(k) {
  cl  <- cutree(hclust_res, k = k)
  sil <- silhouette(cl, dist_mat)
  mean(sil[, "sil_width"])
})

best_k_sil <- which.max(sil_widths) + 1
cat("Silhouette-optimal k:", best_k_sil, "\n\n")

clusters   <- cutree(hclust_res, k = N_CLUSTERS)
cluster_pal <- setNames(
  RColorBrewer::brewer.pal(max(N_CLUSTERS, 3), "Set2")[1:N_CLUSTERS],
  paste0("C", seq_len(N_CLUSTERS))
)

cluster_df <- data.frame(
  OTU_Label = otu_labels,
  OTU_ID    = otu_ids,
  Cluster   = paste0("C", clusters),
  stringsAsFactors = FALSE
) %>%
  left_join(
    model_data$tax_mat[otu_ids, tax_cols_use, drop = FALSE] %>%
      as.data.frame() %>%
      rownames_to_column("OTU_ID"),
    by = "OTU_ID"
  )

write.csv(cluster_df, file.path(TABLE_DIR, "S3B_cluster_taxonomy_NB.csv"), row.names = FALSE)

cluster_order <- order(clusters)
res_cor_plot  <- getResidualCor(best_model)

annot_heatmap <- data.frame(
  Cluster  = factor(cluster_df$Cluster[cluster_order],
                    levels = paste0("C", seq_len(N_CLUSTERS))),
  Division = replace_na(cluster_df$Division[cluster_order], "Unknown"),
  row.names = otu_labels[cluster_order]
)

n_div <- length(unique(annot_heatmap$Division))
annot_colours <- list(
  Cluster  = cluster_pal,
  Division = setNames(
    colorRampPalette(brewer.pal(9, "Set1"))(n_div),
    unique(annot_heatmap$Division)
  )
)

png(file.path(FIG_DIR, "S3B_residual_heatmap_NB.png"),
    width = 16, height = 14, units = "in", res = 300)
pheatmap(
  res_cor_plot[cluster_order, cluster_order],
  color             = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks            = seq(-1, 1, length.out = 102),
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_row          = cumsum(table(sort(clusters))),
  gaps_col          = cumsum(table(sort(clusters))),
  annotation_row    = annot_heatmap,
  annotation_col    = annot_heatmap,
  annotation_colors = annot_colours,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  treeheight_row    = 0,
  treeheight_col    = 0,
  main              = paste0(
    "Residual OTU co-occurrence (NB abundance model) | ", trace_reduction_pct,
    "% of co-occurrence variation explained by design | k = ", N_CLUSTERS, " clusters"
  ),
  fontsize     = 9,
  border_color = NA,
  legend_breaks  = c(-1, -0.5, 0, 0.5, 1),
  legend_labels  = c("-1", "-0.5", "0", "0.5", "1")
)
dev.off()
cat("Residual heatmap saved.\n\n")


# ==============================================================================
# 3-C: NETWORK ANALYSIS OF STRONG RESIDUAL ASSOCIATIONS
# ==============================================================================

cat("--- 3-C: Network analysis ---\n\n")

edges <- pair_df %>%
  filter(abs(res_cor) > R_STRONG) %>%
  select(OTU_j, OTU_k, res_cor, env_cor, surprise, quadrant,
         same_class, Class_j, Class_k) %>%
  mutate(
    association = if_else(res_cor > 0, "positive", "negative"),
    abs_r       = abs(res_cor)
  )

cat("Edges in network (|r| >", R_STRONG, "):", nrow(edges), "\n\n")

node_otus <- unique(c(edges$OTU_j, edges$OTU_k))

node_df <- tax_lookup %>%
  filter(OTU_Label %in% node_otus) %>%
  mutate(
    Class_plot    = if_else(is.na(Class) | Class == "", "Unclassified", Class),
    Division_plot = if_else(is.na(Division) | Division == "", "Unclassified", Division)
  ) %>%
  left_join(cluster_df %>% select(OTU_Label, Cluster), by = "OTU_Label")

g <- graph_from_data_frame(
  d        = edges %>% select(from = OTU_j, to = OTU_k, weight = abs_r, association, res_cor),
  directed = FALSE,
  vertices = node_df %>% select(name = OTU_Label, everything())
)

V(g)$node_degree      <- as.numeric(igraph::degree(g))
V(g)$node_strength    <- as.numeric(igraph::strength(g, weights = E(g)$weight))
V(g)$node_betweenness <- as.numeric(igraph::betweenness(g, weights = 1 / E(g)$weight, normalized = TRUE))
V(g)$pos_degree <- sapply(V(g)$name, function(otu) {
  sum(edges$association == "positive" & (edges$OTU_j == otu | edges$OTU_k == otu))
})
V(g)$neg_degree <- sapply(V(g)$name, function(otu) {
  sum(edges$association == "negative" & (edges$OTU_j == otu | edges$OTU_k == otu))
})

cat("Network summary:\n")
cat("  Nodes:", vcount(g), "| Edges:", ecount(g),
    "| Density:", round(edge_density(g), 4),
    "| Components:", components(g)$no, "\n\n")

tax_priority <- tax_lookup %>%
  select(OTU_Label, Species, Genus, Family, Order, Class) %>%
  mutate(
    lowest_tax = case_when(
      !is.na(Species) & Species != "" ~ paste0(Species, " (species)"),
      !is.na(Genus)   & Genus   != "" ~ paste0(Genus,   " (genus)"),
      !is.na(Family)  & Family  != "" ~ paste0(Family,  " (family)"),
      !is.na(Order)   & Order   != "" ~ paste0(Order,   " (order)"),
      !is.na(Class)   & Class   != "" ~ paste0(Class,   " (class)"),
      TRUE                            ~ "Unclassified"
    )
  ) %>%
  select(OTU_Label, lowest_tax, Class)

node_metrics <- as_data_frame(g, what = "vertices") %>%
  rename(OTU_Label = name) %>%
  select(OTU_Label, Class_plot, Division_plot, Cluster,
         node_degree, node_strength, node_betweenness,
         pos_degree, neg_degree) %>%
  left_join(tax_priority, by = "OTU_Label") %>%
  mutate(
    signed_dominance = (pos_degree - neg_degree) / (pos_degree + neg_degree),
    across(c(node_strength, node_betweenness), ~ round(.x, 4))
  ) %>%
  select(OTU_Label, lowest_tax, Class, Cluster,
         node_degree, node_strength, node_betweenness,
         pos_degree, neg_degree, signed_dominance) %>%
  arrange(desc(node_betweenness))

write.csv(node_metrics, file.path(TABLE_DIR, "S3C_node_metrics_NB.csv"), row.names = FALSE)
write.csv(as_data_frame(g, what = "edges"),
          file.path(TABLE_DIR, "S3C_edge_list_NB.csv"), row.names = FALSE)


# ==============================================================================
# 3-D: HUB / KEYSTONE TAXA IDENTIFICATION
# ==============================================================================

cat("--- 3-D: Important taxa ---\n\n")

cat("Node metric distributions (to guide threshold choice):\n")
print(quantile(node_metrics$node_betweenness, probs = c(0.5, 0.75, 0.9, 0.95, 0.99)))
print(quantile(node_metrics$node_strength,    probs = c(0.5, 0.75, 0.9, 0.95, 0.99)))
cat("\n")
cat("NOTE: BETWEEN_THRESH / STRENGTH_THRESH below were tuned for the original\n")
cat("PA-model network density. Since NB residual co-occurrence networks can\n")
cat("differ in density, check the quantiles above before trusting these cutoffs.\n\n")

BETWEEN_THRESH  <- 0.05
STRENGTH_THRESH <- 42

metrics <- node_metrics %>%
  mutate(
    high_betweenness = node_betweenness > BETWEEN_THRESH,
    high_strength    = node_strength    > STRENGTH_THRESH
  )

high_between_tbl <- metrics %>%
  filter(high_betweenness) %>%
  arrange(desc(node_betweenness)) %>%
  mutate(category = "High betweenness")

high_strength_tbl <- metrics %>%
  filter(high_strength & !high_betweenness) %>%
  arrange(desc(node_strength)) %>%
  mutate(category = "High strength")

neg_dom_tbl <- node_metrics %>%
  filter((pos_degree + neg_degree) > 0) %>%
  mutate(
    pct_neg       = neg_degree / (pos_degree + neg_degree),
    neg_dominated = pct_neg > 0.5,
    category      = "Negative-dominated"
  ) %>%
  filter(neg_dominated) %>%
  arrange(desc(neg_degree))

important_taxa_table <- bind_rows(high_between_tbl, high_strength_tbl, neg_dom_tbl) %>%
  mutate(
    pct_pos = round(pos_degree / (pos_degree + neg_degree) * 100, 1),
    pct_neg = round(neg_degree / (pos_degree + neg_degree) * 100, 1)
  ) %>%
  select(
    Category = category, Taxon = lowest_tax, Class,
    Betweenness = node_betweenness, Degree = node_degree,
    `Node Strength` = node_strength,
    `Positive edges (%)` = pct_pos, `Negative edges (%)` = pct_neg
  )

hub_labels <- metrics %>%
  filter(high_betweenness | high_strength) %>%
  arrange(desc(node_betweenness), desc(node_strength)) %>%
  mutate(
    hub_number = row_number(),
    flags = case_when(
      high_betweenness & high_strength ~ "betweenness + strength",
      high_betweenness                 ~ "betweenness",
      high_strength                    ~ "strength"
    )
  ) %>%
  select(OTU_Label, flags, hub_number, node_betweenness, node_strength,
         node_degree, Class, pos_degree, neg_degree)

write.csv(important_taxa_table, file.path(TABLE_DIR, "S3D_important_taxa_NB.csv"), row.names = FALSE)

# ── Network plot ──────────────────────────────────────────────────────────────

tg <- as_tbl_graph(g) %>%
  left_join(
    metrics %>% select(name = OTU_Label, high_betweenness, high_strength,
                       node_betweenness, node_strength),
    by = "name"
  ) %>%
  left_join(hub_labels %>% select(name = OTU_Label, hub_number), by = "name") %>%
  mutate(
    node_label = if_else(!is.na(hub_number), as.character(hub_number), ""),
    is_hub     = !is.na(hub_number)
  )

all_classes_net <- unique(as_data_frame(tg, what = "vertices")$Class_plot)
all_classes_net <- all_classes_net[!is.na(all_classes_net)]
n_cls <- length(all_classes_net)

class_cols <- setNames(
  colorRampPalette(
    c(brewer.pal(12, "Paired"), brewer.pal(9, "Set1"), brewer.pal(12, "Set3"))
  )(n_cls),
  all_classes_net
)

set.seed(42)
p_net_numbered <- ggraph(tg, layout = "stress") +
  geom_edge_link(aes(colour = association, width = weight, alpha = weight), show.legend = TRUE) +
  geom_node_point(data = . %>% filter(!is_hub), aes(colour = Class), size = 5, alpha = 0.75) +
  geom_node_point(data = . %>% filter(is_hub), aes(colour = Class), size = 14, alpha = 1) +
  geom_node_point(data = . %>% filter(is_hub), colour = "white", size = 10, alpha = 1) +
  geom_node_text(data = . %>% filter(is_hub), aes(label = node_label, colour = Class),
                 size = 4, fontface = "bold") +
  scale_edge_colour_manual(
    values = c(positive = "#4575b4", negative = "#d73027"),
    name   = "Association",
    labels = c(positive = "Positive (co-occurrence)", negative = "Negative (avoidance)")
  ) +
  scale_edge_width(range = c(0.2, 2.0), guide = "none") +
  scale_edge_alpha(range = c(0.2, 0.8), guide = "none") +
  scale_colour_manual(values = class_cols, name = "Taxonomic class",
                      guide = guide_legend(override.aes = list(shape = 16, size = 4, alpha = 1))) +
  theme_graph(base_family = "sans", base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9),
        legend.position = "right")

ggsave(file.path(FIG_DIR, "S3C_network_numbered_stress_NB.png"),
       p_net_numbered, width = 14, height = 12, dpi = 200)
print(p_net_numbered)
cat("Numbered network plot (NB model) saved.\n\n")


# ==============================================================================
# 3-D SUPPLEMENT: DIAGNOSTIC — CONDITION DISTRIBUTION OF TOP RESIDUAL PAIRS
#
# This part still uses raw PRESENCE/ABSENCE co-occurrence (Jaccard, phi) as
# an informal sanity check, regardless of which family fit the model — the
# question here is "do these OTU pairs co-occur across conditions in the raw
# data", which is a presence/absence question by nature. Presence/absence is
# derived from the abundance count_mat (count > 0) rather than from a
# separate pa_mat, since this script never built a binary response matrix.
# ==============================================================================

cat("--- DIAGNOSTIC: Condition distribution of top residual pairs ---\n\n")

N_DIAG_PAIRS <- 50

diag_pairs <- pair_df %>%
  filter(strong_res) %>%
  arrange(desc(abs(res_cor))) %>%
  slice_head(n = N_DIAG_PAIRS)

pa_full <- (model_data$otu_filt[, model_data$top_otus] > 0) * 1L
colnames(pa_full) <- model_data$otu_id_to_label[model_data$top_otus]

meta_diag <- model_data$metadata
stopifnot(identical(rownames(pa_full), rownames(meta_diag)))

condition_vars <- c("site", "season", "substrate")

diag_results <- lapply(condition_vars, function(cvar) {
  levels_cvar <- levels(meta_diag[[cvar]])
  pair_rows <- lapply(seq_len(nrow(diag_pairs)), function(i) {
    otu_j <- diag_pairs$OTU_j[i]
    otu_k <- diag_pairs$OTU_k[i]
    if (!otu_j %in% colnames(pa_full) || !otu_k %in% colnames(pa_full)) return(NULL)
    level_rows <- lapply(levels_cvar, function(lv) {
      idx <- which(meta_diag[[cvar]] == lv)
      if (length(idx) < 3) return(NULL)
      pj <- pa_full[idx, otu_j]; pk <- pa_full[idx, otu_k]
      n_samples      <- length(idx)
      n_both_present <- sum(pj == 1 & pk == 1)
      n_either       <- sum(pj == 1 | pk == 1)
      jaccard <- if (n_either > 0) n_both_present / n_either else NA_real_
      phi <- if (sd(pj) > 0 && sd(pk) > 0) cor(pj, pk) else NA_real_
      data.frame(
        condition_var = cvar, condition_level = lv, OTU_j = otu_j, OTU_k = otu_k,
        res_cor = diag_pairs$res_cor[i], n_samples = n_samples,
        n_both = n_both_present, n_either = n_either,
        jaccard = round(jaccard, 3), phi = round(phi, 3), stringsAsFactors = FALSE
      )
    })
    do.call(rbind, level_rows)
  })
  do.call(rbind, pair_rows)
})

diag_df <- do.call(rbind, diag_results)

pair_summary <- diag_df %>%
  group_by(condition_var, OTU_j, OTU_k, res_cor) %>%
  summarise(
    n_levels_total    = n(),
    n_levels_cooccur  = sum(jaccard > 0, na.rm = TRUE),
    n_levels_positive = sum(phi > 0, na.rm = TRUE),
    n_levels_strong   = sum(phi > 0.2, na.rm = TRUE),
    mean_jaccard      = round(mean(jaccard, na.rm = TRUE), 3),
    mean_phi          = round(mean(phi, na.rm = TRUE), 3),
    min_phi           = round(min(phi, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  mutate(
    coverage = round(n_levels_cooccur / n_levels_total, 2),
    distribution = case_when(
      coverage == 1.00 ~ "All levels",
      coverage >= 0.50 ~ "Most levels",
      coverage >  0    ~ "Few levels",
      TRUE             ~ "Absent"
    )
  ) %>%
  arrange(condition_var, desc(mean_phi))

write.csv(pair_summary, file.path(TABLE_DIR, "S3_pair_distribution_summary_NB.csv"), row.names = FALSE)

cat("=== DISTRIBUTION FLAGS ===\n\n")
for (cv in condition_vars) {
  cat("Condition:", cv, "\n")
  tbl <- pair_summary %>% filter(condition_var == cv) %>%
    count(distribution) %>% mutate(pct = round(n / sum(n) * 100, 1))
  print(tbl)
  cat("\n")
}

cat("=== NB SENSITIVITY CHECK SCRIPT COMPLETE ===\n")
cat("Compare AIC/BIC, quadrant breakdowns, network structure, and hub taxa\n")
cat("against the original binomial PA model outputs to assess sensitivity.\n")
