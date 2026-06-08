################################################################################
# GLLVM PRESENCE-ABSENCE ANALYSIS PIPELINE — SEASONAL FOULING COMMUNITY
#
# PRIMARY MODEL: best_inc_f_v2  (PA, M2 all-2-way interactions, fouling)
# NULL MODEL:    uncons_inc_f_v2 (PA, M7 unconstrained, fouling)
################################################################################



library(gllvm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(patchwork)
library(knitr)
library(scales)
library(vegan)
library(RColorBrewer)
library(tibble)
library(ggsci)
library(knitr)
library(kableExtra)
library(ggrepel)
library(flextable)
library(officer)

cat("=== GLLVM ANALYSIS SCRIPT ===\n\n")

setwd("~/Github/MAPLE_Seasonal_Plastisphere/Scripts/03.StatisticalAnalysis")

SAVE_DIR_BASE  <- "~/Github/MAPLE_Seasonal_Plastisphere/Processed_data/gllvm_models"

OUT_DIR   <- "~/Github/MAPLE_Seasonal_Plastisphere/Results/gllvm_results"
TABLE_DIR <- file.path(OUT_DIR, "Tables")
FIG_DIR   <- file.path(OUT_DIR, "Figures")
for (d in c(OUT_DIR, FIG_DIR, TABLE_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# SECTION 0: LOAD DATA AND MODELS
# ==============================================================================


model_data <- readRDS(file.path(SAVE_DIR_BASE, "data_for_gllvm_FINAL.rds"))
metadata   <- model_data$metadata

patch_model <- function(fit, num_lv, include_X = TRUE) {
  fit$call$y           <- fit$y
  fit$call$num.lv      <- num_lv
  fit$call$studyDesign <- fit$studyDesign
  if (!is.null(fit$offset)) fit$call$offset <- fit$offset
  if (include_X && !is.null(fit$X)) fit$call$X <- fit$X
  fit
}

best_model <- patch_model(
  readRDS(file.path(SAVE_DIR_BASE, "GLLVM_final_incM2.rds")),
  num_lv = 2, include_X = TRUE
)

null_model <- patch_model(
  readRDS(file.path(SAVE_DIR_BASE, "GLLVM_null_incM7.rds")),
  num_lv = 2, include_X = FALSE
)

stopifnot(!is.null(best_model), !is.null(null_model))
cat("Models loaded successfully.\n\n")


# ==============================================================================
# SECTION 1: COEFFICIENT EXTRACTION
# Extract all OTU × predictor coefficients with taxonomy, significance
# (Benjamini-Hochberg corrected), and reliability flags.
# ==============================================================================

SEASON_LEVELS    <- c("Winter", "Spring", "Summer", "Fall", "Winter2")
SUBSTRATE_LEVELS <- c("Glass", "PE", "Weathered_PE", "PET", "Weathered_PET")
SITE_LEVELS      <- c("SELVA", "TBS")
SUBSTRATE_LABELS <- c("Glass", "PE", "W-PE", "PET", "W-PET")



extract_coefficients <- function(model, d, track) {
  
  Xcoef    <- coef(model)$Xcoef
  se_Xcoef <- model$sd$Xcoef
  
  # Relative abundance weights — used to prioritise ecologically
  # important OTUs in summaries (rare OTUs count less than dominant ones)
  top_otus  <- names(sort(colSums(d$otu_filt), decreasing = TRUE)[1:200])
  otu_sub   <- d$otu_filt[, top_otus]
  rel_abund <- colMeans(otu_sub / rowSums(otu_sub))
  names(rel_abund) <- d$otu_id_to_label[top_otus]
  
  # z-statistics and two-tailed p-values from the normal approximation
  # (standard for GLLVM variational approximation; Niku et al. 2019 MEE)
  z_mat    <- Xcoef / se_Xcoef
  p_mat    <- 2 * pnorm(abs(z_mat), lower.tail = FALSE)
  
  # BH correction applied globally across all OTU × term combinations
  p_BH_vec <- p.adjust(as.vector(p_mat), method = "BH")
  p_BH_mat <- matrix(p_BH_vec, nrow = nrow(p_mat), ncol = ncol(p_mat),
                     dimnames = dimnames(p_mat))
  
  cat("  Total OTU × term tests:", length(p_mat), "\n")
  cat("  Significant after BH (α = 0.05):",
      sum(p_BH_mat < 0.05, na.rm = TRUE), "\n\n")
  
  res <- do.call(rbind, lapply(colnames(Xcoef), function(term) {
    data.frame(
      OTU_Label = rownames(Xcoef),
      Term      = term,
      Estimate  = Xcoef[, term],
      SE        = se_Xcoef[, term],
      z         = z_mat[, term],
      p_raw     = p_mat[, term],
      p_BH      = p_BH_mat[, term],
      stringsAsFactors = FALSE
    )
  }))
  
  res <- res %>%
    mutate(
      lower    = Estimate - 1.96 * SE,
      upper    = Estimate + 1.96 * SE,
      # sig_BH: statistically significant after multiple testing correction
      sig_BH   = p_BH < 0.05,
      # CI_excl_zero: the 95% CI does not overlap zero — stronger evidence
      # than p-value alone because it accounts for both direction and precision
      CI_excl_zero     = lower > 0 | upper < 0,
      # reliable: both significant AND CI excludes zero AND |coef| ≤ 10
      # We use this as the primary filter for biological interpretation.
      # Large |coef| values (> 15) arise from near-complete separation:
      # a taxon essentially absent or always present in one condition.
      # Direction remains valid but magnitude cannot be trusted.
      extreme_estimate = abs(Estimate) > 15,
      reliable         = sig_BH & CI_excl_zero & !extreme_estimate,
      # For directional analyses we also keep extreme-but-significant ones,
      # flagged separately: sig_any = significant regardless of magnitude
      sig_any          = sig_BH & CI_excl_zero,
      # Predictor classification for grouped summaries
      predictor_type = case_when(
        grepl("site", Term) & grepl("season", Term) & grepl("substrate", Term) ~ "3-way",
        grepl("site", Term) & grepl("season",    Term) ~ "site × season",
        grepl("site", Term) & grepl("substrate", Term) ~ "site × substrate",
        grepl("season", Term) & grepl("substrate", Term) ~ "season × substrate",
        grepl("site",      Term) ~ "site",
        grepl("season",    Term) ~ "season",
        grepl("substrate", Term) ~ "substrate",
        TRUE ~ "other"
      )
    ) %>%
    # Remove numerical artefacts (failed convergence for individual terms)
    filter(SE > 1e-10, SE < 1e4)
  
  # Attach taxonomy
  label_to_id <- setNames(names(d$otu_id_to_label), d$otu_id_to_label)
  tax_cols    <- intersect(
    c("Domain","Supergroup","Division","Subdivision","Class","Order","Family","Genus","Species"),
    colnames(d$tax_mat)
  )
  tax_df <- d$tax_mat[label_to_id[unique(res$OTU_Label)], tax_cols, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("OTU_ID") %>%
    mutate(OTU_Label = d$otu_id_to_label[OTU_ID])
  
  res <- res %>%
    left_join(select(tax_df, -OTU_ID), by = "OTU_Label") %>%
    mutate(
      mean_rel_abund = rel_abund[OTU_Label],
      abund_weight   = mean_rel_abund / max(mean_rel_abund, na.rm = TRUE)
    )
  
  cat("  Extreme estimates (|coef| > ", 15, "):", sum(res$extreme_estimate, na.rm = TRUE), "\n")
  cat("  Reliable (sig + CI excl. 0 + not extreme):", sum(res$reliable, na.rm = TRUE), "\n\n")
  res
}

results <- extract_coefficients(best_model, model_data, "fouling_PA")

write.csv(results, file.path(TABLE_DIR, "all_coefficients_PA_fouling.csv"), row.names = FALSE)
cat("Full coefficient table saved.\n\n")


# ==============================================================================
#  DOES SITE, SEASON AND SUBSTRATE STRUCTURE THE COMMUNITY?
# ==============================================================================

cat("=== Q1: Community structuring — variance partitioning ===\n\n")


### FROM GLLVM PACKAGE

vp_gllvm = VP(best_model)
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
doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "vp_gllvm.docx")

# ==============================================================================
# Substrate effects focus: SEASON × SUBSTRATE ANALYSIS
# ==============================================================================

# ── Helper: build OTU display label ──────────────────────────────────────────
# Priority: Species > Genus + "sp." + ASV > (keep existing label if Genus NA)
make_otu_label <- function(species, genus, otu_label) {
  case_when(
    !is.na(species) & species != "" ~ species,
    !is.na(genus)   & genus   != "" ~ paste0(genus, " sp. ", otu_label),
    TRUE                            ~ otu_label   # fallback: keep label as-is
  )
}

#1. compute total effects per site

compute_total_effects_bysite <- function(results) {
  
  season_levels    <- c("Winter", "Spring", "Summer", "Fall", "Winter2")
  substrate_levels <- c("Glass", "PE", "Weathered_PE", "PET", "Weathered_PET")
  site_levels      <- c("SELVA", "TBS")
  
  otu_meta <- results %>%
    distinct(OTU_Label, Class, Family, Genus, Species,
             mean_rel_abund, abund_weight)
  
  two_way <- function(a, b) c(paste0(a, ":", b), paste0(b, ":", a))
  all_terms <- unique(results$Term)
  
  rows <- lapply(site_levels, function(site) {
    lapply(season_levels, function(seas) {
      lapply(substrate_levels, function(sub) {
        
        # ── Reference condition: SELVA, Winter, Glass ─────────────────────
        if (site == "SELVA" & seas == "Winter" & sub == "Glass") {
          return(otu_meta %>%
                   mutate(
                     total_estimate = 0,
                     total_SE       = 0,
                     total_lower    = 0,
                     total_upper    = 0,
                     extreme        = FALSE,
                     any_sig        = FALSE,
                     site           = site,
                     season         = seas,
                     substrate      = sub
                   ))
        }
        
        # ── Build list of terms to sum ─────────────────────────────────────
        # Site terms (if TBS)
        site_term         <- if (site == "TBS") "siteTBS"
        site_season_term  <- if (site == "TBS" & seas != "Winter")
          two_way("siteTBS", paste0("season", seas))
        site_sub_term     <- if (site == "TBS" & sub != "Glass")
          two_way("siteTBS", paste0("substrate", sub))
        
        # Season and substrate terms
        season_term       <- if (seas != "Winter") paste0("season", seas)
        sub_term          <- if (sub  != "Glass")  paste0("substrate", sub)
        seas_sub_term     <- if (seas != "Winter" & sub != "Glass")
          two_way(paste0("season", seas), paste0("substrate", sub))
        
        terms_to_sum <- c(
          site_term,
          season_term,
          sub_term,
          site_season_term,
          site_sub_term,
          seas_sub_term
        )
        terms_to_sum <- terms_to_sum[terms_to_sum %in% all_terms]
        
        if (length(terms_to_sum) == 0) {
          warning("No terms found for: ", site, " / ", seas, " × ", sub,
                  " — returning zeros.")
          return(otu_meta %>%
                   mutate(
                     total_estimate = 0,
                     total_SE       = 0,
                     total_lower    = 0,
                     total_upper    = 0,
                     extreme        = FALSE,
                     any_sig        = FALSE,
                     site           = site,
                     season         = seas,
                     substrate      = sub
                   ))
        }
        
        results %>%
          filter(Term %in% terms_to_sum) %>%
          group_by(OTU_Label) %>%
          summarise(
            total_estimate = sum(Estimate,  na.rm = TRUE),
            total_SE       = sqrt(sum(SE^2, na.rm = TRUE)),
            .groups        = "drop"
          ) %>%
          mutate(
            total_lower = total_estimate - 1.96 * total_SE,
            total_upper = total_estimate + 1.96 * total_SE,
            extreme     = abs(total_estimate) > 15,
            any_sig     = total_lower > 0 | total_upper < 0,
            site        = site,
            season      = seas,
            substrate   = sub
          ) %>%
          left_join(otu_meta, by = "OTU_Label")
      })
    })
  })
  
  bind_rows(unlist(rows, recursive = FALSE)) %>%
    mutate(
      site      = factor(site,      levels = site_levels),
      season    = factor(season,    levels = season_levels),
      substrate = factor(substrate, levels = substrate_levels),
      display_label = make_otu_label(Species, Genus, OTU_Label)
    )
}

# Run it
total_effects_bysite <- compute_total_effects_bysite(results)

# Sanity check
cat("Rows:", nrow(total_effects_bysite), "\n")
cat("OTUs with any_sig TRUE:", sum(total_effects_bysite$any_sig, na.rm = TRUE), "\n")
cat("Extreme:", sum(total_effects_bysite$extreme, na.rm = TRUE), "\n")

# Save
write.csv(total_effects_bysite,
          file.path(TABLE_DIR, "Q3_total_effects_bysite.csv"),
          row.names = FALSE)





# ==============================================================================
# SUBSTRATE RESPONSE MATRIX
# One row per OTU × season × site
# Columns: PE, wPE, PET, wPET — direction vs glass reference
# ▲ = significantly elevated, ▼ = significantly depleted, · = ns
# ==============================================================================

cat("=== Substrate response matrix ===\n\n")

SUBSTRATE_PLASTIC_LEVELS <- c("PE", "Weathered_PE", "PET", "Weathered_PET")
SUBSTRATE_PLASTIC_LABELS <- c("PE", "wPE", "PET", "wPET")

# ── Step 1: Assign direction per OTU × substrate × season × site ──────────────
# Uses total_effects_bysite computed earlier.
# any_sig = TRUE means the 95% CI of the total effect excludes zero.
# extreme is already flagged in total_effects_bysite (|total_estimate| > 15).

direction_long <- total_effects_bysite %>%
  filter(
    substrate %in% SUBSTRATE_PLASTIC_LEVELS,
    !extreme                              # exclude unreliable large estimates
  ) %>%
  mutate(
    direction = case_when(
      any_sig & total_estimate > 0 ~ "higher_occurrence",
      any_sig & total_estimate < 0 ~ "lower_occurrence",
      TRUE                          ~ "ns"
    ),
    # readable substrate label
    substrate_label = factor(
      SUBSTRATE_PLASTIC_LABELS[match(as.character(substrate),
                                     SUBSTRATE_PLASTIC_LEVELS)],
      levels = SUBSTRATE_PLASTIC_LABELS
    )
  ) %>%
  select(site, season, OTU_Label, display_label, Class,
         substrate_label, direction, total_estimate, abund_weight)

# ── Step 2: Pivot to wide matrix ───────────────────────────────────────────────
# Each row = OTU × season × site
# Each column = one plastic substrate

response_matrix <- direction_long %>%
  pivot_wider(
    id_cols     = c(site, season, OTU_Label, display_label, Class, abund_weight),
    names_from  = substrate_label,
    values_from = direction,
    values_fill = "ns"
  )

# ── Step 3: Keep only rows with at least one significant result ────────────────
plastic_cols <- SUBSTRATE_PLASTIC_LABELS  # c("PE", "wPE", "PET", "wPET")

response_matrix <- response_matrix %>%
  filter(
    rowSums(across(all_of(plastic_cols), ~ . != "ns")) > 0
  ) %>%
  arrange(site, season, Class, display_label)

cat("Total OTU × season rows with ≥1 significant result:\n")
print(table(response_matrix$site, response_matrix$season))
cat("\n")

# ── Step 4: Split by site and save ────────────────────────────────────────────
response_matrix_SELVA <- response_matrix %>% filter(site == "SELVA")
response_matrix_TBS   <- response_matrix %>% filter(site == "TBS")

write.csv(response_matrix,
          file.path(TABLE_DIR, "substrate_response_matrix_all.csv"),
          row.names = FALSE)
write.csv(response_matrix_SELVA,
          file.path(TABLE_DIR, "substrate_response_matrix_SELVA.csv"),
          row.names = FALSE)
write.csv(response_matrix_TBS,
          file.path(TABLE_DIR, "substrate_response_matrix_TBS.csv"),
          row.names = FALSE)

cat("Response matrices saved.\n\n")

# ── Step 5: Quick summary — how many OTUs per direction pattern ────────────────
# Collapse the 4 substrate columns into a single pattern string for inspection
# e.g. "elevated/ns/ns/depleted" — useful for seeing what patterns dominate

response_matrix <- response_matrix %>%
  mutate(
    pattern = paste(PE, wPE, PET, wPET, sep = "/")
    # patterns will now read e.g. "higher_occurrence/ns/ns/lower_occurrence"
  )

cat("Most common patterns (top 20):\n")
response_matrix %>%
  count(site, pattern, sort = TRUE) %>%
  group_by(site) %>%
  slice_max(n, n = 20) %>%
  print(n = 40)


# ==============================================================================
# SUBSTRATE RESPONSE CLASSIFICATION
# Classification at OTU × site × substrate level across seasons
# OTUs can belong to multiple groups
# ==============================================================================

cat("=== Substrate response classification ===\n\n")

SUBSTRATE_PLASTIC_LEVELS <- c("PE", "Weathered_PE", "PET", "Weathered_PET")
SUBSTRATE_PLASTIC_LABELS <- c("PE", "wPE", "PET", "wPET")

# ── Step 1: Direction per OTU × site × substrate × season ─────────────────────
# Starting from total_effects_bysite (already computed)
# any_sig + direction already established in direction_long

direction_long <- total_effects_bysite %>%
  filter(
    substrate %in% SUBSTRATE_PLASTIC_LEVELS,
    !extreme
  ) %>%
  mutate(
    direction = case_when(
      any_sig & total_estimate > 0 ~ "higher_occurrence",
      any_sig & total_estimate < 0 ~ "lower_occurrence",
      TRUE                          ~ "ns"
    ),
    substrate_label = factor(
      SUBSTRATE_PLASTIC_LABELS[match(as.character(substrate),
                                     SUBSTRATE_PLASTIC_LEVELS)],
      levels = SUBSTRATE_PLASTIC_LABELS
    )
  ) %>%
  select(site, season, OTU_Label, display_label, Class,
         substrate_label, direction, total_estimate, abund_weight)


# ── Step 2: Summarise per OTU × site × substrate across seasons ───────────────
# For each OTU × site × substrate combination, record:
#   - whether it was ever elevated across seasons
#   - whether it was ever depleted across seasons
#   - how many seasons it was significant at all

substrate_summary <- direction_long %>%
  group_by(site, OTU_Label, display_label, Class,
           abund_weight, substrate_label) %>%
  summarise(
    ever_higher          = any(direction == "higher_occurrence"),
    ever_lower           = any(direction == "lower_occurrence"),
    n_higher             = sum(direction == "higher_occurrence"),
    n_lower              = sum(direction == "lower_occurrence"),
    n_sig                = sum(direction != "ns"),
    n_conditions_tested  = n(),
    prop_sig             = round(n_sig / n_conditions_tested, 2),
    mean_estimate        = round(mean(total_estimate[direction != "ns"],
                                      na.rm = TRUE), 2),
    seasons_higher       = paste(sort(season[direction == "higher_occurrence"]),
                                 collapse = ", "),
    seasons_lower        = paste(sort(season[direction == "lower_occurrence"]),
                                 collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    # Mixed: higher occurrence in some seasons, lower in others — same substrate
    is_mixed    = ever_higher & ever_lower,
    # Consistently higher occurrence: never lower
    is_higher   = ever_higher & !ever_lower,
    # Consistently lower occurrence: never higher
    is_lower    = ever_lower  & !ever_higher,
    consistency = case_when(
      n_sig >= 3 ~ "consistent",
      n_sig == 2 ~ "moderate",
      n_sig == 1 ~ "occasional",
      TRUE       ~ "none"
    )
  )


# ── Step 3: Count enriched/avoided substrates per OTU × site ──────────────────
# Used to classify plastic generalists and avoiders

otu_site_counts <- substrate_summary %>%
  group_by(site, OTU_Label, display_label, Class, abund_weight) %>%
  summarise(
    n_higher_substrates = sum(is_higher),
    n_lower_substrates  = sum(is_lower),
    higher_substrates   = paste(substrate_label[is_higher], collapse = ", "),
    lower_substrates    = paste(substrate_label[is_lower],  collapse = ", "),
    .groups = "drop"
  )


# ── Step 4: Assign group membership ───────────────────────────────────────────
# Groups are NOT mutually exclusive — an OTU can appear in multiple groups.
# Classification is at OTU × site × substrate level.

# Plastic generalist — higher occurrence on all 4 plastic substrates
# Most robust: signal required across all substrate types
group_plastic_generalist <- substrate_summary %>%
  filter(is_higher) %>%
  left_join(select(otu_site_counts, site, OTU_Label,
                   n_higher_substrates, higher_substrates),
            by = c("site", "OTU_Label")) %>%
  filter(n_higher_substrates == 4) %>%
  mutate(group = "Plastic generalist") %>%
  select(site, OTU_Label, display_label, Class,
         abund_weight, substrate_label, group, consistency,
         n_higher, n_lower, n_sig,
         seasons_higher, seasons_lower, mean_estimate)

# Plastic avoider — lower occurrence on 3+ plastic substrates
group_plastic_avoider <- substrate_summary %>%
  filter(is_lower) %>%
  left_join(select(otu_site_counts, site, OTU_Label,
                   n_lower_substrates, lower_substrates),
            by = c("site", "OTU_Label")) %>%
  filter(n_lower_substrates == 4) %>%
  mutate(group = "Plastic avoider") %>%
  select(site, OTU_Label, display_label, Class,
         abund_weight, substrate_label, group, consistency,
         n_higher, n_lower, n_sig,
         seasons_higher, seasons_lower, mean_estimate)

# Multi-substrate — higher occurrence on 2-3 substrates
# Intermediate evidence, not pan-plastic
group_multi_higher <- substrate_summary %>%
  filter(is_higher) %>%
  left_join(select(otu_site_counts, site, OTU_Label,
                   n_higher_substrates, higher_substrates),
            by = c("site", "OTU_Label")) %>%
  filter(n_higher_substrates %in% c(2, 3)) %>%
  mutate(group = paste0("Multi-substrate (",
                        n_higher_substrates, " substrates)")) %>%
  select(site, OTU_Label, display_label, Class,
         abund_weight, substrate_label, group, consistency,
         n_higher, n_lower, n_sig,
         seasons_higher, seasons_lower, mean_estimate)


# ── Step 5: Combine all groups ─────────────────────────────────────────────────

all_groups <- bind_rows(
  group_plastic_generalist,
  group_plastic_avoider,
  group_multi_higher
)

# Summary counts
cat("Group counts per site:\n")
all_groups %>%
  distinct(site, OTU_Label, group) %>%
  count(site, group) %>%
  arrange(site, group) %>%
  print()

cat("\nConsistency breakdown per group:\n")
all_groups %>%
  distinct(site, OTU_Label, substrate_label, group, consistency) %>%
  count(group, consistency) %>%
  arrange(group, consistency) %>%
  print()



# ── Step 7: Save ──────────────────────────────────────────────────────────────

write.csv(all_groups,
          file.path(TABLE_DIR, "substrate_response_groups.csv"),
          row.names = FALSE)

write.csv(substrate_summary,
          file.path(TABLE_DIR, "substrate_summary_per_OTU_site_substrate.csv"),
          row.names = FALSE)

cat("Classification tables saved.\n\n")





## RESULTS

# ── Plastic generalists ───────────────────────────────────────────────────────

cat("=== Plastic generalists ===\n\n")

# Count per site
cat("Count per site:\n")
group_plastic_generalist %>%
  distinct(site, OTU_Label) %>%
  count(site) %>%
  arrange(site) %>%
  print()

# Taxa list with consistency
cat("\nTaxa and consistency:\n")
group_plastic_generalist %>%
  distinct(site, OTU_Label, display_label, Class) %>%
  arrange(site, Class, display_label) %>%
  print()

# How many seasons of higher occurrence per OTU × substrate
cat("\nSeasons of higher occurrence per OTU × substrate:\n")
group_plastic_generalist %>%
  select(site, OTU_Label, display_label, substrate_label,
         n_higher, n_sig, consistency, seasons_higher) %>%
  arrange(site, display_label, substrate_label) %>%
  print(n = Inf)


# ── Plastic avoiders ──────────────────────────────────────────────────────────

cat("=== Plastic avoiders ===\n\n")

# Count per site
cat("Count per site:\n")
group_plastic_avoider %>%
  distinct(site, OTU_Label) %>%
  count(site) %>%
  arrange(site) %>%
  print()

# Taxa list with consistency
cat("\nTaxa and consistency:\n")
group_plastic_avoider %>%
  distinct(site, OTU_Label, display_label, Class, consistency) %>%
  arrange(site, Class, display_label) %>%
  print()

# Which substrates avoided per OTU
cat("\nSubstrates with lower occurrence per OTU:\n")
group_plastic_avoider %>%
  select(site, OTU_Label, display_label, substrate_label,
         n_lower, n_sig, consistency, seasons_lower) %>%
  arrange(site, display_label, substrate_label) %>%
  print(n = Inf)


# ── Multi-substrate higher occurrence ─────────────────────────────────────────

cat("=== Multi-substrate higher occurrence (2-3 substrates) ===\n\n")

# Count per site and group
cat("Count per site and group:\n")
group_multi_higher %>%
  distinct(site, OTU_Label, group) %>%
  count(site, group) %>%
  arrange(site, group) %>%
  print()

# Taxa list
cat("\nTaxa:\n")
group_multi_higher %>%
  distinct(site, OTU_Label, display_label, Class, group, consistency) %>%
  arrange(site, group, Class) %>%
  print()

# Which substrates per OTU
cat("\nSubstrates with higher occurrence per OTU:\n")
group_multi_higher %>%
  select(site, OTU_Label, display_label, substrate_label,
         n_higher, n_sig, consistency, seasons_higher) %>%
  arrange(site, display_label, substrate_label) %>%
  print(n = Inf)


# ── Consistency breakdown across all groups ───────────────────────────────────

cat("=== Consistency breakdown across all groups ===\n\n")

all_groups %>%
  distinct(site, OTU_Label, substrate_label, group, consistency) %>%
  count(group, consistency) %>%
  arrange(group, consistency) %>%
  print()


# ── Cross-site responses ──────────────────────────────────────────────────────

cat("=== Cross-site responses ===\n\n")

# OTUs appearing in both sites within the same group
cross_site <- all_groups %>%
  distinct(site, OTU_Label, display_label, Class, group) %>%
  group_by(OTU_Label, display_label, Class, group) %>%
  filter(n_distinct(site) > 1) %>%
  ungroup() %>%
  arrange(group, Class)

# Count per group
cat("OTUs with consistent group membership at both sites:\n")
cross_site %>%
  distinct(OTU_Label, group) %>%
  count(group) %>%
  arrange(group) %>%
  print()

# Full list with consistency per site
cat("\nFull list:\n")
all_groups %>%
  distinct(site, OTU_Label, display_label, Class, group, consistency) %>%
  semi_join(cross_site, by = c("OTU_Label", "group")) %>%
  arrange(group, Class, OTU_Label, site) %>%
  print(n = Inf)



#Plastic generalists table

library(dplyr)
library(flextable)
library(officer)

# ── Step 1: Pivot to one column per substrate ─────────────────────────────────
generalist_wide <- group_plastic_generalist %>%
  mutate(substrate_label = as.character(substrate_label)) %>%
  select(site, OTU_Label, display_label, Class,
         substrate_label, seasons_higher) %>%
  pivot_wider(
    names_from  = substrate_label,
    values_from = seasons_higher,
    values_fill = "—"
  ) %>%
  arrange(site, Class, display_label) %>%
  select(Taxon = display_label, Class, Site = site, PE, wPE, PET, wPET)

# ── Step 2: Build flextable ───────────────────────────────────────────────────
# A4 with 2.5 cm margins = ~16.5 cm = 6.5 inches usable width
ft <- flextable(generalist_wide) %>%
  italic(j = "Taxon") %>%
  bold(part = "header") %>%
  bg(part = "header", bg = "#D5E8F0") %>%
  bg(i = seq(2, nrow(generalist_wide), by = 2),
     bg = "#F5F5F5", part = "body") %>%
  width(j = "Taxon", width = 1.6) %>%
  width(j = "Class", width = 1.4) %>%
  width(j = "Site",  width = 0.5) %>%
  width(j = "PE",    width = 0.75) %>%
  width(j = "wPE",   width = 0.75) %>%
  width(j = "PET",   width = 0.75) %>%
  width(j = "wPET",  width = 0.75) %>%
  align(align = "left", part = "all") %>%
  set_table_properties(layout = "fixed") %>%
  fontsize(size = 9, part = "header") %>%
  fontsize(size = 9, part = "header") %>%
  font(fontname = "Arial", part = "all") %>%
  border_outer(part = "all", border = fp_border(color = "#AAAAAA", width = 0.5)) %>%
  border_inner(part = "all", border = fp_border(color = "#CCCCCC", width = 0.5)) %>%
  add_header_row(
    values    = c("", "", "", "Seasons with higher occurrence"),
    colwidths = c(1, 1, 1, 4)
  ) %>%
  bold(part = "header") %>%
  align(i = 1, align = "center", part = "header") %>%
  bg(i = 1, part = "header", bg = "#D5E8F0")

# ── Step 3: Export ────────────────────────────────────────────────────────────
doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = file.path(TABLE_DIR, "plastic_generalists_table.docx"))









