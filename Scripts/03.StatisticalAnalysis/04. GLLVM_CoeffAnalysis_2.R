################################################################################
# GLLVM PRESENCE-ABSENCE ANALYSIS PIPELINE — SEASONAL FOULING COMMUNITY
#
# MODELS: Two per-site models (SELVA, TBS) with formula ~ season/substrate
#         Equivalent to ~ season + season:substrate, giving direct single
#         coefficients for each plastic-vs-glass contrast within each season.
#         No post-hoc coefficient summing or covariance propagation needed.
#
# NULL MODELS: Two per-site unconstrained models (no fixed effects)
#              Used for variance partitioning per site.
#
# KEY CHANGE FROM PREVIOUS VERSION:
#   Old: one global model ~ (site + season + substrate)^2
#        plastic-vs-glass effect in a given season required summing
#        substratePE + seasonX:substratePE, ignoring their covariance
#   New: separate per-site models ~ season/substrate
#        plastic-vs-glass effect in season X = single coefficient
#        seasonX:substratePE with its own SE and p-value directly
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
library(kableExtra)
library(ggrepel)
library(flextable)
library(officer)

cat("=== GLLVM PER-SITE ANALYSIS SCRIPT ===\n\n")

setwd("~/Github/MAPLE_Seasonal_Plastisphere/Scripts/03.StatisticalAnalysis")

SAVE_DIR_BASE <- "~/Github/MAPLE_Seasonal_Plastisphere/Processed_data/gllvm_models"
OUT_DIR       <- "~/Github/MAPLE_Seasonal_Plastisphere/Results/gllvm_results"
TABLE_DIR     <- file.path(OUT_DIR, "Tables")
FIG_DIR       <- file.path(OUT_DIR, "Figures")
for (d in c(OUT_DIR, FIG_DIR, TABLE_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# CONSTANTS
# ==============================================================================

SEED                     <- 123
NUM_LV                   <- 2
SEASON_LEVELS            <- c("Winter", "Spring", "Summer", "Fall", "Winter2")
SUBSTRATE_LEVELS         <- c("Glass", "PE", "Weathered_PE", "PET", "Weathered_PET")
SITE_LEVELS              <- c("SELVA", "TBS")
SUBSTRATE_PLASTIC_LEVELS <- c("PE", "Weathered_PE", "PET", "Weathered_PET")
SUBSTRATE_PLASTIC_LABELS <- c("PE", "wPE", "PET", "wPET")


# ==============================================================================
# SECTION 1: LOAD DATA AND SPLIT BY SITE
# ==============================================================================

cat("=== SECTION 1: LOAD DATA ===\n\n")

model_data <- readRDS(file.path(SAVE_DIR_BASE, "data_for_gllvm_FINAL.rds"))

# Convenience handles
otu_filt        <- model_data$otu_filt
tax_mat         <- model_data$tax_mat
metadata        <- model_data$metadata
otu_id_to_label <- model_data$otu_id_to_label
pa_mat_full     <- model_data$pa_mat

cat("Full PA matrix:", nrow(pa_mat_full), "x", ncol(pa_mat_full), "\n\n")

# ── Split by site ──────────────────────────────────────────────────────────────
# The OTU set is identical across both sites (same columns) so results are
# directly comparable. Only rows (samples) differ.

cat("--- Splitting by site ---\n")

pa_mat  <- list()
meta    <- list()
sd_site <- list()

for (s in SITE_LEVELS) {
  idx         <- which(metadata$site == s)
  pa_mat[[s]] <- pa_mat_full[idx, ]
  meta[[s]]   <- metadata[idx, ]
  meta[[s]]$season    <- droplevels(meta[[s]]$season)
  meta[[s]]$substrate <- droplevels(meta[[s]]$substrate)
  sd_site[[s]] <- data.frame(
    clean_sample_names = meta[[s]]$clean_sample_names
  )
  cat(" ", s, ":", nrow(pa_mat[[s]]), "samples |",
      nlevels(meta[[s]]$season), "seasons |",
      nlevels(meta[[s]]$substrate), "substrates\n")
}
cat("\n")


# ==============================================================================
# SECTION 2: MODEL FITTING
#
# Formula: ~ season/substrate
#   Expands to: ~ season + season:substrate
#   Gives one coefficient per plastic substrate per season, directly
#   representing plastic-vs-glass within that season.
#   Reference: Winter:Glass (absorbed into per-OTU intercept)
#
# All other settings identical to original fitting script:
#   family = binomial, num.lv = 2, VA, probit, same random effect structure
# ==============================================================================

cat("=== SECTION 2: MODEL FITTING ===\n\n")

fit_best <- list()
fit_null <- list()

for (s in SITE_LEVELS) {
  
  cat("--- Fitting BEST model for site:", s, "---\n")
  set.seed(SEED)
  fit_best[[s]] <- gllvm(
    y           = pa_mat[[s]],
    X           = meta[[s]][, c("season", "substrate"), drop = FALSE],
    formula     = ~ season/substrate,
    family      = "binomial",
    num.lv      = NUM_LV,
    studyDesign = sd_site[[s]],
    row.eff     = ~(1 | clean_sample_names),
    method      = "VA",
    link        = "probit",
    control     = list(reltol = 1e-6),
    seed        = SEED
  )
  saveRDS(fit_best[[s]],
          file.path(SAVE_DIR_BASE, paste0("GLLVM_best_persite_", s, ".rds")))
  cat("  AIC:   ", round(AIC(fit_best[[s]]), 1), "\n")
  cat("  logLik:", round(as.numeric(logLik(fit_best[[s]])), 1), "\n\n")
  
  cat("--- Fitting NULL model for site:", s, "---\n")
  set.seed(SEED)
  fit_null[[s]] <- gllvm(
    y           = pa_mat[[s]],
    family      = "binomial",
    num.lv      = NUM_LV,
    studyDesign = sd_site[[s]],
    row.eff     = ~(1 | clean_sample_names),
    method      = "VA",
    link        = "probit",
    control     = list(reltol = 1e-6),
    seed        = SEED
  )
  saveRDS(fit_null[[s]],
          file.path(SAVE_DIR_BASE, paste0("GLLVM_null_persite_", s, ".rds")))
  cat("  AIC:   ", round(AIC(fit_null[[s]]), 1), "\n")
  cat("  logLik:", round(as.numeric(logLik(fit_null[[s]])), 1), "\n\n")
}

cat("All models fitted and saved.\n\n")


# ==============================================================================
# SECTION 3: COEFFICIENT EXTRACTION
#
# With ~ season/substrate each coefficient is already a direct contrast:
#   seasonWinter:substratePE   = PE vs Glass in Winter (reference season)
#   seasonSummer:substratePE   = PE vs Glass in Summer
#   seasonSpring:substrateWPE  = W-PE vs Glass in Spring
#
# No summing needed. BH correction applied across all OTU × term tests
# within each site separately.
#
# term_season and term_substrate parse the coefficient name so downstream
# filtering works the same way as before.
# is_substrate_contrast flags the plastic-vs-glass terms specifically.
# ==============================================================================

cat("=== SECTION 3: COEFFICIENT EXTRACTION ===\n\n")

make_otu_label <- function(species, genus, otu_label) {
  case_when(
    !is.na(species) & species != "" ~ species,
    !is.na(genus)   & genus   != "" ~ paste0(genus, " sp. ", otu_label),
    TRUE                            ~ otu_label
  )
}

extract_coefficients_persite <- function(model, site_name, d) {
  
  cat("  Extracting coefficients for site:", site_name, "\n")
  
  Xcoef    <- coef(model)$Xcoef
  se_Xcoef <- model$sd$Xcoef
  
  # Relative abundance weights from full dataset for cross-site comparability
  top_otus_site <- names(sort(colSums(d$otu_filt), decreasing = TRUE)[1:200])
  otu_sub       <- d$otu_filt[, top_otus_site]
  rel_abund     <- colMeans(otu_sub / rowSums(otu_sub))
  names(rel_abund) <- d$otu_id_to_label[top_otus_site]
  
  z_mat    <- Xcoef / se_Xcoef
  p_mat    <- 2 * pnorm(abs(z_mat), lower.tail = FALSE)
  
  p_BH_vec <- p.adjust(as.vector(p_mat), method = "BH")
  p_BH_mat <- matrix(p_BH_vec, nrow = nrow(p_mat), ncol = ncol(p_mat),
                     dimnames = dimnames(p_mat))
  
  cat("    Total OTU × term tests:", length(p_mat), "\n")
  cat("    Significant after BH (α = 0.05):",
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
      site             = site_name,
      lower            = Estimate - 1.96 * SE,
      upper            = Estimate + 1.96 * SE,
      sig_BH           = p_BH < 0.05,
      CI_excl_zero     = lower > 0 | upper < 0,
      #extreme_estimate = abs(Estimate) > 20,
      #reliable         = sig_BH & CI_excl_zero & !extreme_estimate,
      sig_any          = sig_BH & CI_excl_zero,
      
      # Parse season and substrate from term name.
      # Terms look like:
      #   "seasonSummer"                — season main effect on glass
      #   "seasonWinter:substratePE"    — PE vs glass in Winter
      #   "seasonSummer:substratePE"    — PE vs glass in Summer
      term_season = case_when(
        grepl("season([^:]+):substrate", Term) ~
          sub("season([^:]+):substrate.*", "\\1", Term),
        grepl("^season", Term) ~
          sub("season(.*)", "\\1", Term),
        TRUE ~ NA_character_
      ),
      term_substrate = case_when(
        grepl(":substrate", Term) ~
          sub(".*:substrate(.*)", "\\1", Term),
        TRUE ~ NA_character_
      ),
      # TRUE only for plastic-vs-glass contrast terms
      is_substrate_contrast = !is.na(term_substrate)
    ) %>%
    filter(SE > 1e-10, SE < 1e4)
  
  # Attach taxonomy
  label_to_id <- setNames(names(d$otu_id_to_label), d$otu_id_to_label)
  tax_cols    <- intersect(
    c("Domain","Supergroup","Division","Subdivision","Class",
      "Order","Family","Genus","Species"),
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
      abund_weight   = mean_rel_abund / max(mean_rel_abund, na.rm = TRUE),
      display_label  = make_otu_label(Species, Genus, OTU_Label)
    )
  
  # cat("    Extreme estimates (|coef| > 15):",
  #     sum(res$extreme_estimate, na.rm = TRUE), "\n")
  cat("    Reliable (sig + CI excl. 0):",
      sum(res$sig_any, na.rm = TRUE), "\n\n")
  res
}

results_all <- bind_rows(lapply(SITE_LEVELS, function(s) {
  extract_coefficients_persite(fit_best[[s]], s, model_data)
}))

write.csv(results_all,
          file.path(TABLE_DIR, "all_coefficients_PA_persite.csv"),
          row.names = FALSE)
cat("Full coefficient table saved.\n\n")



# ==============================================================================
# CATERPILLAR PLOTS — substrate contrast coefficients per site
# ==============================================================================

SEASON_COLOURS <- c(
  Winter  = "#4575b4",
  Spring  = "#74c476",
  Summer  = "#fd8d3c",
  Fall    = "#d73027",
  Winter2 = "#756bb1"
)

plot_caterpillar <- function(site_name, results) {
  
  dat <- results %>%
    filter(
      site == site_name,
      is_substrate_contrast,
      term_substrate %in% SUBSTRATE_PLASTIC_LEVELS,
      sig_any                              # only significant results
    ) %>%
    mutate(
      season = factor(term_season, levels = SEASON_LEVELS),
      substrate_label = factor(
        SUBSTRATE_PLASTIC_LABELS[match(term_substrate, SUBSTRATE_PLASTIC_LEVELS)],
        levels = SUBSTRATE_PLASTIC_LABELS
      )
    )
  
  if (nrow(dat) == 0) {
    cat("  No significant substrate contrasts for site:", site_name, "\n")
    return(NULL)
  }
  
  # Sort OTUs by median estimate across substrates for readability
  otu_order <- dat %>%
    group_by(display_label) %>%
    summarise(med = median(Estimate), .groups = "drop") %>%
    arrange(med) %>%
    pull(display_label)
  
  dat <- dat %>%
    mutate(display_label = factor(display_label, levels = otu_order))
  
  ggplot(dat, aes(x = Estimate, y = display_label)) +
    
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    
    geom_linerange(
      aes(xmin = lower, xmax = upper),
      linewidth = 0.5, colour = "grey40"
    ) +
    
    geom_point(
      aes(size   = abund_weight,
          colour = season)
    ) +
    
    scale_colour_manual(values = SEASON_COLOURS, name = "Season") +
    scale_size_continuous(name = "Abundance weight", range = c(1, 4)) +
    
    facet_wrap(~ substrate_label, nrow = 1, scales = "free_x") +
    
    labs(
      title    = paste("Substrate contrast coefficients —", site_name),
      subtitle = "Significant plastic vs glass contrasts within each season (BH-corrected, CI excludes zero)",
      x        = "Coefficient estimate (95% CI)",
      y        = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.y = element_line(colour = "grey92"),
      panel.grid.minor   = element_blank(),
      strip.background   = element_rect(fill = "grey90"),
      legend.position    = "right",
      axis.text.y        = element_text(size = 7, face = "italic")
    )
}

# Generate and save per site
for (s in SITE_LEVELS) {
  p <- plot_caterpillar(s, results_all)
  
  if (!is.null(p)) {
    n_otus <- length(unique(filter(results_all,
                                   site == s,
                                   sig_any,
                                   is_substrate_contrast)$display_label))
    ggsave(
      filename = file.path(FIG_DIR, paste0("caterpillar_", s, ".pdf")),
      plot     = p,
      width    = 14,
      height   = max(4, n_otus * 0.35),
      units    = "in"
    )
    print(p)
    cat("Caterpillar plot saved for site:", s, "\n")
  }
}


# Table of significant substrate contrasts (matches caterpillar plot)
sig_substrate_table <- results_all %>%
  filter(
    is_substrate_contrast,
    term_substrate %in% SUBSTRATE_PLASTIC_LEVELS,
    sig_any
  ) %>%
  mutate(
    season = factor(term_season, levels = SEASON_LEVELS),
    substrate_label = factor(
      SUBSTRATE_PLASTIC_LABELS[match(term_substrate, SUBSTRATE_PLASTIC_LEVELS)],
      levels = SUBSTRATE_PLASTIC_LABELS
    ),
    direction = ifelse(Estimate > 0, "higher_occurrence", "lower_occurrence"),
    Estimate  = round(Estimate, 3),
    SE        = round(SE, 3),
    lower     = round(lower, 3),
    upper     = round(upper, 3),
    p_BH      = signif(p_BH, 3)
  ) %>%
  arrange(site, display_label, substrate_label, season) %>%
  select(
    Site          = site,
    Taxon         = display_label,
    Class,
    Substrate     = substrate_label,
    Season        = season,
    Estimate,
    SE,
    Lower_95CI    = lower,
    Upper_95CI    = upper,
    p_BH,
    Direction     = direction
  )

write.csv(sig_substrate_table,
          file.path(TABLE_DIR, "significant_substrate_contrasts.csv"),
          row.names = FALSE)

cat("Significant substrate contrasts:\n")
print(sig_substrate_table)
cat("\nTotal rows:", nrow(sig_substrate_table), "\n")
cat("Per site:\n")
sig_substrate_table %>% count(Site) %>% print()



# Summary of direction per OTU × site
direction_summary <- sig_substrate_table %>%
  group_by(Site, Taxon, Class, Direction) %>%
  summarise(
    n_significant     = n(),
    substrates        = paste(sort(unique(as.character(Substrate))), collapse = ", "),
    seasons           = paste(sort(unique(as.character(Season))),    collapse = ", "),
    mean_estimate     = round(mean(Estimate), 3),
    .groups = "drop"
  ) %>%
  arrange(Site, Direction, desc(n_significant))

cat("=== Higher occurrence on plastic ===\n")
direction_summary %>%
  filter(Direction == "higher_occurrence") %>%
  print(n = Inf)

cat("\n=== Lower occurrence on plastic ===\n")
direction_summary %>%
  filter(Direction == "lower_occurrence") %>%
  print(n = Inf)











# ==============================================================================
# SECTION 4: VARIANCE PARTITIONING
# VP() called separately per site on best and null models
# ==============================================================================

cat("=== SECTION 4: VARIANCE PARTITIONING ===\n\n")

vp_list <- lapply(SITE_LEVELS, function(s) {
  cat("  VP for site:", s, "\n")
  VP(fit_best[[s]])
})
names(vp_list) <- SITE_LEVELS

vp_display <- bind_rows(lapply(SITE_LEVELS, function(s) {
  vp_values <- colMeans(vp_list[[s]]$PropExplainedVarSp) * 100
  data.frame(
    Site                   = s,
    Component              = names(vp_values),
    Variance_explained_pct = round(as.numeric(vp_values), 1),
    stringsAsFactors       = FALSE
  )
}))

cat("\nVariance partitioning results:\n")
print(vp_display)

write.csv(vp_display,
          file.path(TABLE_DIR, "variance_partitioning_persite.csv"),
          row.names = FALSE)

kbl(vp_display,
    booktabs  = TRUE,
    align     = c("l", "l", "r"),
    col.names = c("Site", "Component", "Variance explained (%)")) %>%
  kable_classic(full_width = FALSE, html_font = "Arial") %>%
  collapse_rows(columns = 1, valign = "top")


# ==============================================================================
# SECTION 5: SUBSTRATE RESPONSE ANALYSIS
#
# direction_long is built directly from is_substrate_contrast == TRUE terms.
# Each row is already one plastic-vs-glass contrast — no summing needed.
# Everything from the response matrix onward is identical to the original.
# ==============================================================================

cat("=== SECTION 5: SUBSTRATE RESPONSE ANALYSIS ===\n\n")

# ── Step 1: Direction per OTU × site × season × substrate ─────────────────────

direction_long <- results_all %>%
  filter(
    is_substrate_contrast,
    term_substrate %in% SUBSTRATE_PLASTIC_LEVELS,
    !extreme_estimate
  ) %>%
  mutate(
    season = factor(term_season, levels = SEASON_LEVELS),
    substrate_label = factor(
      SUBSTRATE_PLASTIC_LABELS[match(term_substrate, SUBSTRATE_PLASTIC_LEVELS)],
      levels = SUBSTRATE_PLASTIC_LABELS
    ),
    direction = case_when(
      sig_any & Estimate > 0 ~ "higher_occurrence",
      sig_any & Estimate < 0 ~ "lower_occurrence",
      TRUE                    ~ "ns"
    )
  ) %>%
  select(site, season, OTU_Label, display_label, Class,
         substrate_label, direction, Estimate, abund_weight)

cat("Direction long — rows:", nrow(direction_long), "\n")
cat("Significant (any direction):", sum(direction_long$direction != "ns"), "\n\n")


# ── Step 2: Substrate response matrix ─────────────────────────────────────────

response_matrix <- direction_long %>%
  pivot_wider(
    id_cols     = c(site, season, OTU_Label, display_label, Class, abund_weight),
    names_from  = substrate_label,
    values_from = direction,
    values_fill = "ns"
  ) %>%
  filter(
    rowSums(across(all_of(SUBSTRATE_PLASTIC_LABELS), ~ . != "ns")) > 0
  ) %>%
  mutate(
    pattern = paste(PE, wPE, PET, wPET, sep = "/")
  ) %>%
  arrange(site, season, Class, display_label)

cat("OTU × season rows with ≥1 significant result:\n")
print(table(response_matrix$site, response_matrix$season))

write.csv(response_matrix,
          file.path(TABLE_DIR, "substrate_response_matrix_persite.csv"),
          row.names = FALSE)
write.csv(filter(response_matrix, site == "SELVA"),
          file.path(TABLE_DIR, "substrate_response_matrix_SELVA.csv"),
          row.names = FALSE)
write.csv(filter(response_matrix, site == "TBS"),
          file.path(TABLE_DIR, "substrate_response_matrix_TBS.csv"),
          row.names = FALSE)

cat("\nMost common patterns (top 15 per site):\n")
response_matrix %>%
  count(site, pattern, sort = TRUE) %>%
  group_by(site) %>%
  slice_max(n, n = 15) %>%
  print(n = 30)


# ── Step 3: Substrate summary per OTU × site × substrate across seasons ────────

substrate_summary <- direction_long %>%
  group_by(site, OTU_Label, display_label, Class,
           abund_weight, substrate_label) %>%
  summarise(
    ever_higher         = any(direction == "higher_occurrence"),
    ever_lower          = any(direction == "lower_occurrence"),
    n_higher            = sum(direction == "higher_occurrence"),
    n_lower             = sum(direction == "lower_occurrence"),
    n_sig               = sum(direction != "ns"),
    n_conditions_tested = n(),
    prop_sig            = round(n_sig / n_conditions_tested, 2),
    mean_estimate       = round(mean(Estimate[direction != "ns"], na.rm = TRUE), 2),
    seasons_higher      = paste(sort(as.character(season[direction == "higher_occurrence"])),
                                collapse = ", "),
    seasons_lower       = paste(sort(as.character(season[direction == "lower_occurrence"])),
                                collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    is_mixed  = ever_higher & ever_lower,
    is_higher = ever_higher & !ever_lower,
    is_lower  = ever_lower  & !ever_higher,
    consistency = case_when(
      n_sig >= 3 ~ "consistent",
      n_sig == 2 ~ "moderate",
      n_sig == 1 ~ "occasional",
      TRUE       ~ "none"
    )
  )


# ── Step 4: Count enriched/avoided substrates per OTU × site ──────────────────

otu_site_counts <- substrate_summary %>%
  group_by(site, OTU_Label, display_label, Class, abund_weight) %>%
  summarise(
    n_higher_substrates = sum(is_higher),
    n_lower_substrates  = sum(is_lower),
    higher_substrates   = paste(substrate_label[is_higher], collapse = ", "),
    lower_substrates    = paste(substrate_label[is_lower],  collapse = ", "),
    .groups = "drop"
  )


# ── Step 5: Group classification ──────────────────────────────────────────────

# Plastic generalist — higher occurrence on all 4 plastic substrates
group_plastic_generalist <- substrate_summary %>%
  filter(is_higher) %>%
  left_join(select(otu_site_counts, site, OTU_Label,
                   n_higher_substrates, higher_substrates),
            by = c("site", "OTU_Label")) %>%
  filter(n_higher_substrates == 4) %>%
  mutate(group = "Plastic generalist") %>%
  select(site, OTU_Label, display_label, Class, abund_weight,
         substrate_label, group, consistency,
         n_higher, n_lower, n_sig,
         seasons_higher, seasons_lower, mean_estimate)

# Plastic avoider — lower occurrence on all 4 plastic substrates
group_plastic_avoider <- substrate_summary %>%
  filter(is_lower) %>%
  left_join(select(otu_site_counts, site, OTU_Label,
                   n_lower_substrates, lower_substrates),
            by = c("site", "OTU_Label")) %>%
  filter(n_lower_substrates == 4) %>%
  mutate(group = "Plastic avoider") %>%
  select(site, OTU_Label, display_label, Class, abund_weight,
         substrate_label, group, consistency,
         n_higher, n_lower, n_sig,
         seasons_higher, seasons_lower, mean_estimate)

# Multi-substrate — higher occurrence on 2-3 substrates
group_multi_higher <- substrate_summary %>%
  filter(is_higher) %>%
  left_join(select(otu_site_counts, site, OTU_Label,
                   n_higher_substrates, higher_substrates),
            by = c("site", "OTU_Label")) %>%
  filter(n_higher_substrates %in% c(2, 3)) %>%
  mutate(group = paste0("Multi-substrate (", n_higher_substrates, " substrates)")) %>%
  select(site, OTU_Label, display_label, Class, abund_weight,
         substrate_label, group, consistency,
         n_higher, n_lower, n_sig,
         seasons_higher, seasons_lower, mean_estimate)

all_groups <- bind_rows(
  group_plastic_generalist,
  group_plastic_avoider,
  group_multi_higher
)

write.csv(all_groups,
          file.path(TABLE_DIR, "substrate_response_groups_persite.csv"),
          row.names = FALSE)
write.csv(substrate_summary,
          file.path(TABLE_DIR, "substrate_summary_per_OTU_site_substrate.csv"),
          row.names = FALSE)
cat("Classification tables saved.\n\n")


# ── Step 6: Results summaries ─────────────────────────────────────────────────

cat("=== Plastic generalists ===\n\n")

cat("Count per site:\n")
group_plastic_generalist %>%
  distinct(site, OTU_Label) %>%
  count(site) %>%
  print()

cat("\nTaxa and consistency:\n")
group_plastic_generalist %>%
  distinct(site, OTU_Label, display_label, Class) %>%
  arrange(site, Class, display_label) %>%
  print()

cat("\nSeasons of higher occurrence per OTU × substrate:\n")
group_plastic_generalist %>%
  select(site, OTU_Label, display_label, substrate_label,
         n_higher, n_sig, consistency, seasons_higher) %>%
  arrange(site, display_label, substrate_label) %>%
  print(n = Inf)


cat("=== Plastic avoiders ===\n\n")

cat("Count per site:\n")
group_plastic_avoider %>%
  distinct(site, OTU_Label) %>%
  count(site) %>%
  print()

cat("\nTaxa and consistency:\n")
group_plastic_avoider %>%
  distinct(site, OTU_Label, display_label, Class, consistency) %>%
  arrange(site, Class, display_label) %>%
  print()

cat("\nSubstrates with lower occurrence per OTU:\n")
group_plastic_avoider %>%
  select(site, OTU_Label, display_label, substrate_label,
         n_lower, n_sig, consistency, seasons_lower) %>%
  arrange(site, display_label, substrate_label) %>%
  print(n = Inf)


cat("=== Multi-substrate higher occurrence (2-3 substrates) ===\n\n")

cat("Count per site and group:\n")
group_multi_higher %>%
  distinct(site, OTU_Label, group) %>%
  count(site, group) %>%
  arrange(site, group) %>%
  print()

cat("\nTaxa:\n")
group_multi_higher %>%
  distinct(site, OTU_Label, display_label, Class, group, consistency) %>%
  arrange(site, group, Class) %>%
  print()

cat("\nSubstrates with higher occurrence per OTU:\n")
group_multi_higher %>%
  select(site, OTU_Label, display_label, substrate_label,
         n_higher, n_sig, consistency, seasons_higher) %>%
  arrange(site, display_label, substrate_label) %>%
  print(n = Inf)


cat("=== Consistency breakdown across all groups ===\n\n")
all_groups %>%
  distinct(site, OTU_Label, substrate_label, group, consistency) %>%
  count(group, consistency) %>%
  arrange(group, consistency) %>%
  print()


cat("=== Cross-site responses ===\n\n")

cross_site <- all_groups %>%
  distinct(site, OTU_Label, display_label, Class, group) %>%
  group_by(OTU_Label, display_label, Class, group) %>%
  filter(n_distinct(site) > 1) %>%
  ungroup() %>%
  arrange(group, Class)

cat("OTUs with consistent group membership at both sites:\n")
cross_site %>%
  distinct(OTU_Label, group) %>%
  count(group) %>%
  print()

cat("\nFull list:\n")
all_groups %>%
  distinct(site, OTU_Label, display_label, Class, group, consistency) %>%
  semi_join(cross_site, by = c("OTU_Label", "group")) %>%
  arrange(group, Class, OTU_Label, site) %>%
  print(n = Inf)


# ── Step 7: Plastic generalists table (Word export) ───────────────────────────

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

ft <- flextable(generalist_wide) %>%
  italic(j = "Taxon") %>%
  bold(part = "header") %>%
  bg(part = "header", bg = "#D5E8F0") %>%
  bg(i = seq(2, nrow(generalist_wide), by = 2), bg = "#F5F5F5", part = "body") %>%
  width(j = "Taxon", width = 1.6) %>%
  width(j = "Class", width = 1.4) %>%
  width(j = "Site",  width = 0.5) %>%
  width(j = c("PE", "wPE", "PET", "wPET"), width = 0.75) %>%
  align(align = "left", part = "all") %>%
  set_table_properties(layout = "fixed") %>%
  fontsize(size = 9, part = "all") %>%
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

doc <- read_docx() %>% body_add_flextable(ft)
print(doc, target = file.path(TABLE_DIR, "plastic_generalists_table_persite.docx"))

cat("\nScript complete.\n")