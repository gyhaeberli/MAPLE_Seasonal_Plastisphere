################################################################################
# GLLVM PRESENCE-ABSENCE PIPELINE — SEASONAL PLASTISPHERE  (FINAL)
#
# Per-site GLLVMs, formula ~ season/substrate: each coefficient is one plastic
# vs glass within a season (the season main effect cancels in the glass-
# referenced contrast, so no coefficient summing / covariance propagation).
#
# Separation guard (the key correctness fix):
#   - SE < SE_DEGENERATE_FLOOR (or NA)  -> separation artifact (SE collapsed to
#     ~0, giving fake p ~ 0 and a zero-width CI). EXCLUDED.
#   - |coef| > COEF_MAX                 -> absurd-value sanity cap only.
#   Separation rows are dropped from the BH family and from `reliable`. Large
#   coefficients with HEALTHY SEs (e.g. probit ~10, SE ~0.7) are genuine strong
#   effects and are KEPT — this is what restored the Syndiniales wPET-Summer hit.
#
# Reporting:
#   - Each reliable contrast is annotated with raw biological-sample prevalence
#     ("k/n plastic vs k/n glass" for that substrate x season).
#   - Section 6 is a model-free occurrence cross-check (ANY plastic vs glass),
#     run BOTH per-season and SEASON-POOLED (the pooled one has real power),
#     each reporting raw prevalence fractions.
#
# Replicate non-independence (4 PCR reps / sample) is handled in the model by the
# random row effect (1 | clean_sample_names) — confirmed working (row-effect
# SD ~0.7-0.9). Prevalence tables collapse reps to biological samples (present
# if detected in >= 1 replicate).
################################################################################

library(gllvm)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

cat("=== GLLVM PER-SITE ANALYSIS ===\n\n")

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

SE_DEGENERATE_FLOOR  <- 0.01   # SE below this (or NA) = separation artifact
COEF_MAX             <- 30     # |coef| above this = absurd sanity cap only
MIN_EFFECT           <- 0      # optional min |coef|; 0 = off (SE floor suffices)
MIN_SAMPLES_PER_CELL <- 3      # min biological samples per side for a Fisher test
DETECT_RULE          <- "any"  # taxon present in a sample if in >=1 PCR replicate


# ==============================================================================
# SECTION 1: LOAD DATA AND SPLIT BY SITE
# ==============================================================================

cat("=== SECTION 1: LOAD DATA ===\n\n")

model_data      <- readRDS(file.path(SAVE_DIR_BASE, "data_for_gllvm_FINAL.rds"))
otu_filt        <- model_data$otu_filt
tax_mat         <- model_data$tax_mat
metadata        <- model_data$metadata
otu_id_to_label <- model_data$otu_id_to_label
pa_mat_full     <- model_data$pa_mat

cat("Full PA matrix:", nrow(pa_mat_full), "x", ncol(pa_mat_full), "\n\n")

# Pin reference levels: Glass (substrate) and Winter (season) — every contrast
# depends on these being the baselines.
metadata$substrate <- factor(metadata$substrate, levels = SUBSTRATE_LEVELS)
metadata$season    <- factor(metadata$season,    levels = SEASON_LEVELS)
stopifnot(levels(metadata$substrate)[1] == "Glass",
          levels(metadata$season)[1]    == "Winter")

pa_mat  <- list(); meta <- list(); sd_site <- list()
for (s in SITE_LEVELS) {
  idx         <- which(metadata$site == s)
  pa_mat[[s]] <- pa_mat_full[idx, ]
  meta[[s]]   <- metadata[idx, ]
  meta[[s]]$season    <- droplevels(meta[[s]]$season)
  meta[[s]]$substrate <- droplevels(meta[[s]]$substrate)
  sd_site[[s]] <- data.frame(clean_sample_names = meta[[s]]$clean_sample_names)
  cat(" ", s, ":", nrow(pa_mat[[s]]), "rows |",
      nlevels(meta[[s]]$season), "seasons |", nlevels(meta[[s]]$substrate), "substrates\n")
}
cat("\n")


# ==============================================================================
# SECTION 1b: BIOLOGICAL-SAMPLE PREVALENCE (model-free, computed ONCE)
#
# Collapse the 4 PCR replicates per sample to one presence/absence per taxon,
# then count present (k) and total (n) biological samples per
# site x season x substrate x taxon. Reused for (a) annotating the GLLVM table
# with raw prevalence and (b) the Section 6 occurrence tests.
# ==============================================================================

cat("=== SECTION 1b: BIOLOGICAL-SAMPLE PREVALENCE ===\n\n")

# check.names = FALSE keeps the original (non-syntactic) OTU labels intact.
pa_df <- as.data.frame(pa_mat_full, check.names = FALSE)
pa_df$clean_sample_names <- metadata$clean_sample_names
pa_df$site      <- as.character(metadata$site)
pa_df$season    <- as.character(metadata$season)
pa_df$substrate <- as.character(metadata$substrate)

sample_pa <- pa_df %>%
  pivot_longer(cols = -c(clean_sample_names, site, season, substrate),
               names_to = "OTU_Label", values_to = "present_rep") %>%
  group_by(site, season, substrate, clean_sample_names, OTU_Label) %>%
  summarise(present = if (DETECT_RULE == "any") as.integer(any(present_rep > 0))
                      else as.integer(mean(present_rep > 0) >= 0.5),
            .groups = "drop")

# n = biological samples in the cell, k = how many had the taxon
cell_prev <- sample_pa %>%
  group_by(site, season, substrate, OTU_Label) %>%
  summarise(n = n(), k = sum(present), .groups = "drop")

glass_prev <- cell_prev %>%
  filter(substrate == "Glass") %>%
  select(site, season, OTU_Label, n_glass = n, k_glass = k)

# Per-specific-plastic prevalence vs glass (for annotating the GLLVM table)
prev_perplastic <- cell_prev %>%
  filter(substrate %in% SUBSTRATE_PLASTIC_LEVELS) %>%
  rename(n_plastic = n, k_plastic = k) %>%
  left_join(glass_prev, by = c("site", "season", "OTU_Label")) %>%
  mutate(Prevalence = sprintf("%d/%d vs %d/%d",
                              k_plastic, n_plastic,
                              coalesce(k_glass, 0L), coalesce(n_glass, 0L)))

cat("Biological samples per site x season x substrate (median):",
    median(cell_prev$n), "\n\n")


# ==============================================================================
# SECTION 2: MODEL FITTING
#
# ~ season/substrate == ~ season + season:substrate. One BEST model per site.
# (No separate null model is fitted — VP() works directly on the fitted model,
# and the null was unused in the previous pipeline.)
# ==============================================================================

cat("=== SECTION 2: MODEL FITTING ===\n\n")

fit_best <- list()

for (s in SITE_LEVELS) {
  cat("--- Fitting model for site:", s, "---\n")
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
  saveRDS(fit_best[[s]], file.path(SAVE_DIR_BASE, paste0("GLLVM_best_persite_", s, ".rds")))
  cat("  AIC:", round(AIC(fit_best[[s]]), 1),
      "| logLik:", round(as.numeric(logLik(fit_best[[s]])), 1),
      "| row-effect SD:", round(as.numeric(fit_best[[s]]$params$sigma), 3), "\n\n")
}


#load models
fit_best[["SELVA"]] <- readRDS(file.path(SAVE_DIR_BASE, "GLLVM_best_persite_SELVA.rds"))
fit_best[["TBS"]]   <- readRDS(file.path(SAVE_DIR_BASE, "GLLVM_best_persite_TBS.rds"))

# ==============================================================================
# SECTION 3: COEFFICIENT EXTRACTION (+ separation guard) AND RELIABLE TABLE
# ==============================================================================

cat("=== SECTION 3: COEFFICIENT EXTRACTION ===\n\n")

make_otu_label <- function(species, genus, otu_label) {
  case_when(!is.na(species) & species != "" ~ species,
            !is.na(genus)   & genus   != "" ~ paste0(genus, " sp. ", otu_label),
            TRUE                            ~ otu_label)
}

extract_coefficients_persite <- function(model, site_name, d) {
  cat("  Site:", site_name, "\n")
  Xcoef    <- coef(model)$Xcoef
  se_Xcoef <- model$sd$Xcoef

  # abundance weights (for caterpillar point size)
  top_otus  <- names(sort(colSums(d$otu_filt), decreasing = TRUE)[1:200])
  rel_abund <- colMeans(d$otu_filt[, top_otus] / rowSums(d$otu_filt[, top_otus]))
  names(rel_abund) <- d$otu_id_to_label[top_otus]

  res <- do.call(rbind, lapply(colnames(Xcoef), function(term)
    data.frame(OTU_Label = rownames(Xcoef), Term = term,
               Estimate = Xcoef[, term], SE = se_Xcoef[, term],
               stringsAsFactors = FALSE))) %>%
    mutate(
      z     = Estimate / SE,
      p_raw = 2 * pnorm(abs(z), lower.tail = FALSE),
      # separation flags
      is_degenerate_SE = is.na(SE) | SE < SE_DEGENERATE_FLOOR,
      is_extreme_coef  = abs(Estimate) > COEF_MAX,
      is_separation    = is_degenerate_SE | is_extreme_coef,
      lower = Estimate - 1.96 * SE,
      upper = Estimate + 1.96 * SE,
      site  = site_name,
      term_season = case_when(
        grepl("season([^:]+):substrate", Term) ~ sub("season([^:]+):substrate.*", "\\1", Term),
        grepl("^season", Term)                 ~ sub("season(.*)", "\\1", Term),
        TRUE ~ NA_character_),
      term_substrate = case_when(
        grepl(":substrate", Term) ~ sub(".*:substrate(.*)", "\\1", Term),
        TRUE ~ NA_character_),
      is_substrate_contrast = !is.na(term_substrate)
    ) %>%
    filter(is.finite(Estimate), is.finite(SE) | is.na(SE))

  # BH over trustworthy (non-separation) tests only
  res$p_BH <- NA_real_
  ok <- !res$is_separation & is.finite(res$p_raw)
  res$p_BH[ok] <- p.adjust(res$p_raw[ok], method = "BH")
  res <- res %>%
    mutate(sig_BH   = !is.na(p_BH) & p_BH < 0.05,
           reliable = sig_BH & !is_separation & (abs(Estimate) >= MIN_EFFECT))

  cat("    separation flagged:", sum(res$is_separation),
      "(SE-degenerate", sum(res$is_degenerate_SE), "| extreme", sum(res$is_extreme_coef), ")",
      "| reliable substrate contrasts:",
      sum(res$reliable & res$is_substrate_contrast, na.rm = TRUE), "\n")

  # taxonomy
  label_to_id <- setNames(names(d$otu_id_to_label), d$otu_id_to_label)
  tax_cols <- intersect(c("Domain","Supergroup","Division","Subdivision",
                          "Class","Order","Family","Genus","Species"), colnames(d$tax_mat))
  tax_df <- d$tax_mat[label_to_id[unique(res$OTU_Label)], tax_cols, drop = FALSE] %>%
    as.data.frame() %>% rownames_to_column("OTU_ID") %>%
    mutate(OTU_Label = d$otu_id_to_label[OTU_ID])

  res %>%
    left_join(select(tax_df, -OTU_ID), by = "OTU_Label") %>%
    mutate(mean_rel_abund = rel_abund[OTU_Label],
           abund_weight   = mean_rel_abund / max(mean_rel_abund, na.rm = TRUE),
           display_label  = make_otu_label(Species, Genus, OTU_Label))
}

results_all <- bind_rows(lapply(SITE_LEVELS, function(s)
  extract_coefficients_persite(fit_best[[s]], s, model_data)))
write.csv(results_all, file.path(TABLE_DIR, "all_coefficients_PA_persite.csv"), row.names = FALSE)
cat("\n")

# ── Reliable contrast table, annotated with raw biological-sample prevalence ──
sig_substrate_table <- results_all %>%
  filter(is_substrate_contrast, term_substrate %in% SUBSTRATE_PLASTIC_LEVELS, reliable) %>%
  left_join(select(prev_perplastic, site, season, substrate, OTU_Label, Prevalence),
            by = c("site", "term_season" = "season", "term_substrate" = "substrate", "OTU_Label")) %>%
  mutate(
    Season    = factor(term_season, levels = SEASON_LEVELS),
    Substrate = factor(SUBSTRATE_PLASTIC_LABELS[match(term_substrate, SUBSTRATE_PLASTIC_LEVELS)],
                       levels = SUBSTRATE_PLASTIC_LABELS),
    Direction = ifelse(Estimate > 0, "higher", "lower"),
    Estimate = round(Estimate, 3), SE = round(SE, 3),
    Lower_95CI = round(lower, 3), Upper_95CI = round(upper, 3),
    p_BH = signif(p_BH, 3)) %>%
  arrange(site, display_label, Substrate, Season) %>%
  select(Site = site, Taxon = display_label, Class, Substrate, Season,
         Estimate, SE, Lower_95CI, Upper_95CI, p_BH, Direction,
         Prevalence)   # "k/n plastic vs k/n glass" (biological samples, that season)

write.csv(sig_substrate_table, file.path(TABLE_DIR, "significant_substrate_contrasts.csv"), row.names = FALSE)
cat("Reliable substrate contrasts (with raw prevalence):", nrow(sig_substrate_table), "\n")
print(as_tibble(sig_substrate_table), n = Inf)
cat("\n")


# ==============================================================================
# CATERPILLAR PLOT — reliable contrasts per site
# ==============================================================================

SEASON_COLOURS <- c(Winter="#4575b4", Spring="#74c476", Summer="#fd8d3c",
                    Fall="#d73027", Winter2="#756bb1")

for (s in SITE_LEVELS) {
  dat <- results_all %>%
    filter(site == s, is_substrate_contrast,
           term_substrate %in% SUBSTRATE_PLASTIC_LEVELS, reliable) %>%
    mutate(season = factor(term_season, levels = SEASON_LEVELS),
           substrate_label = factor(SUBSTRATE_PLASTIC_LABELS[match(term_substrate, SUBSTRATE_PLASTIC_LEVELS)],
                                    levels = SUBSTRATE_PLASTIC_LABELS))
  if (nrow(dat) == 0) { cat("No reliable contrasts to plot for", s, "\n"); next }
  ord <- dat %>% group_by(display_label) %>% summarise(m = median(Estimate)) %>%
    arrange(m) %>% pull(display_label)
  dat$display_label <- factor(dat$display_label, levels = ord)
  p <- ggplot(dat, aes(Estimate, display_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = lower, xmax = upper), colour = "grey40") +
    geom_point(aes(size = abund_weight, colour = season)) +
    scale_colour_manual(values = SEASON_COLOURS, name = "Season") +
    scale_size_continuous(name = "Abundance", range = c(1, 4)) +
    facet_wrap(~ substrate_label, nrow = 1, scales = "free_x") +
    labs(title = paste("Plastic vs glass (reliable) -", s),
         x = "Coefficient (95% CI)", y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 7, face = "italic"),
          panel.grid.minor = element_blank())
  ggsave(file.path(FIG_DIR, paste0("caterpillar_", s, ".pdf")),
         p, width = 14, height = max(4, length(ord) * 0.35), units = "in")
  cat("Caterpillar saved for", s, "\n")
}
cat("\n")


# ==============================================================================
# SECTION 4: VARIANCE PARTITIONING  (clean component labels)
#
# CAVEAT: the "Season:Substrate" component is the 20-column interaction block,
# not a clean substrate main effect — do not compare directly to the global
# PERMANOVA partition.
# ==============================================================================

cat("=== SECTION 4: VARIANCE PARTITIONING ===\n\n")

relabel_vp <- function(x)
  ifelse(grepl("substrate", x), "Season:Substrate",
  ifelse(grepl("Row random", x), "Sample (random row)",
  ifelse(grepl("^LV", x), x,
  ifelse(grepl("season", x), "Season", x))))

vp_display <- bind_rows(lapply(SITE_LEVELS, function(s) {
  v <- colMeans(VP(fit_best[[s]])$PropExplainedVarSp) * 100
  data.frame(Site = s, Component = relabel_vp(names(v)),
             Variance_explained_pct = round(as.numeric(v), 1))
}))
print(vp_display)
write.csv(vp_display, file.path(TABLE_DIR, "variance_partitioning_persite.csv"), row.names = FALSE)
cat("\n")


# ==============================================================================
# SECTION 5: PER-TAXON PLASTIC RESPONSE (rolled up from reliable contrasts)
# ==============================================================================

cat("=== SECTION 5: PER-TAXON PLASTIC RESPONSE ===\n\n")

direction_long <- results_all %>%
  filter(is_substrate_contrast, term_substrate %in% SUBSTRATE_PLASTIC_LEVELS) %>%
  mutate(reliable = coalesce(reliable, FALSE),
         season = factor(term_season, levels = SEASON_LEVELS),
         substrate_label = SUBSTRATE_PLASTIC_LABELS[match(term_substrate, SUBSTRATE_PLASTIC_LEVELS)],
         direction = case_when(reliable & Estimate > 0 ~ "higher",
                               reliable & Estimate < 0 ~ "lower",
                               TRUE ~ "ns"))

plastic_response <- direction_long %>%
  group_by(site, OTU_Label, display_label, Class, substrate_label) %>%
  summarise(ever_higher = any(direction == "higher"),
            ever_lower  = any(direction == "lower"), .groups = "drop") %>%
  mutate(is_higher = ever_higher & !ever_lower,
         is_lower  = ever_lower  & !ever_higher) %>%
  group_by(site, OTU_Label, display_label, Class) %>%
  summarise(n_higher_substrates = sum(is_higher),
            n_lower_substrates  = sum(is_lower),
            higher_substrates = paste(sort(unique(substrate_label[is_higher])), collapse = ", "),
            lower_substrates  = paste(sort(unique(substrate_label[is_lower])),  collapse = ", "),
            .groups = "drop") %>%
  filter(n_higher_substrates > 0 | n_lower_substrates > 0) %>%
  mutate(response = case_when(
    n_higher_substrates == 4 ~ "Plastic generalist (4/4 higher)",
    n_lower_substrates  == 4 ~ "Plastic avoider (4/4 lower)",
    n_higher_substrates >= 2 ~ paste0("Multi-plastic higher (", n_higher_substrates, ")"),
    n_higher_substrates == 1 & n_lower_substrates == 0 ~ "Single-plastic higher",
    n_lower_substrates  >= 1 & n_higher_substrates == 0 ~ "Plastic lower",
    TRUE ~ "Mixed")) %>%
  arrange(site, desc(n_higher_substrates), Class, display_label)

write.csv(plastic_response, file.path(TABLE_DIR, "plastic_response_per_taxon.csv"), row.names = FALSE)
cat("Per-taxon response summary:\n"); print(plastic_response, n = Inf)
cat("\nResponse-type counts:\n"); print(count(plastic_response, site, response))
cat("\n")



#with season
plastic_response_season <- direction_long %>%
  group_by(site, OTU_Label, display_label, Class, substrate_label) %>%
  summarise(ever_higher = any(direction == "higher"),
            ever_lower  = any(direction == "lower"),
            higher_seasons = paste(sort(unique(as.character(season[direction == "higher"]))), collapse = ", "),
            lower_seasons  = paste(sort(unique(as.character(season[direction == "lower"]))),  collapse = ", "),
            .groups = "drop") %>%
  mutate(is_higher = ever_higher & !ever_lower,
         is_lower  = ever_lower  & !ever_higher) %>%
  group_by(site, OTU_Label, display_label, Class) %>%
  summarise(n_higher_substrates = sum(is_higher),
            n_lower_substrates  = sum(is_lower),
            higher_substrates = paste(sort(unique(substrate_label[is_higher])), collapse = ", "),
            lower_substrates  = paste(sort(unique(substrate_label[is_lower])),  collapse = ", "),
            higher_seasons = paste(sort(unique(unlist(strsplit(higher_seasons[is_higher], ", ")))), collapse = ", "),
            lower_seasons  = paste(sort(unique(unlist(strsplit(lower_seasons[is_lower],  ", ")))), collapse = ", "),
            .groups = "drop") %>%
  filter(n_higher_substrates > 0 | n_lower_substrates > 0) %>%
  mutate(response = case_when(
    n_higher_substrates == 4 ~ "Plastic generalist (4/4 higher)",
    n_lower_substrates  == 4 ~ "Plastic avoider (4/4 lower)",
    n_higher_substrates >= 2 ~ paste0("Multi-plastic higher (", n_higher_substrates, ")"),
    n_higher_substrates == 1 & n_lower_substrates == 0 ~ "Single-plastic higher",
    n_lower_substrates  >= 1 & n_higher_substrates == 0 ~ "Plastic lower",
    TRUE ~ "Mixed")) %>%
  arrange(site, desc(n_higher_substrates), Class, display_label)

print(plastic_response_season)



# ==============================================================================
# PUBLICATION TABLE — LaTeX-style Word export
# ==============================================================================

library(officer)
library(flextable)

# ── 1. Build the two sections from sig_substrate_table ────────────────────────

substrate_levels <- c("PE", "wPE", "PET", "wPET")

df <- sig_substrate_table %>%
  mutate(
    substrate = as.character(Substrate),
    season    = as.character(Season),
    value     = paste0(season, " (", Prevalence, ")")
  ) %>%
  select(display_label = Taxon, Class, substrate, Direction, value)

make_wide <- function(direction) {
  df %>%
    filter(Direction == direction) %>%
    group_by(display_label, Class, substrate) %>%
    summarise(value = paste(value, collapse = "; "), .groups = "drop") %>%
    mutate(substrate = factor(substrate, levels = substrate_levels)) %>%
    pivot_wider(
      names_from   = substrate,
      values_from  = value,
      values_fill  = "—",
      names_expand = TRUE
    ) %>%
    arrange(Class, display_label)
}

higher_wide <- make_wide("higher")
lower_wide  <- make_wide("lower")

# ── 2. Tag rows with section, then bind ───────────────────────────────────────

higher_wide$section <- "Elevated on plastic"
lower_wide$section  <- "Reduced on plastic"

combined <- bind_rows(higher_wide, lower_wide)

# Row indices for each section (used for border placement)
n_higher   <- nrow(higher_wide)
n_lower    <- nrow(lower_wide)
last_row   <- nrow(combined)

# ── 3. Define borders ─────────────────────────────────────────────────────────

thin_rule  <- fp_border(color = "#444444", width = 0.75)
thick_rule <- fp_border(color = "#000000", width = 2.5)
no_border  <- fp_border(width = 0)

# ── 4. Build flextable (section column used only for grouping, then hidden) ───

ft <- flextable(combined %>% select(-section)) %>%
  
  # Column labels
  set_header_labels(
    display_label = "Taxon",
    Class         = "Class",
    PE            = "PE",
    wPE           = "wPE",
    PET           = "PET",
    wPET          = "wPET"
  ) %>%
  
  # ── Section header rows ────────────────────────────────────────────────────
  # Prepend a merged label row for each section using add_header_row is not
  # ideal here; instead we style row 1 and row (n_higher+1) as section banners
  # by inserting blank rows above them via compose().
  # Simpler and more robust: use bg + bold on the first row of each group,
  # merged across all columns with a labelless style.
  
  # Section banners: row 1 = "Elevated", row n_higher+1 = "Reduced"
  # We inject these as real data rows so merge works cleanly.
  
  # ── Typography ────────────────────────────────────────────────────────────
  font(fontname = "Arial", part = "all") %>%
  fontsize(size = 9,  part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  italic(j = "display_label", part = "body") %>%
  bold(part = "header") %>%
  
  # ── Column widths (inches) — landscape A4 content width ~10.1 in ─────────────
  width(j = "display_label",                width = 2.2) %>%
  width(j = "Class",                        width = 1.4) %>%
  width(j = c("PE", "wPE", "PET", "wPET"), width = 1.6) %>%
  
  # ── Outer frame: thin top, thin bottom ───────────────────────────────────
  hline_top(part = "header", border = thin_rule) %>%
  hline_top(part = "body",   border = thin_rule) %>%
  hline(i = last_row, part = "body", border = thin_rule) %>%
  
  # ── Inner body lines: no horizontal rules except the section divider ──────
  border_remove() %>%                             # wipe all inner borders first
  hline_top(part = "header", border = thin_rule) %>%  # restore outer frame
  hline_top(part = "body",   border = thin_rule) %>%
  hline(i = last_row, part = "body", border = thin_rule) %>%
  
  # ── Section divider: thick double-hline after last higher row ─────────────
  hline(i = n_higher, part = "body", border = thick_rule) %>%
  
  # ── Vertical lines: none (LaTeX booktabs style) ───────────────────────────
  
  # ── Row backgrounds: plain white throughout ───────────────────────────────
  bg(bg = "white", part = "all") %>%
  
  # ── Alignment ─────────────────────────────────────────────────────────────
  align(j = c("PE", "wPE", "PET", "wPET"), align = "left",   part = "body") %>%
  align(j = c("display_label", "Class"),   align = "left",   part = "body") %>%
  align(part = "header",                   align = "center") %>%
  
  # ── Wrap & layout ─────────────────────────────────────────────────────────
  set_table_properties(layout = "fixed") %>%
  valign(valign = "top", part = "body") %>%
  
  # ── Section label: bold italic in first column of each section ────────────
  bold(i  = seq_len(n_higher),              j = "display_label", part = "body") %>%
  italic(i = seq_len(n_higher),             j = "display_label", part = "body") %>%
  bold(i  = seq(n_higher + 1, last_row),    j = "display_label", part = "body") %>%
  italic(i = seq(n_higher + 1, last_row),   j = "display_label", part = "body")

# ── 5. Add section label annotations above each block via a header footnote ──
# The cleanest LaTeX-equivalent in Word is a bold paragraph before the table.

# ── 6. Write to Word ──────────────────────────────────────────────────────────
out_path <- file.path(TABLE_DIR, "publication_table.docx")

read_docx() %>%
  body_add_par("Table X. Reliable plastic vs glass contrasts", style = "heading 1") %>%
  body_add_par(
    paste0("Each cell lists Season (k_plastic/n_plastic vs k_glass/n_glass biological samples). ",
           "Top section: taxa with significantly higher occurrence on plastic. ",
           "Bottom section: taxa with significantly lower occurrence on plastic. ",
           "Thick rule separates sections. Corrected p < 0.05 (BH)."),
    style = "Normal"
  ) %>%
  body_add_par("Elevated on plastic", style = "Normal") %>%
  body_add_flextable(ft) %>%
  body_end_section_landscape() %>%    # <-- sets this section to landscape
  print(target = out_path)
cat("Publication table saved to:", out_path, "\n")



# ==============================================================================
# SECTION 6: OCCURRENCE CROSS-CHECK — ANY PLASTIC vs GLASS (model-free)
#
# Pools the four plastics into one "plastic" group and tests plastic vs glass at
# the biological-sample level with Fisher's exact test, reporting raw prevalence.
# Run two ways:
#   (a) PER SEASON  — keeps seasonal context, but n_glass ~ 3 per season caps the
#       achievable p, so it is descriptive only (won't clear BH).
#   (b) SEASON-POOLED — pools seasons too (n_glass ~ 15, n_plastic ~ 60), which
#       has real inferential power; this is the test to lean on.
# ==============================================================================

cat("=== SECTION 6: OCCURRENCE (ANY PLASTIC vs GLASS) ===\n\n")

# Helper: take a counts table with n_plastic/k_plastic/n_glass/k_glass (+ keys),
# run the guarded Fisher test, BH within site, and add prevalence + direction.
run_occ <- function(df) {
  df <- df %>%
    mutate(enough       = n_plastic >= MIN_SAMPLES_PER_CELL & n_glass >= MIN_SAMPLES_PER_CELL,
           prop_plastic = k_plastic / n_plastic,
           prop_glass   = k_glass   / n_glass,
           Prevalence   = sprintf("%d/%d vs %d/%d", k_plastic, n_plastic, k_glass, n_glass),
           fisher_p     = mapply(function(k1, n1, k2, n2, e)
                                   if (isTRUE(e)) fisher.test(matrix(c(k1, n1 - k1, k2, n2 - k2),
                                                                     nrow = 2, byrow = TRUE))$p.value
                                   else NA_real_,
                                 k_plastic, n_plastic, k_glass, n_glass, enough))
  df$fisher_BH <- NA_real_
  for (s in SITE_LEVELS) {
    ix <- which(df$site == s & df$enough & !is.na(df$fisher_p))
    if (length(ix)) df$fisher_BH[ix] <- p.adjust(df$fisher_p[ix], method = "BH")
  }
  df %>% mutate(reliable_occ = enough & !is.na(fisher_BH) & fisher_BH < 0.05,
                occ_dir = case_when(reliable_occ & prop_plastic > prop_glass ~ "higher",
                                    reliable_occ & prop_plastic < prop_glass ~ "lower",
                                    TRUE ~ "ns"))
}

tax_lookup <- results_all %>% distinct(OTU_Label, display_label, Class)
pooled_plastic <- cell_prev %>% mutate(grp = if_else(substrate == "Glass", "glass", "plastic"))

# ── (a) PER SEASON ────────────────────────────────────────────────────────────
occ_season <- pooled_plastic %>%
  group_by(site, season, OTU_Label, grp) %>%
  summarise(n = sum(n), k = sum(k), .groups = "drop") %>%
  pivot_wider(names_from = grp, values_from = c(n, k), values_fill = 0) %>%
  filter(n_glass > 0, n_plastic > 0) %>%
  run_occ() %>%
  left_join(tax_lookup, by = "OTU_Label")

write.csv(occ_season, file.path(TABLE_DIR, "occurrence_anyplastic_per_season.csv"), row.names = FALSE)
cat("(a) PER-SEASON any-plastic-vs-glass:\n")
cat("    tested:", sum(occ_season$enough),
    "| reliable (BH<0.05):", sum(occ_season$reliable_occ, na.rm = TRUE),
    "  (descriptive only — n_glass ~3/season caps power)\n")
occ_season %>% filter(reliable_occ) %>%
  transmute(site, Taxon = display_label, Class, Season = season,
            Prevalence, prop_plastic = round(prop_plastic, 2),
            prop_glass = round(prop_glass, 2), fisher_BH = signif(fisher_BH, 3),
            direction = occ_dir) %>%
  arrange(site, desc(prop_plastic - prop_glass)) %>% print(n = Inf)
cat("\n")

# ── (b) SEASON-POOLED (the powered test) ──────────────────────────────────────
gllvm_taxon <- results_all %>%
  filter(is_substrate_contrast, term_substrate %in% SUBSTRATE_PLASTIC_LEVELS) %>%
  group_by(site, OTU_Label) %>%
  summarise(gllvm_any_reliable   = any(coalesce(reliable, FALSE)),
            gllvm_any_separation = any(is_separation), .groups = "drop")

occ_pooled <- pooled_plastic %>%
  group_by(site, OTU_Label, grp) %>%
  summarise(n = sum(n), k = sum(k), .groups = "drop") %>%
  pivot_wider(names_from = grp, values_from = c(n, k), values_fill = 0) %>%
  filter(n_glass > 0, n_plastic > 0) %>%
  run_occ() %>%
  left_join(tax_lookup, by = "OTU_Label") %>%
  left_join(gllvm_taxon, by = c("site", "OTU_Label")) %>%
  mutate(group = case_when(occ_dir == "higher" ~ "Plastic specialist (occurrence)",
                           occ_dir == "lower"  ~ "Plastic avoider (occurrence)",
                           TRUE                ~ "No reliable preference"))

write.csv(occ_pooled, file.path(TABLE_DIR, "occurrence_anyplastic_season_pooled.csv"), row.names = FALSE)
cat("(b) SEASON-POOLED any-plastic-vs-glass (powered):\n")
cat("    tested:", sum(occ_pooled$enough),
    "| reliable (BH<0.05):", sum(occ_pooled$reliable_occ, na.rm = TRUE), "\n\n")

cat("Reliable season-pooled occurrence preferences (with raw prevalence):\n")
occ_pooled %>% filter(reliable_occ) %>%
  transmute(site, Taxon = display_label, Class, group,
            Prevalence, prop_plastic = round(prop_plastic, 2),
            prop_glass = round(prop_glass, 2), fisher_BH = signif(fisher_BH, 3),
            gllvm_agrees = coalesce(gllvm_any_reliable, FALSE)) %>%
  arrange(site, group, desc(prop_plastic - prop_glass)) %>% print(n = Inf)

cat("\nSeason-pooled classification counts:\n")
print(count(occ_pooled, site, group))

cat("\nScript complete.\n")
