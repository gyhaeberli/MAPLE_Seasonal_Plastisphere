################################################################################
# PERMANOVA ANALYSIS
# Object: ps_all$pr2_collap_rep50$fP.lulu$ST
################################################################################

setwd("~/GitHub/MAPLE_Seasonal_Plastisphere/Scripts/03.StatisticalAnalysis")

dir.create("~/Github/MAPLE_Seasonal_Plastisphere/Results/PERMANOVA/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("~/Github/MAPLE_Seasonal_Plastisphere/Results/PERMANOVA/Tables",  recursive = TRUE, showWarnings = FALSE)

res_dir <- "~/Github/MAPLE_Seasonal_Plastisphere/Results/PERMANOVA/"
fig_dir <- "~/Github/MAPLE_Seasonal_Plastisphere/Results/PERMANOVA/Figures/"
tab_dir <- "~/Github/MAPLE_Seasonal_Plastisphere/Results/PERMANOVA/Tables/"

library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(pairwiseAdonis)
library(tidyr)
library(RColorBrewer)
library(knitr)
library(kableExtra)

dataset_name <- "pr2_collap_rep50_fP.lulu_ST"

# ==============================================================================
# LOAD DATA
# ==============================================================================

ps_all <- readRDS("~/GitHub/MAPLE_Seasonal_Plastisphere/Processed_data/Phyloseq_objects/FINAL_PHYLOSEQ_OBJECTS.rds")

ps <- ps_all$pr2_collap_rep50$fP.lulu$ST

# ==============================================================================
# STEP 0: DATA TIDYING
# ==============================================================================

cat("=== STEP 0: Data Tidying ===\n")

n_samples_original <- nsamples(ps)
n_taxa_original    <- ntaxa(ps)

cat("Original:", n_taxa_original, "taxa,", n_samples_original, "samples\n")

ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)

cat("After tidying:", ntaxa(ps), "taxa,", nsamples(ps), "samples\n\n")

# ==============================================================================
# STEP 1: EXTRACT OTU TABLE
# ==============================================================================

cat("=== STEP 1: Extracting OTU Table ===\n")

mat <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) {
  mat <- t(mat)
}

original_sample_names <- rownames(mat)
cat("Matrix dimensions:", dim(mat), "\n")
cat("Samples with 0 reads:", sum(rowSums(mat) == 0), "\n")
cat("Taxa with 0 reads:",    sum(colSums(mat) == 0), "\n\n")

# ==============================================================================
# STEP 2: rclr TRANSFORMATION
# ==============================================================================


cat("=== STEP 2: rclr Transformation ===\n")

al.rclr <- decostand(mat, method = "rclr", MARGIN = 1, na.rm = TRUE)
rownames(al.rclr) <- original_sample_names

cat("Row names intact:", identical(rownames(al.rclr), original_sample_names), "\n\n")

# ==============================================================================
# STEP 3: AITCHISON DISTANCE MATRIX
# ==============================================================================


cat("=== STEP 3: Aitchison Distance Matrix ===\n")

aitchison_dist <- dist(al.rclr, method = "euclidean")
cat("Distance matrix size:", attr(aitchison_dist, "Size"), "x", attr(aitchison_dist, "Size"), "\n\n")

# ==============================================================================
# STEP 4: PCA
# ==============================================================================


cat("=== STEP 4: PCA (rda on rclr data) ===\n")

my.rda        <- rda(al.rclr)
var_explained <- round(100 * my.rda$CA$eig[1:2] / my.rda$tot.chi, 1)

cat("Total inertia:", my.rda$tot.chi, "\n")
cat("PC1:", var_explained[1], "%\n")
cat("PC2:", var_explained[2], "%\n\n")

# ==============================================================================
# STEP 5: EXTRACT PC SCORES
# ==============================================================================

cat("=== STEP 5: PC Scores ===\n")

pca_scores              <- as.data.frame(scores(my.rda, display = "sites", choices = 1:2))
colnames(pca_scores)    <- c("PC1", "PC2")
rownames(pca_scores)    <- original_sample_names
pca_scores$sample_id    <- rownames(pca_scores)

cat("PC scores extracted for", nrow(pca_scores), "samples\n\n")

# ==============================================================================
# STEP 6: MERGE WITH METADATA
# ==============================================================================

cat("=== STEP 6: Metadata Merge ===\n")

sample_meta          <- as.data.frame(sample_data(ps))
class(sample_meta)   <- "data.frame"
sample_meta$sample_id <- rownames(sample_meta)

common_ids <- intersect(pca_scores$sample_id, sample_meta$sample_id)
cat("Common IDs:", length(common_ids), "\n")

pca_data          <- merge(pca_scores, sample_meta, by = "sample_id")
pca_data$season   <- factor(pca_data$season,
                            levels = c("Winter", "Spring", "Summer", "Fall", "Winter2"))

has_substrate <- "substrate" %in% colnames(sample_meta) &&
  length(unique(na.omit(sample_meta$substrate))) > 1

cat("Substrate present:", has_substrate, "\n")
cat("Merged data:", nrow(pca_data), "rows x", ncol(pca_data), "cols\n\n")

# ==============================================================================
# STEP 7: CENTROIDS (per unique biological sample, collapsed across PCR replicates)
# ==============================================================================

cat("=== STEP 7: Centroids ===\n")

# base_sample_name strips the trailing replicate digit (e.g. "Sample1" from "Sample11")
pca_data$base_sample_name <- sub("[0-9]$", "", pca_data$clean_sample_names)

pca_centroids <- pca_data %>%
  group_by(base_sample_name, site, season, substrate) %>%
  summarise(
    PC1          = mean(PC1, na.rm = TRUE),
    PC2          = mean(PC2, na.rm = TRUE),
    n_replicates = n(),
    .groups      = "drop"
  )

cat("Centroids:", nrow(pca_centroids), "(from", nrow(pca_data), "replicates)\n")
cat("Replicates per sample: min =", min(pca_centroids$n_replicates),
    ", max =", max(pca_centroids$n_replicates), "\n\n")

# ==============================================================================
# STEP 8A: PCA PLOT — ALL REPLICATES, FACETED BY SITE × SEASON
# ==============================================================================

cat("=== STEP 8A: Replicate PCA Plot ===\n")

p_rep <- ggplot(pca_data, aes(x = PC1, y = PC2, color = substrate)) +
  geom_point(size = 2.5, alpha = 0.7) +
  facet_grid(season ~ site, scales = "free") +
  labs(
    x        = paste0("PC1 (", var_explained[1], "%)"),
    y        = paste0("PC2 (", var_explained[2], "%)"),
    title    = "PCA on Aitchison Distances (rclr-transformed data)",
    subtitle = "All Biological Replicates — Faceted by Site and Season",
    color    = "Substrate"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle    = element_text(hjust = 0.5, size = 10)
  )

print(p_rep)

ggsave(
  filename = paste0(fig_dir, "aitchison_pca_replicates_site_season_", dataset_name, ".png"),
  plot     = p_rep,
  width    = 14, height = 10, dpi = 300
)

# ==============================================================================
# STEP 8B: PCA PLOT — INDIVIDUAL SITES
# ==============================================================================

cat("=== STEP 8B: Per-Site Replicate Plots ===\n")

sites_vec <- unique(pca_data$site)

# --- SELVA ---
site_data_SELVA <- pca_data[pca_data$site == "SELVA", ]

p_site_SELVA <- ggplot(site_data_SELVA, aes(x = PC1, y = PC2, color = substrate)) +
  geom_point(size = 3, alpha = 0.7) +
  facet_grid(season ~ site, scales = "free") +
  labs(
    x        = paste0("PC1 (", var_explained[1], "%)"),
    y        = paste0("PC2 (", var_explained[2], "%)"),
    title    = "PCA — Site: SELVA",
    subtitle = "Faceted by Season",
    color    = "Substrate"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = paste0(fig_dir, "aitchison_pca_replicates_", dataset_name, "_selva.png"),
  plot     = p_site_SELVA,
  width    = 12, height = 8, dpi = 300
)
cat("  Saved: selva plot\n")

# --- TBS ---
site_data_TBS <- pca_data[pca_data$site == "TBS", ]

p_site_TBS <- ggplot(site_data_TBS, aes(x = PC1, y = PC2, color = substrate)) +
  geom_point(size = 3, alpha = 0.7) +
  facet_grid(season ~ site, scales = "free") +
  labs(
    x        = paste0("PC1 (", var_explained[1], "%)"),
    y        = paste0("PC2 (", var_explained[2], "%)"),
    title    = "PCA — Site: TBS",
    subtitle = "Faceted by Season",
    color    = "Substrate"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = paste0(fig_dir, "aitchison_pca_replicates_", dataset_name, "_tbs.png"),
  plot     = p_site_TBS,
  width    = 12, height = 8, dpi = 300
)
cat("  Saved: TBS plot\n\n")

# ==============================================================================
# STEP 8C: SPIDER PLOT
# ==============================================================================

cat("=== STEP 8C: Spider Plot ===\n")

pca_data_with_centroids <- pca_data %>%
  left_join(
    pca_centroids %>% select(base_sample_name, site, season, substrate,
                             PC1_centroid = PC1, PC2_centroid = PC2),
    by = c("base_sample_name", "site", "season", "substrate")
  )

p_spider <- ggplot(pca_data_with_centroids,
                   aes(x = PC1, y = PC2, color = substrate)) +
  geom_segment(aes(xend = PC1_centroid, yend = PC2_centroid),
               alpha = 0.3, linewidth = 0.5) +
  geom_point(size = 2, alpha = 0.8) +
  geom_point(data = pca_centroids,
             aes(x = PC1, y = PC2, color = substrate),
             size = 3, alpha = 0.8, shape = 17) +
  facet_grid(season ~ site, scales = "free") +
  labs(
    x        = paste0("PC1 (", var_explained[1], "%)"),
    y        = paste0("PC2 (", var_explained[2], "%)"),
    #title    = "PCA on Aitchison Distances (rclr-transformed data)",
    #subtitle = "Spider Plot: Biological Replicates connected to Centroids",
    color    = "Substrate"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle    = element_text(hjust = 0.5, size = 10)
  )

print(p_spider)

ggsave(
  filename = paste0(fig_dir, "aitchison_pca_spider_site_season_", dataset_name, ".png"),
  plot     = p_spider,
  width    = 14, height = 10, dpi = 300
)
cat("  Saved: spider plot\n\n")


#other plots
# Define a season color palette (colorblind-friendly)

season_colors <- c(
  Winter  = "#4E79A7",
  Spring  = "#59A14F",
  Summer  = "#E15759",
  Fall    = "#F28E2B",
  Winter2 = "#9B6AC8"
)

# Define substrate shapes
substrate_shapes <- c(
  Glass         = 16 ,   
  PE            = 0,   
  PET           = 2,   
  Weathered_PE  = 15,    
  Weathered_PET = 17
)

p_spider2 <- ggplot(pca_data_with_centroids,
                    aes(x = PC1, y = PC2, color = season)) +
  # Ellipses per season (dashed, season-colored)
  stat_ellipse(aes(group = season, color = season),
               linetype = "dashed", linewidth = 0.7, level = 0.95) +
  # Replicate points — color = season, shape = substrate
  geom_point(aes(shape = substrate), size = 3, alpha = 0.8) +
  facet_wrap(~ site) +
  scale_color_manual(values = season_colors) +
  scale_shape_manual(values = substrate_shapes) +
  labs(
    x     = paste0("PC1 (", var_explained[1], "%)"),
    y     = paste0("PC2 (", var_explained[2], "%)"),
    color = "Season",
    shape = "Substrate"
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4, shape = 16, linetype = "blank")),
    shape = guide_legend(override.aes = list(size = 4, color = "gray30"))
  ) +
  theme_bw() +
  theme(
    strip.text        = element_text(size = 16, face = "bold"),
    strip.background  = element_rect(fill = "lightgray"),
    legend.position   = "right",
    legend.box        = "vertical",
    legend.text = element_text(size = 14),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.background  = element_blank(),
    plot.background   = element_blank()
  )
print(p_spider2)

ggsave(
  filename = paste0(fig_dir, "1.spider_facet_site_color_season_shape_substrate_", dataset_name, ".png"),
  plot     = p_spider2,
  width    = 14, height = 7, dpi = 300
)

ggsave(
  filename = paste0(fig_dir, "1.spider_facet_site_color_season_shape_substrate_", dataset_name, ".pdf"),
  plot     = p_spider2,
  width    = 14, height = 7
)

# ==============================================================================
# STEP 9: PREPARE METADATA FOR PERMANOVA
# ==============================================================================


cat("=== STEP 9: Metadata Preparation for PERMANOVA ===\n")

rownames(pca_data) <- pca_data$sample_id
metadata <- pca_data %>% select(-sample_id, -PC1, -PC2)

# Reorder to match distance matrix row order
dist_labels <- attr(aitchison_dist, "Labels")
metadata    <- metadata[match(dist_labels, rownames(metadata)), ]

metadata$site      <- as.factor(metadata$site)
metadata$season    <- factor(metadata$season,
                             levels = c("Winter", "Spring", "Summer", "Fall", "Winter2"))
metadata$substrate <- as.factor(metadata$substrate)

cat("Samples in PERMANOVA metadata:", nrow(metadata), "\n")
cat("Site levels:     ", paste(levels(metadata$site),      collapse = ", "), "\n")
cat("Season levels:   ", paste(levels(metadata$season),    collapse = ", "), "\n")
cat("Substrate levels:", paste(levels(metadata$substrate), collapse = ", "), "\n\n")

# ==============================================================================
# STEP 10: PERMANOVA MODELS
# ==============================================================================

cat("=== STEP 10: PERMANOVA Models ===\n\n")

# --- Model 1: Full factorial (omnibus, no by= argument) ---
# Tests the overall model as a whole. The R² and p-value refer to the combined
# effect of all terms together.
cat("--- Model 1: Full Factorial (omnibus) ---\n")
set.seed(123)
perm_full <- adonis2(aitchison_dist ~ site * season * substrate,
                     data = metadata, permutations = 999, method = "euclidean")
print(perm_full)

# --- Model 2: Full factorial by terms (Type III / marginal) ---
# by = "terms" reports each term tested *after* all other terms are already in
# the model (marginal effects).

cat("\n--- Model 2: Full Factorial by Terms (marginal) ---\n")
set.seed(123)
perm_full_byterms <- adonis2(aitchison_dist ~ site * season * substrate,
                             data = metadata, permutations = 999,
                             method = "euclidean", by = "terms")
print(perm_full_byterms)

# --- Model 3: Additive model (no interactions) ---
# Tests main effects only (site, season, substrate), each adjusted for the
# others.

set.seed(123)
perm_additive <- adonis2(aitchison_dist ~ site + season + substrate,
                         data = metadata, permutations = 999, method = "euclidean")
print(perm_additive)

# --- Model 4: Two-way interactions only ---
# Includes all pairwise interactions but not the three-way interaction.


cat("\n--- Model 4: Two-way Interactions ---\n")
set.seed(123)
perm_twoway <- adonis2(aitchison_dist ~ site * season + site * substrate + season * substrate,
                       data = metadata, permutations = 999, method = "euclidean", by = "terms")
print(perm_twoway)

# --- Model 5: Sequential (Type I SS) ---
# by = "terms" with sequential ordering: each term is tested *given* the terms
# before it (not all others). Order matters here.

cat("\n--- Model 5: Sequential (Type I) ---\n")
set.seed(123)
perm_sequential <- adonis2(aitchison_dist ~ site + season + substrate,
                           data = metadata, permutations = 999,
                           method = "euclidean", by = "terms")
print(perm_sequential)

# --- Model 6: Substrate alone ---
cat("\n--- Model 6: Substrate only ---\n")
set.seed(123)
perm_substrate <- adonis2(aitchison_dist ~ substrate,
                          data = metadata, permutations = 999, method = "euclidean")
print(perm_substrate)

# --- Model 7: Site alone ---
cat("\n--- Model 7: Site only ---\n")
set.seed(123)
perm_site <- adonis2(aitchison_dist ~ site,
                     data = metadata, permutations = 999, method = "euclidean")
print(perm_site)

# --- Model 8: Season alone ---
cat("\n--- Model 8: Season only ---\n")
set.seed(123)
perm_season <- adonis2(aitchison_dist ~ season,
                       data = metadata, permutations = 999, method = "euclidean")
print(perm_season)

# ==============================================================================
# STEP 11: HOMOGENEITY OF DISPERSION (PERMDISP)
# ==============================================================================

#
# betadisper() calculates each sample's distance to its group centroid.
# permutest() then tests whether those distances differ significantly across groups.
#
# If PERMDISP is NOT significant: dispersion is homogeneous, and a significant
# PERMANOVA result confidently reflects a difference in community composition.
#
# If PERMDISP IS significant: some of the PERMANOVA signal may be due to
# differences in spread rather than location 

cat("\n=== STEP 11: PERMDISP (Homogeneity of Dispersions) ===\n\n")

cat("--- Site ---\n")
set.seed(123)
disp_site      <- betadisper(aitchison_dist, metadata$site)
perm_disp_site <- permutest(disp_site, permutations = 999)
print(perm_disp_site)
cat("Average distances to centroid by site:\n")
print(disp_site$group.distances)

cat("\n--- Season ---\n")
set.seed(123)
disp_season      <- betadisper(aitchison_dist, metadata$season)
perm_disp_season <- permutest(disp_season, permutations = 999)
print(perm_disp_season)
cat("Average distances to centroid by season:\n")
print(disp_season$group.distances)

cat("\n--- Substrate ---\n")
set.seed(123)
disp_substrate      <- betadisper(aitchison_dist, metadata$substrate)
perm_disp_substrate <- permutest(disp_substrate, permutations = 999)
print(perm_disp_substrate)
cat("Average distances to centroid by substrate:\n")
print(disp_substrate$group.distances)


# Reporting results in a table




##############  RESULTS TABLE ################################################


# ── Pull values from your models object ──────────────────────────────────────

perm <- as.data.frame(perm_full_byterms)
disp <- list(
  site      = perm_disp_site$tab,
  season    = perm_disp_season$tab,
  substrate = perm_disp_substrate$tab
)

# ── Helper: significance stars ────────────────────────────────────────────────

sig_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p <= 0.001, "***",
                ifelse(p <= 0.01,  "**",
                       ifelse(p <= 0.05,  "*",
                              ifelse(p <= 0.1,   ".", "ns")))))
}

# ── Build PERMANOVA rows ──────────────────────────────────────────────────────

row_names <- c(
  "Site",
  "Season",
  "Substrate",
  "Site × Season",
  "Site × Substrate",
  "Season × Substrate",
  "Site × Season × Substrate",
  "Residual",
  "Total"
)

perm_table <- data.frame(
  Parameter = row_names,
  df        = perm$Df,
  R2        = round(perm$R2, 3),
  F         = c(round(perm$F[1:7], 3), NA, NA),
  p         = c(perm$`Pr(>F)`[1:7], NA, NA),
  stringsAsFactors = FALSE
)

perm_table$Sig  <- sig_stars(perm_table$p)
perm_table$p    <- ifelse(is.na(perm_table$p), "", formatC(perm_table$p, digits = 3, format = "f"))
perm_table$F    <- ifelse(is.na(perm_table$F), "", formatC(perm_table$F, digits = 3, format = "f"))
perm_table$R2   <- ifelse(perm_table$R2 == 0, "", formatC(perm_table$R2, digits = 3, format = "f"))

# ── Build PERMDISP rows ───────────────────────────────────────────────────────

disp_table <- data.frame(
  Parameter = c("Site", "Season", "Substrate"),
  df        = c(disp$site$Df[1],      disp$season$Df[1],      disp$substrate$Df[1]),
  R2        = "",
  F         = round(c(disp$site$F[1], disp$season$F[1],       disp$substrate$F[1]), 3),
  p         = c(disp$site$`Pr(>F)`[1], disp$season$`Pr(>F)`[1], disp$substrate$`Pr(>F)`[1]),
  stringsAsFactors = FALSE
)

disp_table$Sig <- sig_stars(disp_table$p)
disp_table$p   <- formatC(disp_table$p, digits = 3, format = "f")
disp_table$F   <- formatC(disp_table$F, digits = 3, format = "f")

# ── Combine with section headers ──────────────────────────────────────────────

# header_row <- data.frame(
#   Parameter = NA, df = NA, R2 = NA, F = NA, p = NA, Sig = NA,
#   stringsAsFactors = FALSE
# )

full_table <- bind_rows(
  perm_table,
  disp_table
)

# ── Render with kableExtra ────────────────────────────────────────────────────

# Which rows have significant p-values?
sig_rows <- which(full_table$Sig %in% c("*", "**", "***")) 

colnames(full_table) <- c("Parameter", "df", "R²", "F", "p-value", "")

kbl(full_table,
    #caption  = "Table 1. PERMANOVA (adonis2) and PERMDISP results for eukaryotic community composition (Aitchison distance, 999 permutations).",
    booktabs = TRUE,
    align    = c("l", "r", "r", "r", "r", "l"),
    na       = "") %>%
  kable_classic(full_width = FALSE, html_font = "Arial") %>%
  pack_rows("PERMANOVA", 1, 9) %>%
  pack_rows("Homogeneity of dispersion (PERMDISP)", 10, 12) %>%
  row_spec(c(8, 9), italic = TRUE, color = "gray") %>%
  column_spec(5, bold = ifelse(seq_len(nrow(full_table)) %in% sig_rows, TRUE, FALSE))


# mean dispersion values

disp_distances <- bind_rows(
  data.frame(
    Factor  = "Site",
    Group   = names(disp_site$group.distances),
    `Mean distance to centroid` = round(disp_site$group.distances, 2),
    check.names = FALSE
  ),
  data.frame(
    Factor  = "Season",
    Group   = names(disp_season$group.distances),
    `Mean distance to centroid` = round(disp_season$group.distances, 2),
    check.names = FALSE
  ),
  data.frame(
    Factor  = "Substrate",
    Group   = names(disp_substrate$group.distances),
    `Mean distance to centroid` = round(disp_substrate$group.distances, 2),
    check.names = FALSE
  )
)

kbl(disp_distances,
    booktabs    = TRUE,
    align       = c("l", "l", "r"),
    na          = "",
    row.names   = FALSE) %>%
  kable_classic(full_width = FALSE, html_font = "Arial") %>%
  pack_rows("Site (p = 0.001 ***)",      1, 2) %>%
  pack_rows("Season (p = 0.001 ***)",    3, 7) %>%
  pack_rows("Substrate (p = 0.856 ns)", 8, 12)

library(flextable)
library(officer)
library(dplyr)

# PERMANOVA and PERMDISP
print(full_table)
print(disp_distances)

# Fix column names FIRST (before any dplyr)
names(full_table)[names(full_table) == "" | is.na(names(full_table))] <- "empty"
colnames(full_table) <- make.names(names(full_table), unique = TRUE)

# Now assign clean names explicitly
colnames(full_table)[1:5] <- c("Parameter", "df", "R2", "F", "p_value")

# Remove row names safely
full_table <- as.data.frame(full_table)
rownames(full_table) <- NULL

# Convert numeric columns (handles blanks safely)
full_table$df <- as.numeric(full_table$df)
full_table$F  <- as.numeric(full_table$F)
full_table$R2 <- suppressWarnings(as.numeric(full_table$R2))

ft <- flextable(full_table) %>%
  autofit()
doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "PermanovaResults.docx")

#CENTROID Distances


ft <- flextable(disp_distances) %>%
  autofit()

doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "CentroidDist.docx")

# ==============================================================================
# STEP 12: PAIRWISE COMPARISONS (POST-HOC)
# ==============================================================================

# Multiple testing correction: running many tests inflates the chance of false
# positives. Two corrections are applied:
#   - Bonferroni: strict; multiply each p-value by the number of tests.
#     Appropriate when you want to control the family-wise error rate tightly.
#   - BH (Benjamini-Hochberg): less strict; controls the false discovery rate
#     (expected proportion of false positives among significant results).
#     Preferred when you have many comparisons and want more power.

cat("\n=== STEP 12: Pairwise Comparisons ===\n\n")

sig_stars <- function(p) {
  ifelse(p <= 0.001, "***",
         ifelse(p <= 0.01,  "**",
                ifelse(p <= 0.05,  "*",
                       ifelse(p <= 0.1,   ".", " "))))
}

# --- Pairwise: Site ---
cat("--- Pairwise: Site ---\n")

site_levels    <- levels(metadata$site)
pairwise_site  <- data.frame()

for (i in 1:(length(site_levels) - 1)) {
  for (j in (i + 1):length(site_levels)) {
    sub_data <- metadata[metadata$site %in% c(site_levels[i], site_levels[j]), ]
    sub_dist <- as.dist(as.matrix(aitchison_dist)[
      metadata$site %in% c(site_levels[i], site_levels[j]),
      metadata$site %in% c(site_levels[i], site_levels[j])])
    set.seed(123)
    res    <- adonis2(sub_dist ~ site, data = sub_data, permutations = 999)
    res_df <- as.data.frame(res)
    pairwise_site <- rbind(pairwise_site, data.frame(
      Comparison = paste(site_levels[i], "vs", site_levels[j]),
      R2         = round(res_df$R2[1], 4),
      F_value    = round(res_df[["F"]][1], 3),
      raw_p      = round(res_df$`Pr(>F)`[1], 4),
      stringsAsFactors = FALSE
    ))
  }
}

pairwise_site$p_Bonferroni <- round(p.adjust(pairwise_site$raw_p, "bonferroni"), 4)
pairwise_site$p_BH         <- round(p.adjust(pairwise_site$raw_p, "BH"), 4)
pairwise_site$Sig_raw      <- sig_stars(pairwise_site$raw_p)
pairwise_site$Sig_Bonf     <- sig_stars(pairwise_site$p_Bonferroni)
pairwise_site$Sig_BH       <- sig_stars(pairwise_site$p_BH)

print(pairwise_site)

# --- Pairwise: Season ---
cat("\n--- Pairwise: Season ---\n")

season_levels    <- levels(metadata$season)
pairwise_season  <- data.frame()

for (i in 1:(length(season_levels) - 1)) {
  for (j in (i + 1):length(season_levels)) {
    sub_data <- metadata[metadata$season %in% c(season_levels[i], season_levels[j]), ]
    sub_dist <- as.dist(as.matrix(aitchison_dist)[
      metadata$season %in% c(season_levels[i], season_levels[j]),
      metadata$season %in% c(season_levels[i], season_levels[j])])
    set.seed(123)
    res    <- adonis2(sub_dist ~ season, data = sub_data, permutations = 999)
    res_df <- as.data.frame(res)
    pairwise_season <- rbind(pairwise_season, data.frame(
      Comparison = paste(season_levels[i], "vs", season_levels[j]),
      R2         = round(res_df$R2[1], 4),
      F_value    = round(res_df[["F"]][1], 3),
      raw_p      = round(res_df$`Pr(>F)`[1], 4),
      stringsAsFactors = FALSE
    ))
  }
}

pairwise_season$p_Bonferroni <- round(p.adjust(pairwise_season$raw_p, "bonferroni"), 4)
pairwise_season$p_BH         <- round(p.adjust(pairwise_season$raw_p, "BH"), 4)
pairwise_season$Sig_raw      <- sig_stars(pairwise_season$raw_p)
pairwise_season$Sig_Bonf     <- sig_stars(pairwise_season$p_Bonferroni)
pairwise_season$Sig_BH       <- sig_stars(pairwise_season$p_BH)

print(pairwise_season)

# --- Pairwise: Substrate ---
cat("\n--- Pairwise: Substrate ---\n")

substrate_levels    <- levels(metadata$substrate)
pairwise_substrate  <- data.frame()

for (i in 1:(length(substrate_levels) - 1)) {
  for (j in (i + 1):length(substrate_levels)) {
    sub_data <- metadata[metadata$substrate %in% c(substrate_levels[i], substrate_levels[j]), ]
    sub_dist <- as.dist(as.matrix(aitchison_dist)[
      metadata$substrate %in% c(substrate_levels[i], substrate_levels[j]),
      metadata$substrate %in% c(substrate_levels[i], substrate_levels[j])])
    set.seed(123)
    res    <- adonis2(sub_dist ~ substrate, data = sub_data, permutations = 999)
    res_df <- as.data.frame(res)
    pairwise_substrate <- rbind(pairwise_substrate, data.frame(
      Comparison = paste(substrate_levels[i], "vs", substrate_levels[j]),
      R2         = round(res_df$R2[1], 4),
      F_value    = round(res_df[["F"]][1], 3),
      raw_p      = round(res_df$`Pr(>F)`[1], 4),
      stringsAsFactors = FALSE
    ))
  }
}

pairwise_substrate$p_Bonferroni <- round(p.adjust(pairwise_substrate$raw_p, "bonferroni"), 4)
pairwise_substrate$p_BH         <- round(p.adjust(pairwise_substrate$raw_p, "BH"), 4)
pairwise_substrate$Sig_raw      <- sig_stars(pairwise_substrate$raw_p)
pairwise_substrate$Sig_Bonf     <- sig_stars(pairwise_substrate$p_Bonferroni)
pairwise_substrate$Sig_BH       <- sig_stars(pairwise_substrate$p_BH)

print(pairwise_substrate)


## TABLE with all

# Helper to format pairwise df for display
format_pairwise <- function(df) {
  df %>%
    mutate(
      R2      = formatC(R2,      digits = 3, format = "f"),
      F_value = formatC(F_value, digits = 3, format = "f"),
      raw_p   = formatC(raw_p,   digits = 3, format = "f"),
      p_BH    = formatC(p_BH,    digits = 3, format = "f")
    ) %>%
    select(
      Comparison,
      `R²`       = R2,
      `F`        = F_value,
      `p (raw)`  = raw_p,
      ` `        = Sig_raw,
      `p (BH)`   = p_BH,
      `  `       = Sig_BH
    )
}

site_fmt      <- format_pairwise(pairwise_site)
season_fmt    <- format_pairwise(pairwise_season)
substrate_fmt <- format_pairwise(pairwise_substrate)

# Row indices for pack_rows
n_site      <- nrow(site_fmt)
n_season    <- nrow(season_fmt)
n_substrate <- nrow(substrate_fmt)

full_pairwise <- bind_rows(site_fmt, season_fmt, substrate_fmt)

# Significant rows for bold p (BH column is col 6)
sig_rows_bh <- which(
  c(pairwise_site$p_BH, pairwise_season$p_BH, pairwise_substrate$p_BH) <= 0.05
)

kbl(full_pairwise,
    booktabs  = TRUE,
    align     = c("l", "r", "r", "r", "l", "r", "l"),
    na        = "",
    row.names = FALSE) %>%
    #caption   = "Table S2. Pairwise PERMANOVA comparisons (Aitchison distance, 999 permutations). p-values adjusted using Benjamini-Hochberg (BH) correction.") %>%
  kable_classic(full_width = TRUE, html_font = "Arial") %>%
  pack_rows("Site",      1,                        n_site) %>%
  pack_rows("Season",    n_site + 1,               n_site + n_season) %>%
  pack_rows("Substrate", n_site + n_season + 1,    n_site + n_season + n_substrate) %>%
  column_spec(6, bold = ifelse(seq_len(nrow(full_pairwise)) %in% sig_rows_bh, TRUE, FALSE)) 


ft <- flextable(full_pairwise) %>%
  autofit()
doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "PairwisePerm.docx")


# ==============================================================================
# STEP 14: INTERACTION ANALYSIS — PLOTS
# ==============================================================================

cat("\n=== STEP 14: Interaction Plots ===\n\n")

metadata_with_pca <- cbind(metadata, pca_data[match(rownames(metadata), pca_data$sample_id),
                                              c("PC1", "PC2")])

# Plot 1: Substrate within each Season
p_int1 <- ggplot(metadata_with_pca, aes(x = substrate, y = PC1, fill = substrate)) +
  geom_boxplot(alpha = 0.6, outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_wrap(~ season, scales = "free_y") +
  labs(
    title    = paste("Substrate Effect on PC1 within Each Season —", dataset_name),
    subtitle = "Does substrate matter differently across seasons?",
    y = paste0("PC1 (", var_explained[1], "%)"),
    x = "Substrate"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "lightblue"),
        strip.text = element_text(face = "bold"))

ggsave(paste0(fig_dir, "interaction_substrate_season_PC1_", dataset_name, ".png"),
       p_int1, width = 12, height = 8, dpi = 300)

# Plot 2: Full interaction — Substrate within Site × Season
p_int2 <- ggplot(metadata_with_pca, aes(x = substrate, y = PC1, fill = substrate)) +
  geom_boxplot(alpha = 0.6, outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
  facet_grid(site ~ season, scales = "free_y") +
  labs(
    title    = paste("Substrate Effect within Site × Season —", dataset_name),
    subtitle = "Full interaction visualization",
    y = paste0("PC1 (", var_explained[1], "%)"),
    x = "Substrate"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold", size = 9))

ggsave(paste0(fig_dir, "interaction_full_PC1_", dataset_name, ".png"),
       p_int2, width = 14, height = 8, dpi = 300)

# Plot 3: Interaction line plot (mean ± SE by season × substrate)
# --- Why line plots for interactions? ---
# Parallel lines = no interaction (substrate effect is the same across seasons).
# Non-parallel (crossing) lines = interaction (substrate effect changes by season).

mean_data <- metadata_with_pca %>%
  group_by(season, substrate) %>%
  summarise(
    mean_PC1 = mean(PC1, na.rm = TRUE),
    se_PC1   = sd(PC1,   na.rm = TRUE) / sqrt(n()),
    .groups  = "drop"
  )

p_int3 <- ggplot(mean_data, aes(x = season, y = mean_PC1,
                                color = substrate, group = substrate)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_PC1 - se_PC1, ymax = mean_PC1 + se_PC1),
                width = 0.2) +
  labs(
    title    = paste("Substrate × Season Interaction Pattern —", dataset_name),
    subtitle = "Non-parallel lines indicate interaction",
    y        = paste0("Mean PC1 (", var_explained[1], "%) ± SE"),
    x        = "Season",
    color    = "Substrate"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(fig_dir, "interaction_lines_substrate_season_", dataset_name, ".png"),
       p_int3, width = 10, height = 7, dpi = 300)

cat("  Interaction plots saved.\n\n")

# ==============================================================================
# STEP 15: DISPERSION PLOTS
# ==============================================================================

cat("=== STEP 15: Dispersion Plots ===\n\n")

disp_df <- data.frame(
  sample_id      = names(disp_site$distances),
  site           = metadata$site,
  season         = metadata$season,
  substrate      = metadata$substrate,
  disp_site      = disp_site$distances,
  disp_season    = disp_season$distances,
  disp_substrate = disp_substrate$distances
)

p_disp_site <- ggplot(disp_df, aes(x = site, y = disp_site, fill = site)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
  labs(title    = "Dispersion by Site",
       subtitle = paste0("PERMDISP p = ", round(perm_disp_site$tab$`Pr(>F)`[1], 3)),
       y = "Distance to Centroid", x = "Site") +
  theme_bw() + theme(legend.position = "none")

p_disp_season <- ggplot(disp_df, aes(x = season, y = disp_season, fill = season)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
  labs(title    = "Dispersion by Season",
       subtitle = paste0("PERMDISP p = ", round(perm_disp_season$tab$`Pr(>F)`[1], 3)),
       y = "Distance to Centroid", x = "Season") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p_disp_substrate <- ggplot(disp_df, aes(x = substrate, y = disp_substrate, fill = substrate)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
  labs(title    = "Dispersion by Substrate",
       subtitle = paste0("PERMDISP p = ", round(perm_disp_substrate$tab$`Pr(>F)`[1], 3)),
       y = "Distance to Centroid", x = "Substrate") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

combined_disp <- gridExtra::grid.arrange(p_disp_site, p_disp_season, p_disp_substrate,
                                         ncol = 3)
ggsave(paste0(fig_dir, "dispersion_analysis_", dataset_name, ".png"),
       combined_disp, width = 14, height = 5, dpi = 300)

cat("  Dispersion plots saved.\n\n")

# Dispersion summary
cat("Dispersion summary — mean distance to centroid:\n")
cat("  Site:\n");      print(tapply(disp_site$distances,      metadata$site,      mean))
cat("  Season:\n");    print(tapply(disp_season$distances,    metadata$season,    mean))
cat("  Substrate:\n"); print(tapply(disp_substrate$distances, metadata$substrate, mean))

# ==============================================================================
# STEP 16: SUBSTRATE EFFECT WITHIN EACH SEASON
# ==============================================================================

cat("\n=== STEP 16: Substrate × Season Conditional Tests ===\n\n")

season_results_df <- data.frame()

# Winter
idx_winter  <- metadata$season == "Winter"
sub_meta_winter <- metadata[idx_winter, ]
sub_dist_winter <- as.dist(as.matrix(aitchison_dist)[idx_winter, idx_winter])
set.seed(123)
res_winter  <- adonis2(sub_dist_winter ~ substrate, data = sub_meta_winter, permutations = 999)
print(res_winter)
season_results_df <- rbind(season_results_df, data.frame(
  Season = "Winter", N = sum(idx_winter),
  R2 = round(res_winter$R2[1], 4), F_stat = round(res_winter$F[1], 3),
  raw_p = res_winter$`Pr(>F)`[1], stringsAsFactors = FALSE))

# Spring
idx_spring  <- metadata$season == "Spring"
sub_meta_spring <- metadata[idx_spring, ]
sub_dist_spring <- as.dist(as.matrix(aitchison_dist)[idx_spring, idx_spring])
set.seed(123)
res_spring  <- adonis2(sub_dist_spring ~ substrate, data = sub_meta_spring, permutations = 999)
print(res_spring)
season_results_df <- rbind(season_results_df, data.frame(
  Season = "Spring", N = sum(idx_spring),
  R2 = round(res_spring$R2[1], 4), F_stat = round(res_spring$F[1], 3),
  raw_p = res_spring$`Pr(>F)`[1], stringsAsFactors = FALSE))

# Summer
idx_summer  <- metadata$season == "Summer"
sub_meta_summer <- metadata[idx_summer, ]
sub_dist_summer <- as.dist(as.matrix(aitchison_dist)[idx_summer, idx_summer])
set.seed(123)
res_summer  <- adonis2(sub_dist_summer ~ substrate, data = sub_meta_summer, permutations = 999)
print(res_summer)
season_results_df <- rbind(season_results_df, data.frame(
  Season = "Summer", N = sum(idx_summer),
  R2 = round(res_summer$R2[1], 4), F_stat = round(res_summer$F[1], 3),
  raw_p = res_summer$`Pr(>F)`[1], stringsAsFactors = FALSE))

# Fall
idx_fall    <- metadata$season == "Fall"
sub_meta_fall <- metadata[idx_fall, ]
sub_dist_fall <- as.dist(as.matrix(aitchison_dist)[idx_fall, idx_fall])
set.seed(123)
res_fall    <- adonis2(sub_dist_fall ~ substrate, data = sub_meta_fall, permutations = 999)
print(res_fall)
season_results_df <- rbind(season_results_df, data.frame(
  Season = "Fall", N = sum(idx_fall),
  R2 = round(res_fall$R2[1], 4), F_stat = round(res_fall$F[1], 3),
  raw_p = res_fall$`Pr(>F)`[1], stringsAsFactors = FALSE))

# Winter2
idx_winter2 <- metadata$season == "Winter2"
sub_meta_winter2 <- metadata[idx_winter2, ]
sub_dist_winter2 <- as.dist(as.matrix(aitchison_dist)[idx_winter2, idx_winter2])
set.seed(123)
res_winter2 <- adonis2(sub_dist_winter2 ~ substrate, data = sub_meta_winter2, permutations = 999)
print(res_winter2)
season_results_df <- rbind(season_results_df, data.frame(
  Season = "Winter2", N = sum(idx_winter2),
  R2 = round(res_winter2$R2[1], 4), F_stat = round(res_winter2$F[1], 3),
  raw_p = res_winter2$`Pr(>F)`[1], stringsAsFactors = FALSE))

# Apply corrections
season_results_df$p_Bonferroni <- round(p.adjust(season_results_df$raw_p, "bonferroni"), 4)
season_results_df$p_BH         <- round(p.adjust(season_results_df$raw_p, "BH"), 4)
season_results_df$Sig_raw      <- sig_stars(season_results_df$raw_p)
season_results_df$Sig_Bonf     <- sig_stars(season_results_df$p_Bonferroni)
season_results_df$Sig_BH       <- sig_stars(season_results_df$p_BH)

cat("\nSubstrate effects within each season:\n")
print(season_results_df, row.names = FALSE)

kable(season_results_df) %>% kable_styling(latex_options = "striped")

# ==============================================================================
# STEP 17: SUBSTRATE EFFECT WITHIN EACH SITE × SEASON
# ==============================================================================

cat("\n=== STEP 17: Substrate Effect within Site × Season Contexts ===\n\n")

site_season_results <- data.frame()

# SELVA × Winter
idx_sw <- metadata$site == "SELVA" & metadata$season == "Winter"
if (sum(idx_sw) >= 5 && length(unique(metadata$substrate[idx_sw])) > 1) {
  sub_m <- metadata[idx_sw, ]; sub_d <- as.dist(as.matrix(aitchison_dist)[idx_sw, idx_sw])
  set.seed(123); r <- adonis2(sub_d ~ substrate, data = sub_m, permutations = 999)
  site_season_results <- rbind(site_season_results, data.frame(
    Site = "SELVA", Season = "Winter", N = sum(idx_sw),
    R2 = round(r$R2[1], 4), F_stat = round(r$F[1], 3), raw_p = r$`Pr(>F)`[1]))
}

# TBS × Winter
idx_tw <- metadata$site == "TBS" & metadata$season == "Winter"
if (sum(idx_tw) >= 5 && length(unique(metadata$substrate[idx_tw])) > 1) {
  sub_m <- metadata[idx_tw, ]; sub_d <- as.dist(as.matrix(aitchison_dist)[idx_tw, idx_tw])
  set.seed(123); r <- adonis2(sub_d ~ substrate, data = sub_m, permutations = 999)
  site_season_results <- rbind(site_season_results, data.frame(
    Site = "TBS", Season = "Winter", N = sum(idx_tw),
    R2 = round(r$R2[1], 4), F_stat = round(r$F[1], 3), raw_p = r$`Pr(>F)`[1]))
}

# SELVA × Spring
idx_ss <- metadata$site == "SELVA" & metadata$season == "Spring"
if (sum(idx_ss) >= 5 && length(unique(metadata$substrate[idx_ss])) > 1) {
  sub_m <- metadata[idx_ss, ]; sub_d <- as.dist(as.matrix(aitchison_dist)[idx_ss, idx_ss])
  set.seed(123); r <- adonis2(sub_d ~ substrate, data = sub_m, permutations = 999)
  site_season_results <- rbind(site_season_results, data.frame(
    Site = "SELVA", Season = "Spring", N = sum(idx_ss),
    R2 = round(r$R2[1], 4), F_stat = round(r$F[1], 3), raw_p = r$`Pr(>F)`[1]))
}

# TBS × Spring
idx_ts <- metadata$site == "TBS" & metadata$season == "Spring"
if (sum(idx_ts) >= 5 && length(unique(metadata$substrate[idx_ts])) > 1) {
  sub_m <- metadata[idx_ts, ]; sub_d <- as.dist(as.matrix(aitchison_dist)[idx_ts, idx_ts])
  set.seed(123); r <- adonis2(sub_d ~ substrate, data = sub_m, permutations = 999)
  site_season_results <- rbind(site_season_results, data.frame(
    Site = "TBS", Season = "Spring", N = sum(idx_ts),
    R2 = round(r$R2[1], 4), F_stat = round(r$F[1], 3), raw_p = r$`Pr(>F)`[1]))
}

# SELVA × Summer
idx_ssu <- metadata$site == "SELVA" & metadata$season == "Summer"
if (sum(idx_ssu) >= 5 && length(unique(metadata$substrate[idx_ssu])) > 1) {
  sub_m <- metadata[idx_ssu, ]; sub_d <- as.dist(as.matrix(aitchison_dist)[idx_ssu, idx_ssu])
  set.seed(123); r <- adonis2(sub_d ~ substrate, data = sub_m, permutations = 999)
  site_season_results <- rbind(site_season_results, data.frame(
    Site = "SELVA", Season = "Summer", N = sum(idx_ssu),
    R2 = round(r$R2[1], 4), F_stat = round(r$F[1], 3), raw_p = r$`Pr(>F)`[1]))
}

# TBS × Summer
idx_tsu <- metadata$site == "TBS" & metadata$season == "Summer"
if (sum(idx_tsu) >= 5 && length(unique(metadata$substrate[idx_tsu])) > 1) {
  sub_m <- metadata[idx_tsu, ]; sub_d <- as.dist(as.matrix(aitchison_dist)[idx_tsu, idx_tsu])
  set.seed(123); r <- adonis2(sub_d ~ substrate, data = sub_m, permutations = 999)
  site_season_results <- rbind(site_season_results, data.frame(
    Site = "TBS", Season = "Summer", N = sum(idx_tsu),
    R2 = round(r$R2[1], 4), F_stat = round(r$F[1], 3), raw_p = r$`Pr(>F)`[1]))
}

# SELVA × Fall
idx_sf <- metadata$site == "SELVA" & metadata$season == "Fall"
if (sum(idx_sf) >= 5 && length(unique(metadata$substrate[idx_sf])) > 1) {
  sub_m <- metadata[idx_sf, ]; sub_d <- as.dist(as.matrix(aitchison_dist)[idx_sf, idx_sf])
  set.seed(123); r <- adonis2(sub_d ~ substrate, data = sub_m, permutations = 999)
  site_season_results <- rbind(site_season_results, data.frame(
    Site = "SELVA", Season = "Fall", N = sum(idx_sf),
    R2 = round(r$R2[1], 4), F_stat = round(r$F[1], 3), raw_p = r$`Pr(>F)`[1]))
}

# TBS × Fall
idx_tf <- metadata$site == "TBS" & metadata$season == "Fall"
if (sum(idx_tf) >= 5 && length(unique(metadata$substrate[idx_tf])) > 1) {
  sub_m <- metadata[idx_tf, ]; sub_d <- as.dist(as.matrix(aitchison_dist)[idx_tf, idx_tf])
  set.seed(123); r <- adonis2(sub_d ~ substrate, data = sub_m, permutations = 999)
  site_season_results <- rbind(site_season_results, data.frame(
    Site = "TBS", Season = "Fall", N = sum(idx_tf),
    R2 = round(r$R2[1], 4), F_stat = round(r$F[1], 3), raw_p = r$`Pr(>F)`[1]))
}

# SELVA × Winter2
idx_sw2 <- metadata$site == "SELVA" & metadata$season == "Winter2"
if (sum(idx_sw2) >= 5 && length(unique(metadata$substrate[idx_sw2])) > 1) {
  sub_m <- metadata[idx_sw2, ]; sub_d <- as.dist(as.matrix(aitchison_dist)[idx_sw2, idx_sw2])
  set.seed(123); r <- adonis2(sub_d ~ substrate, data = sub_m, permutations = 999)
  site_season_results <- rbind(site_season_results, data.frame(
    Site = "SELVA", Season = "Winter2", N = sum(idx_sw2),
    R2 = round(r$R2[1], 4), F_stat = round(r$F[1], 3), raw_p = r$`Pr(>F)`[1]))
}

# TBS × Winter2
idx_tw2 <- metadata$site == "TBS" & metadata$season == "Winter2"
if (sum(idx_tw2) >= 5 && length(unique(metadata$substrate[idx_tw2])) > 1) {
  sub_m <- metadata[idx_tw2, ]; sub_d <- as.dist(as.matrix(aitchison_dist)[idx_tw2, idx_tw2])
  set.seed(123); r <- adonis2(sub_d ~ substrate, data = sub_m, permutations = 999)
  site_season_results <- rbind(site_season_results, data.frame(
    Site = "TBS", Season = "Winter2", N = sum(idx_tw2),
    R2 = round(r$R2[1], 4), F_stat = round(r$F[1], 3), raw_p = r$`Pr(>F)`[1]))
}

# Apply corrections
site_season_results$p_Bonferroni <- round(p.adjust(site_season_results$raw_p, "bonferroni"), 4)
site_season_results$p_BH         <- round(p.adjust(site_season_results$raw_p, "BH"), 4)
site_season_results$Sig_raw      <- sig_stars(site_season_results$raw_p)
site_season_results$Sig_Bonf     <- sig_stars(site_season_results$p_Bonferroni)
site_season_results$Sig_BH       <- sig_stars(site_season_results$p_BH)

cat("Substrate effects within each Site × Season:\n")
print(site_season_results, row.names = FALSE)

kable(site_season_results) %>% kable_styling(latex_options = "striped")

# ==============================================================================
# STEP 18: PAIRWISE SUBSTRATE IN KEY CONTEXTS (Fall & Winter2)
# ==============================================================================

cat("\n=== STEP 18: Pairwise Substrate Comparisons — Fall & Winter2 ===\n\n")

run_substrate_pairwise_in_context <- function(aitchison_dist, metadata,
                                              site_val, season_val) {
  idx       <- metadata$site == site_val & metadata$season == season_val
  sub_meta  <- metadata[idx, ]
  sub_dist  <- as.dist(as.matrix(aitchison_dist)[idx, idx])
  subs      <- levels(droplevels(sub_meta$substrate))
  out       <- data.frame()
  for (i in 1:(length(subs) - 1)) {
    for (j in (i + 1):length(subs)) {
      pair_idx  <- sub_meta$substrate %in% c(subs[i], subs[j])
      pair_meta <- sub_meta[pair_idx, ]
      pair_dist <- as.dist(as.matrix(sub_dist)[pair_idx, pair_idx])
      set.seed(123)
      r   <- adonis2(pair_dist ~ substrate, data = pair_meta, permutations = 999)
      rdf <- as.data.frame(r)
      out <- rbind(out, data.frame(
        Comparison = paste(subs[i], "vs", subs[j]),
        R2         = round(rdf$R2[1], 4),
        F_value    = round(rdf[["F"]][1], 3),
        raw_p      = round(rdf$`Pr(>F)`[1], 4),
        stringsAsFactors = FALSE))
    }
  }
  out$p_Bonferroni <- round(p.adjust(out$raw_p, "bonferroni"), 4)
  out$p_BH         <- round(p.adjust(out$raw_p, "BH"), 4)
  out$Sig_raw      <- sig_stars(out$raw_p)
  out$Sig_Bonf     <- sig_stars(out$p_Bonferroni)
  out$Sig_BH       <- sig_stars(out$p_BH)
  return(out)
}

cat("SELVA — Fall:\n")
pw_SELVA_Fall    <- run_substrate_pairwise_in_context(aitchison_dist, metadata, "SELVA", "Fall")
print(pw_SELVA_Fall)

cat("\nTBS — Fall:\n")
pw_TBS_Fall      <- run_substrate_pairwise_in_context(aitchison_dist, metadata, "TBS", "Fall")
print(pw_TBS_Fall)

cat("\nSELVA — Winter2:\n")
pw_SELVA_Winter2 <- run_substrate_pairwise_in_context(aitchison_dist, metadata, "SELVA", "Winter2")
print(pw_SELVA_Winter2)

cat("\nTBS — Winter2:\n")
pw_TBS_Winter2   <- run_substrate_pairwise_in_context(aitchison_dist, metadata, "TBS", "Winter2")
print(pw_TBS_Winter2)

# Combined table
all_context_pairwise <- bind_rows(
  pw_SELVA_Fall    %>% mutate(Context = "SELVA-Fall"),
  pw_TBS_Fall      %>% mutate(Context = "TBS-Fall"),
  pw_SELVA_Winter2 %>% mutate(Context = "SELVA-Winter2"),
  pw_TBS_Winter2   %>% mutate(Context = "TBS-Winter2")
) %>% select(Context, Comparison, R2, F_value, raw_p, Sig_raw,
             p_Bonferroni, Sig_Bonf, p_BH, Sig_BH)

write.csv(all_context_pairwise,
          paste0(tab_dir, "pairwise_substrates_contexts_", dataset_name, ".csv"),
          row.names = FALSE)

cat("\nContext pairwise table saved.\n")

# ==============================================================================
# STEP 19: GROUPED SUBSTRATE COMPARISONS
# ==============================================================================

cat("\n=== STEP 19: Grouped Substrate Comparisons ===\n\n")

metadata$substrate_grouped <- as.character(metadata$substrate)

# --- Glass vs New Plastic (PE + PET) ---
cat("--- Glass vs New Plastic ---\n")

idx_gnp  <- metadata$substrate_grouped %in% c("Glass", "PE", "PET")
meta_gnp <- metadata[idx_gnp, ]
dist_gnp <- as.dist(as.matrix(aitchison_dist)[idx_gnp, idx_gnp])
meta_gnp$substrate_group <- ifelse(meta_gnp$substrate_grouped == "Glass", "Glass", "New_Plastic")
meta_gnp$substrate_group <- as.factor(meta_gnp$substrate_group)
cat("Sample sizes:\n"); print(table(meta_gnp$substrate_group))
set.seed(123)
perm_gnp <- adonis2(dist_gnp ~ substrate_group, data = meta_gnp, permutations = 999)
print(perm_gnp)
set.seed(123)
disp_gnp      <- betadisper(dist_gnp, meta_gnp$substrate_group)
perm_disp_gnp <- permutest(disp_gnp, permutations = 999)
print(perm_disp_gnp)

# --- Glass vs Weathered Plastic ---
cat("\n--- Glass vs Weathered Plastic ---\n")

idx_gwp  <- metadata$substrate_grouped %in% c("Glass", "Weathered_PE", "Weathered_PET")
meta_gwp <- metadata[idx_gwp, ]
dist_gwp <- as.dist(as.matrix(aitchison_dist)[idx_gwp, idx_gwp])
meta_gwp$substrate_group <- ifelse(meta_gwp$substrate_grouped == "Glass",
                                   "Glass", "Weathered_Plastic")
meta_gwp$substrate_group <- as.factor(meta_gwp$substrate_group)
cat("Sample sizes:\n"); print(table(meta_gwp$substrate_group))
set.seed(123)
perm_gwp <- adonis2(dist_gwp ~ substrate_group, data = meta_gwp, permutations = 999)
print(perm_gwp)
set.seed(123)
disp_gwp      <- betadisper(dist_gwp, meta_gwp$substrate_group)
perm_disp_gwp <- permutest(disp_gwp, permutations = 999)
print(perm_disp_gwp)

# --- Glass vs All Plastic ---
cat("\n--- Glass vs All Plastic ---\n")

idx_gap  <- metadata$substrate_grouped %in% c("Glass", "PE", "PET", "Weathered_PE", "Weathered_PET")
meta_gap <- metadata[idx_gap, ]
dist_gap <- as.dist(as.matrix(aitchison_dist)[idx_gap, idx_gap])
meta_gap$substrate_group <- ifelse(meta_gap$substrate_grouped == "Glass",
                                   "Glass", "Plastic")
meta_gap$substrate_group <- as.factor(meta_gap$substrate_group)
cat("Sample sizes:\n"); print(table(meta_gap$substrate_group))
set.seed(123)
perm_gap <- adonis2(dist_gap ~ substrate_group, data = meta_gap, permutations = 999)
print(perm_gap)
set.seed(123)
disp_gap      <- betadisper(dist_gap, meta_gap$substrate_group)
perm_disp_gap <- permutest(disp_gap, permutations = 999)
print(perm_disp_gap)

# --- New Plastic vs Weathered Plastic ---
cat("\n--- New Plastic vs Weathered Plastic ---\n")

idx_nwp  <- metadata$substrate_grouped %in% c("PE", "PET", "Weathered_PE", "Weathered_PET")
meta_nwp <- metadata[idx_nwp, ]
dist_nwp <- as.dist(as.matrix(aitchison_dist)[idx_nwp, idx_nwp])
meta_nwp$substrate_group <- ifelse(meta_nwp$substrate_grouped %in% c("PE", "PET"),
                                   "New_Plastic", "Weathered_Plastic")
meta_nwp$substrate_group <- as.factor(meta_nwp$substrate_group)
cat("Sample sizes:\n"); print(table(meta_nwp$substrate_group))
set.seed(123)
perm_nwp <- adonis2(dist_nwp ~ substrate_group, data = meta_nwp, permutations = 999)
print(perm_nwp)
set.seed(123)
disp_nwp      <- betadisper(dist_nwp, meta_nwp$substrate_group)
perm_disp_nwp <- permutest(disp_nwp, permutations = 999)
print(perm_disp_nwp)

# ==============================================================================
# STEP 20: SAVE ALL OBJECTS
# ==============================================================================

cat("\n=== STEP 20: Saving All Objects ===\n")

models <- list(
  dataset_name         = dataset_name,
  rda                  = my.rda,
  aitchison_dist       = aitchison_dist,
  rclr_data            = al.rclr,
  pca_data             = pca_data,
  pca_centroids        = pca_centroids,
  metadata             = metadata,
  var_explained        = var_explained,
  perm_full            = perm_full,
  perm_full_byterms    = perm_full_byterms,
  perm_additive        = perm_additive,
  perm_twoway          = perm_twoway,
  perm_sequential      = perm_sequential,
  perm_site            = perm_site,
  perm_season          = perm_season,
  perm_substrate       = perm_substrate,
  disp_site            = disp_site,
  disp_season          = disp_season,
  disp_substrate       = disp_substrate,
  perm_disp_site       = perm_disp_site,
  perm_disp_season     = perm_disp_season,
  perm_disp_substrate  = perm_disp_substrate,
  pairwise_site        = pairwise_site,
  pairwise_season      = pairwise_season,
  pairwise_substrate   = pairwise_substrate,
  season_results       = season_results_df,
  site_season_results  = site_season_results,
  context_pairwise     = list(
    SELVA_Fall    = pw_SELVA_Fall,
    TBS_Fall      = pw_TBS_Fall,
    SELVA_Winter2 = pw_SELVA_Winter2,
    TBS_Winter2   = pw_TBS_Winter2
  ),
  grouped_substrate = list(
    glass_vs_new_plastic       = list(perm = perm_gnp, disp = disp_gnp, perm_disp = perm_disp_gnp),
    glass_vs_weathered_plastic = list(perm = perm_gwp, disp = disp_gwp, perm_disp = perm_disp_gwp),
    glass_vs_all_plastic       = list(perm = perm_gap, disp = disp_gap, perm_disp = perm_disp_gap),
    new_vs_weathered_plastic   = list(perm = perm_nwp, disp = disp_nwp, perm_disp = perm_disp_nwp)
  )
)

saveRDS(models, file = paste0(tab_dir, "models_", dataset_name, ".rds"))
cat("  Models saved:", paste0(tab_dir, "models_", dataset_name, ".rds"), "\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

