################################################################################
# COMBINED ST + WF ANALYSIS
# Joint ordination and source tracking
################################################################################

# STRUCTURE:
#   PART 0:   Data preparation
#   PART 0.5: Data exploration & dataset characteristics  ← NEW
#   PART 1:   Relative abundance visualisation (combined biofilm + water)
#   PART 2:   Build combined phyloseq
#   PART 3:   Joint ordination (rclr PCA + Aitchison distance)
#   PART 4:   Global PERMANOVA (Water vs Biofilm)
#   PART 5:   Procrustes analysis
#   PART 6:   Source tracking (FEAST)
################################################################################

library(phyloseq)
library(vegan)
library(ggplot2)
library(ggVennDiagram)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(knitr)
library(kableExtra)
library(FEAST)
library(glmmTMB)
library(emmeans)
library(car)
library(lme4)
library(lmerTest)
library(pals)
library(colorspace)
library(reshape2)

################################################################################
########################## PATH CONFIGURATION ##################################
################################################################################

GITHUB_DIR   <- "~/GitHub/MAPLE_Seasonal_Plastisphere"
SCRIPTS_DIR  <- file.path(GITHUB_DIR, "Scripts", "03.StatisticalAnalysis")
PROC_DIR     <- file.path(GITHUB_DIR, "Processed_data")

PS_DIR       <- file.path(PROC_DIR, "Phyloseq_objects")
FEAST_DIR    <- file.path(PROC_DIR, "FEAST_results")
FEAST_CONC   <- file.path(FEAST_DIR, "concurrent")
FEAST_LAG    <- file.path(FEAST_DIR, "lagged")
FEAST_SUB    <- file.path(FEAST_DIR, "per_substrate")
FEAST_COLL   <- file.path(FEAST_DIR, "collapsed")

## Create directories if they don't exist
invisible(lapply(
  c(PS_DIR, FEAST_DIR, FEAST_CONC, FEAST_LAG, FEAST_SUB, FEAST_COLL),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

setwd(SCRIPTS_DIR)

# ==============================================================================
# PART 0: DATA PREPARATION
# ==============================================================================

ps_all <- readRDS(file.path(PS_DIR, "FINAL_PHYLOSEQ_OBJECTS.rds"))

# 1. Pull both raw objects from the same source
ps_st_raw <- ps_all$pr2_collap_rep50$fP.lulu$ST
ps_wf_raw <- ps_all$pr2_collap_rep50$fP.lulu$Water




# ==============================================================================
# Goal: test whether water communities differ between deployment and retrieval
# ==============================================================================

# --- Build retrieval and deployment subsets ---
wf_meta <- as.data.frame(sample_data(ps_wf_raw))
class(wf_meta) <- "data.frame"


ps_wf_retrieval  <- prune_samples(!is.na(wf_meta$st_role) & wf_meta$st_role == "Retrieval",  ps_wf_raw)
ps_wf_deployment <- prune_samples(!is.na(wf_meta$st_role) & wf_meta$st_role == "Deployment", ps_wf_raw)

ps_wf_retrieval  <- prune_taxa(taxa_sums(ps_wf_retrieval)  > 0, ps_wf_retrieval)
ps_wf_deployment <- prune_taxa(taxa_sums(ps_wf_deployment) > 0, ps_wf_deployment)

cat("Retrieval samples:", nsamples(ps_wf_retrieval),  "| taxa:", ntaxa(ps_wf_retrieval),  "\n")
cat("Deployment samples:", nsamples(ps_wf_deployment), "| taxa:", ntaxa(ps_wf_deployment), "\n")



# --- Merge into one object for joint ordination ---
sample_data(ps_wf_retrieval)$st_role  <- "Retrieval"
sample_data(ps_wf_deployment)$st_role <- "Deployment"

ps_wf_both <- merge_phyloseq(ps_wf_retrieval, ps_wf_deployment)
ps_wf_both <- prune_taxa(taxa_sums(ps_wf_both) > 0, ps_wf_both)

meta_wf_both <- as.data.frame(sample_data(ps_wf_both))
class(meta_wf_both) <- "data.frame"

cat("Combined:", nsamples(ps_wf_both), "samples |", ntaxa(ps_wf_both), "taxa\n")


# --- rclr transformation + Aitchison distance ---
mat_wf <- as(otu_table(ps_wf_both), "matrix")
if (taxa_are_rows(ps_wf_both)) mat_wf <- t(mat_wf)

rclr_wf <- decostand(mat_wf, method = "rclr", MARGIN = 1, na.rm = TRUE)
rownames(rclr_wf) <- rownames(mat_wf)  # preserve names immediately after decostand

aitchison_wf <- dist(rclr_wf, method = "euclidean")

# --- PCA ---
rda_wf  <- rda(rclr_wf)
var_exp <- round(100 * rda_wf$CA$eig[1:2] / rda_wf$tot.chi, 1)
cat("PC1:", var_exp[1], "%  PC2:", var_exp[2], "%\n")

pca_wf <- as.data.frame(scores(rda_wf, display = "sites", choices = 1:2))
colnames(pca_wf) <- c("PC1", "PC2")

# Add metadata by direct rowname matching
pca_wf$site    <- meta_wf_both[rownames(pca_wf), "site"]
pca_wf$season  <- meta_wf_both[rownames(pca_wf), "season"]
pca_wf$st_role <- meta_wf_both[rownames(pca_wf), "st_role"]

pca_wf$season  <- factor(pca_wf$season,
                         levels = c("Winter","Spring","Summer","Fall","Winter2"))
pca_wf$st_role <- factor(pca_wf$st_role,
                         levels = c("Deployment","Retrieval"))

# Verify before plotting
cat("PCA rows:", nrow(pca_wf), "\n")
table(pca_wf$site,    useNA = "always")
table(pca_wf$st_role, useNA = "always")

# --- PCA plot ---
p_wf_role <- ggplot(pca_wf, aes(x = PC1, y = PC2,
                                color = season, shape = st_role)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_shape_manual(values = c("Deployment" = 16, "Retrieval" = 17)) +
  facet_wrap(~ site) +
  labs(
    title    = "Aitchison PCA: Retrieval vs Deployment water samples",
    subtitle = "Circle = Deployment | Triangle = Retrieval | Colour = Season",
    x        = paste0("PC1 (", var_exp[1], "%)"),
    y        = paste0("PC2 (", var_exp[2], "%)"),
    color    = "Season",
    shape    = "Role"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgray"),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5)
  )
print(p_wf_role)

# --- PERMANOVA ---
perm_meta <- meta_wf_both[rownames(as.matrix(aitchison_wf)), ]

# Model 1: role alone
perm_role <- adonis2(aitchison_wf ~ st_role,
                     data         = perm_meta,
                     permutations = 999,
                     by           = "terms")
cat("\nPERMANOVA Model 1 — st_role only:\n")
print(perm_role)

# Model 2: site + season + role
# Asks: is role significant BEYOND seasonal and spatial structure?
perm_role_full <- adonis2(aitchison_wf ~ site + season + st_role,
                          data         = perm_meta,
                          permutations = 999,
                          by           = "terms")
cat("\nPERMANOVA Model 2 — site + season + st_role:\n")
print(perm_role_full)

# --- Dispersion test ---
# Checks whether groups differ in spread rather than location.
# Significant result = interpret PERMANOVA with caution.
bd <- betadisper(aitchison_wf, perm_meta$st_role)
cat("\nDispersion test (betadisper):\n")
print(permutest(bd, permutations = 999))

# ==============================================================================
# CONCLUSION: Deployment and retrieval water communities are not meaningfully
# different (role R² = 0.003, p = 0.027 after season + site; betadisper p = 0.374).
# Season dominates water community structure (R² = 0.94).
# We merge both roles into a single seasonal water source for downstream analysis.
# ==============================================================================

# --- Build merged water phyloseq (all ST water samples) ---
ps_wf_collap <- prune_samples(!is.na(wf_meta$st_role), ps_wf_raw)
ps_wf_collap <- prune_taxa(taxa_sums(ps_wf_collap) > 0, ps_wf_collap)

cat("Merged water (all ST):", nsamples(ps_wf_collap), 
    "samples |", ntaxa(ps_wf_collap), "taxa\n")

# Sanity check: seasons and roles present
meta_wf_collap <- as.data.frame(sample_data(ps_wf_collap))
class(meta_wf_collap) <- "data.frame"
cat("\nSamples per season:\n")
print(table(meta_wf_collap$season, useNA = "always"))
cat("\nSamples per role:\n")
print(table(meta_wf_collap$st_role, useNA = "always"))
cat("\nSamples per site x season:\n")
print(table(meta_wf_collap$site, meta_wf_collap$season, useNA = "always"))

# --- Tag sample types and merge with biofilm to share ASV namespace ---
sample_data(ps_st_raw)$sample_type    <- "Biofilm"
sample_data(ps_wf_collap)$sample_type <- "Water"

ps_combined <- merge_phyloseq(ps_st_raw, ps_wf_collap)
ps_combined <- prune_taxa(taxa_sums(ps_combined) > 0, ps_combined)
ps_combined <- prune_samples(sample_sums(ps_combined) > 0, ps_combined)

# Set factor levels for downstream analyses
sample_data(ps_combined)$sample_type <- factor(
  sample_data(ps_combined)$sample_type, levels = c("Water", "Biofilm"))
sample_data(ps_combined)$season <- factor(
  sample_data(ps_combined)$season,
  levels = c("Winter", "Spring", "Summer", "Fall", "Winter2"))
sample_data(ps_combined)$site <- as.factor(sample_data(ps_combined)$site)

# --- Split back by sample type ---
ps_st_collap <- subset_samples(ps_combined, sample_type == "Biofilm")
ps_st_collap <- prune_taxa(taxa_sums(ps_st_collap) > 0, ps_st_collap)

ps_wf_collap <- subset_samples(ps_combined, sample_type == "Water")
ps_wf_collap <- prune_taxa(taxa_sums(ps_wf_collap) > 0, ps_wf_collap)

cat("\nFinal objects:\n")
cat("Biofilm:", ntaxa(ps_st_collap), "taxa,", nsamples(ps_st_collap), "samples\n")
cat("Water:  ", ntaxa(ps_wf_collap), "taxa,", nsamples(ps_wf_collap), "samples\n")

# --- Save ---
# ST water samples
saveRDS(ps_wf_collap,
        file.path(PS_DIR, "ps_wf_collap_rep50_allST.rds"))
# ST water + biofilm
saveRDS(ps_combined,
        file.path(PS_DIR, "ps_st_wf_allST_combined.rds"))




# ==============================================================================
# PART 0.5: DATA EXPLORATION & DATASET CHARACTERISTICS
# ==============================================================================
# Goal: understand the structure, size, and balance of your dataset BEFORE
# running any statistics
# ==============================================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PART 0.5: DATASET CHARACTERISTICS\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")


#General characterisitics

#OTUs
ntaxa(ps_st_collap) 
ntaxa(ps_wf_collap)

#Reads
sum(sample_sums(ps_st_collap))
sum(sample_sums(ps_wf_collap))

#unassigned at Division level
table(tax_table(ps_st_collap)[, "Division"], useNA = "always")
table(tax_table(ps_wf_collap)[, "Division"], useNA = "always")


#number of divisions that are not NA
length(unique(na.omit(tax_table(ps_st_collap)[, "Division"])))
length(unique(na.omit(tax_table(ps_wf_collap)[, "Division"])))

#% of unassigned at Division
(table(tax_table(ps_st_collap)[, "Division"], useNA = "always")/ntaxa(ps_st_collap)) *100
(table(tax_table(ps_wf_collap)[, "Division"], useNA = "always")/ntaxa(ps_wf_collap)) *100


# ------------------------------------------------------------------------------
# 0.5A: Sample inventory — balanced design check
# ------------------------------------------------------------------------------
# A balanced design means every combination of your grouping variables has the
# same number of replicates. Imbalance is not fatal but affects some statistics.
# ------------------------------------------------------------------------------

cat("--- 0.5A: Sample inventory ---\n\n")

meta_st <- as.data.frame(sample_data(ps_st_collap))
class(meta_st) <- "data.frame"

meta_wf <- as.data.frame(sample_data(ps_wf_collap))
class(meta_wf) <- "data.frame"

# Biofilm: how many replicates per site × season × substrate?
inv_st <- meta_st %>%
  count(site, season, substrate, name = "n_replicates") %>%
  arrange(site, season, substrate)

cat("Biofilm sample inventory (site × season × substrate):\n")
print(kable(inv_st, format = "simple"))

# Check balance: are all cells equal?
n_expected_st <- median(inv_st$n_replicates)
unbalanced_st <- inv_st %>% filter(n_replicates != n_expected_st)
if (nrow(unbalanced_st) == 0) {
  cat("\n✓ Biofilm design is BALANCED — all cells have", n_expected_st, "replicates.\n\n")
} else {
  cat("\n⚠ Biofilm design is UNBALANCED. These cells deviate from expected (", n_expected_st, "):\n")
  print(unbalanced_st)
  cat("\n")
}

# Water: how many replicates per site × season?
inv_wf <- meta_wf %>%
  count(site, season, name = "n_replicates") %>%
  arrange(site, season)

cat("Water sample inventory (site × season):\n")
print(kable(inv_wf, format = "simple"))

n_expected_wf <- median(inv_wf$n_replicates)
unbalanced_wf <- inv_wf %>% filter(n_replicates != n_expected_wf)
if (nrow(unbalanced_wf) == 0) {
  cat("\n✓ Water design is BALANCED — all cells have", n_expected_wf, "replicates.\n\n")
} else {
  cat("\n⚠ Water design is UNBALANCED. These cells deviate from expected (", n_expected_wf, "):\n")
  print(unbalanced_wf)
  cat("\n")
}

# ------------------------------------------------------------------------------
# 0.5B: Read depth summary
# ------------------------------------------------------------------------------
# "Read depth" = total number of sequences per sample.
# Very low-depth samples (outliers) can distort distance-based analyses.
# Check for any samples that look suspicious before proceeding.
# ------------------------------------------------------------------------------

cat("--- 0.5B: Read depth (library size) ---\n\n")

#general summary
summary(sample_sums(ps_st_collap))
summary(sample_sums(ps_wf_collap))

reads_st <- data.frame(
  sample      = sample_names(ps_st_collap),
  sample_type = "Biofilm",
  reads       = sample_sums(ps_st_collap),
  site        = meta_st$site,
  season      = meta_st$season,
  substrate   = meta_st$substrate
)

reads_wf <- data.frame(
  sample      = sample_names(ps_wf_collap),
  sample_type = "Water",
  reads       = sample_sums(ps_wf_collap),
  site        = meta_wf$site,
  season      = meta_wf$season,
  substrate   = NA
)

reads_all <- bind_rows(reads_st, reads_wf)

# Overall summary by sample type
reads_summary <- reads_all %>%
  group_by(sample_type) %>%
  summarise(
    n_samples  = n(),
    min_reads  = min(reads),
    median_reads = median(reads),
    mean_reads = round(mean(reads)),
    max_reads  = max(reads),
    sd_reads   = round(sd(reads)),
    .groups = "drop"
  )

cat("Read depth summary by sample type:\n")
print(kable(reads_summary, format = "simple"))
cat("\n")

# Flag potential outliers: samples with reads < 10% of median within type
reads_all <- reads_all %>%
  group_by(sample_type) %>%
  mutate(low_depth_flag = reads < 0.10 * median(reads)) %>%
  ungroup()

low_depth <- reads_all %>% filter(low_depth_flag)
if (nrow(low_depth) == 0) {
  cat("✓ No extreme low-depth samples detected.\n\n")
} else {
  cat("⚠ Low-depth samples (< 10% of median for their type):\n")
  print(low_depth %>% 
          select(sample, sample_type, site, season, substrate, reads) %>%
          arrange(reads))
  cat("\n")
}

# Visual: read depth distribution
p_reads <- ggplot(reads_all, aes(x = reads, fill = sample_type)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~ sample_type, scales = "free") +
  scale_fill_manual(values = c("Biofilm" = "#D6604D", "Water" = "#4393C3")) +
  labs(
    title = "Read depth distribution per sample type",
    x     = "Reads per sample",
    y     = "Count",
    fill  = "Sample type"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")
print(p_reads)



# ------------------------------------------------------------------------------
# 0.5C: OTU richness (observed taxa) per sample
# ------------------------------------------------------------------------------
# "OTU richness" = number of taxa with at least 1 read in a sample.
# Low-richness samples may indicate failed extractions or sequencing issues.
# This is checked BEFORE rarefaction — raw counts only.
# ------------------------------------------------------------------------------
cat("--- 0.5C: OTU richness (observed taxa per sample) ---\n\n")

summary(rowSums(t(otu_table(ps_st_collap)) > 0))
summary(rowSums(t(otu_table(ps_wf_collap)) > 0))

# Count number of OTUs with > 0 reads per sample
otu_counts_st <- data.frame(
  sample      = sample_names(ps_st_collap),
  sample_type = "Biofilm",
  n_otus      = rowSums(t(otu_table(ps_st_collap)) > 0),  # t() transposes
  site        = meta_st$site,
  season      = meta_st$season,
  substrate   = meta_st$substrate
)

otu_counts_wf <- data.frame(
  sample      = sample_names(ps_wf_collap),
  sample_type = "Water",
  n_otus      = rowSums(t(otu_table(ps_wf_collap)) > 0),
  site        = meta_wf$site,
  season      = meta_wf$season,
  substrate   = NA
)

# Note: if your otu_table has taxa as columns, rowSums works directly.
# If taxa are rows (phyloseq default), you need to transpose first:
# rowSums(t(otu_table(ps_st_collap)) > 0)
# Check with: taxa_are_rows(ps_st_collap)

otu_all <- bind_rows(otu_counts_st, otu_counts_wf)

# Summary table
otu_summary <- otu_all %>%
  group_by(sample_type) %>%
  summarise(
    n_samples    = n(),
    min_otus     = min(n_otus),
    median_otus  = median(n_otus),
    mean_otus    = round(mean(n_otus)),
    max_otus     = max(n_otus),
    sd_otus      = round(sd(n_otus)),
    .groups = "drop"
  )

cat("OTU richness summary by sample type:\n")
print(kable(otu_summary, format = "simple"))
cat("\n")

# Flag low-richness samples: < 10% of median within type
otu_all <- otu_all %>%
  group_by(sample_type) %>%
  mutate(low_richness_flag = n_otus < 0.10 * median(n_otus)) %>%
  ungroup()

low_richness <- otu_all %>% filter(low_richness_flag)

if (nrow(low_richness) == 0) {
  cat("✓ No extreme low-richness samples detected.\n\n")
} else {
  cat("⚠ Low-richness samples (< 10% of median for their type):\n")
  print(low_richness %>%
          select(sample, sample_type, site, season, substrate, n_otus) %>%
          arrange(n_otus))
  cat("\n")
}

# Visual: OTU richness distribution
p_otus <- ggplot(otu_all, aes(x = n_otus, fill = sample_type)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~ sample_type, scales = "free") +
  scale_fill_manual(values = c("Biofilm" = "#D6604D", "Water" = "#4393C3")) +
  labs(
    title = "OTU richness distribution per sample type",
    x     = "Number of OTUs per sample",
    y     = "Count",
    fill  = "Sample type"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

print(p_otus)


# ------------------------------------------------------------------------------
# 0.5C: Taxa summary — total, shared, unique per substrate + seawater
# ------------------------------------------------------------------------------
cat("--- 0.5C: Taxa summary ---\n\n")

# 1. Build ASV sets per substrate group (union = ANY sample within group)
substrates <- c("PE", "PET", "Glass", "Weathered_PE", "Weathered_PET")
meta_st <- as.data.frame(sample_data(ps_st_collap))

substrate_taxa <- lapply(substrates, function(sub) {
  keep_samples <- rownames(meta_st)[meta_st$substrate == sub]
  ps_sub <- prune_samples(keep_samples, ps_st_collap)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  taxa_names(ps_sub)
})
names(substrate_taxa) <- substrates

# Add seawater
substrate_taxa[["Seawater"]] <- taxa_names(
  prune_taxa(taxa_sums(ps_wf_collap) > 0, ps_wf_collap)
)

# 2. High-level overlap table (biofilm vs seawater)
biofilm_asvs  <- Reduce(union, substrate_taxa[1:5])
seawater_asvs <- substrate_taxa[["Seawater"]]
all_asvs      <- union(biofilm_asvs, seawater_asvs)

taxa_overlap <- data.frame(
  Category  = c("Total unique (Biofilm + Seawater)", "Shared", "Biofilm only", "Seawater only"),
  n_taxa    = c(
    length(all_asvs),
    length(intersect(biofilm_asvs, seawater_asvs)),
    length(setdiff(biofilm_asvs, seawater_asvs)),
    length(setdiff(seawater_asvs, biofilm_asvs))
  )
)
taxa_overlap$pct_of_total <- round(100 * taxa_overlap$n_taxa / length(all_asvs), 1)

cat("High-level taxa overlap (Biofilm vs Seawater):\n")
print(kable(taxa_overlap, format = "simple"))
cat("\n")

# 3. Pairwise shared ASV matrix across all 6 groups
groups    <- names(substrate_taxa)
count_mat <- outer(groups, groups, Vectorize(function(a, b) {
  length(intersect(substrate_taxa[[a]], substrate_taxa[[b]]))
}))
rownames(count_mat) <- groups
colnames(count_mat) <- groups

pct_mat <- round(100 * count_mat / diag(count_mat), 1)

display_mat <- matrix(
  paste0(count_mat, " (", pct_mat, "%)"),
  nrow = length(groups),
  dimnames = list(groups, groups)
)

# Extra row: ASVs shared across ALL 6 groups
core_asvs <- length(Reduce(intersect, substrate_taxa))
core_row  <- paste0(core_asvs, " (",
                    round(100 * core_asvs / sapply(groups, function(g) length(substrate_taxa[[g]])), 1),
                    "%)")
names(core_row) <- groups
display_mat <- rbind(display_mat, "Shared across all" = core_row)

cat("Pairwise ASV overlap (n shared, % of row group):\n")
print(kable(display_mat, format = "simple"))
cat("\n")

# 4. Simple 2-circle Venn: Biofilm vs Seawater
venn_list <- list(Biofilm = biofilm_asvs, Seawater = seawater_asvs)

venn_taxa <- ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "#E8F4F8", high = "#2196F3") +
  scale_color_manual(values = c("black", "black")) +
  labs(title = "ASV overlap: Biofilm vs Seawater") +
  theme(legend.position = "none")
print(venn_taxa)



# Prevalence: how many samples does each taxon appear in?

# --- Group-level prevalence ---

# 1. Extract OTU tables and sample data
otu_st <- as(otu_table(ps_st_collap), "matrix")
otu_wf <- as(otu_table(ps_wf_collap), "matrix")

# Make sure taxa are rows
if (!taxa_are_rows(ps_st_collap)) otu_st <- t(otu_st)
if (!taxa_are_rows(ps_wf_collap)) otu_wf <- t(otu_wf)

sd_st <- sample_data(ps_st_collap) %>% as("data.frame") %>%
  mutate(group = paste(site, season, substrate, sep = "_"))

sd_wf <- sample_data(ps_wf_collap) %>% as("data.frame") %>%
  mutate(group = paste(site, season, substrate, sep = "_"))

# 2. For each taxon, check presence/absence per group
# A taxon is "present in a group" if it appears in >= 1 sample of that group
# group = seaosn x site x substrate combination

group_prev <- function(otu_mat, sample_df) {
  groups <- unique(sample_df$group)
  sapply(rownames(otu_mat), function(taxon) {
    sum(sapply(groups, function(g) {
      samples_in_group <- rownames(sample_df)[sample_df$group == g]
      any(otu_mat[taxon, samples_in_group] > 0)
    }))
  })
}

prev_st_group <- group_prev(otu_st, sd_st)
prev_wf_group <- group_prev(otu_wf, sd_wf)

n_groups_st <- length(unique(sd_st$group))
n_groups_wf <- length(unique(sd_wf$group))

# 3. Summary table
prev_group_summary <- data.frame(
  sample_type            = c("Biofilm", "Water"),
  n_groups               = c(n_groups_st, n_groups_wf),
  n_taxa_total           = c(length(prev_st_group), length(prev_wf_group)),
  n_in_one_group_only    = c(sum(prev_st_group == 1), sum(prev_wf_group == 1)),
  pct_one_group_only     = round(100 * c(sum(prev_st_group == 1) / length(prev_st_group),
                                         sum(prev_wf_group == 1) / length(prev_wf_group)), 1),
  n_in_all_groups        = c(sum(prev_st_group == n_groups_st),
                             sum(prev_wf_group == n_groups_wf)),
  median_group_prevalence = c(median(prev_st_group), median(prev_wf_group))
)

cat("Group-level prevalence (site × season × substrate):\n")
print(kable(prev_group_summary, format = "simple"))


# ------------------------------------------------------------------------------
# 0.5D: Relative abundance summary at Phylum level
# ------------------------------------------------------------------------------
# Quick sanity check: what are the dominant groups in each sample type?
# Helps catch any obvious contamination or labelling errors.
# ------------------------------------------------------------------------------

cat("--- 0.5D: Top 20 phyla (divison) by mean relative abundance ---\n\n")

summarise_top_phyla <- function(ps, label, n = 20) {
  
  # Determine phylum rank name (case-insensitive)
  div_rank <- rank_names(ps)[grep("^Division$", rank_names(ps), ignore.case = TRUE)]
  
  ps_div <- tax_glom(ps, taxrank = div_rank, NArm = FALSE)
  ps_rel    <- transform_sample_counts(ps_div, function(x) x / sum(x) * 100)
  df        <- psmelt(ps_rel)
  
  top <- df %>%
    rename(Division = all_of(div_rank)) %>%
    mutate(Division = ifelse(is.na(Division) | Division == "", "Unclassified", Division)) %>%
    group_by(Division) %>%
    summarise(mean_rel_abund = round(mean(Abundance), 2),
              sd_rel_abund   = round(sd(Abundance), 2),
              .groups = "drop") %>%
    arrange(desc(mean_rel_abund)) %>%
    slice_head(n = n) %>%
    mutate(sample_type = label)
  
  return(top)
}

phyla_st <- summarise_top_phyla(ps_st_collap, "Biofilm")
phyla_wf <- summarise_top_phyla(ps_wf_collap, "Water")

cat("Biofilm — top phyla:\n")
print(kable(phyla_st %>% select(-sample_type), format = "simple"))
cat("\nWater — top phyla:\n")
print(kable(phyla_wf %>% select(-sample_type), format = "simple"))
cat("\n")

# ------------------------------------------------------------------------------
# 0.5E: Per-season shared taxa 
# ------------------------------------------------------------------------------

cat("--- 0.5E: Shared taxa by site × season ---\n\n")

# Build combined object temporarily for this check only
ps_temp_combined <- merge_phyloseq(ps_st_collap, ps_wf_collap)

meta_temp <- as.data.frame(sample_data(ps_temp_combined))
class(meta_temp) <- "data.frame"

seasons_present <- unique(na.omit(as.character(meta_temp$season)))
seasons_present <- seasons_present[seasons_present %in% 
                                     c("Winter","Spring","Summer","Fall","Winter2")]
seasons_present <- factor(seasons_present, 
                          levels = c("Winter","Spring","Summer","Fall","Winter2"))
seasons_present <- levels(seasons_present)[levels(seasons_present) %in% seasons_present]

sites_present   <- levels(as.factor(meta_temp$site))

shared_taxa_table <- expand.grid(site = sites_present, season = seasons_present,
                                 stringsAsFactors = FALSE)
shared_taxa_table[, c("n_ST","n_WF","n_shared","pct_ST_shared")] <- NA

for (i in seq_len(nrow(shared_taxa_table))) {
  s    <- shared_taxa_table$site[i]
  seas <- shared_taxa_table$season[i]
  
  ps_sub_st <- prune_samples(
    sample_data(ps_temp_combined)$site        == s &
      sample_data(ps_temp_combined)$season      == seas &
      sample_data(ps_temp_combined)$sample_type == "Biofilm",
    ps_temp_combined)
  
  ps_sub_wf <- prune_samples(
    sample_data(ps_temp_combined)$site        == s &
      sample_data(ps_temp_combined)$season      == seas &
      sample_data(ps_temp_combined)$sample_type == "Water",
    ps_temp_combined)
  
  if (nsamples(ps_sub_st) == 0 || nsamples(ps_sub_wf) == 0) next
  
  t_st <- taxa_names(prune_taxa(taxa_sums(ps_sub_st) > 0, ps_sub_st))
  t_wf <- taxa_names(prune_taxa(taxa_sums(ps_sub_wf) > 0, ps_sub_wf))
  
  shared_taxa_table$n_ST[i]          <- length(t_st)
  shared_taxa_table$n_WF[i]          <- length(t_wf)
  shared_taxa_table$n_shared[i]      <- length(intersect(t_st, t_wf))
  shared_taxa_table$pct_ST_shared[i] <- round(100 * length(intersect(t_st, t_wf)) /
                                                length(t_st), 1)
}

cat("Shared taxa by site × season:\n")
print(kable(shared_taxa_table, format = "simple"))

rm(ps_temp_combined)  # clean up; will be rebuilt properly in Part 2

## Heatmap of shared taxa

# ------------------------------------------------------------------------------
# Helper: build pairwise % shared matrix (relative to smaller group)
# ------------------------------------------------------------------------------
make_pct_matrix <- function(taxa_list) {
  groups <- names(taxa_list)
  mat <- outer(groups, groups, Vectorize(function(a, b) {
    n_shared  <- length(intersect(taxa_list[[a]], taxa_list[[b]]))
    n_smaller <- min(length(taxa_list[[a]]), length(taxa_list[[b]]))
    round(100 * n_shared / n_smaller, 1)
  }))
  rownames(mat) <- groups
  colnames(mat) <- groups
  mat
}

# ------------------------------------------------------------------------------
# HEATMAP 1: Site × Season
# ------------------------------------------------------------------------------
meta_st <- as.data.frame(sample_data(ps_st_collap))
meta_wf <- as.data.frame(sample_data(ps_wf_collap))

season_levels <- c("Winter", "Spring", "Summer", "Fall", "Winter2")

# Build ASV sets per site × season for biofilm
site_season_combos <- unique(meta_st[, c("site", "season")])

ss_taxa <- apply(site_season_combos, 1, function(row) {
  keep   <- rownames(meta_st)[meta_st$site == row["site"] &
                                meta_st$season == row["season"]]
  ps_sub <- prune_samples(keep, ps_st_collap)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  taxa_names(ps_sub)
})
names(ss_taxa) <- paste(site_season_combos$site,
                        site_season_combos$season, sep = "_")

# Add seawater per site × season
wf_combos <- unique(meta_wf[, c("site", "season")])
for (i in seq_len(nrow(wf_combos))) {
  s    <- wf_combos$site[i]
  seas <- wf_combos$season[i]
  keep <- rownames(meta_wf)[meta_wf$site == s & meta_wf$season == seas]
  ps_sub <- prune_samples(keep, ps_wf_collap)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  ss_taxa[[paste0("WF_", s, "_", seas)]] <- taxa_names(ps_sub)
}

# Define ordered group levels: site, then sample type, then season
group_order <- c(
  paste("SELVA", season_levels, sep = "_"),
  paste("TBS",   season_levels, sep = "_"),
  paste0("WF_SELVA_", season_levels),
  paste0("WF_TBS_",   season_levels)
)
group_order <- group_order[group_order %in% names(ss_taxa)]

mat1 <- make_pct_matrix(ss_taxa)
df1  <- melt(mat1, varnames = c("Group1", "Group2"), value.name = "pct_shared")
df1$Group1 <- factor(df1$Group1, levels = rev(group_order))
df1$Group2 <- factor(df1$Group2, levels = group_order)

p_heatmap_ss <- ggplot(df1, aes(x = Group2, y = Group1, fill = pct_shared)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#EFF3FF", high = "#2171B5",
                      name = "% shared\n(smaller group)") +
  labs(title = "ASV sharing: Site × Season", x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        plot.title   = element_text(hjust = 0.5, face = "bold"))

print(p_heatmap_ss)

# ------------------------------------------------------------------------------
# HEATMAP 2: Substrate type
# ------------------------------------------------------------------------------
substrate_order <- c("Seawater", "Glass", "PE", "Weathered_PE", 
                     "PET", "Weathered_PET")
substrate_order <- substrate_order[substrate_order %in% names(substrate_taxa)]

mat2 <- make_pct_matrix(substrate_taxa)
df2  <- melt(mat2, varnames = c("Group1", "Group2"), value.name = "pct_shared")
df2$Group1 <- factor(df2$Group1, levels = rev(substrate_order))
df2$Group2 <- factor(df2$Group2, levels = substrate_order)

p_heatmap_sub <- ggplot(df2, aes(x = Group2, y = Group1, fill = pct_shared)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#EFF3FF", high = "#2171B5",
                      name = "% shared\n(smaller group)") +
  labs(title = "ASV sharing: Substrate type", x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        plot.title   = element_text(hjust = 0.5, face = "bold"))

print(p_heatmap_sub)

# ==============================================================================
# PART 1: RELATIVE ABUNDANCE VISUALISATION
# ==============================================================================
# Water samples are shown as an extra "substrate" column in the same plot,
# separated from biofilm substrates by a vertical line.
# Deployment (reps 1-3) and Retrieval (reps 4-6) water samples are shown
# separately within the Seawater group, divided by a dashed line.
# ==============================================================================

# ------------------------------------------------------------------------------
# Helper: aggregate to chosen rank and compute relative abundance
# Returns a long data frame ready for plotting
# ------------------------------------------------------------------------------

prepare_taxbar_base <- function(ps, top_n = 30, rank = "Family") {
  
  ps_rank <- tax_glom(ps, taxrank = rank, NArm = FALSE)
  ps_rel  <- transform_sample_counts(ps_rank, function(x) x / sum(x) * 100)
  df      <- psmelt(ps_rel)
  
  df <- df %>%
    rename_with(~ "Taxon", .cols = all_of(rank)) %>%
    mutate(
      Taxon  = case_when(is.na(Taxon) | Taxon == "" ~ "Unclassified", TRUE ~ Taxon),
      season = factor(season, levels = c("Winter","Spring","Summer","Fall","Winter2"))
    )
  
  return(df)
}

# ------------------------------------------------------------------------------
# Combined plot: Water as extra "substrate" columns
# Deployment water = reps 1-3, Retrieval water = reps 4-6
# ------------------------------------------------------------------------------

prepare_taxbar_combined <- function(ps_st, ps_wf, top_n = 30, rank = "Family") {
  
  df_st <- prepare_taxbar_base(ps_st, top_n, rank = rank)
  df_wf <- prepare_taxbar_base(ps_wf, top_n, rank = rank)
  
  # Label water substrate
  df_wf$substrate <- "Seawater"
  
  # Renumber water bio_replicate so deployment = 1-3, retrieval = 4-6
  # This prevents overlap on x-axis when both roles are present
  df_wf <- df_wf %>%
    mutate(bio_replicate = case_when(
      st_role == "Deployment" ~ bio_replicate,        # keep 1-3
      st_role == "Retrieval"  ~ bio_replicate + 3L,   # shift to 4-6
      TRUE                    ~ bio_replicate
    ))
  
  # Combine — only keep columns present in both
  keep_cols <- intersect(colnames(df_st), colnames(df_wf))
  df_all    <- bind_rows(df_st[, keep_cols], df_wf[, keep_cols])
  
  df_all <- df_all %>%
    mutate(Taxon = gsub("^X(?=[0-9])", "", Taxon, perl = TRUE))
  
  # Select top N taxa by mean abundance across all samples
  top_taxa <- df_all %>%
    group_by(Taxon) %>%
    summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n) %>%
    pull(Taxon)
  
  df_all <- df_all %>%
    mutate(Taxon_plot = if_else(Taxon %in% top_taxa, Taxon, "Other"))
  
  lvls <- c(setdiff(top_taxa, c("Other","Unclassified")), "Other","Unclassified")
  df_all$Taxon_plot <- factor(df_all$Taxon_plot, levels = rev(lvls))
  
  # --- Build x-axis levels ---
  # Biofilm: 5 substrates x 3 reps = 15 bars
  # Water:   6 reps (3 deployment + 3 retrieval) = 6 bars
  biofilm_subs   <- c("Glass","PE","Weathered_PE","PET","Weathered_PET")
  substrate_levels <- c(biofilm_subs, "Seawater")
  n_water_reps   <- max(df_wf$bio_replicate, na.rm = TRUE)  # should be 6
  
  df_all <- df_all %>%
    mutate(
      substrate = factor(substrate, levels = substrate_levels),
      x_label   = paste0(substrate, "_", bio_replicate),
      x_display = case_when(
        # Biofilm: label at middle rep (rep 2)
        substrate != "Seawater" & bio_replicate == 2 ~ as.character(substrate),
        # Water: label deployment at rep 2, retrieval at rep 5
        substrate == "Seawater" & bio_replicate == 2 ~ "Seawater\nDeployment",
        substrate == "Seawater" & bio_replicate == 5 ~ "Seawater\nRetrieval",
        TRUE ~ ""
      )
    )
  
  correct_levels <- c(
    as.vector(sapply(biofilm_subs, function(s) paste0(s, "_", 1:3))),
    paste0("Seawater_", 1:n_water_reps)
  )
  df_all$x_label <- factor(df_all$x_label, levels = correct_levels)
  
  label_map <- setNames(df_all$x_display, df_all$x_label)
  label_map <- label_map[!duplicated(names(label_map))]
  
  return(list(
    df           = df_all,
    label_map    = label_map,
    rank         = rank,
    n_water_reps = n_water_reps
  ))
}

# ------------------------------------------------------------------------------
# Plot function
# ------------------------------------------------------------------------------

make_taxbar_combined <- function(df, label_map, rank = "Family",
                                 n_water_reps = 6,
                                 plot_title = NULL) {
  
  if (is.null(plot_title)) 
    plot_title <- paste0(rank, "-level taxonomy: Biofilm vs Water")
  
  n_colors <- length(levels(df$Taxon_plot)) - 2
  pal <- setNames(
    c(pals::polychrome(n_colors), "#D3D3D3", "#FFFFFF"),
    c(levels(df$Taxon_plot)[!levels(df$Taxon_plot) %in% c("Other","Unclassified")],
      "Other", "Unclassified")
  )
  
  # Vertical separators:
  # Between biofilm substrates (every 3 bars): 3.5, 6.5, 9.5, 12.5
  # Between biofilm and water: 15.5
  # Between deployment and retrieval water: 18.5
  biofilm_separators  <- seq(3.5, by = 3, length.out = 4)
  water_separator     <- 15.5
  water_mid_separator <- 15.5 + 3  # after 3 deployment reps
  
  ggplot(df, aes(x = x_label, y = Abundance, fill = Taxon_plot)) +
    geom_bar(stat = "identity", position = "stack", width = 0.85, color = NA) +
    # Between biofilm substrates
    geom_vline(xintercept = biofilm_separators,
               color = "grey40", linewidth = 0.4, linetype = "dashed") +
    # Between biofilm and water
    geom_vline(xintercept = water_separator,
               color = "grey10", linewidth = 1.0) +
    # Between deployment and retrieval water
    geom_vline(xintercept = water_mid_separator,
               color = "grey40", linewidth = 0.4, linetype = "dashed") +
    facet_grid(season ~ site, scales = "free_x", space = "free_x") +
    scale_x_discrete(labels = label_map) +
    scale_fill_manual(values = pal, name = rank) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 101),
                       labels = scales::label_percent(scale = 1)) +
    labs(
      x = "Substrate / Sample type",
      y = "Relative Abundance (%)"
    ) +
    theme_bw(base_size = 15) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      strip.background = element_rect(fill = "#f0f0f0"),
      strip.text       = element_text(face = "bold", size = 9),
      panel.spacing.y  = unit(1.8, "lines"),
      panel.spacing    = unit(0.3, "lines"),
      legend.position  = "right",
      legend.key.size  = unit(0.4, "cm"),
      legend.text      = element_text(size = 9),
      legend.title     = element_text(size = 11, face = "bold")
    )
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

df_combined_taxbar <- prepare_taxbar_combined(
  ps_st_collap, ps_wf_collap, 
  top_n = 15, 
  rank  = "Division"
)

p_taxbar_combined <- make_taxbar_combined(
  df           = df_combined_taxbar$df,
  label_map    = df_combined_taxbar$label_map,
  rank         = df_combined_taxbar$rank,
  n_water_reps = df_combined_taxbar$n_water_reps
)

print(p_taxbar_combined)


# ==============================================================================
# PART 2: BUILD COMBINED PHYLOSEQ (collapsed)
# ==============================================================================
# merge_phyloseq fills absent taxa with zeros across datasets automatically.
# Always prune after merging to remove any ghost taxa/samples.
# ==============================================================================

cat("\n--- Building combined phyloseq object ---\n")

saveRDS(ps_combined,
        file.path(PS_DIR, "ps_st_wf_seasonal_combined.rds"))


# ==============================================================================
# PART 3: JOINT ORDINATION (rclr PCA + Aitchison distance)
# ==============================================================================
# rclr = robust centred log-ratio transformation. Unlike standard clr,
# it handles zeros without pseudocount addition, which makes it appropriate
# for sparse amplicon data.
# Aitchison distance = Euclidean distance on rclr-transformed data.
# PCA on rclr-transformed data is called Aitchison PCA or robust PCA.
# ==============================================================================

cat("\n--- Joint ordination (Aitchison PCA) ---\n")

meta_combined <- as.data.frame(sample_data(ps_combined))
class(meta_combined) <- "data.frame"

# --- rclr transformation ---
mat_combined <- as(otu_table(ps_combined), "matrix")
if (taxa_are_rows(ps_combined)) mat_combined <- t(mat_combined)
# Now rows = samples, columns = taxa

original_names  <- rownames(mat_combined)
rclr_combined   <- decostand(mat_combined, method = "rclr", MARGIN = 1, na.rm = TRUE)
rownames(rclr_combined) <- original_names

aitchison_combined <- dist(rclr_combined, method = "euclidean")

# --- PCA ---
rda_combined  <- rda(rclr_combined)
var_explained <- round(100 * rda_combined$CA$eig[1:2] / rda_combined$tot.chi, 1)
cat("PC1:", var_explained[1], "%  PC2:", var_explained[2], "%\n")

pca_scores <- as.data.frame(scores(rda_combined, display = "sites", choices = 1:2))
colnames(pca_scores) <- c("PC1","PC2")
pca_scores$sample_id <- original_names

meta_combined$sample_id <- rownames(meta_combined)
pca_data_combined <- merge(pca_scores, meta_combined, by = "sample_id")

pca_data_combined$season <- factor(pca_data_combined$season,
                                   levels = c("Winter","Spring","Summer","Fall","Winter2"))
pca_data_combined$sample_type <- factor(pca_data_combined$sample_type,
                                        levels = c("Water","Biofilm"))

# --- PCA plot 1: coloured by season, shaped by sample type, faceted by site ---
p_combined_site <- ggplot(pca_data_combined,
                          aes(x = PC1, y = PC2, color = season, shape = sample_type)) +
  geom_point(size = 2.5, alpha = 0.75) +
  scale_shape_manual(values = c("Water" = 14, "Biofilm" = 17)) +
  facet_wrap(~ site) +
  # labs(
  #   x        = paste0("PC1 (", var_explained[1], "%)"),
  #   y        = paste0("PC2 (", var_explained[2], "%)"),
  #   title    = "Combined Aitchison PCA: Biofilm vs Water",
  #   subtitle = "Circles = Water | Triangles = Biofilm | Colour = Season",
  #   color    = "Season", shape = "Sample type"
  # ) +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "lightgray"),
        plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
print(p_combined_site)

# --- PCA plot 2: coloured by sample type, faceted by season ---
p_combined_season <- ggplot(pca_data_combined,
                            aes(x = PC1, y = PC2, color = sample_type, shape = site)) +
  geom_point(size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c("Water" = "#4393C3", "Biofilm" = "#D6604D")) +
  facet_wrap(~ season, nrow = 1) +
  labs(
    x     = paste0("PC1 (", var_explained[1], "%)"),
    y     = paste0("PC2 (", var_explained[2], "%)"),
    title = "Combined Aitchison PCA — faceted by season",
    color = "Sample type", shape = "Site"
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_rect(fill = "lightgray"),
        plot.title = element_text(hjust = 0.5, face = "bold"))
print(p_combined_season)


# ==============================================================================
# PART 4: GLOBAL PERMANOVA — Water vs Biofilm
# ==============================================================================
# PERMANOVA (adonis2) partitions the variance in your distance matrix.
# by = "terms" tests each term sequentially — order matters for this method.
# The richer model (site + season + sample_type) asks: is sample_type still
# significant after accounting for spatial and temporal structure?
# ==============================================================================

cat("\n--- PERMANOVA: Water vs Biofilm ---\n")

meta_perm <- meta_combined[rownames(as.matrix(aitchison_combined)), ]

# Simple test: sample_type alone
perm_global <- adonis2(aitchison_combined ~ sample_type,
                       data = meta_perm, permutations = 999, by = "terms")
cat("\nModel 1 — sample_type only:\n")
print(perm_global)

# Richer model: account for site and season as covariates first
perm_full_type <- adonis2(aitchison_combined ~ site + season + sample_type,
                          data = meta_perm, permutations = 999, by = "terms")
cat("\nModel 2 — site + season + sample_type:\n")
print(perm_full_type)


# ==============================================================================
# PART 5: PROCRUSTES ANALYSIS
# ==============================================================================
# Tests whether the ordination structure of biofilm communities mirrors that
# of water communities (i.e., do sites/seasons that differ in water also
# differ similarly in biofilm?).
# We use centroids (site × season averages) to reduce noise from replicates.
# The Procrustes correlation and p-value (from protest) summarise the match.
# ==============================================================================

cat("\nProcrustes analysis: Biofilm vs Water ordinations\n")

rclr_df <- as.data.frame(rclr_combined, check.names = FALSE)
rclr_df$site        <- meta_combined[rownames(rclr_combined), "site"]
rclr_df$season      <- meta_combined[rownames(rclr_combined), "season"]
rclr_df$sample_type <- meta_combined[rownames(rclr_combined), "sample_type"]

taxa_cols <- setdiff(colnames(rclr_df), c("site","season","sample_type"))

# Compute centroids per sample_type × site × season
rclr_centroids <- rclr_df %>%
  group_by(site, season, sample_type) %>%
  summarise(across(all_of(taxa_cols), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")

centroid_labels     <- paste(rclr_centroids$sample_type,
                             rclr_centroids$site,
                             rclr_centroids$season, sep = "_")
rclr_centroid_mat   <- as.matrix(rclr_centroids[, taxa_cols])
rownames(rclr_centroid_mat) <- centroid_labels

rclr_wf_mat <- rclr_centroid_mat[grepl("^Water_",   rownames(rclr_centroid_mat)), ]
rclr_st_mat <- rclr_centroid_mat[grepl("^Biofilm_", rownames(rclr_centroid_mat)), ]
rownames(rclr_wf_mat) <- sub("^Water_",   "", rownames(rclr_wf_mat))
rownames(rclr_st_mat) <- sub("^Biofilm_", "", rownames(rclr_st_mat))

# Match order — Procrustes requires same row ordering
common_labels <- intersect(rownames(rclr_st_mat), rownames(rclr_wf_mat))
rclr_st_mat   <- rclr_st_mat[common_labels, ]
rclr_wf_mat   <- rclr_wf_mat[common_labels, ]

scores_wf <- scores(rda(rclr_wf_mat), display = "sites", choices = 1:2)
scores_st <- scores(rda(rclr_st_mat), display = "sites", choices = 1:2)

proc      <- procrustes(X = scores_st, Y = scores_wf, symmetric = TRUE)
proc_test <- protest(X = scores_st, Y = scores_wf, permutations = 999, symmetric = TRUE)
print(proc_test)

proc_scores <- data.frame(
  xW = proc$Yrot[, 1], yW = proc$Yrot[, 2],
  xB = proc$X[, 1],    yB = proc$X[, 2],
  label = rownames(proc$X)
)
proc_scores$site   <- sub("_.*", "",    proc_scores$label)
proc_scores$season <- sub("^[^_]*_", "", proc_scores$label)

p_procrustes <- ggplot(proc_scores) +
  geom_segment(aes(x = xW, y = yW, xend = xB, yend = yB),
               color = "grey60", linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(x = xW, y = yW, color = "Water"),   size = 3, shape = 16) +
  geom_point(aes(x = xB, y = yB, color = "Biofilm"), size = 3, shape = 17) +
  geom_label(aes(x = (xW + xB) / 2, y = (yW + yB) / 2,
                 label = paste(site, season, sep = "\n"), fill = site),
             size = 2.5, alpha = 0.6, linewidth = 0.2, color = "black") +
  scale_color_manual(values = c("Water" = "#4393C3", "Biofilm" = "#D6604D")) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    title    = "Procrustes: Biofilm vs Water centroids",
    subtitle = paste0("Correlation = ", round(sqrt(1 - proc_test$ss), 3),
                      "  p = ", proc_test$signif),
    x = "Dimension 1", y = "Dimension 2",
    color = "Sample type", fill = "Site"
  ) +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
print(p_procrustes)




# ==============================================================================
# PART 6: SOURCE TRACKING (FEAST)
# ==============================================================================
# STRUCTURE:
#   6.0  Setup        — merge phyloseq, build shared count matrix + metadata
#   6.1  Core function — run_feast_context() helper
#   6.2  Option A     — Concurrent (same season, all substrates collapsed)
#   6.3  Option A     — Lagged    (preceding season, all substrates collapsed)
#   6.4  Option B     — Concurrent, per substrate
#   6.5  Option C     — Fully collapsed (all replicates pooled → one-sample t-test)
#   6.6  Summarise    — summarise_feast_results() for Options A & B
#   6.7  Stats        — Beta regression models (Options A & B)
#   6.8  Export       — summary tables + save
#   6.9  Plots        — line plots and bar plots
# ==============================================================================

cat("\n--- Source Tracking with FEAST ---\n")


# ==============================================================================
# 6.0  SETUP
# ==============================================================================

ps_combined_raw <- ps_combined  # already built in Part 0

count_mat <- as(otu_table(ps_combined_raw), "matrix")
if (taxa_are_rows(ps_combined_raw)) count_mat <- t(count_mat)

meta_raw <- as.data.frame(sample_data(ps_combined_raw))
class(meta_raw) <- "data.frame"
meta_raw$season      <- factor(meta_raw$season,
                               levels = c("Winter","Spring","Summer","Fall","Winter2"))
meta_raw$site        <- as.factor(meta_raw$site)
meta_raw$sample_type <- as.character(meta_raw$sample_type)

season_order <- c("Winter","Spring","Summer","Fall","Winter2")

cat("FEAST input:", ntaxa(ps_combined_raw), "taxa,",
    nsamples(ps_combined_raw), "samples\n")


# ==============================================================================
# 6.1  CORE FUNCTION: run_feast_context()
# ==============================================================================
# Runs one FEAST call for a single context (site x season, optionally
# filtered to one substrate). Each call has:
#   Sources = WF samples at source_site x source_season
#   Sinks   = ST samples at sink_site x sink_season (x substrate if specified)
# Returns a named list with results + metadata, or NULL if skipped.
# ==============================================================================

run_feast_context <- function(count_mat, meta_raw,
                              sink_site, sink_season,
                              source_season,
                              substrate_filter = NULL,
                              dir_path, outfile_prefix) {
  
  # --- Identify sink samples ---
  sink_idx <- meta_raw$sample_type == "Biofilm" &
    as.character(meta_raw$site)   == sink_site   &
    as.character(meta_raw$season) == sink_season
  
  if (!is.null(substrate_filter))
    sink_idx <- sink_idx & as.character(meta_raw$substrate) == substrate_filter
  
  sink_samples <- rownames(meta_raw)[sink_idx]
  
  # --- Identify source samples ---
  source_idx <- meta_raw$sample_type == "Water" &
    as.character(meta_raw$site)   == sink_site    &
    as.character(meta_raw$season) == source_season
  
  source_samples <- rownames(meta_raw)[source_idx]
  
  # --- Skip if not enough samples ---
  if (length(sink_samples)   == 0) { cat("    Skipping: no sink samples\n");   return(NULL) }
  if (length(source_samples) == 0) { cat("    Skipping: no source samples\n"); return(NULL) }
  
  # --- Build count matrix for this context ---
  all_samples <- c(source_samples, sink_samples)
  C <- count_mat[all_samples, ]
  
  cat("    Before filtering: Samples:", nrow(C), "| Taxa:", ncol(C), "\n")
  
  C <- C[rowSums(C) > 0, ]
  C <- C[, colSums(C) > 0]
  
  # Drop samples below minimum read depth
  min_reads        <- 100
  low_read_samples <- rownames(C)[rowSums(C) < min_reads]
  if (length(low_read_samples) > 0) {
    cat("    Dropping low-read samples:", paste(low_read_samples, collapse = ", "), "\n")
    C <- C[rowSums(C) >= min_reads, ]
  }
  
  # Re-derive sample lists after filtering
  source_samples <- intersect(source_samples, rownames(C))
  sink_samples   <- intersect(sink_samples,   rownames(C))
  all_samples    <- c(source_samples, sink_samples)
  
  cat("    After filtering: Samples:", nrow(C), "| Taxa:", ncol(C),
      "| Sources:", length(source_samples), "| Sinks:", length(sink_samples), "\n")
  
  if (length(sink_samples) == 0 || length(source_samples) == 0) {
    cat("    Skipping: samples dropped to zero after filtering\n"); return(NULL)
  }
  if (ncol(C) < 2) {
    cat("    Skipping: fewer than 2 taxa after filtering\n"); return(NULL)
  }
  
  # --- Build FEAST metadata ---
  feast_meta <- data.frame(
    SourceSink = c(rep("Source", length(source_samples)),
                   rep("Sink",   length(sink_samples))),
    id         = c(rep(NA, length(source_samples)),
                   seq_len(length(sink_samples))),
    Env        = c(rep(paste0("Water_", sink_site, "_", source_season),
                       length(source_samples)),
                   rep("Biofilm", length(sink_samples))),
    stringsAsFactors = FALSE,
    row.names        = all_samples
  )
  
  # --- Run FEAST ---
  outfile <- paste0(outfile_prefix, "_", sink_site, "_", sink_season,
                    if (!is.null(substrate_filter)) paste0("_", substrate_filter) else "")
  
  set.seed(123)
  tryCatch({
    result <- FEAST(
      C                      = C,
      metadata               = feast_meta,
      different_sources_flag = 0,
      dir_path               = dir_path,
      outfile                = outfile
    )
    cat("    Done:", outfile, "\n")
    return(list(
      result         = result,
      sink_samples   = sink_samples,
      source_samples = source_samples,
      sink_site      = sink_site,
      sink_season    = sink_season,
      source_season  = source_season,
      substrate      = substrate_filter
    ))
  }, error = function(e) {
    cat("    Error in", outfile, ":", conditionMessage(e), "\n")
    return(NULL)
  })
}


# ==============================================================================
# 6.2  OPTION A — CONCURRENT (same season, all substrates collapsed)
# ==============================================================================
# Each FEAST call: water source = same site x same season
# Biofilm sinks  = all substrates pooled within that context
# ==============================================================================

cat("\n--- FEAST Option A: Concurrent (all substrates, same season) ---\n")

feast_results_concurrent <- list()

for (site_val in levels(meta_raw$site)) {
  for (season_val in season_order) {
    cat("  ", site_val, "-", season_val, "\n")
    key <- paste(site_val, season_val, sep = "_")
    feast_results_concurrent[[key]] <- run_feast_context(
      count_mat        = count_mat,
      meta_raw         = meta_raw,
      sink_site        = site_val,
      sink_season      = season_val,
      source_season    = season_val,
      substrate_filter = NULL,
      dir_path         = FEAST_CONC,
      outfile_prefix   = "FEAST_concurrent"
    )
  }
}


# ==============================================================================
# 6.3  OPTION A — LAGGED (preceding season, all substrates collapsed)
# ==============================================================================
# Asks whether the water community from the *previous* season predicts
# the current biofilm community. Winter has no preceding season so is skipped.
# ==============================================================================

cat("\n--- FEAST Option A: Lagged (all substrates, preceding season) ---\n")

feast_results_lagged <- list()

for (site_val in levels(meta_raw$site)) {
  for (season_idx in 2:length(season_order)) {
    
    sink_season   <- season_order[season_idx]
    source_season <- season_order[season_idx - 1]
    
    cat("  ", site_val, "- sink:", sink_season, "| source:", source_season, "\n")
    key <- paste(site_val, sink_season, "from", source_season, sep = "_")
    
    feast_results_lagged[[key]] <- run_feast_context(
      count_mat        = count_mat,
      meta_raw         = meta_raw,
      sink_site        = site_val,
      sink_season      = sink_season,
      source_season    = source_season,
      substrate_filter = NULL,
      dir_path         = FEAST_LAG,
      outfile_prefix   = "FEAST_lagged"
    )
  }
}


# ==============================================================================
# 6.4  OPTION B — CONCURRENT, PER SUBSTRATE
# ==============================================================================
# Same as Option A (concurrent) but each substrate gets its own FEAST call.
# Allows you to ask whether water contribution differs by substrate type.
# ==============================================================================

cat("\n--- FEAST Option B: Concurrent, per substrate ---\n")

feast_results_per_substrate <- list()

substrates <- unique(meta_raw$substrate[meta_raw$sample_type == "Biofilm"])
substrates <- substrates[!is.na(substrates)]

for (site_val in levels(meta_raw$site)) {
  for (season_val in season_order) {
    for (sub_val in substrates) {
      cat("  ", site_val, "-", season_val, "-", sub_val, "\n")
      key <- paste(site_val, season_val, sub_val, sep = "_")
      feast_results_per_substrate[[key]] <- run_feast_context(
        count_mat        = count_mat,
        meta_raw         = meta_raw,
        sink_site        = site_val,
        sink_season      = season_val,
        source_season    = season_val,
        substrate_filter = sub_val,
        dir_path         = FEAST_SUB,
        outfile_prefix   = "FEAST_substrate"
      )
    }
  }
}


# ==============================================================================
# 6.5  OPTION C — FULLY COLLAPSED FEAST + ONE-SAMPLE T-TEST
# ==============================================================================
# Motivation: the simplest defensible test of "does water significantly
# contribute to biofilm?", mirroring what the original FEAST paper (Shenhav
# et al. 2019) does.
#
# Approach:
#   Sinks   = ALL biofilm replicates pooled across all sites, seasons,
#             and substrates into a single count matrix
#   Sources = ALL water replicates pooled across the same scope
#
# This produces one water-attribution proportion per biofilm replicate.
# A one-sample two-sided t-test then asks:
#   H0: mean(water_prop) == 0   (water contributes nothing)
#   H1: mean(water_prop) != 0
#
# Limitation: pooling ignores site/season/substrate structure, so this
# cannot distinguish *where* or *when* water contribution is higher.
# It answers only the global question: "Is the contribution > 0?"
# ==============================================================================

cat("\n--- FEAST Option C: Fully collapsed (all replicates pooled) ---\n")

# --- Build pooled count matrix ---
all_biofilm_samples <- rownames(meta_raw)[meta_raw$sample_type == "Biofilm"]
all_water_samples   <- rownames(meta_raw)[meta_raw$sample_type == "Water"]

# Restrict water to seasonal samples only (same scope as Options A/B)
# They are already filtered in meta_raw via ps_wf_collap, so no extra step needed.

C_collapsed <- count_mat[c(all_water_samples, all_biofilm_samples), ]
C_collapsed <- C_collapsed[rowSums(C_collapsed) > 0, ]
C_collapsed <- C_collapsed[, colSums(C_collapsed) > 0]

# Drop low-read samples (same threshold as run_feast_context)
min_reads  <- 100
keep       <- rowSums(C_collapsed) >= min_reads
dropped    <- rownames(C_collapsed)[!keep]
if (length(dropped) > 0)
  cat("  Dropping low-read samples:", paste(dropped, collapse = ", "), "\n")
C_collapsed <- C_collapsed[keep, ]

# Re-derive sample lists after any drops
all_water_samples   <- intersect(all_water_samples,   rownames(C_collapsed))
all_biofilm_samples <- intersect(all_biofilm_samples, rownames(C_collapsed))
all_collapsed       <- c(all_water_samples, all_biofilm_samples)

cat("  Pooled sources (water):", length(all_water_samples),
    "| Pooled sinks (biofilm):", length(all_biofilm_samples), "\n")
cat("  Taxa after filtering:", ncol(C_collapsed), "\n")

# --- Build FEAST metadata ---
feast_meta_collapsed <- data.frame(
  SourceSink = c(rep("Source", length(all_water_samples)),
                 rep("Sink",   length(all_biofilm_samples))),
  id         = c(rep(NA, length(all_water_samples)),
                 seq_len(length(all_biofilm_samples))),
  Env        = c(rep("Water_all", length(all_water_samples)),
                 rep("Biofilm",   length(all_biofilm_samples))),
  stringsAsFactors = FALSE,
  row.names        = all_collapsed
)

# --- Run FEAST ---
set.seed(123)
feast_result_collapsed <- tryCatch({
  FEAST(
    C                      = C_collapsed,
    metadata               = feast_meta_collapsed,
    different_sources_flag = 0,
    dir_path               = FEAST_COLL,
    outfile                = "FEAST_collapsed_all"
  )
}, error = function(e) {
  cat("  Error in collapsed FEAST:", conditionMessage(e), "\n")
  NULL
})


# --- Extract per-replicate water proportions ---
if (!is.null(feast_result_collapsed)) {
  
  n_sinks_collapsed <- length(all_biofilm_samples)
  
  prop_elements_collapsed <- Filter(function(x) {
    is.numeric(x) && is.vector(x) && length(x) == n_sinks_collapsed
  }, feast_result_collapsed)
  
  if (length(prop_elements_collapsed) > 0) {
    
    prop_mat_collapsed           <- do.call(cbind, prop_elements_collapsed)
    colnames(prop_mat_collapsed) <- names(prop_elements_collapsed)
    rownames(prop_mat_collapsed) <- all_biofilm_samples
    
    water_cols_collapsed <- grep("Water", colnames(prop_mat_collapsed), value = TRUE)
    water_prop_collapsed <- if (length(water_cols_collapsed) > 0) {
      rowSums(prop_mat_collapsed[, water_cols_collapsed, drop = FALSE], na.rm = TRUE)
    } else {
      rep(NA, n_sinks_collapsed)
    }
    
    # Attach metadata for context
    sink_meta_collapsed <- meta_raw[all_biofilm_samples,
                                    c("site","season","substrate","bio_replicate")]
    
    replicate_df_collapsed <- data.frame(
      sample     = all_biofilm_samples,
      site       = sink_meta_collapsed$site,
      season     = sink_meta_collapsed$season,
      substrate  = sink_meta_collapsed$substrate,
      bio_rep    = sink_meta_collapsed$bio_replicate,
      water_prop = water_prop_collapsed,
      stringsAsFactors = FALSE
    )
    
    # --- One-sample two-sided t-test ---
    # H0: mean water proportion == 0
    # We use the raw (untransformed) proportions here because the t-test
    # is asking about the location of the distribution, not modelling structure.
    # Check normality assumption visually first.
    
    cat("\n  Distribution of water proportions (collapsed FEAST):\n")
    print(summary(water_prop_collapsed))
    
    hist(water_prop_collapsed,
         main = "Collapsed FEAST: water attribution per biofilm replicate",
         xlab = "Estimated proportion from water column",
         col  = "#4393C3", border = "white", breaks = 20)
    
    ttest_collapsed <- t.test(
      water_prop_collapsed,
      mu          = 0,
      alternative = "two.sided"
    )
    
    cat("\n  One-sample t-test (H0: mean water proportion = 0):\n")
    print(ttest_collapsed)
    
    cat("\n  Summary:\n")
    cat("  Mean water attribution:", round(mean(water_prop_collapsed, na.rm = TRUE) * 100, 1), "%\n")
    cat("  SD:                    ", round(sd(water_prop_collapsed,   na.rm = TRUE) * 100, 1), "%\n")
    cat("  t =", round(ttest_collapsed$statistic, 3),
        "| df =", ttest_collapsed$parameter,
        "| p =", signif(ttest_collapsed$p.value, 3), "\n")
    
  } else {
    cat("  No proportion elements found in collapsed FEAST result.\n")
    replicate_df_collapsed <- NULL
    ttest_collapsed        <- NULL
  }
  
} else {
  replicate_df_collapsed <- NULL
  ttest_collapsed        <- NULL
}


# ==============================================================================
# 6.6  SUMMARISE OPTIONS A & B
# ==============================================================================
# Collapses per-replicate results into per-context means for reporting.
# Returns a data frame with one row per site x season (x substrate) context.
# ==============================================================================

summarise_feast_results <- function(feast_results_list, label) {
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("FEAST Summary:", label, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n", sep = "")
  
  summary_rows <- list()
  
  for (key in names(feast_results_list)) {
    
    res <- feast_results_list[[key]]
    if (is.null(res)) next
    
    result_mat <- res$result
    n_sinks    <- length(res$sink_samples)
    
    prop_elements <- Filter(function(x) {
      is.numeric(x) && is.vector(x) && length(x) == n_sinks
    }, result_mat)
    
    if (length(prop_elements) == 0) {
      cat("    No proportion elements found for key:", key, "— skipping\n"); next
    }
    
    prop_mat           <- do.call(cbind, prop_elements)
    colnames(prop_mat) <- names(prop_elements)
    rownames(prop_mat) <- res$sink_samples
    
    water_cols  <- grep("Water",   colnames(prop_mat), value = TRUE)
    unknown_col <- grep("Unknown", colnames(prop_mat), value = TRUE)
    
    water_per_sink <- if (length(water_cols) > 0)
      rowSums(prop_mat[, water_cols, drop = FALSE], na.rm = TRUE) else rep(NA, n_sinks)
    
    unknown_per_sink <- if (length(unknown_col) > 0)
      prop_mat[, unknown_col] else rep(NA, n_sinks)
    
    summary_rows[[key]] <- data.frame(
      Context       = key,
      site          = res$sink_site,
      sink_season   = res$sink_season,
      source_season = res$source_season,
      substrate     = ifelse(is.null(res$substrate), "All", res$substrate),
      n_sinks       = n_sinks,
      n_sources     = length(res$source_samples),
      mean_water    = round(mean(water_per_sink,   na.rm = TRUE), 4),
      sd_water      = round(sd(water_per_sink,     na.rm = TRUE), 4),
      se_water      = round(sd(water_per_sink,     na.rm = TRUE) /
                              sqrt(sum(!is.na(water_per_sink))), 4),
      mean_unknown  = round(mean(unknown_per_sink, na.rm = TRUE), 4),
      sd_unknown    = round(sd(unknown_per_sink,   na.rm = TRUE), 4),
      se_unknown    = round(sd(unknown_per_sink,   na.rm = TRUE) /
                              sqrt(sum(!is.na(unknown_per_sink))), 4),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(summary_rows) == 0) { cat("No results to summarise.\n"); return(NULL) }
  
  summary_df <- bind_rows(summary_rows)
  cat("Results:\n")
  print(summary_df, row.names = FALSE)
  return(summary_df)
}

summary_concurrent    <- summarise_feast_results(feast_results_concurrent,
                                                 "Concurrent (same season)")
summary_lagged        <- summarise_feast_results(feast_results_lagged,
                                                 "Lagged (preceding season)")
summary_per_substrate <- summarise_feast_results(feast_results_per_substrate,
                                                 "Per substrate (concurrent)")


# ==============================================================================
# 6.7  STATISTICAL TESTS — BETA REGRESSION (Options A & B)
# ==============================================================================
# Tests whether water attribution differs significantly by season, site,
# or substrate. Uses beta regression because the outcome is a proportion
# (bounded 0–1). Exact 0s and 1s are nudged to 0.001/0.999 as beta
# regression is undefined at the boundary.
#
# Note on sample size: with only 3 replicates per context, these models
# have limited power. Results should be interpreted as exploratory.
# ==============================================================================

library(glmmTMB)
library(emmeans)
library(car)

# --- Helper: extract per-replicate proportions from a FEAST results list ---
extract_replicate_df <- function(feast_results_list, meta_raw) {
  rows <- list()
  for (key in names(feast_results_list)) {
    res <- feast_results_list[[key]]
    if (is.null(res)) next
    n_sinks <- length(res$sink_samples)
    prop_elements <- Filter(function(x) is.numeric(x) && is.vector(x) &&
                              length(x) == n_sinks, res$result)
    if (length(prop_elements) == 0) next
    prop_mat           <- do.call(cbind, prop_elements)
    colnames(prop_mat) <- names(prop_elements)
    rownames(prop_mat) <- res$sink_samples
    water_cols <- grep("Water", colnames(prop_mat), value = TRUE)
    water_prop <- if (length(water_cols) > 0)
      rowSums(prop_mat[, water_cols, drop = FALSE], na.rm = TRUE) else rep(NA, n_sinks)
    sink_meta <- meta_raw[res$sink_samples, c("site","season","substrate","bio_replicate")]
    rows[[key]] <- data.frame(
      sample      = res$sink_samples,
      site        = sink_meta$site,
      sink_season = res$sink_season,
      substrate   = if (!is.null(res$substrate)) res$substrate else sink_meta$substrate,
      bio_rep     = sink_meta$bio_replicate,
      water_prop  = water_prop,
      stringsAsFactors = FALSE
    )
  }
  df <- bind_rows(rows) %>%
    mutate(
      water_prop_adj = case_when(
        water_prop == 0 ~ 0.001,
        water_prop == 1 ~ 0.999,
        TRUE            ~ water_prop
      ),
      sink_season = factor(sink_season, levels = c("Winter","Spring","Summer","Fall","Winter2")),
      substrate   = factor(substrate),
      site        = factor(site)
    )
  return(df)
}

replicate_df           <- extract_replicate_df(feast_results_concurrent,    meta_raw)
replicate_df_substrate <- extract_replicate_df(feast_results_per_substrate, meta_raw)

# --- Model 1: Concurrent — season + site effects ---
# Note: the random effect (1|site:sink_season) accounts for non-independence
# of replicates within the same context. With only 3 replicates per context
# this is a conservative approach; interpret p-values cautiously.
m_concurrent <- glmmTMB(
  water_prop_adj ~ sink_season + site + (1 | site:sink_season),
  data   = replicate_df %>% filter(!is.na(water_prop_adj)),
  family = beta_family(link = "logit")
)

cat("\n--- Model 1: Concurrent (season + site) ---\n")
print(car::Anova(m_concurrent, type = "III"))

cat("\nPost-hoc: season\n")
print(pairs(emmeans(m_concurrent, ~ sink_season), adjust = "tukey"))

cat("\nPost-hoc: site\n")
print(pairs(emmeans(m_concurrent, ~ site), adjust = "tukey"))


# --- Model 2: Per substrate — substrate + season + site effects ---
m_substrate <- glmmTMB(
  water_prop_adj ~ substrate + sink_season + site + (1 | site:sink_season),
  data   = replicate_df_substrate %>% filter(!is.na(water_prop_adj)),
  family = beta_family(link = "logit")
)

cat("\n--- Model 2: Per substrate (substrate + season + site) ---\n")
print(car::Anova(m_substrate, type = "III"))

cat("\nPost-hoc: substrate\n")
print(pairs(emmeans(m_substrate, ~ substrate), adjust = "tukey"))


# ==============================================================================
# 6.8  EXPORT
# ==============================================================================

# --- ANOVA tables ---
anova_concurrent <- car::Anova(m_concurrent, type = "III") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Term") %>%
  mutate(Model = "Concurrent (no substrate)")

anova_substrate <- car::Anova(m_substrate, type = "III") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Term") %>%
  mutate(Model = "Per substrate")

anova_combined <- bind_rows(anova_concurrent, anova_substrate) %>%
  select(Model, Term, Chisq, Df, `Pr(>Chisq)`) %>%
  mutate(Chisq = round(Chisq, 3), `Pr(>Chisq)` = round(`Pr(>Chisq)`, 4))

print(anova_combined)
write.csv(anova_combined,
          file.path(FEAST_DIR, "Supplementary_FEAST_stats.csv"),
          row.names = FALSE)

# --- Collapsed t-test result ---
if (!is.null(ttest_collapsed)) {
  ttest_export <- data.frame(
    test           = "One-sample two-sided t-test",
    hypothesis     = "Mean water attribution != 0",
    mean_water_pct = round(mean(water_prop_collapsed, na.rm = TRUE) * 100, 2),
    sd_water_pct   = round(sd(water_prop_collapsed,   na.rm = TRUE) * 100, 2),
    t_statistic    = round(ttest_collapsed$statistic, 4),
    df             = ttest_collapsed$parameter,
    p_value        = signif(ttest_collapsed$p.value, 4),
    ci_lower_pct   = round(ttest_collapsed$conf.int[1] * 100, 2),
    ci_upper_pct   = round(ttest_collapsed$conf.int[2] * 100, 2)
  )
  write.csv(ttest_export,
            file.path(FEAST_DIR, "Supplementary_FEAST_ttest_collapsed.csv"),
            row.names = FALSE)
}

# --- Full summary table (Options A, A-lagged, B) ---
supp_table <- bind_rows(
  summary_concurrent    %>% mutate(Mode = "Concurrent"),
  summary_lagged        %>% mutate(Mode = "Lagged"),
  summary_per_substrate %>% mutate(Mode = "Per substrate")
) %>%
  mutate(
    `Water % (mean +/- SE)`   = paste0(round(mean_water   * 100, 1), " +/- ",
                                       round(se_water     * 100, 1)),
    `Unknown % (mean +/- SE)` = paste0(round(mean_unknown * 100, 1), " +/- ",
                                       round(se_unknown   * 100, 1)),
    sink_season = factor(sink_season, levels = c("Winter","Spring","Summer","Fall","Winter2"))
  ) %>%
  arrange(Mode, site, sink_season, substrate) %>%
  select(
    Mode,
    Site            = site,
    `Sink season`   = sink_season,
    `Source season` = source_season,
    Substrate       = substrate,
    `N sinks`       = n_sinks,
    `N sources`     = n_sources,
    `Water % (mean +/- SE)`,
    `Unknown % (mean +/- SE)`
  )

write.csv(supp_table,
          file.path(FEAST_DIR, "Supplementary_FEAST_table.csv"),
          row.names = FALSE)

# --- Save all R objects ---
saveRDS(
  list(
    feast_results_concurrent    = feast_results_concurrent,
    feast_results_lagged        = feast_results_lagged,
    feast_results_per_substrate = feast_results_per_substrate,
    feast_result_collapsed      = feast_result_collapsed,
    replicate_df_collapsed      = replicate_df_collapsed,
    ttest_collapsed             = ttest_collapsed,
    summary_concurrent          = summary_concurrent,
    summary_lagged              = summary_lagged,
    summary_per_substrate       = summary_per_substrate
  ),
  file.path(FEAST_DIR, "FEAST_all_results.rds")
)



cat("\n=== FEAST COMPLETE ===\n")


save.image("Combined_ST_WF_analysis_collapsed.RData")

summarise_for_paper <- function(summary_df, label) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat(label, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  # Overall summary
  cat("--- OVERALL ---\n")
  cat("Water attribution:\n")
  cat("  Mean:   ", round(mean(summary_df$mean_water, na.rm = TRUE) * 100, 1), "%\n")
  cat("  Median: ", round(median(summary_df$mean_water, na.rm = TRUE) * 100, 1), "%\n")
  cat("  Min:    ", round(min(summary_df$mean_water, na.rm = TRUE) * 100, 1), "%\n")
  cat("  Max:    ", round(max(summary_df$mean_water, na.rm = TRUE) * 100, 1), "%\n\n")
  
  cat("Unknown attribution:\n")
  cat("  Mean:   ", round(mean(summary_df$mean_unknown, na.rm = TRUE) * 100, 1), "%\n")
  cat("  Median: ", round(median(summary_df$mean_unknown, na.rm = TRUE) * 100, 1), "%\n")
  cat("  Min:    ", round(min(summary_df$mean_unknown, na.rm = TRUE) * 100, 1), "%\n")
  cat("  Max:    ", round(max(summary_df$mean_unknown, na.rm = TRUE) * 100, 1), "%\n\n")
  
  # By site
  cat("--- BY SITE ---\n")
  site_summary <- aggregate(cbind(mean_water, mean_unknown) ~ site,
                            data = summary_df,
                            FUN = function(x) round(mean(x) * 100, 1))
  print(site_summary)
  cat("\n")
  
  # By season
  cat("--- BY SEASON ---\n")
  season_summary <- aggregate(cbind(mean_water, mean_unknown) ~ sink_season,
                              data = summary_df,
                              FUN = function(x) round(mean(x) * 100, 1))
  print(season_summary)
  cat("\n")
  
  # By substrate (only if present and meaningful)
  if ("substrate" %in% colnames(summary_df) &&
      length(unique(na.omit(summary_df$substrate))) > 1) {
    cat("--- BY SUBSTRATE ---\n")
    sub_summary <- aggregate(cbind(mean_water, mean_unknown) ~ substrate,
                             data = summary_df,
                             FUN = function(x) round(mean(x) * 100, 1))
    print(sub_summary)
    cat("\n")
  }
  
  # By site x season
  cat("--- BY SITE × SEASON ---\n")
  context_summary <- aggregate(cbind(mean_water, mean_unknown) ~ site + sink_season,
                               data = summary_df,
                               FUN = function(x) round(mean(x) * 100, 1))
  print(context_summary)
  cat("\n")
}

summarise_for_paper(summary_concurrent,    "Concurrent (same season)")
summarise_for_paper(summary_lagged,        "Lagged (preceding season)")
summarise_for_paper(summary_per_substrate, "Per substrate (concurrent)")
